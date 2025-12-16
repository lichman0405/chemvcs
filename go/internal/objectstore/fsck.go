package objectstore

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
)

// FsckOptions controls fsck behavior.
type FsckOptions struct {
	// Full performs comprehensive checks including pack verification
	Full bool
	// Verbose enables detailed output
	Verbose bool
}

// FsckIssue represents a problem found during fsck.
type FsckIssue struct {
	Severity string // "error", "warning", "info"
	Object   string // Hash of affected object (if applicable)
	Message  string
}

// FsckReport contains the results of an fsck run.
type FsckReport struct {
	ObjectsChecked int
	PacksChecked   int
	Issues         []FsckIssue
	HasErrors      bool
}

// Fsck performs integrity checks on the object store.
func (s *Store) Fsck(refsLister RefsLister, opts FsckOptions) (*FsckReport, error) {
	report := &FsckReport{
		Issues: []FsckIssue{},
	}

	if opts.Verbose {
		fmt.Println("Starting filesystem check...")
	}

	// Check 1: Verify all loose objects
	if err := s.checkLooseObjects(report, opts.Verbose); err != nil {
		return nil, fmt.Errorf("failed to check loose objects: %w", err)
	}

	// Check 2: Verify packfiles (if Full mode)
	if opts.Full {
		if err := s.checkPackfiles(report, opts.Verbose); err != nil {
			return nil, fmt.Errorf("failed to check packfiles: %w", err)
		}
	}

	// Check 3: Verify reference integrity
	if refsLister != nil {
		if err := s.checkReferences(refsLister, report, opts.Verbose); err != nil {
			return nil, fmt.Errorf("failed to check references: %w", err)
		}
	}

	if opts.Verbose {
		fmt.Printf("\nFilesystem check completed:\n")
		fmt.Printf("  Objects checked: %d\n", report.ObjectsChecked)
		if opts.Full {
			fmt.Printf("  Packs checked: %d\n", report.PacksChecked)
		}
		fmt.Printf("  Issues found: %d\n", len(report.Issues))

		// Summarize by severity
		errors, warnings, infos := 0, 0, 0
		for _, issue := range report.Issues {
			switch issue.Severity {
			case "error":
				errors++
			case "warning":
				warnings++
			case "info":
				infos++
			}
		}
		if errors > 0 {
			fmt.Printf("    Errors: %d\n", errors)
		}
		if warnings > 0 {
			fmt.Printf("    Warnings: %d\n", warnings)
		}
		if infos > 0 {
			fmt.Printf("    Info: %d\n", infos)
		}
	}

	return report, nil
}

// checkLooseObjects verifies hash integrity of all loose objects.
func (s *Store) checkLooseObjects(report *FsckReport, verbose bool) error {
	if verbose {
		fmt.Println("Checking loose objects...")
	}

	checked := 0
	for i := 0; i < 256; i++ {
		shard := fmt.Sprintf("%02x", i)
		shardPath := filepath.Join(s.root, shard)

		if _, err := os.Stat(shardPath); os.IsNotExist(err) {
			continue
		}

		entries, err := os.ReadDir(shardPath)
		if err != nil {
			return fmt.Errorf("failed to read shard %s: %w", shard, err)
		}

		for _, entry := range entries {
			if entry.IsDir() {
				continue
			}

			hash := entry.Name()
			path := filepath.Join(shardPath, hash)

			// Read object data
			data, err := os.ReadFile(path)
			if err != nil {
				report.Issues = append(report.Issues, FsckIssue{
					Severity: "error",
					Object:   hash,
					Message:  fmt.Sprintf("Cannot read object: %v", err),
				})
				report.HasErrors = true
				continue
			}

			// Verify hash
			computed := sha256.Sum256(data)
			computedHash := hex.EncodeToString(computed[:])

			if computedHash != hash {
				report.Issues = append(report.Issues, FsckIssue{
					Severity: "error",
					Object:   hash,
					Message:  fmt.Sprintf("Hash mismatch: expected %s, got %s", hash[:8], computedHash[:8]),
				})
				report.HasErrors = true
				continue
			}

			// Try to parse as JSON (optional validation)
			var temp interface{}
			if err := json.Unmarshal(data, &temp); err != nil {
				// Might be a binary blob - that's okay
			}

			checked++
			if verbose && checked%1000 == 0 {
				fmt.Printf("Checked %d loose objects...\n", checked)
			}
		}
	}

	report.ObjectsChecked += checked
	if verbose {
		fmt.Printf("Checked %d loose objects\n", checked)
	}

	return nil
}

// checkPackfiles verifies integrity of all packfiles.
func (s *Store) checkPackfiles(report *FsckReport, verbose bool) error {
	packDir := filepath.Join(s.root, "pack")
	if _, err := os.Stat(packDir); os.IsNotExist(err) {
		if verbose {
			fmt.Println("No pack directory found")
		}
		return nil
	}

	entries, err := os.ReadDir(packDir)
	if err != nil {
		return fmt.Errorf("failed to read pack directory: %w", err)
	}

	if verbose {
		fmt.Println("Checking packfiles...")
	}

	// Find all .pack files
	packs := []string{}
	for _, entry := range entries {
		if !entry.IsDir() && filepath.Ext(entry.Name()) == ".pack" {
			packs = append(packs, entry.Name())
		}
	}

	for _, packFile := range packs {
		packPath := filepath.Join(packDir, packFile)
		idxPath := packPath[:len(packPath)-5] + ".idx"

		// Check 1: Index file exists
		if _, err := os.Stat(idxPath); os.IsNotExist(err) {
			report.Issues = append(report.Issues, FsckIssue{
				Severity: "error",
				Object:   packFile,
				Message:  "Missing index file",
			})
			report.HasErrors = true
			continue
		}

		// Check 2: Verify pack file checksum
		if err := s.verifyPackChecksum(packPath); err != nil {
			report.Issues = append(report.Issues, FsckIssue{
				Severity: "error",
				Object:   packFile,
				Message:  fmt.Sprintf("Pack checksum verification failed: %v", err),
			})
			report.HasErrors = true
		}

		// Check 3: Verify each object in pack
		pr, err := OpenPackReader(s.root, packFile[:len(packFile)-5])
		if err != nil {
			report.Issues = append(report.Issues, FsckIssue{
				Severity: "error",
				Object:   packFile,
				Message:  fmt.Sprintf("Cannot open pack: %v", err),
			})
			report.HasErrors = true
			continue
		}

		objectsInPack := 0
		for _, entry := range pr.index.entries {
			hash := entry.hash
			data, objType, err := pr.GetObject(hash)
			if err != nil {
				report.Issues = append(report.Issues, FsckIssue{
					Severity: "error",
					Object:   hash,
					Message:  fmt.Sprintf("Cannot read from pack %s: %v", packFile, err),
				})
				report.HasErrors = true
				continue
			}

			// Verify hash
			computed := sha256.Sum256(data)
			computedHash := hex.EncodeToString(computed[:])

			if computedHash != hash {
				report.Issues = append(report.Issues, FsckIssue{
					Severity: "error",
					Object:   hash,
					Message:  fmt.Sprintf("Hash mismatch in pack %s", packFile),
				})
				report.HasErrors = true
			}

			// Basic type validation
			if objType > 2 {
				report.Issues = append(report.Issues, FsckIssue{
					Severity: "warning",
					Object:   hash,
					Message:  fmt.Sprintf("Unknown object type %d in pack %s", objType, packFile),
				})
			}

			objectsInPack++
		}

		report.ObjectsChecked += objectsInPack
		report.PacksChecked++

		if verbose {
			fmt.Printf("Checked pack %s: %d objects\n", packFile, objectsInPack)
		}
	}

	return nil
}

// verifyPackChecksum verifies the trailing checksum of a packfile.
func (s *Store) verifyPackChecksum(packPath string) error {
	data, err := os.ReadFile(packPath)
	if err != nil {
		return err
	}

	if len(data) < 32 {
		return fmt.Errorf("pack file too short")
	}

	// Last 32 bytes are the checksum
	storedChecksum := data[len(data)-32:]
	content := data[:len(data)-32]

	computed := sha256.Sum256(content)
	computedChecksum := computed[:]

	for i := 0; i < 32; i++ {
		if storedChecksum[i] != computedChecksum[i] {
			return fmt.Errorf("checksum mismatch")
		}
	}

	return nil
}

// checkReferences verifies that all refs point to valid objects.
func (s *Store) checkReferences(refsLister RefsLister, report *FsckReport, verbose bool) error {
	if verbose {
		fmt.Println("Checking references...")
	}

	refs, err := refsLister.ListAllRefs()
	if err != nil {
		return fmt.Errorf("failed to list refs: %w", err)
	}

	for _, hash := range refs {
		if hash == "" {
			continue
		}

		// Check if object exists
		exists := false
		if _, err := s.GetBlob(hash); err == nil {
			exists = true
		} else if _, err := s.GetObject(hash); err == nil {
			exists = true
		} else if _, err := s.GetSnapshot(hash); err == nil {
			exists = true
		}

		if !exists {
			report.Issues = append(report.Issues, FsckIssue{
				Severity: "error",
				Object:   hash,
				Message:  "Ref points to missing object",
			})
			report.HasErrors = true
		}
	}

	if verbose {
		fmt.Printf("Checked %d references\n", len(refs))
	}

	return nil
}

// AddIssue is a helper for adding issues to a report.
func (r *FsckReport) AddIssue(severity, object, message string) {
	r.Issues = append(r.Issues, FsckIssue{
		Severity: severity,
		Object:   object,
		Message:  message,
	})
	if severity == "error" {
		r.HasErrors = true
	}
}
