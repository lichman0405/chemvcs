package objectstore

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"
)

// GCOptions controls garbage collection behavior.
type GCOptions struct {
	// Prune objects unreachable for this duration (0 = immediate)
	PruneAge time.Duration
	// DryRun only reports what would be deleted
	DryRun bool
	// Verbose enables detailed output
	Verbose bool
}

// GCStats contains statistics about a GC run.
type GCStats struct {
	ReachableObjects   int
	UnreachableObjects int
	RemovedObjects     int
	ReclaimedBytes     int64
	Duration           time.Duration
}

// GarbageCollect performs mark-and-sweep garbage collection.
// It identifies unreachable objects and optionally removes them.
func (s *Store) GarbageCollect(refsLister RefsLister, opts GCOptions) (*GCStats, error) {
	start := time.Now()
	stats := &GCStats{}

	if opts.Verbose {
		fmt.Println("Garbage collection starting...")
	}

	// Phase 1: Mark - find all reachable objects
	reachable, err := s.markReachable(refsLister, opts.Verbose)
	if err != nil {
		return nil, fmt.Errorf("mark phase failed: %w", err)
	}
	stats.ReachableObjects = len(reachable)

	if opts.Verbose {
		fmt.Printf("Found %d reachable objects\n", stats.ReachableObjects)
	}

	// Phase 2: Sweep - identify unreachable objects
	unreachable, err := s.findUnreachable(reachable)
	if err != nil {
		return nil, fmt.Errorf("sweep phase failed: %w", err)
	}
	stats.UnreachableObjects = len(unreachable)

	if opts.Verbose {
		fmt.Printf("Found %d unreachable objects\n", stats.UnreachableObjects)
	}

	// Phase 3: Prune - remove old unreachable objects
	if !opts.DryRun && len(unreachable) > 0 {
		removed, reclaimed, err := s.pruneObjects(unreachable, opts.PruneAge, opts.Verbose)
		if err != nil {
			return nil, fmt.Errorf("prune phase failed: %w", err)
		}
		stats.RemovedObjects = removed
		stats.ReclaimedBytes = reclaimed
	}

	stats.Duration = time.Since(start)

	if opts.Verbose {
		if opts.DryRun {
			fmt.Printf("Dry run: would remove %d objects\n", stats.UnreachableObjects)
		} else {
			fmt.Printf("Removed %d objects (%.2f MB)\n", stats.RemovedObjects, float64(stats.ReclaimedBytes)/(1024*1024))
		}
		fmt.Printf("GC completed in %v\n", stats.Duration)
	}

	return stats, nil
}

// RefsLister provides access to all refs in the repository.
type RefsLister interface {
	ListAllRefs() ([]string, error)
}

// markReachable performs BFS traversal from all refs to mark reachable objects.
func (s *Store) markReachable(refsLister RefsLister, verbose bool) (map[string]bool, error) {
	reachable := make(map[string]bool)
	queue := []string{}

	// Collect all refs as starting points
	refs, err := refsLister.ListAllRefs()
	if err != nil {
		return nil, fmt.Errorf("failed to list refs: %w", err)
	}

	for _, ref := range refs {
		if ref != "" {
			queue = append(queue, ref)
		}
	}

	if verbose {
		fmt.Printf("Starting from %d refs\n", len(queue))
	}

	// BFS traversal
	visited := 0
	for len(queue) > 0 {
		hash := queue[0]
		queue = queue[1:]

		if reachable[hash] {
			continue
		}
		reachable[hash] = true
		visited++

		if verbose && visited%1000 == 0 {
			fmt.Printf("Marked %d objects...\n", visited)
		}

		// Try to read as object (tree) and collect refs
		obj, err := s.GetObject(hash)
		if err == nil {
			for _, ref := range obj.Refs {
				if ref.ID != "" && !reachable[ref.ID] {
					queue = append(queue, ref.ID)
				}
			}
			continue
		}

		// Try as snapshot
		snap, err := s.GetSnapshot(hash)
		if err == nil {
			if snap.Root != "" && !reachable[snap.Root] {
				queue = append(queue, snap.Root)
			}
			for _, parent := range snap.Parents {
				if parent != "" && !reachable[parent] {
					queue = append(queue, parent)
				}
			}
			continue
		}

		// It's a blob or inaccessible - no refs to follow
	}

	return reachable, nil
}

// findUnreachable identifies all loose objects not in the reachable set.
func (s *Store) findUnreachable(reachable map[string]bool) ([]string, error) {
	var unreachable []string

	// Walk all shard directories
	for i := 0; i < 256; i++ {
		shard := fmt.Sprintf("%02x", i)
		shardPath := filepath.Join(s.root, shard)

		if _, err := os.Stat(shardPath); os.IsNotExist(err) {
			continue
		}

		entries, err := os.ReadDir(shardPath)
		if err != nil {
			return nil, fmt.Errorf("failed to read shard %s: %w", shard, err)
		}

		for _, entry := range entries {
			if entry.IsDir() {
				continue
			}

			hash := entry.Name()
			if !reachable[hash] {
				unreachable = append(unreachable, hash)
			}
		}
	}

	return unreachable, nil
}

// pruneObjects removes unreachable objects older than the given age.
func (s *Store) pruneObjects(unreachable []string, pruneAge time.Duration, verbose bool) (int, int64, error) {
	removed := 0
	var reclaimed int64

	cutoff := time.Now().Add(-pruneAge)

	for _, hash := range unreachable {
		path := s.objectPath(hash)

		// Check file age
		info, err := os.Stat(path)
		if err != nil {
			if os.IsNotExist(err) {
				continue // Already gone
			}
			return removed, reclaimed, fmt.Errorf("failed to stat %s: %w", hash, err)
		}

		// Skip if too recent (grace period)
		if info.ModTime().After(cutoff) {
			if verbose {
				fmt.Printf("Skipping recent object %s (age: %v)\n", hash[:8], time.Since(info.ModTime()))
			}
			continue
		}

		// Remove the object
		size := info.Size()
		if err := os.Remove(path); err != nil {
			return removed, reclaimed, fmt.Errorf("failed to remove %s: %w", hash, err)
		}

		removed++
		reclaimed += size

		if verbose && removed%100 == 0 {
			fmt.Printf("Removed %d objects (%.2f MB)...\n", removed, float64(reclaimed)/(1024*1024))
		}
	}

	// Clean up empty shard directories
	for i := 0; i < 256; i++ {
		shard := fmt.Sprintf("%02x", i)
		shardPath := filepath.Join(s.root, shard)

		entries, err := os.ReadDir(shardPath)
		if err != nil || len(entries) > 0 {
			continue
		}

		// Directory is empty, remove it
		os.Remove(shardPath)
	}

	return removed, reclaimed, nil
}

// PackLooseObjects packs loose objects into a packfile.
// Returns the name of the created pack (without extension).
func (s *Store) PackLooseObjects(refsLister RefsLister, opts PackOptions) (string, *PackStats, error) {
	start := time.Now()
	stats := &PackStats{}

	if opts.Verbose {
		fmt.Println("Packing loose objects...")
	}

	// Determine which objects to pack
	var toPack []string

	if opts.All {
		// Pack all loose objects
		for i := 0; i < 256; i++ {
			shard := fmt.Sprintf("%02x", i)
			shardPath := filepath.Join(s.root, shard)

			if _, err := os.Stat(shardPath); os.IsNotExist(err) {
				continue
			}

			entries, err := os.ReadDir(shardPath)
			if err != nil {
				return "", nil, fmt.Errorf("failed to read shard %s: %w", shard, err)
			}

			for _, entry := range entries {
				if !entry.IsDir() {
					toPack = append(toPack, entry.Name())
				}
			}
		}
	} else {
		// Pack only reachable objects
		reachable, err := s.markReachable(refsLister, opts.Verbose)
		if err != nil {
			return "", nil, fmt.Errorf("failed to mark reachable objects: %w", err)
		}

		for hash := range reachable {
			// Only pack if it's a loose object (not already in a pack)
			if exists, _ := s.exists(s.objectPath(hash)); exists {
				toPack = append(toPack, hash)
			}
		}
	}

	if len(toPack) == 0 {
		if opts.Verbose {
			fmt.Println("No loose objects to pack")
		}
		return "", stats, nil
	}

	stats.ObjectsConsidered = len(toPack)

	if opts.Verbose {
		fmt.Printf("Found %d loose objects to pack\n", len(toPack))
	}

	// Create packfile
	packName := fmt.Sprintf("pack-%d-%d", time.Now().Unix(), os.Getpid())
	pw, err := NewPackWriter(s.root, packName)
	if err != nil {
		return "", nil, fmt.Errorf("failed to create pack writer: %w", err)
	}

	// Add objects to pack
	for i, hash := range toPack {
		data, objType, err := s.readObjectForPack(hash)
		if err != nil {
			if opts.Verbose {
				fmt.Printf("Warning: failed to read object %s: %v\n", hash[:8], err)
			}
			continue
		}

		if err := pw.AddObject(hash, objType, data); err != nil {
			return "", nil, fmt.Errorf("failed to add object %s to pack: %w", hash, err)
		}

		stats.ObjectsPacked++
		stats.BytesBeforePack += int64(len(data))

		if opts.Verbose && (i+1)%1000 == 0 {
			fmt.Printf("Packed %d/%d objects...\n", i+1, len(toPack))
		}
	}

	// Finalize pack
	if err := pw.Finalize(); err != nil {
		return "", nil, fmt.Errorf("failed to finalize pack: %w", err)
	}

	// Get pack file size
	packPath := filepath.Join(s.root, "pack", packName+".pack")
	if info, err := os.Stat(packPath); err == nil {
		stats.BytesAfterPack = info.Size()
	}

	stats.Duration = time.Since(start)

	if opts.Verbose {
		ratio := 100.0 * (1.0 - float64(stats.BytesAfterPack)/float64(stats.BytesBeforePack))
		fmt.Printf("Pack created: %s\n", packName)
		fmt.Printf("Packed %d objects\n", stats.ObjectsPacked)
		fmt.Printf("Size: %.2f MB → %.2f MB (%.1f%% savings)\n",
			float64(stats.BytesBeforePack)/(1024*1024),
			float64(stats.BytesAfterPack)/(1024*1024),
			ratio)
		fmt.Printf("Completed in %v\n", stats.Duration)
	}

	// Reload packs to include the new one
	if err := s.ReloadPacks(); err != nil {
		fmt.Fprintf(os.Stderr, "Warning: failed to reload packs: %v\n", err)
	}

	// Optionally remove packed loose objects
	if !opts.KeepLoose {
		removed := 0
		for _, hash := range toPack {
			path := s.objectPath(hash)
			if err := os.Remove(path); err == nil {
				removed++
			}
		}
		if opts.Verbose {
			fmt.Printf("Removed %d loose objects\n", removed)
		}
	}

	return packName, stats, nil
}

// PackOptions controls pack behavior.
type PackOptions struct {
	// All packs all loose objects (not just reachable)
	All bool
	// KeepLoose keeps loose objects after packing
	KeepLoose bool
	// Verbose enables detailed output
	Verbose bool
}

// PackStats contains statistics about a pack operation.
type PackStats struct {
	ObjectsConsidered int
	ObjectsPacked     int
	BytesBeforePack   int64
	BytesAfterPack    int64
	Duration          time.Duration
}

// readObjectForPack reads an object and determines its type for packing.
func (s *Store) readObjectForPack(hash string) (data []byte, objType byte, err error) {
	path := s.objectPath(hash)
	data, err = os.ReadFile(path)
	if err != nil {
		return nil, 0, err
	}

	// Try to parse as JSON object/snapshot to determine type
	// This is a heuristic - we check for known fields
	dataStr := string(data)
	if strings.Contains(dataStr, `"version"`) && strings.Contains(dataStr, `"type"`) {
		objType = objTypeTree
	} else if strings.Contains(dataStr, `"tree"`) && strings.Contains(dataStr, `"author"`) {
		objType = objTypeSnapshot
	} else {
		objType = objTypeBlob
	}

	return data, objType, nil
}
