package objectstore

import (
	"bytes"
	"compress/zlib"
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
)

const (
	// Pack file magic header
	packMagic   = "PACK"
	packVersion = uint32(1)

	// Object type identifiers
	objTypeBlob     = byte(0)
	objTypeTree     = byte(1)
	objTypeSnapshot = byte(2)

	// Index file magic
	idxMagic   = "PIDX"
	idxVersion = uint32(1)
)

// PackWriter creates a packfile containing multiple objects with compression.
type PackWriter struct {
	packPath     string
	idxPath      string
	packFile     *os.File
	entries      []packEntry
	objCount     uint32
	bytesWritten int64
}

type packEntry struct {
	hash   string
	offset int64
	size   int64
}

// NewPackWriter creates a new pack writer.
// packName should be unique (e.g., "pack-<timestamp>-<random>").
func NewPackWriter(objectsRoot, packName string) (*PackWriter, error) {
	packDir := filepath.Join(objectsRoot, "pack")
	if err := os.MkdirAll(packDir, 0755); err != nil {
		return nil, fmt.Errorf("failed to create pack directory: %w", err)
	}

	packPath := filepath.Join(packDir, packName+".pack")
	idxPath := filepath.Join(packDir, packName+".idx")

	packFile, err := os.Create(packPath)
	if err != nil {
		return nil, fmt.Errorf("failed to create pack file: %w", err)
	}

	pw := &PackWriter{
		packPath: packPath,
		idxPath:  idxPath,
		packFile: packFile,
		entries:  []packEntry{},
	}

	// Write pack header
	if err := pw.writeHeader(); err != nil {
		packFile.Close()
		os.Remove(packPath)
		return nil, err
	}

	return pw, nil
}

// writeHeader writes the pack file header.
func (pw *PackWriter) writeHeader() error {
	// Magic + version + placeholder for object count (will be updated in Finalize)
	header := []byte(packMagic)
	if _, err := pw.packFile.Write(header); err != nil {
		return fmt.Errorf("failed to write magic: %w", err)
	}
	pw.bytesWritten += int64(len(header))

	if err := binary.Write(pw.packFile, binary.BigEndian, packVersion); err != nil {
		return fmt.Errorf("failed to write version: %w", err)
	}
	pw.bytesWritten += 4

	if err := binary.Write(pw.packFile, binary.BigEndian, uint32(0)); err != nil {
		return fmt.Errorf("failed to write object count placeholder: %w", err)
	}
	pw.bytesWritten += 4

	return nil
}

// AddObject adds an object to the pack.
// objType should be one of objTypeBlob, objTypeTree, objTypeSnapshot.
// data is the raw uncompressed object data.
func (pw *PackWriter) AddObject(hash string, objType byte, data []byte) error {
	offset := pw.bytesWritten

	// Write object type
	if _, err := pw.packFile.Write([]byte{objType}); err != nil {
		return fmt.Errorf("failed to write object type: %w", err)
	}
	pw.bytesWritten++

	// Write uncompressed size (varint)
	sizeVarint := encodeVarint(uint64(len(data)))
	if _, err := pw.packFile.Write(sizeVarint); err != nil {
		return fmt.Errorf("failed to write size: %w", err)
	}
	pw.bytesWritten += int64(len(sizeVarint))

	// Compress and write data
	var compBuf bytes.Buffer
	zw := zlib.NewWriter(&compBuf)
	if _, err := zw.Write(data); err != nil {
		return fmt.Errorf("failed to compress data: %w", err)
	}
	if err := zw.Close(); err != nil {
		return fmt.Errorf("failed to finalize compression: %w", err)
	}

	compData := compBuf.Bytes()
	if _, err := pw.packFile.Write(compData); err != nil {
		return fmt.Errorf("failed to write compressed data: %w", err)
	}
	pw.bytesWritten += int64(len(compData))

	// Record entry for index
	pw.entries = append(pw.entries, packEntry{
		hash:   hash,
		offset: offset,
		size:   pw.bytesWritten - offset,
	})
	pw.objCount++

	return nil
}

// Finalize writes the pack checksum and generates the index file.
func (pw *PackWriter) Finalize() error {
	// Compute checksum of entire pack content
	if _, err := pw.packFile.Seek(0, io.SeekStart); err != nil {
		return fmt.Errorf("failed to seek to start: %w", err)
	}
	hasher := sha256.New()
	if _, err := io.Copy(hasher, pw.packFile); err != nil {
		return fmt.Errorf("failed to compute checksum: %w", err)
	}
	checksum := hasher.Sum(nil)

	// Write checksum to end of pack
	if _, err := pw.packFile.Seek(0, io.SeekEnd); err != nil {
		return fmt.Errorf("failed to seek to end: %w", err)
	}
	if _, err := pw.packFile.Write(checksum); err != nil {
		return fmt.Errorf("failed to write checksum: %w", err)
	}

	// Update object count in header
	if _, err := pw.packFile.Seek(8, io.SeekStart); err != nil {
		return fmt.Errorf("failed to seek to header: %w", err)
	}
	if err := binary.Write(pw.packFile, binary.BigEndian, pw.objCount); err != nil {
		return fmt.Errorf("failed to update object count: %w", err)
	}

	if err := pw.packFile.Sync(); err != nil {
		return fmt.Errorf("failed to sync pack file: %w", err)
	}
	if err := pw.packFile.Close(); err != nil {
		return fmt.Errorf("failed to close pack file: %w", err)
	}

	// Write index file
	if err := pw.writeIndex(checksum); err != nil {
		os.Remove(pw.packPath)
		return err
	}

	return nil
}

// writeIndex generates the pack index file.
func (pw *PackWriter) writeIndex(packChecksum []byte) error {
	idxFile, err := os.Create(pw.idxPath)
	if err != nil {
		return fmt.Errorf("failed to create index file: %w", err)
	}
	defer idxFile.Close()

	// Sort entries by hash for binary search
	sort.Slice(pw.entries, func(i, j int) bool {
		return pw.entries[i].hash < pw.entries[j].hash
	})

	// Write index header
	if _, err := idxFile.Write([]byte(idxMagic)); err != nil {
		return fmt.Errorf("failed to write index magic: %w", err)
	}
	if err := binary.Write(idxFile, binary.BigEndian, idxVersion); err != nil {
		return fmt.Errorf("failed to write index version: %w", err)
	}
	if err := binary.Write(idxFile, binary.BigEndian, pw.objCount); err != nil {
		return fmt.Errorf("failed to write object count: %w", err)
	}

	// Write fanout table (256 entries)
	fanout := make([]uint32, 256)
	for i := range pw.entries {
		firstByte := pw.entries[i].hash[0]
		if firstByte >= '0' && firstByte <= '9' {
			firstByte = firstByte - '0'
		} else if firstByte >= 'a' && firstByte <= 'f' {
			firstByte = firstByte - 'a' + 10
		} else if firstByte >= 'A' && firstByte <= 'F' {
			firstByte = firstByte - 'A' + 10
		}
		for j := int(firstByte); j < 256; j++ {
			fanout[j]++
		}
	}
	for i := 0; i < 256; i++ {
		if err := binary.Write(idxFile, binary.BigEndian, fanout[i]); err != nil {
			return fmt.Errorf("failed to write fanout: %w", err)
		}
	}

	// Write hash list
	for _, entry := range pw.entries {
		hashBytes, err := hex.DecodeString(entry.hash)
		if err != nil {
			return fmt.Errorf("invalid hash %s: %w", entry.hash, err)
		}
		if _, err := idxFile.Write(hashBytes); err != nil {
			return fmt.Errorf("failed to write hash: %w", err)
		}
	}

	// Write offset list
	for _, entry := range pw.entries {
		if err := binary.Write(idxFile, binary.BigEndian, entry.offset); err != nil {
			return fmt.Errorf("failed to write offset: %w", err)
		}
	}

	// Write pack checksum
	if _, err := idxFile.Write(packChecksum); err != nil {
		return fmt.Errorf("failed to write pack checksum: %w", err)
	}

	return idxFile.Sync()
}

// PackReader reads objects from a packfile.
type PackReader struct {
	packPath string
	idxPath  string
	index    *packIndex
}

type packIndex struct {
	fanout  [256]uint32
	entries []indexEntry
}

type indexEntry struct {
	hash   string
	offset int64
}

// OpenPackReader opens an existing packfile for reading.
func OpenPackReader(objectsRoot, packName string) (*PackReader, error) {
	packDir := filepath.Join(objectsRoot, "pack")
	packPath := filepath.Join(packDir, packName+".pack")
	idxPath := filepath.Join(packDir, packName+".idx")

	idx, err := readIndex(idxPath)
	if err != nil {
		return nil, err
	}

	return &PackReader{
		packPath: packPath,
		idxPath:  idxPath,
		index:    idx,
	}, nil
}

// readIndex reads and parses the pack index file.
func readIndex(idxPath string) (*packIndex, error) {
	idxFile, err := os.Open(idxPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open index file: %w", err)
	}
	defer idxFile.Close()

	// Read and verify header
	magic := make([]byte, 4)
	if _, err := idxFile.Read(magic); err != nil {
		return nil, fmt.Errorf("failed to read index magic: %w", err)
	}
	if string(magic) != idxMagic {
		return nil, fmt.Errorf("invalid index magic: %s", string(magic))
	}

	var version, objCount uint32
	if err := binary.Read(idxFile, binary.BigEndian, &version); err != nil {
		return nil, fmt.Errorf("failed to read index version: %w", err)
	}
	if version != idxVersion {
		return nil, fmt.Errorf("unsupported index version: %d", version)
	}

	if err := binary.Read(idxFile, binary.BigEndian, &objCount); err != nil {
		return nil, fmt.Errorf("failed to read object count: %w", err)
	}

	idx := &packIndex{}

	// Read fanout table
	for i := 0; i < 256; i++ {
		if err := binary.Read(idxFile, binary.BigEndian, &idx.fanout[i]); err != nil {
			return nil, fmt.Errorf("failed to read fanout: %w", err)
		}
	}

	// Read hash list
	hashes := make([]string, objCount)
	for i := uint32(0); i < objCount; i++ {
		hashBytes := make([]byte, 32)
		if _, err := idxFile.Read(hashBytes); err != nil {
			return nil, fmt.Errorf("failed to read hash: %w", err)
		}
		hashes[i] = hex.EncodeToString(hashBytes)
	}

	// Read offset list
	offsets := make([]int64, objCount)
	for i := uint32(0); i < objCount; i++ {
		if err := binary.Read(idxFile, binary.BigEndian, &offsets[i]); err != nil {
			return nil, fmt.Errorf("failed to read offset: %w", err)
		}
	}

	// Build index entries
	idx.entries = make([]indexEntry, objCount)
	for i := uint32(0); i < objCount; i++ {
		idx.entries[i] = indexEntry{
			hash:   hashes[i],
			offset: offsets[i],
		}
	}

	return idx, nil
}

// GetObject retrieves an object from the pack by its hash.
// Returns the raw uncompressed data and object type.
func (pr *PackReader) GetObject(hash string) (data []byte, objType byte, err error) {
	// Binary search in index
	idx := sort.Search(len(pr.index.entries), func(i int) bool {
		return pr.index.entries[i].hash >= hash
	})

	if idx >= len(pr.index.entries) || pr.index.entries[idx].hash != hash {
		return nil, 0, fmt.Errorf("object not found in pack: %s", hash)
	}

	offset := pr.index.entries[idx].offset

	// Open pack file and seek to offset
	packFile, err := os.Open(pr.packPath)
	if err != nil {
		return nil, 0, fmt.Errorf("failed to open pack file: %w", err)
	}
	defer packFile.Close()

	if _, err := packFile.Seek(offset, io.SeekStart); err != nil {
		return nil, 0, fmt.Errorf("failed to seek to offset: %w", err)
	}

	// Read object type
	objTypeBuf := make([]byte, 1)
	if _, err := packFile.Read(objTypeBuf); err != nil {
		return nil, 0, fmt.Errorf("failed to read object type: %w", err)
	}
	objType = objTypeBuf[0]

	// Read size (varint)
	size, err := readVarint(packFile)
	if err != nil {
		return nil, 0, fmt.Errorf("failed to read size: %w", err)
	}

	// Decompress data
	zr, err := zlib.NewReader(packFile)
	if err != nil {
		return nil, 0, fmt.Errorf("failed to create decompressor: %w", err)
	}
	defer zr.Close()

	data = make([]byte, size)
	if _, err := io.ReadFull(zr, data); err != nil {
		return nil, 0, fmt.Errorf("failed to decompress data: %w", err)
	}

	return data, objType, nil
}

// HasObject checks if the pack contains an object with the given hash.
func (pr *PackReader) HasObject(hash string) bool {
	idx := sort.Search(len(pr.index.entries), func(i int) bool {
		return pr.index.entries[i].hash >= hash
	})
	return idx < len(pr.index.entries) && pr.index.entries[idx].hash == hash
}

// ListHashes returns all object hashes in this pack.
func (pr *PackReader) ListHashes() []string {
	hashes := make([]string, len(pr.index.entries))
	for i, entry := range pr.index.entries {
		hashes[i] = entry.hash
	}
	return hashes
}

// encodeVarint encodes an unsigned integer as a variable-length integer.
func encodeVarint(n uint64) []byte {
	buf := make([]byte, 0, 10)
	for n >= 0x80 {
		buf = append(buf, byte(n)|0x80)
		n >>= 7
	}
	buf = append(buf, byte(n))
	return buf
}

// readVarint reads a variable-length integer from a reader.
func readVarint(r io.Reader) (uint64, error) {
	var n uint64
	var shift uint
	for {
		b := make([]byte, 1)
		if _, err := r.Read(b); err != nil {
			return 0, err
		}
		n |= uint64(b[0]&0x7f) << shift
		if b[0]&0x80 == 0 {
			break
		}
		shift += 7
		if shift >= 64 {
			return 0, fmt.Errorf("varint overflow")
		}
	}
	return n, nil
}
