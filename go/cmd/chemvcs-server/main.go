package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/lishi/chemvcs/internal/server"
)

const version = "0.1.0-dev"

func main() {
	var (
		port     = flag.Int("port", 8080, "Port to listen on")
		repoRoot = flag.String("repo-root", "./repos", "Root directory for repositories")
		showVer  = flag.Bool("version", false, "Show version information")
	)

	flag.Parse()

	if *showVer {
		fmt.Printf("chemvcs-server version %s\n", version)
		os.Exit(0)
	}

	// Create server configuration
	config := server.Config{
		RepoRoot: *repoRoot,
		Port:     *port,
	}

	// Create server
	srv, err := server.NewServer(config)
	if err != nil {
		log.Fatalf("Failed to create server: %v", err)
	}

	// Start server
	log.Printf("Starting ChemVCS server...")
	if err := srv.Start(*port); err != nil {
		log.Fatalf("Server error: %v", err)
	}
}
