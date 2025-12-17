package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/lishi/chemvcs/internal/server"
)

const version = "0.7.0"

func main() {
	var (
		port      = flag.Int("port", 8080, "Port to listen on")
		repoRoot  = flag.String("repo-root", "./repos", "Root directory for repositories")
		authTok   = flag.String("auth-token", "", "Bearer token for API access (or set CHEMVCS_SERVER_AUTH_TOKEN)")
		authRepos = flag.String("auth-repos", "", "Comma-separated allowed repos for auth-token (e.g. owner/repo,owner2/repo2) or '*' (or set CHEMVCS_SERVER_AUTH_REPOS)")
		adminTok  = flag.String("admin-token", "", "Admin bearer token (bypasses repo scoping and allows listing repos) (or set CHEMVCS_SERVER_ADMIN_TOKEN)")
		showVer   = flag.Bool("version", false, "Show version information")
	)

	flag.Parse()

	if *showVer {
		fmt.Printf("chemvcs-server version %s\n", version)
		os.Exit(0)
	}

	// Create server configuration
	if *authTok == "" {
		*authTok = os.Getenv("CHEMVCS_SERVER_AUTH_TOKEN")
	}
	if *authRepos == "" {
		*authRepos = os.Getenv("CHEMVCS_SERVER_AUTH_REPOS")
	}
	if *adminTok == "" {
		*adminTok = os.Getenv("CHEMVCS_SERVER_ADMIN_TOKEN")
	}

	var repos []string
	for _, p := range strings.Split(*authRepos, ",") {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		repos = append(repos, p)
	}

	config := server.Config{
		RepoRoot:   *repoRoot,
		Port:       *port,
		AuthToken:  *authTok,
		AuthRepos:  repos,
		AdminToken: *adminTok,
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
