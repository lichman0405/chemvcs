"""Custom exceptions for ChemVCS Python layer."""


class ChemVCSError(Exception):
    """Base exception for ChemVCS errors."""
    pass


class RepositoryNotFoundError(ChemVCSError):
    """Raised when a ChemVCS repository cannot be found."""
    pass


class ObjectNotFoundError(ChemVCSError):
    """Raised when a requested object does not exist."""
    pass


class ParseError(ChemVCSError):
    """Raised when parsing a file format fails."""
    pass


class ValidationError(ChemVCSError):
    """Raised when data validation fails."""
    pass


class CLIError(ChemVCSError):
    """Raised when chemvcs CLI command fails."""
    pass
