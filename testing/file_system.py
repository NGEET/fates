"""Module for filesystem abstraction for I/O"""

from pathlib import Path


class FileSystem:
    """Minimal filesystem abstractions for I/O"""

    def __init__(self, root: Path | str | None = None):
        self.root = Path(root) if root else Path(__file__).resolve().parent

    def exists(self, path: Path) -> bool:
        """Returns whether an input path exists

        Args:
            path (Path): input file path

        Returns:
            bool: path exists
        """
        return path.exists()

    def makedirs(self, path: Path):
        """Creates directory specified in path

        Args:
            path (Path): input path
        """
        path.mkdir(parents=True, exist_ok=True)

    def write_file(self, path: Path, content: str):
        """Writes content to a file

        Args:
            path (Path): path to file
            content (str): content to write
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    def append_file(self, path: Path, content: str):
        """Appends content to a file

        Args:
            path (Path): path to file
            content (str): content to append
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as f:
            f.write(content)

    def read_text(self, path: Path) -> str:
        """Reads and returns text from a file

        Args:
            path (Path): full path to file

        Returns:
            str: file content
        """
        return path.read_text(encoding="utf-8")

    def read_lines(self, path: Path) -> list[str]:
        """Reads and returns lines of a file

        Args:
            path (Path): path to file

        Returns:
            list[str]: lines of file
        """
        return path.read_text(encoding="utf-8").splitlines(keepends=True)

    def write_lines(self, path: Path, content: list[str]):
        """Writes lines to a file

        Args:
            path (Path): path to file
            content (list[str]): content to write
        """
        with path.open("w", encoding="utf-8") as f:
            f.writelines(content)
