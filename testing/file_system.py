"""Module for filesystem abstraction for test generators"""

from pathlib import Path


class FileSystem:
    """Minimal filesystem abstractions for test generators"""

    def exists(self, path: str) -> bool:
        """Returns whether an input path exists

        Args:
            path (str): input file path

        Returns:
            bool: path exists
        """
        return Path(path).exists()

    def makedirs(self, path: str):
        """Creates directory specified in path

        Args:
            path (str): input path
        """
        Path(path).mkdir(parents=True, exist_ok=True)

    def write_file(self, path: str, content: str):
        """Writes content to a file

        Args:
            path (str): path to file
            content (str): content to write
        """
        full = Path(path)
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text(content, encoding="utf-8")

    def append_file(self, path: str, content: str):
        """Appends content to a file

        Args:
            path (str): path to file
            content (str): content to append
        """
        full = Path(path)
        full.parent.mkdir(parents=True, exist_ok=True)
        with full.open("a", encoding="utf-8") as f:
            f.write(content)

    def read_text(self, path: str) -> str:
        """Reads and returns text from a file

        Args:
            path (str): full path to file

        Returns:
            str: file content
        """
        full = Path(path)
        return full.read_text(encoding="utf-8")

    def read_lines(self, path: str) -> list[str]:
        """Reads and returns lines of a file

        Args:
            path (str): path to file

        Returns:
            list[str]: lines of file
        """
        full = Path(path)
        with full.open("r", encoding="utf-8") as f:
            return f.readlines()

    def write_lines(self, path: str, content: list[str]):
        """Writes lines to a file

        Args:
            path (str): path to file
            content (list[str]): content to write
        """
        full = Path(path)
        with full.open("w", encoding="utf-8") as f:
            return f.writelines(content)
