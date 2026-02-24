"""Main entry point for petsurfer-km group-level analysis"""

from __future__ import annotations

# Use pysqlite3 as a drop-in replacement for sqlite3 if available.
# This is needed for FreeSurfer's fspython, which lacks the _sqlite3 C extension.
# pysqlite3-binary only has wheels for Linux x86-64
#
# Be sure to run `fspython -m pip install pysqlite3-binary` if running from
# inside fspython
try:
    __import__('pysqlite3')
    import sys
    sys.modules['sqlite3'] = sys.modules.pop('pysqlite3')
except ImportError:
    pass

import argparse
import logging
import os
import shutil
import sys
from pathlib import Path

def main(argv: list[str] | None = None) -> None:
    """
    Main entry point for petsurfer-km group-level analysis

    Args:
        argv: Command-line arguments (defaults to sys.argv[1:])
    """
    print("Hello world!")
    sys.exit(0)


if __name__ == "__main__":
    main()
