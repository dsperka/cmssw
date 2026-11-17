#!/usr/bin/env python3
"""Drop files from the page cache, so that the next job reads them from the disk.

No root needed: POSIX_FADV_DONTNEED on files we can read.
usage: python3 drop_page_cache.py file1 [file2 ...]
"""
import os
import sys

for path in sys.argv[1:]:
    fd = os.open(path, os.O_RDONLY)
    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    os.close(fd)
