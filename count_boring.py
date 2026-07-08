#!/usr/bin/env python3

import pyximport
pyximport.install(language_level=3) # automatic compliation of cython modules
import enumdags

# Settings
MAXSIZE = 7
NUMPROCESSES = 7 # number of parallel processes

# Count total and boring DAGs of size at most MAXSIZE
if __name__ == '__main__':
    for n in range(1, MAXSIZE+1):
        enumdags.countdags(n, NUMPROCESSES)
