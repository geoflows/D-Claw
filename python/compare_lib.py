#!/usr/bin/env python
# encoding: utf-8
r"""
Run diff on all files in two directories, primarily used for comparing
libraries.

This is not quite done yet!  Use at your own risk...
"""

import getopt
import glob
import os
import sys
import time

help_message = """
Run diff on all files in two directories, primarily used for comparing
libraries.

Arguments -
  compare_lib dir1 dir2 [-p --pattern match pattern] [-d --diff diff program]
"""


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[3:], "hp:d:", ["help", "pattern=", "diff="])
        except getopt.error as msg:
            raise Usage(msg)

        # Defaults
        pattern = "*.f"
        diff = "diff"

        # These are required
        dir1 = argv[1]
        dir2 = argv[2]

        # option processing
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-p", "--pattern"):
                pattern = value
            if option in ("-d", "--diff"):
                diff = value

    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg), file=sys.stderr)
        print("\t for help use --help", file=sys.stderr)
        return 2

    os.chdir(dir1)
    print("In directory ", dir1)

    for fname in glob.glob(pattern):
        f2 = dir2 + "/" + fname[:-6] + ".f"
        if os.path.isfile(f2):
            print("comparing %s to %s" % (fname, f2))
            os.system("%s -w %s %s" % (diff, fname, f2))
        else:
            print("no file %s " % f2)
        time.sleep(1)


if __name__ == "__main__":
    sys.exit(main())
