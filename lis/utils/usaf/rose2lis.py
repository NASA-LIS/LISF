#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

import sys
import re
import os

if len(sys.argv) != 3:
    print("Usage: ")
    print("rose2lis rose-app.conf out_nml_dir")
else:
    if os.path.isdir(sys.argv[2]) == False:
        os.makedirs(sys.argv[2])

    conf = open(sys.argv[1])
    lines = conf.readlines()
    namelist = []
    filelist = []
    l = 0
    for line in lines:
        if line.startswith("[file:"):
            idx1 = line.index(":")
            idx2 = line.index("]")
            file = line[idx1+1:idx2]
            print("file: %s" % file)
            line2 = lines[l+1]
            if "(" in line2:
                # originally, the name list in parenthsis has been removed.
                #line3 = re.sub(r"\([^)]*\)", "", line2)
                # now, remove "(" itself
                line3 = line2.replace("(", "")
                line3 = line3.replace(")", "")
            else:
                line3 = line2
            line3 = line3.replace("source=", "")
            line3 = line3.replace("namelist:", "")
            if ")" in line3:
                line3 = line3.replace(")", "")
            names = line3.split()
            print(names)
            # now we have file and names, then create name list file
            N = len(names)
            fname = "%s/%s" % (sys.argv[2], file)
            f = open(fname, "w")
            for n in range(0, N):
                k = 0
                key = "[namelist:" + names[n] + "]"
                print("key=%s" % key)
                for row in lines:
                    if row.startswith(key):
                        f.write("&%s\n" % names[n])
                        r = k+1
                        while lines[r] != "\n":
                            row_len = len(lines[r])
                            if lines[r][row_len-2] == ',':
                                f.write("%s\n" % (lines[r].rstrip()))
                                print("%s" % lines[r+1])
                                lines[r+1] = lines[r+1].replace("=", " ")
                                print("%s" % lines[r+1])
                            else:
                                f.write("%s,\n" % (lines[r].rstrip()))
                            r = r + 1
                        f.write("/\n")
                    k = k + 1
            f.close()

        l = l + 1
