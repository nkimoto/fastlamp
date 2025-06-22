#!/usr/bin/env python3

# Convert ProbeID to OtherID using Affymetrix annotation file.
# This code convert to EntrezGene (column number is set by TO_COLUMN_ID).
# If a probe maps several IDs, this code selects the first ID in the annotation file.
# @author A. Terada

"""
Copyright (c) 2013, LAMP development team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the LAMP development team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

__author__ = "Aika TERADA"

import sys
import os
from optparse import OptionParser

# Definitions for map-file
ORIG_COLUMN_ID = 0 # Column index to converted ID
TO_COLUMN_ID = 18 # Column index to new ID
MAP_SEPARATOR = '","' # Separator for the line
MAP_SEPARATOR_COLUMN = "///" # Separator for a column.
# Definitions for converted-file
CONV_ORIG_COLUMN_ID = 0 # Prob ID column number
CONV_SEPARATOR = "\t" # Separator for the input file
OUT_CONV_SEPARATOR = "," # Separator for the output file

def readMapFile(filename):
    old2new = {} # the mapping of prob ID for EntrezGene ID
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                line = line[:-1]
                # If line is the comment line, skip
                if (line.startswith('#')):
                    continue
                else:
                    s = line.split(MAP_SEPARATOR)
                    old_id = s[ORIG_COLUMN_ID].strip('"') # String of an old ID
                    new_ids = list( map(str.strip, s[TO_COLUMN_ID].strip('"').split('///'))) # List of new IDs
                    if not (new_ids[0] == "---"):
                        new_ids = list( reversed( new_ids ) )
                        old2new[old_id] = new_ids[0]
        return old2new
    except IOError as e:
        message = "Error in read %s\n" % filename
        sys.stderr.write(message)
        raise e

def convertID( old2new, converted_file, output_file ):
    try:
        total_size = 0
        converted_size = 0
        with open( converted_file, 'r', encoding='utf-8' ) as fr, open( output_file, 'w', encoding='utf-8' ) as fw:
            for line in fr:
                if not ( line.startswith("!") ):
                    total_size += 1
                    line = line[:-1]
                    s = line.split(CONV_SEPARATOR)
                    old_id = s[CONV_ORIG_COLUMN_ID].strip('"')
                    if old_id in old2new:
                        converted_size += 1
                        new_id = old2new[old_id]
                        fw.write("%s" % new_id)
                        for i in s[1:]:
                            fw.write("%s%s" % (OUT_CONV_SEPARATOR, i))
                        fw.write("\n")
        return converted_size, total_size
    except IOError as e:
        message = "Error in convert %s\n" % converted_file
        sys.stderr.write(message)
        raise e

def run( map_file, converted_file, output_file ):
    print("read files ...", file=sys.stderr)
    old2new = readMapFile(map_file)
    print("convert ID...", file=sys.stderr)
    converted_size, total_size = convertID( old2new, converted_file, output_file )
    print("%s/%s IDs were converted" % (converted_size, total_size), file=sys.stderr)

if __name__ == "__main__":
    usage = "usage: %prog map-file converted-file output-file"
    p = OptionParser(usage = usage)
    
    opts, args = p.parse_args()

    if (len(args) < 3):
        message = "Error: input map-file, converted-file and ontput-file.\n"
        sys.stderr.write(message)
        sys.exit(1)
    
    # check the file exist.
    if not os.path.isfile(args[0]):
        message = "IOError: No such file: \'" + args[0] + "\'\n"
        sys.stderr.write(message)
        sys.exit(1)
    if not os.path.isfile(args[1]):
        message = "IOError: No such file: \'" + args[1] + "\'\n"
        sys.stderr.write(message)
        sys.exit(1)
    
    try:
        run(args[0], args[1], args[2])
    except (ValueError, IOError) as e:
        sys.stderr.write("Error: %s\n" % e)
        sys.exit(1)
