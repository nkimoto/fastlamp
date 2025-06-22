#!/usr/bin/env python3

"""Copyright (c) 2013, LAMP development team
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

# Convert the expression to 1 or 0 according to the gene expression profile.
# If the gene profile > threshold, the expression is 1, otherwise 0.
# @author A. Terada

__author__ = "Aika TERADA"

import os
import sys
import argparse

NAME_COLUMN = 0 # gene name column number
EXPRESSION_COLUMN = 1 # gene expression profile number
SEPARATOR = "," # file separator


##
# Read expression file
# exp_file: The expression-file name
# return: the list of tuple:
#    tuple[0]: gene name
#    tuple[2]: expression profile
##
def readExpFile( exp_file ):
    exp_list = []
    try:
        with open( exp_file, 'r', encoding='utf-8' ) as f:
            for line in f:
                s = line.strip().split(SEPARATOR)
                exp_list.append( tuple( [s[NAME_COLUMN], float(s[EXPRESSION_COLUMN])]) )
        return exp_list
    except OSError as e:
        sys.stderr.write(f"Error in read {exp_file}\n")
        raise e

##
# Classify the gene expression profile to 1 or 0.
# exp_list: the list of tuple (obtained by readExpFile)
# threshold: threshold to classify the expression level.
#    tuple[0]: gene name
#    tuple[1]: expression profile
# return: List of tuple
#    tuple[0]: gene name
#    tuple[1]: 1 or 0 according to gene expression profile
##
def deriveFlag( exp_list, threshold ):
    flag_list = []
    flag1 = 0
    for t in exp_list:
        # If exp > threshold the expressions to 1, otherwise 0
        if t[1] > threshold:
            flag_list.append(tuple([t[0], 1]))
            flag1 = flag1 + 1
        else:
            flag_list.append(tuple([t[0], 0]))
    return flag_list, flag1

def output( flag_list, output_file ):
    try:
        with open( output_file, "w", encoding='utf-8' ) as f:
            for t in flag_list:
                f.write(f"{t[0]}{SEPARATOR}{t[1]}\n")
    except OSError as e:
        sys.stderr.write(f"Error during output the result to {output_file}")
        raise e

def run( exp_file, threshold, output_file ):
    exp_list = readExpFile( exp_file )
    flag_list, flag1 = deriveFlag( exp_list, threshold )
    output( flag_list, output_file )
    print(f"{len(flag_list)} genes are included.")
    print(f"{flag1} genes are classified to expression level 1.")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Convert the expression to 1 or 0 according to the gene expression profile.")
    p.add_argument("expression_file", help="expression-file")
    p.add_argument("threshold", type=float, help="threshold to classify the expression level")
    p.add_argument("output_file", help="output-file")
    
    args = p.parse_args()

    # Check the file exist.
    if not os.path.isfile(args.expression_file):
        sys.stderr.write(f"IOError: No such file: '{args.expression_file}'\n")
        sys.exit(1)

    try:
        run( args.expression_file, args.threshold, args.output_file )
    except (OSError, IOError, ValueError) as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)
