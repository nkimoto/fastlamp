#!/usr/bin/env python3

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
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT of THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from functions import functions4fisher, functions4u_test, functions4chi
import sys
import os

pardir = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
sys.path.append(pardir)
sys.path.append(pardir + '/functions')

__author__ = "Takayuki ITOH"

motifRpvalue = []
motifApvalue = []
motifNgenes = []
motifSscore = []
motifName = []

combiRank = []
combiRpvalue = []
combiApvalue = []
combiNgenes = []
combiSscore = []
combiName = []

factor = 0.0
significance = 0.0

def readResult(resfname, csvfname, tabfname):
    isResult = False
    alternative = 1

    # open result file
    with open(resfname, 'r', encoding='utf-8') as resfile:
        # read the data line-by-line
        for line in resfile:
            line = line.strip()
            # read the threshold and correction factor
            if "Correction factor:" in line:
                parts = line.split(" ")
                float(parts[4].split(",")[0])
                float(parts[7])
                continue

            # read the csv filename
            if line.startswith("# item-file") and len(csvfname) == 0:
                parts = line.split(" ")
                csvfname = parts[2]
                continue

            # read the tab filename
            if line.startswith("# value-file") and len(tabfname) == 0:
                parts = line.split(" ")
                tabfname = parts[2]
                continue

            # read the significance level
            if line.startswith("# significance-level:"):
                parts = line.split(" ")
                significance = float(parts[2])
                continue
            
            # read the procedure to compute P-value
            if line.startswith("# P-value computing procedure"):
                parts = line.split(" ")
                pval_procedure = parts[4]
                if len(parts) > 5:
                    alternative_str = parts[5][1:-1]
                    if alternative_str == "greater":
                        alternative = 1
                    elif alternative_str == "less":
                        alternative = -1
                    elif alternative_str == "two.sided":
                        alternative = 0
                continue
            
            # stop inputting the values
            if line.startswith("Time"):
                isResult = False

            # add one set of values
            if isResult:
                line = line.replace(' ', '\t')
                parts = line.split('\t')
                namearray = parts[3].split(',')
                nn = len(namearray)

                # add a motif
                if nn == 1:
                    motifRpvalue.append(float(parts[1]))
                    motifApvalue.append(float(parts[2]))
                    motifName.append(namearray[0])
                    motifNgenes.append(int(parts[5]))
                    motifSscore.append(float(parts[6]))

                # add a combination
                else:
                    combiRank.append(int(parts[0]))
                    combiRpvalue.append(float(parts[1]))
                    combiApvalue.append(float(parts[2]))
                    combiName.append(namearray)
                    combiNgenes.append(int(parts[5]))
                    combiSscore.append(float(parts[6]))

            # start inputting values at the next line
            if line.startswith("Rank"):
                isResult = True
                
    # add motifs missed in the input file
    for i in range(len(combiName)):
        nameset = combiName[i]
        for k in range(len(nameset)):
        
            # search for the motif and specify its ID
            nameid = -1
            for j in range(len(motifName)):
                if nameset[k] in motifName[j]:
                    nameid = j
                    break

            # append the name because it is not described in the input file
            if nameid < 0:
                motifRpvalue.append(-1.0)
                motifApvalue.append(-1.0)
                motifName.append(nameset[k])
                motifNgenes.append(0)
                motifSscore.append(0)

    # retrieve values of notifs not described in the input file
    pval_func = None
    if pval_procedure == "fisher":
        pval_func = functions4fisher.run
    elif pval_procedure == "chi":
        pval_func = functions4chi.run
    else:
        pval_func = functions4u_test.run
    for i in range(len(motifName)):
        if(motifApvalue[i] > 0.0):
            continue
        
        # retrieve values of the current motif
        mnamelist = motifName[i].split(',')
        p_value, down_size = pval_func(csvfname, tabfname, mnamelist, ',', alternative)
        motifRpvalue[i] = p_value
        motifApvalue[i] = p_value * factor
        motifNgenes[i] = down_size

    # return
    return significance
