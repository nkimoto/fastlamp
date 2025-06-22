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
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Run multiple testing correction.
# This script need transaction file, expression-file and significance-level.
# transaction-file: The file includes associations between TFs and genes.
#     Each line indicates a gene.
#     If gene is targeted by the TF, then value is 1, otherwise 0.
# expression-file: Each line indicates a gene. The column1 is gene name.
#     If gene has the feature, the column2 is 1. The other case 0.
# significance-level: The statistical significance threshold.
# @author Terada 10, March, 2014

import sys
import os.path
import time
import datetime
import random
import math
import readFile
import lamp
from optparse import OptionParser
from functools import cmp_to_key

import functions.functionsSuper as fs

set_opts = ("fisher", "u_test", "chi") # methods which used each test

__version__ =  "1.0.1" + " (LAMP ver." + lamp.__version__ + ")"


def version():
    return __version__


def getValuesList( transaction_list ):
    values_list = []
    for t in transaction_list:
        values_list.append( t.value )
    return values_list

##
# Generate the permuted values dataset.
# transaction_list: the original transaction list.
# org_values_list: the values list contained the transaction_list
##
def permute( transaction_list, org_values_list, rand_seed ):
    random.seed( rand_seed )
    random_index_list = list(range( 0, len(org_values_list) ))
    random.shuffle( random_index_list )
    permute_transaction_list = []
    org2shuffled_list = [-1]*len( transaction_list )
    for i in range( 0, len( random_index_list ) ):
        random_index = random_index_list[ i ]
        permute_transaction_list.append( transaction_list[random_index] )
        org2shuffled_list[ transaction_list[random_index].getID() ] = i
    return permute_transaction_list, org2shuffled_list

##
# Calculate the minimum p-value in the permute_transaction_list
# transaction_list: the permuted transaction list
# trans4lcm: the filename to run LCM
# fre_pattern: the instance to obtain the frequent pattern from the min_sup
# func_f: the instance to calculate MASL and the p-value
# max_comb: the limit to maximum combination size
# org2shuffled_list: the mapping from transaction ID from the row dataset to shuffled dataset
##
def calculateMinimumPValue( permute_transaction_list, trans4lcm, fre_pattern, func_f, \
                            max_comb, org2shuffled_list ):
    min_p = 1.0
    min_p_pattern = None # the minimum p-value of the permuted set
    flag = True
    low_sup = fre_pattern.max_support
    i = 0
    freq_time = 0
    while flag:
        starttime = time.time() # time to construct apriori
        fre_pattern.frequentPatterns( trans4lcm, low_sup, max_comb ) # construct frequent patterns
        bound = fre_pattern.getBound( low_sup )
        if bound > 1:
            bound = func_f.funcF( low_sup ) # minimum support value
            fre_pattern.setBound( low_sup, bound )
        cal_list = fre_pattern.getFrequentList( low_sup ) # Itemset calculated its P-value
        endtime = time.time()
        freq_time += endtime - starttime # time to construct apriori
        for cal_item_set, cal_transaction_list in cal_list:
            i = i + 1
            flag_transaction_list = [] # transaction list which has all items in itemset.
            for t in cal_transaction_list:
                shuffled_id = org2shuffled_list[ t ]
                flag_transaction_list.append( shuffled_id )
            p, stat_score = func_f.calPValue( permute_transaction_list, flag_transaction_list )
            if p < min_p:
                min_p = p
                min_p_pattern = cal_item_set
        
        # If the minimum p-value is less than the lower bound of P-value, finish the calculation.
        if (min_p < bound) or (low_sup <= 1):
            flag = False
        # If the minimum p-value is over than MASL, the minimum support is small and repeat the calculation.
        else:
            low_sup = low_sup - 1
    return min_p, min_p_pattern, freq_time

##
# Generate a probability distribution of the minimum P-value using permuted datasets 
# transaction_list: List of itemset and expression value.
# trans4lcm: File name for argument of LCM program. This file is made in this method.
# threshold: The statistical significance threshold.
# set_method: The procedure name for calibration p-value (fisher/u_test).
# lcm_path: LCM path
# max_comb: the limit to maximum combination size
# permute_num: the number of permuted dataset used in FastWY
# outlog: file object to output logs
##
def generateMinPDist(transaction_list, trans4lcm, threshold, set_method, lcm_path, \
                     max_comb, permute_num, outlog, alternative):
    starttime = time.time()

    # Initialize the apriori and functinos using LAMP. 
    fre_pattern, lam_star, max_lambda, correction_term_time, func_f \
                 = lamp.runMultTest( transaction_list, trans4lcm, threshold, set_method, \
                                     lcm_path, max_comb, outlog, alternative )
    
    # calculate the set of minimum p-values using permuted data
    min_p_list = [] # the list stores the minimum p-values
    org_values_list = getValuesList( transaction_list ) # Raw (non-permuted) dataset
    
    # estimate the probability distribution of the minimum p-value using permuted datasets.
    for i in range( 0, permute_num ):
        per_start = time.time()
        permute_transaction_list, org2shuffled_list = permute( transaction_list, org_values_list, i ) # generate the permuted dataset.
        func_f.calTime = 0
        min_p, min_p_pattern, freq_time = calculateMinimumPValue( permute_transaction_list, trans4lcm, fre_pattern,
                                                            func_f, max_comb, org2shuffled_list )
        per_time = time.time() - per_start
        if ( i == 0 ):
            per_time = time.time() - starttime
        min_p_list.append( tuple( [ min_p, min_p_pattern, fre_pattern.getTotal( min_p_pattern ), freq_time, per_time, func_f.calTime ] ) )
        
        outlog.write( "[permute %s] minP %s, minSupport %s, totalTest %s, freqTime %s, totalTime %s, #ofPvalue %s\n" \
                          % (i, min_p_list[i][0], min_p_list[i][1], min_p_list[i][2], \
                             min_p_list[i][3], min_p_list[i][4], min_p_list[i][5]))
        
    return min_p_list, fre_pattern, func_f

##
# Calculate the adjusted significance level
# min_p_list: the list of minimum P-values used by FastWY
# threshold: the statistical significance threshold.
# permute_num: the number of permuted dataset used in FastWY
##
def adjustedThreshold( min_p_list, threshold, permute_num ):
    # calculate the adjusted significance level
    min_p_index = max( int( math.floor(permute_num * threshold) ) - 1, 0 )
    sorted_min_p_list =  sorted( min_p_list, key = cmp_to_key(lambda x,y: (x[0] > y[0]) - (x[0] < y[0])))
    adjusted_threshold = sorted_min_p_list[ min_p_index ][0] # the adjusted significance level.
    
    if min_p_index + 1 >= len(min_p_list):
        return adjusted_threshold, sorted_min_p_list
    while adjusted_threshold == sorted_min_p_list[ min_p_index + 1][0]:
        min_p_index = min_p_index - 1
        adjusted_threshold = sorted_min_p_list[ min_p_index ][0] # the adjusted significance level.
        if min_p_index < 0:
            min_p_index = 0
            break
    
    return adjusted_threshold, sorted_min_p_list

##
# Enumerate significant combinations (P-value <= adjusted threshold)
# transaction_list: the original transaction list.
# trans4lcm: the filename to run LCM
# fre_pattern: the instance to obtain the frequent pattern from the min_sup
# func_f: the instance to calculate MASL and the p-value
# max_comb: the limit to maximum combination size
# adjusted_threshold: adjusted threshold for P-value
# outlog: file object to output logs
##
def enumerateSigComb(transaction_list, trans4lcm, fre_pattern, func_f, \
                     max_comb, adjusted_threshold, outlog):
    # test the raw (non-permuted) dataset
    i = 0
    enrich_lst = []
    flag = True
    low_sup = fre_pattern.max_support
    freq_time = 0
    start_time = time.time()
    while flag:
        outlog.write("support: %d\n" % low_sup)
        freq_start_time = time.time()
        fre_pattern.frequentPatterns( trans4lcm, low_sup, max_comb ) # construct frequent patterns
        bound = fre_pattern.getBound( low_sup )
        if bound > 1:
            bound = func_f.funcF( low_sup ) # minimum support value
            fre_pattern.setBound( low_sup, bound )
        cal_list = fre_pattern.getFrequentList( low_sup ) # Itemset calculated its P-value
        freq_time = freq_time + time.time() - freq_start_time
        func_f.calTime = 0
        for item_set, item_transid_list in cal_list:
            i = i + 1
            outlog.write("--- testing %s: " % i)
            flag_transaction_list = [] # transaction list which has all items in itemset.
            for t in item_transid_list:
                flag_transaction_list.append( t )
            p, stat_score = func_f.calPValue( transaction_list, flag_transaction_list )
            outlog.write( "p " + str(p) + ", stat_score %s\n" % stat_score )
            if ( p <= adjusted_threshold ):
                enrich_lst.append([item_set, p, low_sup, stat_score])
        # If the minimum p-value is less than MASL, finish the calculation.
        if (adjusted_threshold < bound) or (low_sup <= 1):
            flag = False
        # If the minimum p-value is over than MASL, the minimum support is small and repeat the calculation.
        else:
            low_sup = low_sup - 1
        time_enumerate_total = time.time() - start_time
    return enrich_lst, freq_time, time_enumerate_total

##
# Calculate the adjusted p-value of each combination.
# The adjusted p-value is the proportion of the permuted minimum p-values
# that are less than or equal to the original p-value.
# pvalue: an original p-value
# sorted_min_p_list: the sorted list of the minimum p-values
# start_index: the start index to search the minimum p-value
##
def adjustPval( pvalue, sorted_min_p_list, start_index ):
    i = start_index
    while i < len( sorted_min_p_list ):
        min_p = sorted_min_p_list[i][0]
        if pvalue < min_p:
            break
        i = i + 1
    return float(i) / len(sorted_min_p_list), i

def outputResult( transaction_file, flag_file, threshold, permute_num, set_method, max_comb, \
                  columnid2name, enrich_lst, adjusted_threshold, transaction_list, func_f, sorted_min_p_list, alternative ):
    print("# LAMP: Linear time sampler for discovering statistically significant combination of features")
    print("# File version: %s" % __version__)
    print("# Transaction file: %s" % transaction_file)
    print("# Value file: %s" % flag_file)
    print("# Significance level: %s" % threshold)
    print("# The number of permutations: %s" % permute_num)
    if not max_comb == -1:
        print("# Maximum combination size: %s" % max_comb)
    print("# Adjusted significance level: %f" % adjusted_threshold)
    print("# Alternative hypothesis: %s" % alternative)
    print("p-value\tadjusted p-value\tcombination")
    # If the positives and negatives are reversed, the number of positives is calculated.
    if ( alternative < 0 ) and ( set_method in lamp.BINARY_METHODS ):
        for enrich_item in enrich_lst:
            enrich_item[3] = enrich_item[2] - enrich_item[3]
            
    # caluculate the adjusted p-value of each combination.
    start_index = 0
    for enrich_item in enrich_lst:
        item_set = enrich_item[0]
        p = enrich_item[1]
        adj_p, start_index = adjustPval( p, sorted_min_p_list, start_index )
        print("%f\t%f\t" % (p, adj_p), end="")
        out_column = ""
        for i in item_set:
            out_column += columnid2name[i-1] + " "
        print(out_column)
    print("# Statistical test: %s" % set_method)
    print("# Total number of tested combinations: %s" % len(enrich_lst))


def outputMinP( min_p_list ):
    # output minP distribution
    print("p-value\tminimum support\ttotal tests\tfrequent pattern time\ttotal time\t#of p-value calculation")
    for p, s, t, f, a, c in min_p_list:
        print("%f\t%s\t%s\t%s\t%s\t%s" % (p, s, t, f, a, c))


def run(transaction_file, flag_file, threshold, k, set_method, lcm_path, max_comb, log_file, alternative):
    # read 2 files and get transaction list
    sys.stderr.write("Read input files ...\n")
    try:
        transaction_list, columnid2name = readFile.readFiles(transaction_file, flag_file, ',')
        # If the alternative hypothesis is 'less',
        # the positive and negative of observe values are reversed, 
        # and conduct the identical procedure to 'greater'.
        if alternative < 0:
            transaction_list = lamp.reverseValue( transaction_list, set_method )
        max_comb = lamp.convertMaxComb( max_comb, len(columnid2name) )
    except (ValueError, IOError) as e:
        sys.stderr.write("Error: %s\n" % e)
        raise e
    
    # run multiple test
    transaction4lcm53 = transaction_file + ".4lcm53"
    # run
    try:
        with open( log_file, 'w' ) as outlog:

            # Generate the probability distribution of the minimum p-value.
            sys.stderr.write( "Generate the distribution of the minimum p-value ...\n" )
            min_p_list, fre_pattern, func_f = generateMinPDist( transaction_list, transaction4lcm53, threshold,
                                                                set_method, lcm_path, max_comb, k, outlog, alternative )
            
            # Calculate the adjusted significance level
            sys.stderr.write( "Calculate the adjusted significance level ... " )
            adjusted_threshold, sorted_min_p_list = adjustedThreshold( min_p_list, threshold, k )
            sys.stderr.write( "%f\n" % adjusted_threshold )
            
            # Enumerate significant combinations
            sys.stderr.write( "Enumerate significant combinations ...\n" )
            enrich_lst, freq_time, time_enumerate_total = enumerateSigComb( transaction_list, transaction4lcm53,
                                                                            fre_pattern, func_f, max_comb, adjusted_threshold, outlog )
            
    except (IOError, lamp.MASLError, fs.TestMethodError) as e:
        sys.stderr.write("Error: %s\n" % e)
        raise e
    
    # output the result
    outputResult(transaction_file, flag_file, threshold, k, set_method, max_comb, \
                 columnid2name, enrich_lst, adjusted_threshold, transaction_list, func_f, sorted_min_p_list, alternative)

if __name__ == "__main__":
    usage = "usage: %prog [options] transaction_file value_file significance_probability"
    p = OptionParser(usage = usage, version = "%prog " + __version__)
    p.add_option('-p', '--pvalue', dest = "pvalue_procedure",\
                 help = "Choose the p-value calculation procedure from 'fisher' (Fisher's exact test) or 'u_test' (Mann-Whitney's U-test)")
    
    p.add_option('-k', dest = "permute_num", type = "int", default = 100, \
                 help = "The number of permutation. (default: 100)" )

    p.add_option('--lcm', dest = "lcm_path", \
                 default = os.path.dirname(os.path.abspath( __file__ )) + "/lcm53/lcm", \
                 help = "Set LCM program path if you do not use ./lcm53/lcm")

    p.add_option('--max_comb', dest = "max_comb", default = "all", \
                 help = "Set the maximum size of combination to be tested.")

    p.add_option('-e', dest = "log_filename", default = "", help = "The file name to output log.\n")

    p.add_option('--min_p_dist', dest = "min_p_dist", default = "", \
                 help = "The file name to output the distribution of the minimum p-value.\n")

    p.add_option('--alternative', dest = "alternative", default = "greater", \
                 help = "Indicate which alternative hypothesis is used. Select \"greater\", \"less\" or \"two.sided\"\n, and the default is \"greater\".")

    opts, args = p.parse_args()
    
    # check arguments
    if len(args) != 3:
        message = "Error: input [target-file], [expression-file] and [significance-level].\n"
        sys.stderr.write(message)
        raise ValueError(message)
        
    if opts.pvalue_procedure not in set_opts:
        message = "Error: choose \"fisher\", \"chi\" or \"u_test\" by using -p option.\n"
        sys.stderr.write(message)
        raise ValueError(message)
        
    opts.max_comb = opts.max_comb.lower()
    if not opts.max_comb == "all":
        if (opts.max_comb.isdigit()):
            opts.max_comb = int(opts.max_comb)
        else:
            message = "Error: max_comb must be an integer value or all.\n"
            sys.stderr.write(message)
            raise ValueError(message)
        
    # check the file exist.
    if not os.path.isfile(args[0]):
        message = "IOError: No such file: \'" + args[0] + "\'\n"
        sys.stderr.write(message)
        raise IOError(message)
    if not os.path.isfile(args[1]):
        message = "IOError: No such file: \'" + args[1] + "\'\n"
        sys.stderr.write(message)
        raise IOError(message)
    try:
        sig_pro = float(args[2])
        if (sig_pro < 0) or (sig_pro > 1):
            raise ValueError
    except ValueError:
        message = "Error: significance probability must be a float value from 0.0 to 1.0.\n"
        sys.stderr.write(message)
        raise ValueError(message)

    # check the value of alternative hypothesis
    if opts.alternative == "greater":
        opts.alternative = 1
    elif opts.alternative == "less":
        opts.alternative = -1
    elif opts.alternative == "two.sided":
        opts.alternative = 0
    else:
        message = "Error: \"alternative\" should be one of {\"greater\", \"less\", \"two.sided\"}\n"
        sys.stderr.write(message)
        raise ValueError(message)

    # change log file
    d = datetime.datetime.today()
    log_file = "fastwy_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
    if len(opts.log_filename) > 0:
        log_file = opts.log_filename

    try:
        run(args[0], args[1], sig_pro, opts.permute_num, opts.pvalue_procedure, \
            opts.lcm_path, opts.max_comb, log_file, opts.alternative)
    except (ValueError, IOError, lamp.MASLError, fs.TestMethodError) as e:
        sys.stderr.write("Error: %s\n" % e)
        sys.exit(1)
