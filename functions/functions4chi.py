#!/usr/bin/env python

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

# Define fuctions that is used in multiple_test.py
# This source includes calculate P-value and MASL of the chi-square test.
# @author Terada, 16, Apr, 2013
# @editor Terada, 11, Mar, 2015,
#     Implement computation of the 'less' and 'two-sided' Chi-square test. 

from __future__ import division
import sys
import os
import numpy as np
from . import functionsSuper as fs
from . import pvalTable
from numba import njit
import math

pardir = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
sys.path.append(pardir)

# Numba-optimized functions for Chi-square test
@njit(cache=True)
def _numba_chi2_pval(chi):
    """Numba-optimized chi-square p-value calculation for df=1"""
    if chi == 0.0:
        return 1.0
    
    # For df=1, chi2 CDF = 2 * normal_cdf(sqrt(chi)) - 1
    # So chi2 survival function = 2 * (1 - normal_cdf(sqrt(chi)))
    # But we want the upper tail, so we use: 1 - normal_cdf(sqrt(chi))
    z = math.sqrt(chi)
    
    # Standard normal CDF calculation (same as in functionsSuper)
    pi2 = 0.398942280401432677940
    is_value = -1
    y = abs(z)
    c = y * y
    p = 0.0
    z_exp = math.exp(-c * 0.5) * pi2
    
    if y < 2.5:
        for i in range(20, 0, -1):
            p = i * c / (i * 2 + 1 + is_value * p)
            is_value = -is_value
        p = 0.5 - z_exp * y / (1.0 - p)
    else:
        for i in range(20, 0, -1):
            p = i / (y + p)
        p = z_exp / (y + p)
    
    return p

@njit(cache=True)
def _numba_chi2_statistic(observed, expected):
    """Numba-optimized chi-square statistic calculation"""
    chi = 0.0
    yates_correction = 0.5  # Apply Yates correction for 2x2 tables
    
    for i in range(4):  # 2x2 table has 4 cells
        if expected[i] > 0:
            diff = abs(observed[i] - expected[i]) - yates_correction
            if diff < 0:
                diff = 0
            chi += (diff * diff) / expected[i]
    
    return chi

@njit(cache=True)
def _numba_calculate_expected(ovalues, total, total_col1):
    """Numba-optimized expected values calculation"""
    total_col2 = total - total_col1
    total_row1 = ovalues[0] + ovalues[1]  # a + b
    total_row2 = ovalues[2] + ovalues[3]  # c + d
    
    expected = np.zeros(4, dtype=np.float64)
    expected[0] = (total_row1 * total_col1) / total    # E[a]
    expected[1] = (total_row1 * total_col2) / total    # E[b]  
    expected[2] = (total_row2 * total_col1) / total    # E[c]
    expected[3] = (total_row2 * total_col2) / total    # E[d]
    
    return expected

##
# Define class
# This class calculate function f that define in paper. (C(n1, x)/C(n0+n1, x))
# transaction_list: list of transactions
# alternative: alternative hypothesis, 1 -> "greater" or "less", 0 -> two.sided
##
class FunctionOfX(fs.FunctionsSuper):
	def __init__(self, transaction_list, row_size, alternative):
		fs.FunctionsSuper.__init__(self)
		self.__t_size = len(transaction_list) # all transaction size
		self.__f_size = self.sumValue(transaction_list) # transaction size which have flag = 1 (n1)
		self.alternative = alternative # alternative hypothesis. greater or less -> 1, two.sided -> 0.
		self.__pvalTable = pvalTable.PvalTable( row_size ) # P-value table
		self.__chiTable = pvalTable.PvalTable( row_size ) # P-value table
		if self.__f_size == 0:
			sys.stdout.write("Error: There is no up-regulate gene.\n")
			sys.exit()
		# Check the transaction value.
		# If the value is not 1 or 0, raise error.
		# Because fisher's exact test does not handle numerical value.
		for t in transaction_list:
			if not (t.value == 1.0 or t.value == 0.0):
				message = "Error: \"" + t.name + "\" value is " + str(t.value)+".\n" + \
									"       But value is 1 or 0 if you test by fisher's exact test.\n"
				sys.stderr.write(message)
				raise ValueError(message)
				
	
	def getN1(self):
		return self.__f_size

	def getAllSize(self):
		return self.__t_size

	##
	# Calculate probability of occurrence probability about table.
	# a: Top left of table
	# b: Top right of table
	# n1: Sum of top and bottom left (a + c)
	# n0: Sum of top and bottom right (b + d)
	##
	def __probabilityTable(self, ovalues):
		# Always use Numba-optimized version
		# Flatten ovalues for Numba processing
		observed = np.array([ovalues[0][0], ovalues[0][1], ovalues[1][0], ovalues[1][1]], dtype=np.float64)
		expected = _numba_calculate_expected(observed, self.__t_size, self.__f_size)
		
		# Calculate chi-square statistic
		chi = _numba_chi2_statistic(observed, expected)
		
		# Calculate p-value
		p_value = _numba_chi2_pval(chi)
		
		return p_value, chi

	def __chi2pval(self, chi):
		"""High-precision chi-square p-value calculation"""
		if (chi == 0.0):
			return 1.0
		else: # dimension = 1
			return (self.stdNorDistribution(chi**0.5))

	def __calMeans(self, ovalues):
		"""High-precision means calculation"""
		total = self.__t_size
		total_col1 = self.__f_size # the number of all flag 1 transaction (n1)
		total_col2 = total - total_col1 # the number of all flag 0 transactio (n0)
		total_row1 = ovalues[0][0] + ovalues[0][1]
		total_row2 = ovalues[1][0] + ovalues[1][1]
		means = []
		means.append([0]*2)
		means.append([0]*2)
		means[0][0] = float(total_row1 * total_col1) / total
		means[0][1] = float(total_row1 * total_col2) / total
		means[1][0] = float(total_row2 * total_col1) / total
		means[1][1] = float(total_row2 * total_col2) / total
		return means
	
	##
	# calclate MASL
	##
	def funcF(self, x):
		"""High-precision MASL calculation for Chi-square test"""
		p1 = p2 = 1.0
		chi1 = chi2 = 0.0
		total_row1 = self.__f_size
		total = self.__t_size
		# when x < n_u
		if x < total_row1:
			ovalues = [[x, 0], [total_row1 - x, total - total_row1]]
			p1, chi1 = self.__probabilityTable( ovalues )
			ovalues = [[0, x], [total_row1, total - total_row1 - x]]
			p2, chi2 = self.__probabilityTable( ovalues )
		# when x >= n_u
		else:
			ovalues = [[total_row1, x-total_row1], [0, total - x]]
			p1, chi1 = self.__probabilityTable( ovalues )
			ovalues = [[0, x], [total_row1, total - total_row1 - x]]
			p2, chi2 = self.__probabilityTable( ovalues )
		if self.alternative == 0:
			p1 = min( p1 * 2., 1.0 )
			p2 = min( p2 * 2., 1.0 )
		if p1 < p2:
			return p1
		else:
			return p2
	
	##
	# Calculate p-value by using chi-square test.
	# transaction_list: List of transactions
	# flag_itemset_id: Transactions which have items
	# alternative: alternative hypothesis, 1 -> greater, 0 -> two-sided, -1 -> less.
	##
	def calPValue(self, transaction_list, flag_transactions_id):
		ovalues = self.contingencyTable( transaction_list, flag_transactions_id, self.__t_size, self.__f_size )
		total_row1 = sum( ovalues[0] )
		p = self.__pvalTable.getValue( total_row1, ovalues[0][0] )
		chi = self.__chiTable.getValue( total_row1, ovalues[0][0] )
		if p < 0: # calculate P-value and save to the table
			p, chi = self.__probabilityTable(ovalues)
			if (self.alternative > 0): 
				if (ovalues[0][0] < min(self.__f_size, total_row1)/2):
					p = 1. - p
			# when the alternative hypothesis is "two.sided", 
			# the P-value is doubled. 
			else:
				p = min( p * 2., 1.0 )
			self.__pvalTable.putValue( total_row1, ovalues[0][0], p )
			self.__chiTable.putValue( total_row1, ovalues[0][0], chi )
#		sys.stdout.write( "x: %d, a:%d, chi: %s, p: %s\n" % (total_row1, ovalues[0][0], chi, p) )
		return p, ovalues[0][0]

def maxLambda(transaction_list):
	# Count each item size
	item_sizes = {}
	for t in transaction_list:
#		print t.itemset
		for item in t.itemset:
#			print item
			# If item does not exist in item_size, then make mapping to 0
			if item not in item_sizes:
				item_sizes[item] = 0
			item_sizes[item] = item_sizes[item] + 1
	# Get max value in item_sizes
	max_value = 0
	for i in item_sizes.values():
		if i > max_value:
			max_value = i
	return max_value


def run(xls_file, value_file, itemset_str_lst, delimiter, alternative):
	global readFile
	import readFile
	transaction_list, columnid2name = readFile.readFiles(xls_file, value_file, delimiter)
	max_lambda = maxLambda(transaction_list)

	if alternative < 0:
		global lamp
		from lamp import reverseValue
		transaction_list = reverseValue( transaction_list, "chi" )
	func = FunctionOfX(transaction_list, max_lambda, abs(alternative))
	colname2id_dict = readFile.colname2id(columnid2name)

	itemset = set()
	for i in itemset_str_lst:
		item_id = colname2id_dict[i]
		itemset.add(item_id + 1)
		
	flag_transactions_id = []
	for i in range( len(transaction_list) ):
		t = transaction_list[i]
		if len( itemset & t.itemset ) == len(itemset):
			flag_transactions_id.append( i )
	p_value, stat_value = func.calPValue(transaction_list, flag_transactions_id)
	n = len(transaction_list)
	n1 = func.getN1()
	sys.stdout.write("p-value: %s (N: %s, n1: %s, x: %s, chi: %s)\n"
					 % (p_value, n, n1, len(flag_transactions_id), stat_value))
	return p_value, len(flag_transactions_id)

"""
if __name__ == "__main__":
	if (len(sys.argv) < 4):
		sys.stderr.write("Error: functions4chi.py [item-file] [value_file] [itemset]\n")
		sys.exit()
	
	xls_file = sys.argv[1]
	value_file = sys.argv[2]
	itemset_str_lst = sys.argv[3].split(',')
	delimiter = ','
	p_value, down_size = run(xls_file, value_file, itemset_str_lst, delimiter)
"""
