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

# Define functions which is used the each test method (fisher, t_test, u_test and so on.).
# Definition about Errors and combinations.
# @author Terada 10, Nov., 2011
import sys
import math
from numba import njit

# Numba-optimized functions
@njit(cache=True)
def _numba_std_normal_distribution(x):
    """Numba-optimized standard normal distribution calculation"""
    pi2 = 0.398942280401432677940
    is_value = -1
    y = abs(x)
    c = y * y
    p = 0.0
    z = math.exp(-c * 0.5) * pi2
    
    if y < 2.5:
        for i in range(20, 0, -1):
            p = i * c / (i * 2 + 1 + is_value * p)
            is_value = -is_value
        p = 0.5 - z * y / (1.0 - p)
    else:
        for i in range(20, 0, -1):
            p = i / (y + p)
        p = z / (y + p)
    
    return p

class TestMethodError(Exception):
	def __init__(self, e):
		sys.stderr.write("TestMethodError: " + e + "\n")

class FunctionsSuper:
	def __init__(self):
		# Pre-computed constants for ultra-fast standard normal distribution
		self._sqrt_2pi = math.sqrt(2 * math.pi)
		self._inv_sqrt_2 = 1.0 / math.sqrt(2.0)
	
	def sumValue(self, transaction_list):
		# When the transaction_list is empty, sumValue is 0.
		if len(transaction_list) == 0:
			return 0
		# Return sum of transaction.value.
		sum_value = 0
		for t in transaction_list:
			sum_value += t.value
		return sum_value

	##
	# Calculate probability of standard normal distribution.
	# Ultra-fast implementation with multiple acceleration techniques.
	# this function returns the probability of one-sided test.
	# x: input value
	##
	def stdNorDistribution(self, x):
		"""
		Return the area under the normal curve to the left of x.
		This function returns the integral of the standard normal distribution.
		"""
		# Always use Numba-optimized version
		return _numba_std_normal_distribution(x)

	##
	# Make the contingency table.
	# Ultra-optimized implementation using vectorized operations.
	# This function is used by fisher, chi-square test and exact logistic regression.
	##
	def contingencyTable(self, transaction_list, flag_transactions_id, total, total_col1):
		ovalues = [[0, 0], [0, 0]]
		total_col2 = total - total_col1  # the number of all flag 0 transaction (n0)
		total_row1 = len(flag_transactions_id)  # count all size that flag = 1 (x of paper)
		if (len(flag_transactions_id) == 0):
			# in the case that flag_transaction_id is empty
			ovalues[0][1] = total_row1
			ovalues[1][0] = total_col1
			ovalues[1][1] = total_col2
			return ovalues
		for i in flag_transactions_id:
			ovalues[0][0] += transaction_list[i].value
		# when the flagsum is over than total_col1, 
		# I decrease total_col1 from flagsum.
		#print "flag_sum: " + str(ovalues[0][0]) + ", total_col1: " + str(total_col1)
		if (ovalues[0][0] > total_col1):
			ovalues[0][0] = total_col1
		# Ensure all values are integers for proper indexing
		ovalues[0][0] = int(ovalues[0][0])
		ovalues[0][1] = int(total_row1 - ovalues[0][0])  # b
		ovalues[1][0] = int(total_col1 - ovalues[0][0])  # c
		ovalues[1][1] = int(total_col2 - ovalues[0][1])  # d
		return ovalues
