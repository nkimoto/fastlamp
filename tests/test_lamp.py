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

import unittest
import sys
import os
import time
import inspect
import glob

# Add parent directory to path to import lamp module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lamp

from .generate_benchmark_data import generate_data
from .performance_config import PerformanceConfig
import shutil
import subprocess
import tempfile
import argparse

# Global settings
LEGACY_REPO_URL = "https://github.com/a-terada/lamp.git"

# Get parent directory for proper paths
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGS_DIR = os.path.join(PARENT_DIR, "logs")

# Create logs directory
if not os.path.exists(LOGS_DIR):
    os.makedirs(LOGS_DIR)

RESULT_FILE = os.path.join(LOGS_DIR, "test_result.log")


class BaseLampTest(unittest.TestCase):
    """Base class with common test functionality"""

    def setUp(self):
        # Test data files (fixed test data for functional tests)
        # Get parent directory path
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        self.csv_file = os.path.join(parent_dir, "sample/sample_item.csv")
        self.flag_file = os.path.join(parent_dir, "sample/sample_expression_over1.csv")
        self.flag_less_file = os.path.join(
            parent_dir, "sample/sample_expression_less1.csv"
        )
        self.value_file_sample = os.path.join(
            parent_dir, "sample/sample_expression_value.csv"
        )
        self.value_less_file = os.path.join(
            parent_dir, "sample/sample_expression_value_rev.csv"
        )
        self.csv_file2 = os.path.join(parent_dir, "sample/sample_item2.csv")
        self.flag_file2 = os.path.join(
            parent_dir, "sample/sample_expression2_over1.csv"
        )
        self.flag_less_file2 = os.path.join(
            parent_dir, "sample/sample_expression2_less1.csv"
        )

        self.sig_level = 0.05

    def _generate_log_filename(self, method, alternative, description=""):
        """Generate a log filename for each test"""
        # Get the calling test method name
        frame = inspect.currentframe()
        try:
            # Go up the call stack to find the test method
            caller_frame = (
                frame.f_back.f_back
            )  # Skip _generate_log_filename and run_lamp_test
            test_method_name = caller_frame.f_code.co_name
        finally:
            del frame

        # Create class name without 'Test' prefix
        class_name = self.__class__.__name__.replace("Test", "").lower()

        # Build filename components
        filename_parts = [class_name, test_method_name, method, f"alt{alternative}"]

        if description:
            filename_parts.append(description.replace(" ", "_").replace(".", "_"))

        filename = "_".join(filename_parts) + ".log"
        return os.path.join(LOGS_DIR, filename)

    def run_lamp_test(
        self,
        csv_file,
        value_file,
        method,
        arity_lim,
        expected_k,
        expected_lam,
        expected_combinations,
        alternative,
        description="",
    ):
        """Helper method to run LAMP and verify results"""
        lcm_path = os.path.join(os.path.dirname(lamp.__file__ or "."), "lcm53/lcm")

        # Generate unique log filename for this test
        log_file = self._generate_log_filename(method, alternative, description)

        with open(RESULT_FILE, "a+", encoding="utf-8") as fw:
            sys.stdout = fw
            print(
                f"\n=== Starting {self.__class__.__name__}.{inspect.currentframe().f_back.f_code.co_name} ===",
                file=fw,
            )
            print(f"Log file: {log_file}", file=fw)
            print(
                f"Method: {method}, Alternative: {alternative}, Arity limit: {arity_lim}",
                file=fw,
            )

            results, k, lam, columnid2name = lamp.run(
                csv_file,
                value_file,
                self.sig_level,
                method,
                lcm_path,
                arity_lim,
                log_file,
                alternative,
                1,
            )
            print(f"\n--- {method} test results ---", file=fw)
            print(f"k={k}, lam={lam}, combinations={len(results)}", file=fw)
            print(
                f"=== Finished {self.__class__.__name__}.{inspect.currentframe().f_back.f_code.co_name} ===\n",
                file=fw,
            )
        sys.stdout = sys.__stdout__

        # Verify results
        self.assertEqual(lam, expected_lam, "Minimum support check failed")
        self.assertEqual(k, expected_k, "Correction factor check failed")
        self.assertEqual(
            len(results),
            len(expected_combinations),
            "Number of significant combinations mismatch",
        )

        # Verify each combination
        expected_map = {frozenset(c[0]): c for c in expected_combinations}
        for result in results:
            detected_set = frozenset(columnid2name[i - 1] for i in result[0])
            self.assertIn(
                detected_set, expected_map, f"Unexpected combination: {detected_set}"
            )

            expected = expected_map[detected_set]
            self.assertAlmostEqual(
                result[1], expected[1], places=5, msg="P-value mismatch"
            )
            self.assertEqual(result[2], expected[2], msg="Support mismatch")
            self.assertAlmostEqual(
                result[3], expected[3], places=5, msg="Score mismatch"
            )


class TestFisherExactTest(BaseLampTest):
    """Tests for Fisher's exact test"""

    def test_fisher_default(self):
        """Test Fisher's exact test without arity limit"""
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.00699300699301, 5, 5])
        ]
        self.run_lamp_test(
            self.csv_file, self.flag_file, "fisher", -1, 5, 5, expected_combinations, 1
        )

    def test_fisher_arity_limit(self):
        """Test Fisher's exact test with arity limit = 2"""
        expected_combinations = [
            tuple([set(["TF1", "TF2"]), 0.00699300699301, 5, 5]),
            tuple([set(["TF1", "TF3"]), 0.00699300699301, 5, 5]),
            tuple([set(["TF2", "TF3"]), 0.00699300699301, 5, 5]),
        ]
        self.run_lamp_test(
            self.csv_file, self.flag_file, "fisher", 2, 7, 5, expected_combinations, 1
        )

    def test_fisher_greater(self):
        """Test Fisher's exact test with alternative='greater'"""
        self.sig_level = 0.5
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.00699300699301, 5, 5]),
            tuple([set(["TF2"]), 0.034965034965, 6, 5]),
            tuple([set(["TF3"]), 0.034965034965, 6, 5]),
        ]
        self.run_lamp_test(
            self.csv_file, self.flag_file, "fisher", -1, 5, 3, expected_combinations, 1
        )

    def test_fisher_two_sided(self):
        """Test Fisher's exact test with alternative='two.sided'"""
        self.sig_level = 0.5
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.00699300699301, 5, 5]),
            tuple([set(["TF2"]), 0.0405594405594, 6, 5]),
            tuple([set(["TF3"]), 0.0405594405594, 6, 5]),
        ]
        self.run_lamp_test(
            self.csv_file, self.flag_file, "fisher", -1, 5, 3, expected_combinations, 0
        )

    def test_fisher_less(self):
        """Test Fisher's exact test with alternative='less'"""
        self.sig_level = 0.5
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.00699300699301, 5, 0]),
            tuple([set(["TF2"]), 0.034965034965, 6, 1]),
            tuple([set(["TF3"]), 0.034965034965, 6, 1]),
        ]
        self.run_lamp_test(
            self.csv_file,
            self.flag_less_file,
            "fisher",
            -1,
            5,
            3,
            expected_combinations,
            -1,
        )

    def test_fisher_greater_empty(self):
        """Test Fisher's exact test with alternative='greater' - empty results"""
        self.sig_level = 0.5
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file,
            self.flag_less_file,
            "fisher",
            -1,
            5,
            3,
            expected_combinations,
            1,
        )

    def test_fisher_two_sided_less_file(self):
        """Test Fisher's exact test with alternative='two.sided' using less file"""
        self.sig_level = 0.5
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.00699300699301, 5, 0]),
            tuple([set(["TF2"]), 0.0405594405594, 6, 1]),
            tuple([set(["TF3"]), 0.0405594405594, 6, 1]),
        ]
        self.run_lamp_test(
            self.csv_file,
            self.flag_less_file,
            "fisher",
            -1,
            5,
            3,
            expected_combinations,
            0,
        )

    def test_fisher_less_empty(self):
        """Test Fisher's exact test with alternative='less' - empty results"""
        self.sig_level = 0.5
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file, self.flag_file, "fisher", -1, 5, 3, expected_combinations, -1
        )

    def test_fisher_balanced_data_greater(self):
        """Test Fisher's exact test with balanced data - greater"""
        self.sig_level = 0.3
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0104895104895, 5, 5]),
            tuple([set(["TF2"]), 0.0512820512821, 6, 5]),
            tuple([set(["TF3"]), 0.0512820512821, 6, 5]),
        ]
        self.run_lamp_test(
            self.csv_file2,
            self.flag_file2,
            "fisher",
            -1,
            5,
            4,
            expected_combinations,
            1,
        )

    def test_fisher_balanced_data_greater_empty(self):
        """Test Fisher's exact test with balanced data - greater empty"""
        self.sig_level = 0.3
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file2,
            self.flag_less_file2,
            "fisher",
            -1,
            5,
            4,
            expected_combinations,
            1,
        )

    def test_fisher_balanced_data_two_sided(self):
        """Test Fisher's exact test with balanced data - two.sided"""
        self.sig_level = 0.3
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.020979020979, 5, 5])
        ]
        self.run_lamp_test(
            self.csv_file2,
            self.flag_file2,
            "fisher",
            -1,
            5,
            5,
            expected_combinations,
            0,
        )

    def test_fisher_balanced_data_two_sided_less_file(self):
        """Test Fisher's exact test with balanced data - two.sided less file"""
        self.sig_level = 0.3
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.020979020979, 5, 0])
        ]
        self.run_lamp_test(
            self.csv_file2,
            self.flag_less_file2,
            "fisher",
            -1,
            5,
            5,
            expected_combinations,
            0,
        )

    def test_fisher_balanced_data_less_empty(self):
        """Test Fisher's exact test with balanced data - less empty"""
        self.sig_level = 0.3
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file2,
            self.flag_file2,
            "fisher",
            -1,
            5,
            4,
            expected_combinations,
            -1,
        )

    def test_fisher_balanced_data_less(self):
        """Test Fisher's exact test with balanced data - less"""
        self.sig_level = 0.3
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0104895104895, 5, 0]),
            tuple([set(["TF2"]), 0.0512820512821, 6, 1]),
            tuple([set(["TF3"]), 0.0512820512821, 6, 1]),
        ]
        self.run_lamp_test(
            self.csv_file2,
            self.flag_less_file2,
            "fisher",
            -1,
            5,
            4,
            expected_combinations,
            -1,
        )


class TestUTest(BaseLampTest):
    """Tests for Mann-Whitney U-test"""

    def test_utest_default(self):
        """Test U-test without arity limit"""
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.00602414187918, 5, 2.510727])
        ]
        self.run_lamp_test(
            self.csv_file,
            self.value_file_sample,
            "u_test",
            -1,
            5,
            3,
            expected_combinations,
            1,
        )

    def test_utest_arity_limit(self):
        """Test U-test with arity limit = 2"""
        expected_combinations = [
            tuple([set(["TF1", "TF2"]), 0.00602414187918, 5, 2.510727]),
            tuple([set(["TF1", "TF3"]), 0.00602414187918, 5, 2.510727]),
            tuple([set(["TF2", "TF3"]), 0.00602414187918, 5, 2.510727]),
        ]
        self.run_lamp_test(
            self.csv_file,
            self.value_file_sample,
            "u_test",
            2,
            7,
            3,
            expected_combinations,
            1,
        )

    def test_utest_greater_empty(self):
        """Test U-test with alternative='greater' - empty results"""
        self.sig_level = 0.05
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file,
            self.value_less_file,
            "u_test",
            -1,
            5,
            3,
            expected_combinations,
            1,
        )

    def test_utest_two_sided(self):
        """Test U-test with alternative='two.sided'"""
        self.sig_level = 0.1
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.01204828, 5, 2.510727])
        ]
        self.run_lamp_test(
            self.csv_file,
            self.value_file_sample,
            "u_test",
            -1,
            5,
            3,
            expected_combinations,
            0,
        )

    def test_utest_two_sided_less_file(self):
        """Test U-test with alternative='two.sided' using less file"""
        self.sig_level = 0.1
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.01204828, 5, -2.510727])
        ]
        self.run_lamp_test(
            self.csv_file,
            self.value_less_file,
            "u_test",
            -1,
            5,
            3,
            expected_combinations,
            0,
        )

    def test_utest_less_empty(self):
        """Test U-test with alternative='less' - empty results"""
        self.sig_level = 0.05
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file,
            self.value_file_sample,
            "u_test",
            -1,
            5,
            3,
            expected_combinations,
            -1,
        )

    def test_utest_less(self):
        """Test U-test with alternative='less'"""
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.00602414187918, 5, -2.510727])
        ]
        self.run_lamp_test(
            self.csv_file,
            self.value_less_file,
            "u_test",
            -1,
            5,
            3,
            expected_combinations,
            -1,
        )


class TestChiSquareTest(BaseLampTest):
    """Tests for Chi-square test"""

    def test_chi_default(self):
        """Test Chi-square test without arity limit"""
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0086855750272, 5, 5])
        ]
        self.run_lamp_test(
            self.csv_file, self.flag_file, "chi", -1, 5, 5, expected_combinations, 1
        )

    def test_chi_arity_limit(self):
        """Test Chi-square test with arity limit = 2"""
        self.sig_level = 0.1
        expected_combinations = [
            tuple([set(["TF1", "TF2"]), 0.0086855750272, 5, 5]),
            tuple([set(["TF1", "TF3"]), 0.0086855750272, 5, 5]),
            tuple([set(["TF2", "TF3"]), 0.0086855750272, 5, 5]),
        ]
        self.run_lamp_test(
            self.csv_file, self.flag_file, "chi", 2, 7, 5, expected_combinations, 1
        )

    def test_chi_greater(self):
        """Test Chi-square test with alternative='greater'"""
        self.sig_level = 0.5
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0086855750272, 5, 5]),
            tuple([set(["TF2"]), 0.036251012711, 6, 5]),
            tuple([set(["TF3"]), 0.036251012711, 6, 5]),
        ]
        self.run_lamp_test(
            self.csv_file, self.flag_file, "chi", -1, 5, 3, expected_combinations, 1
        )

    def test_chi_greater_empty(self):
        """Test Chi-square test with alternative='greater' - empty results"""
        self.sig_level = 0.5
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file,
            self.flag_less_file,
            "chi",
            -1,
            5,
            3,
            expected_combinations,
            1,
        )

    def test_chi_two_sided(self):
        """Test Chi-square test with alternative='two.sided'"""
        self.sig_level = 0.5
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0173711500544, 5, 5]),
            tuple([set(["TF2"]), 0.0725020254219, 6, 5]),
            tuple([set(["TF3"]), 0.072502025419, 6, 5]),
        ]
        self.run_lamp_test(
            self.csv_file, self.flag_file, "chi", -1, 5, 4, expected_combinations, 0
        )

    def test_chi_two_sided_less_file(self):
        """Test Chi-square test with alternative='two.sided' using less file"""
        self.sig_level = 0.5
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0173711500544, 5, 0]),
            tuple([set(["TF2"]), 0.0725020254219, 6, 1]),
            tuple([set(["TF3"]), 0.072502025419, 6, 1]),
        ]
        self.run_lamp_test(
            self.csv_file,
            self.flag_less_file,
            "chi",
            -1,
            5,
            4,
            expected_combinations,
            0,
        )

    def test_chi_less_empty(self):
        """Test Chi-square test with alternative='less' - empty results"""
        self.sig_level = 0.5
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file, self.flag_file, "chi", -1, 5, 3, expected_combinations, -1
        )

    def test_chi_less(self):
        """Test Chi-square test with alternative='less'"""
        self.sig_level = 0.5
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0086855750272, 5, 0]),
            tuple([set(["TF2"]), 0.036251012711, 6, 1]),
            tuple([set(["TF3"]), 0.036251012711, 6, 1]),
        ]
        self.run_lamp_test(
            self.csv_file,
            self.flag_less_file,
            "chi",
            -1,
            5,
            3,
            expected_combinations,
            -1,
        )

    def test_chi_balanced_data_greater(self):
        """Test Chi-square test with balanced data - greater"""
        self.sig_level = 0.3
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0128374712894, 5, 5]),
            tuple([set(["TF2"]), 0.05259625256, 6, 5]),
            tuple([set(["TF3"]), 0.05259625256, 6, 5]),
        ]
        self.run_lamp_test(
            self.csv_file2, self.flag_file2, "chi", -1, 5, 4, expected_combinations, 1
        )

    def test_chi_balanced_data_greater_empty(self):
        """Test Chi-square test with balanced data - greater empty"""
        self.sig_level = 0.3
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file2,
            self.flag_less_file2,
            "chi",
            -1,
            5,
            4,
            expected_combinations,
            1,
        )

    def test_chi_balanced_data_two_sided(self):
        """Test Chi-square test with balanced data - two.sided"""
        self.sig_level = 0.3
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0256749425788, 5, 5])
        ]
        self.run_lamp_test(
            self.csv_file2, self.flag_file2, "chi", -1, 5, 5, expected_combinations, 0
        )

    def test_chi_balanced_data_two_sided_less_file(self):
        """Test Chi-square test with balanced data - two.sided less file"""
        self.sig_level = 0.3
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0256749425788, 5, 0])
        ]
        self.run_lamp_test(
            self.csv_file2,
            self.flag_less_file2,
            "chi",
            -1,
            5,
            5,
            expected_combinations,
            0,
        )

    def test_chi_balanced_data_less_empty(self):
        """Test Chi-square test with balanced data - less empty"""
        self.sig_level = 0.3
        expected_combinations = []
        self.run_lamp_test(
            self.csv_file2, self.flag_file2, "chi", -1, 5, 4, expected_combinations, -1
        )

    def test_chi_balanced_data_less(self):
        """Test Chi-square test with balanced data - less"""
        self.sig_level = 0.3
        expected_combinations = [
            tuple([set(["TF1", "TF2", "TF3"]), 0.0128374712894, 5, 0]),
            tuple([set(["TF2"]), 0.05259625256, 6, 1]),
            tuple([set(["TF3"]), 0.05259625256, 6, 1]),
        ]
        self.run_lamp_test(
            self.csv_file2,
            self.flag_less_file2,
            "chi",
            -1,
            5,
            4,
            expected_combinations,
            -1,
        )


class LegacyEnvironment:
    """Manages legacy LAMP environment setup and execution"""

    def __init__(self):
        self.directory = None
        self.ready = False

    def setup(self):
        """Set up legacy LAMP environment"""
        # Check if Python 2 is available
        try:
            subprocess.run(["python2", "--version"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Python 2 is not available. Legacy LAMP requires Python 2.")
            return False

        # Create temporary directory for legacy environment
        self.directory = tempfile.mkdtemp(prefix="legacy_lamp_")

        try:
            subprocess.run(
                ["git", "clone", "--depth=1", LEGACY_REPO_URL, self.directory],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                ["make", "-C", os.path.join(self.directory, "lcm53")],
                check=True,
                capture_output=True,
                text=True,
            )

            self.ready = True
            print("Legacy environment set up successfully.")
            return True

        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print("Failed to set up legacy environment.")
            print(f"Error: {e}")
            self.cleanup()
            return False

    def run_benchmark(self, method, item_file, value_file, sig_level):
        """Run legacy LAMP benchmark"""
        if not self.ready or not self.directory:
            return None

        legacy_script_path = os.path.join(self.directory, "lamp.py")
        legacy_lcm_path = os.path.join(self.directory, "lcm53", "lcm")

        # Change to legacy directory to capture the auto-generated log file
        original_cwd = os.getcwd()

        start_time = time.time()
        try:
            os.chdir(self.directory)

            # Find existing lamp_log files before execution
            before_files = set(glob.glob("lamp_log_*.txt"))

            subprocess.run(
                [
                    "python2",
                    legacy_script_path,
                    "-p",
                    method,
                    "--lcm",
                    legacy_lcm_path,
                    item_file,
                    value_file,
                    str(sig_level),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            # Find new lamp_log files after execution
            after_files = set(glob.glob("lamp_log_*.txt"))
            new_files = after_files - before_files

            # Move and rename the auto-generated log file
            if new_files:
                project_root = os.path.dirname(
                    os.path.dirname(os.path.abspath(__file__))
                )
                logs_dir = os.path.join(project_root, "logs")
                if not os.path.exists(logs_dir):
                    os.makedirs(logs_dir)

                for log_file in new_files:
                    legacy_log_name = f"legacy_{method}.log"
                    destination = os.path.join(logs_dir, legacy_log_name)
                    shutil.move(log_file, destination)
                    print(f"      Legacy log: {destination}")
                    break

        except subprocess.CalledProcessError as e:
            print(f"Legacy lamp.py failed to execute.\nStderr: {e.stderr}")
            return None
        finally:
            os.chdir(original_cwd)

        return time.time() - start_time

    def cleanup(self):
        """Clean up legacy environment"""
        if self.directory and os.path.exists(self.directory):
            shutil.rmtree(self.directory)
            self.directory = None
            self.ready = False


class PerformanceTestRunner:
    """Handles individual performance test execution"""

    def __init__(self, config: PerformanceConfig):
        self.config = config

    def run_single_test(self, test_params, log_dir):
        """Run a single performance test with given parameters"""
        samples, items, density, config_name = test_params

        # Generate test data
        generate_data(
            num_samples=samples,
            num_items=items,
            density=density,
            num_significant_items=min(10, items),
        )

        # Test optimized LAMP
        optimized_time = self._run_optimized_lamp(config_name, log_dir)

        return {
            "optimized_time": optimized_time,
            "config": config_name,
            "params": {"samples": samples, "items": items, "density": density},
        }

    def run_single_test_with_legacy(self, test_params, log_dir, legacy_env):
        """Run a single performance test with legacy comparison"""
        result = self.run_single_test(test_params, log_dir)

        # Test legacy LAMP
        if legacy_env and legacy_env.ready:
            parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            legacy_time = legacy_env.run_benchmark(
                self.config.DEFAULT_METHOD,
                os.path.join(parent_dir, "sample/benchmark_item.csv"),
                os.path.join(parent_dir, "sample/benchmark_value.csv"),
                self.config.DEFAULT_SIG_LEVEL,
            )
            result["legacy_time"] = legacy_time
        else:
            result["legacy_time"] = None

        return result

    def _run_optimized_lamp(self, config_name, log_dir):
        """Run optimized LAMP and return execution time"""
        print("    Optimized LAMP...")
        start_time = time.time()
        try:
            parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

            # Generate log filename
            log_file_path = os.path.join(log_dir, f"performance_{config_name}.log")

            enrich_lst, k, lam_star, _ = lamp.run(
                os.path.join(parent_dir, "sample/benchmark_item.csv"),
                os.path.join(parent_dir, "sample/benchmark_value.csv"),
                self.config.DEFAULT_SIG_LEVEL,
                self.config.DEFAULT_METHOD,
                self.config.DEFAULT_LCM_PATH,
                "all",
                log_file_path,
                self.config.DEFAULT_ALTERNATIVE,
                self.config.DEFAULT_JOBS,
            )
            execution_time = time.time() - start_time
            print(f"      Optimized: {execution_time:.4f}s, Results: {len(enrich_lst)}")
            print(f"      Log file: {log_file_path}")
            return execution_time

        except Exception as e:
            print(f"      Optimized Error: {e}")
            return None


class TestPerformance(unittest.TestCase):
    """Performance comparison tests"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.config = PerformanceConfig()
        self.test_runner = PerformanceTestRunner(self.config)
        self.legacy_env = None

    @classmethod
    def setUpClass(cls):
        """Set up legacy LAMP environment for comparison"""
        print("\n" + "=" * 50)
        print("Setting up legacy LAMP environment...")

        cls.legacy_env = LegacyEnvironment()
        cls.legacy_env_ready = cls.legacy_env.setup()

        if not cls.legacy_env_ready:
            print("Python 2 not available. Skipping legacy tests.")

        print("=" * 50)

    @classmethod
    def tearDownClass(cls):
        """Clean up legacy environment and any stray log files"""
        if hasattr(cls, "legacy_env") and cls.legacy_env:
            cls.legacy_env.cleanup()

        # Clean up any auto-generated log files from project root
        cls._cleanup_stray_log_files()

    @classmethod
    def _cleanup_stray_log_files(cls):
        """Clean up any lamp log files that ended up in the project root"""
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # Look for lamp log files in project root
        lamp_log_files = glob.glob(os.path.join(project_root, "lamp_log_*.txt"))

        if lamp_log_files:
            print(
                f"Found {len(lamp_log_files)} stray log files in project root. Moving to logs directory..."
            )

            # Ensure logs directory exists
            logs_dir = os.path.join(project_root, "logs")
            if not os.path.exists(logs_dir):
                os.makedirs(logs_dir)

            # Move each log file
            for log_file in lamp_log_files:
                try:
                    filename = os.path.basename(log_file)
                    destination = os.path.join(logs_dir, filename)
                    shutil.move(log_file, destination)
                    print(f"  Moved {filename} to logs/")
                except Exception as e:
                    print(f"  Failed to move {log_file}: {e}")

        # Also check current working directory (tests directory)
        current_dir_logs = glob.glob("lamp_log_*.txt")
        if current_dir_logs:
            print(
                f"Found {len(current_dir_logs)} log files in tests directory. Moving to logs directory..."
            )
            for log_file in current_dir_logs:
                try:
                    destination = os.path.join(LOGS_DIR, os.path.basename(log_file))
                    shutil.move(log_file, destination)
                    print(f"  Moved {log_file} to logs/")
                except Exception as e:
                    print(f"  Failed to move {log_file}: {e}")

    def setUp(self):
        """Set up instance-level legacy environment if needed"""
        if not hasattr(self.__class__, "legacy_env") or not self.__class__.legacy_env:
            self.legacy_env = LegacyEnvironment()
        else:
            self.legacy_env = self.__class__.legacy_env

    def _create_results_structure(self):
        """Create the standard results data structure"""
        return {
            "sample_scaling": {
                "sizes": [],
                "optimized_times": [],
                "legacy_times": [],
                "configs": [],
            },
            "item_scaling": {
                "counts": [],
                "optimized_times": [],
                "legacy_times": [],
                "configs": [],
            },
            "density_scaling": {
                "densities": [],
                "optimized_times": [],
                "legacy_times": [],
                "configs": [],
            },
        }

    def _ensure_log_directory(self, log_dir):
        """Create log directory if it doesn't exist"""
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
            print(f"Created log directory: {log_dir}")

    def _run_scaling_test_batch(
        self, test_type, test_params_list, results, legacy_env=None
    ):
        """Run a batch of scaling tests and store results"""
        for test_params in test_params_list:
            samples, items, density, config_name = test_params
            print(f"  Testing {config_name}...")

            if legacy_env:
                result = self.test_runner.run_single_test_with_legacy(
                    test_params, self._current_log_dir, legacy_env
                )
                legacy_time = result.get("legacy_time")
                if legacy_time:
                    print(f"      Legacy: {legacy_time:.4f}s")
                else:
                    print("      Legacy: Failed")
            else:
                result = self.test_runner.run_single_test(
                    test_params, self._current_log_dir
                )
                legacy_time = 0

            optimized_time = result.get("optimized_time")
            if optimized_time is not None:
                if test_type == "sample":
                    results["sample_scaling"]["sizes"].append(samples)
                    results["sample_scaling"]["optimized_times"].append(optimized_time)
                    results["sample_scaling"]["legacy_times"].append(legacy_time or 0)
                    results["sample_scaling"]["configs"].append(config_name)
                elif test_type == "item":
                    results["item_scaling"]["counts"].append(items)
                    results["item_scaling"]["optimized_times"].append(optimized_time)
                    results["item_scaling"]["legacy_times"].append(legacy_time or 0)
                    results["item_scaling"]["configs"].append(config_name)
                elif test_type == "density":
                    results["density_scaling"]["densities"].append(density)
                    results["density_scaling"]["optimized_times"].append(optimized_time)
                    results["density_scaling"]["legacy_times"].append(legacy_time or 0)
                    results["density_scaling"]["configs"].append(config_name)

    def run_scalability_tests(self, log_dir="logs"):
        """Run scalability tests with both optimized and legacy LAMP"""
        print("Running scalability performance tests with legacy comparison...")

        self._current_log_dir = log_dir
        self._ensure_log_directory(log_dir)

        # Set up legacy environment for this run
        legacy_env = LegacyEnvironment()
        legacy_env.setup()

        # Results storage
        results = self._create_results_structure()

        # Test 1: Sample size scaling
        print("\n1. Testing sample size scaling...")
        sample_test_params = [
            (
                samples,
                self.config.SCALE_FIXED_ITEMS,
                self.config.SCALE_FIXED_DENSITY,
                f"{samples}s_{self.config.SCALE_FIXED_ITEMS}i_{self.config.SCALE_FIXED_DENSITY}d",
            )
            for samples in self.config.FULL_SAMPLE_SIZES
        ]
        self._run_scaling_test_batch("sample", sample_test_params, results, legacy_env)

        # Test 2: Item count scaling
        print("\n2. Testing item count scaling...")
        item_test_params = [
            (
                self.config.SCALE_FIXED_SAMPLES,
                items,
                self.config.SCALE_FIXED_DENSITY,
                f"{self.config.SCALE_FIXED_SAMPLES}s_{items}i_{self.config.SCALE_FIXED_DENSITY}d",
            )
            for items in self.config.FULL_ITEM_COUNTS
        ]
        self._run_scaling_test_batch("item", item_test_params, results, legacy_env)

        # Test 3: Density scaling
        print("\n3. Testing density scaling...")
        density_test_params = [
            (
                self.config.SCALE_FIXED_SAMPLES,
                self.config.SCALE_FIXED_ITEMS,
                density,
                f"{self.config.SCALE_FIXED_SAMPLES}s_{self.config.SCALE_FIXED_ITEMS}i_{density}d",
            )
            for density in self.config.FULL_DENSITY_VALUES
        ]
        self._run_scaling_test_batch(
            "density", density_test_params, results, legacy_env
        )

        # Clean up legacy environment
        legacy_env.cleanup()
        print("Legacy environment cleaned up.")

        return results

    def run_performance_analysis(self, output_dir="logs"):
        """Run performance analysis with legacy comparison"""
        print("Starting LAMP performance analysis with legacy comparison...")
        print(f"Output directory: {output_dir}")

        self._current_log_dir = output_dir
        self._ensure_log_directory(output_dir)

        # Set up legacy environment
        print("Setting up legacy environment for comparison...")
        legacy_env = LegacyEnvironment()
        if legacy_env.setup():
            print("Legacy environment ready for comparison.")
        else:
            print(
                "Failed to setup legacy environment. Continuing without legacy comparison."
            )
            legacy_env = None

        # Results storage
        results = self._create_results_structure()

        # Sample size scaling (using configuration values)
        print("\n1. Sample size scaling...")
        sample_test_params = [
            (
                samples,
                self.config.SCALE_FIXED_ITEMS,
                self.config.SCALE_FIXED_DENSITY,
                f"{samples}s_{self.config.SCALE_FIXED_ITEMS}i_{self.config.SCALE_FIXED_DENSITY}d",
            )
            for samples in self.config.FULL_SAMPLE_SIZES
        ]
        self._run_scaling_test_batch("sample", sample_test_params, results, legacy_env)

        # Item count scaling
        print("\n2. Item count scaling...")
        item_test_params = [
            (
                self.config.SCALE_FIXED_SAMPLES,
                items,
                self.config.SCALE_FIXED_DENSITY,
                f"{self.config.SCALE_FIXED_SAMPLES}s_{items}i_{self.config.SCALE_FIXED_DENSITY}d",
            )
            for items in self.config.FULL_ITEM_COUNTS
        ]
        self._run_scaling_test_batch("item", item_test_params, results, legacy_env)

        # Density scaling
        print("\n3. Density scaling...")
        density_test_params = [
            (
                self.config.SCALE_FIXED_SAMPLES,
                self.config.SCALE_FIXED_ITEMS,
                density,
                f"{self.config.SCALE_FIXED_SAMPLES}s_{self.config.SCALE_FIXED_ITEMS}i_{density}d",
            )
            for density in self.config.FULL_DENSITY_VALUES
        ]
        self._run_scaling_test_batch(
            "density", density_test_params, results, legacy_env
        )

        # Clean up legacy environment
        if legacy_env:
            legacy_env.cleanup()
            print("Legacy environment cleaned up.")

        # Print summary
        self._print_performance_summary(results)

        # Import graph utilities and create graphs
        try:
            from .performance_graph_utils import (
                create_performance_graphs,
                print_performance_summary,
            )

            # Print detailed summary
            print_performance_summary(results)

            # Create graphs
            graph_path = create_performance_graphs(results, output_dir)

            print("\n🎉 Performance analysis completed!")
            if graph_path:
                print(f"📊 Main graph: {graph_path}")
            print(f"📁 All files saved in: {output_dir}/")

            return results, graph_path

        except ImportError:
            print(
                "Warning: Graph utilities not available. Results will not be plotted."
            )
            print("\n🎉 Performance analysis completed!")
            print(f"📁 All files saved in: {output_dir}/")
            return results, None

    def _print_performance_summary(self, results):
        """Print a comprehensive summary of performance comparison results"""
        print("\n" + "=" * 80)
        print("COMPREHENSIVE PERFORMANCE SUMMARY (Optimized vs Legacy LAMP)")
        print("=" * 80)

        # Sample scaling summary
        if results["sample_scaling"]["sizes"]:
            sizes = results["sample_scaling"]["sizes"]
            optimized_times = results["sample_scaling"]["optimized_times"]
            legacy_times = results["sample_scaling"]["legacy_times"]

            print(f"\n📊 Sample Size Scaling ({min(sizes)} to {max(sizes)} samples):")
            print(
                f"   📈 Optimized: {min(optimized_times):.4f}s to {max(optimized_times):.4f}s"
            )

            if any(t > 0 for t in legacy_times):
                valid_legacy = [t for t in legacy_times if t > 0]
                print(
                    f"   📉 Legacy:    {min(valid_legacy):.4f}s to {max(valid_legacy):.4f}s"
                )

                # Calculate average speedup
                speedups = [
                    legacy_time / optimized_time
                    for legacy_time, optimized_time in zip(
                        legacy_times, optimized_times
                    )
                    if legacy_time > 0 and optimized_time > 0
                ]
                if speedups:
                    avg_speedup = sum(speedups) / len(speedups)
                    max_speedup = max(speedups)
                    min_speedup = min(speedups)
                    print(
                        f"   🚀 Speedup: Avg={avg_speedup:.2f}x, Max={max_speedup:.2f}x, Min={min_speedup:.2f}x"
                    )
            else:
                print("   📉 Legacy:    No valid data")

        # Item scaling summary
        if results["item_scaling"]["counts"]:
            counts = results["item_scaling"]["counts"]
            optimized_times = results["item_scaling"]["optimized_times"]
            legacy_times = results["item_scaling"]["legacy_times"]

            print(f"\n📊 Item Count Scaling ({min(counts)} to {max(counts)} items):")
            print(
                f"   📈 Optimized: {min(optimized_times):.4f}s to {max(optimized_times):.4f}s"
            )

            if any(t > 0 for t in legacy_times):
                valid_legacy = [t for t in legacy_times if t > 0]
                print(
                    f"   📉 Legacy:    {min(valid_legacy):.4f}s to {max(valid_legacy):.4f}s"
                )

                speedups = [
                    legacy_time / optimized_time
                    for legacy_time, optimized_time in zip(
                        legacy_times, optimized_times
                    )
                    if legacy_time > 0 and optimized_time > 0
                ]
                if speedups:
                    avg_speedup = sum(speedups) / len(speedups)
                    max_speedup = max(speedups)
                    min_speedup = min(speedups)
                    print(
                        f"   🚀 Speedup: Avg={avg_speedup:.2f}x, Max={max_speedup:.2f}x, Min={min_speedup:.2f}x"
                    )
            else:
                print("   📉 Legacy:    No valid data")

        # Density scaling summary
        if results["density_scaling"]["densities"]:
            densities = results["density_scaling"]["densities"]
            optimized_times = results["density_scaling"]["optimized_times"]
            legacy_times = results["density_scaling"]["legacy_times"]

            print(
                f"\n📊 Density Scaling ({min(densities)} to {max(densities)} density):"
            )
            print(
                f"   📈 Optimized: {min(optimized_times):.4f}s to {max(optimized_times):.4f}s"
            )

            if any(t > 0 for t in legacy_times):
                valid_legacy = [t for t in legacy_times if t > 0]
                print(
                    f"   📉 Legacy:    {min(valid_legacy):.4f}s to {max(valid_legacy):.4f}s"
                )

                speedups = [
                    legacy_time / optimized_time
                    for legacy_time, optimized_time in zip(
                        legacy_times, optimized_times
                    )
                    if legacy_time > 0 and optimized_time > 0
                ]
                if speedups:
                    avg_speedup = sum(speedups) / len(speedups)
                    max_speedup = max(speedups)
                    min_speedup = min(speedups)
                    print(
                        f"   🚀 Speedup: Avg={avg_speedup:.2f}x, Max={max_speedup:.2f}x, Min={min_speedup:.2f}x"
                    )
            else:
                print("   📉 Legacy:    No valid data")

        print("=" * 80)

    def run_full_performance_analysis(self, log_dir="logs"):
        """Legacy method - redirects to run_performance_analysis"""
        return self.run_performance_analysis(log_dir)

    def benchmark_method(self, method, value_type):
        """Run performance benchmark for a specific method"""
        print(f"\n--- Performance Test: {method.upper()} ---")

        # Generate test data
        generate_data(
            num_samples=1000,
            num_items=100,
            density=0.30,
            num_significant_items=10,
            value_type=value_type,
        )

        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        item_file = os.path.join(parent_dir, "sample/benchmark_item.csv")
        value_file = os.path.join(parent_dir, "sample/benchmark_value.csv")

        # Generate log filenames
        seq_log = os.path.join(LOGS_DIR, f"benchmark_{method}_seq.log")
        par_log = os.path.join(LOGS_DIR, f"benchmark_{method}_par.log")

        # Sequential test
        start_time = time.time()
        results_seq, k_seq, lam_seq, _ = lamp.run(
            item_file,
            value_file,
            self.config.DEFAULT_SIG_LEVEL,
            method,
            self.config.DEFAULT_LCM_PATH,
            "all",
            seq_log,
            self.config.DEFAULT_ALTERNATIVE,
            1,
        )
        seq_time = time.time() - start_time

        # Parallel test
        start_time = time.time()
        num_jobs = os.cpu_count() or 1
        results_par, k_par, lam_par, _ = lamp.run(
            item_file,
            value_file,
            self.config.DEFAULT_SIG_LEVEL,
            method,
            self.config.DEFAULT_LCM_PATH,
            "all",
            par_log,
            self.config.DEFAULT_ALTERNATIVE,
            num_jobs,
        )
        par_time = time.time() - start_time

        # Legacy test
        legacy_time = None
        if self.legacy_env and self.legacy_env.ready:
            legacy_time = self.legacy_env.run_benchmark(
                method, item_file, value_file, self.config.DEFAULT_SIG_LEVEL
            )

        # Report results
        print(f"Sequential: {seq_time:.4f}s (log: {seq_log})")
        print(f"Parallel ({num_jobs} jobs): {par_time:.4f}s (log: {par_log})")
        if legacy_time:
            print(f"Legacy: {legacy_time:.4f}s")
            speedup = (
                legacy_time / seq_time
                if seq_time < legacy_time
                else seq_time / legacy_time
            )
            print(
                f"Speedup vs legacy: {speedup:.2f}x {'faster' if seq_time < legacy_time else 'slower'}"
            )

        # Verify consistency
        self.assertEqual(
            results_seq, results_par, "Sequential and parallel results should match"
        )
        self.assertEqual(k_seq, k_par, "Correction factors should match")
        self.assertEqual(lam_seq, lam_par, "Lambda values should match")

        print(f"Found {len(results_seq)} significant combinations")

    def test_fisher_performance(self):
        self.benchmark_method("fisher", "binary")

    def test_utest_performance(self):
        self.benchmark_method("u_test", "continuous")

    def test_chi_performance(self):
        self.benchmark_method("chi", "binary")

    def test_comprehensive_analysis(self):
        """Run comprehensive performance analysis directly"""
        print("\n--- Comprehensive Performance Analysis ---")

        try:
            # Run the analysis directly
            results, graph_path = self.run_performance_analysis(LOGS_DIR)

            # Verify that results were generated
            self.assertIsNotNone(results, "Results should be generated")

            # Verify that at least some data was collected
            has_data = (
                len(results["sample_scaling"]["sizes"]) > 0
                or len(results["item_scaling"]["counts"]) > 0
                or len(results["density_scaling"]["densities"]) > 0
            )
            self.assertTrue(has_data, "Should have collected some performance data")

            # Verify that the main graph was created (if graph utilities are available)
            if graph_path:
                self.assertTrue(
                    os.path.exists(graph_path), f"Graph should be created: {graph_path}"
                )
                print(f"✅ Analysis completed. Graph saved: {graph_path}")
            else:
                print(
                    "✅ Analysis completed. (No graph generated - utilities unavailable)"
                )

        except Exception as e:
            self.fail(f"Performance analysis failed: {e}")

    @staticmethod
    def run_standalone_analysis(output_dir="logs"):
        """Run standalone performance analysis with legacy comparison"""
        # Create a temporary test instance
        test_instance = TestPerformance()

        print("🚀 LAMP Performance Analysis Tool")
        print("=" * 50)

        try:
            print(
                "Running comprehensive performance analysis with legacy comparison..."
            )
            results, graph_path = test_instance.run_performance_analysis(output_dir)
            return results

        except Exception as e:
            print(f"❌ Analysis failed: {e}")
            raise


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="LAMP Performance Testing and Analysis Tool"
    )

    # Run performance analysis (always with legacy comparison)
    parser.add_argument(
        "--output-dir",
        default="logs",
        help="Directory to save results and graphs (default: logs)",
    )

    args = parser.parse_args()

    # Check if running as standalone analysis or unit tests
    if len(sys.argv) > 1 and any(arg.startswith("--") for arg in sys.argv[1:]):
        # Performance analysis mode - when flags are provided
        try:
            results = TestPerformance.run_standalone_analysis(args.output_dir)

            print("\n🎉 Analysis completed!")
            print(f"📁 All files saved in: {args.output_dir}/")

        except Exception as e:
            print(f"❌ Error: {e}")
            import traceback

            traceback.print_exc()
            sys.exit(1)
    else:
        # Unit test mode - when no flags provided
        unittest.main()
