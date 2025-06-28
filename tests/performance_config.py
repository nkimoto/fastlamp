#!/usr/bin/env python3

"""
Performance test configuration for LAMP testing.

This module contains configuration parameters used by the performance
testing framework in test_lamp.py.
"""

import os


class PerformanceConfig:
    """Configuration for performance tests"""

    # Full analysis configurations
    FULL_SAMPLE_SIZES = [10, 100, 1000, 10000]
    FULL_ITEM_COUNTS = [10, 100, 1000, 10000]
    FULL_DENSITY_VALUES = [0.1, 0.2, 0.3, 0.4]

    # Test parameters
    DEFAULT_SIG_LEVEL = 0.05
    DEFAULT_METHOD = "u_test"
    DEFAULT_LCM_PATH = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lcm53/lcm"
    )
    DEFAULT_ALTERNATIVE = 1
    DEFAULT_JOBS = 1

    # Fixed test parameters
    SCALE_FIXED_ITEMS = 100
    SCALE_FIXED_DENSITY = 0.2
    SCALE_FIXED_SAMPLES = 100
