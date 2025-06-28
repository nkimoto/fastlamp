# LAMP: Limitless-Arity Multiple testing Procedure

LAMP is a Python tool for performing rapid multiple hypothesis testing on combinations of features. It identifies statistically significant combinations from a given dataset using Fisher's exact test, the Mann-Whitney U-test, or the Chi-square test, while controlling the Family-Wise Error Rate (FWER).

It is an implementation of the algorithm described in the following paper:

> Minato, S., Uno, T., Tsuda, K., Terada, A., & Sese, J. (2014). A Fast Method of Statistical Assessment for Combinatorial Hypotheses Based on Frequent Itemset Enumeration. In *Machine Learning and Knowledge Discovery in Databases* (pp. 422-436). Springer, Berlin, Heidelberg.

This repository contains a modernized version of the original LAMP software, migrated to Python 3 with updated dependencies, code style, and **significant performance optimizations**.

## Features

- **Multiple Statistical Tests**: Supports Fisher's exact test, Mann-Whitney U-test, and Chi-square test.
- **FWER Control**: Employs a permutation-based method to control the Family-Wise Error Rate, ensuring the reliability of discovered patterns.
- **Efficient Mining**: Utilizes the LCM (Linear time Closed itemset Miner) for efficient discovery of frequent itemsets.
- **⚡ High Performance**: Multiple optimization techniques for significant speedup over the original implementation.

## Performance Optimizations

FastLAMP incorporates several advanced optimization techniques to achieve significant performance improvements:

### 🚀 Parallel Processing
- **Multiprocessing Support**: P-value calculations are distributed across multiple CPU cores using optimized process pools
- **Intelligent Process Management**: Dynamic adjustment of process count and chunk size based on dataset characteristics
- **Optimized Context**: Uses `fork()` method on Unix systems for faster process creation and memory sharing

### ⚙️ JIT Compilation with Numba
- **Numba-Optimized Functions**: Critical statistical computations are accelerated using Just-In-Time (JIT) compilation
- **Optimized Statistical Calculations**:
  - Hypergeometric distribution for Fisher's exact test
  - Mann-Whitney U-value calculation
  - Standard normal distribution functions
- **Cache-Enabled Compilation**: Functions are compiled once and cached for subsequent runs

### 🔧 High-Performance Mining Backend
- **LCM (Linear time Closed itemset Miner)**: Utilizes a highly optimized C implementation for frequent pattern discovery
- **Closed Pattern Mining**: Efficiently enumerates closed frequent itemsets, reducing redundant computations
- **Memory-Efficient**: Optimized data structures minimize memory footprint during pattern enumeration

### 📊 Algorithmic Optimizations
- **Adaptive Load Balancing**: Smaller chunk sizes for better workload distribution across processes
- **Early Termination**: Smart bounds checking to avoid unnecessary computations
- **Memory-Optimized Data Structures**: Efficient representation of transaction data and patterns

### Performance Benchmarks
Typical performance improvements compared to the original implementation:
- **2-5x speedup** on multi-core systems with parallel processing
- **1.5-3x speedup** from Numba JIT compilation of statistical functions
- **Scalable performance** with increasing dataset size and complexity

## Requirements

- Python 3.12+
- `uv` for package management
- A C compiler (like `gcc`) to build the `lcm` dependency
- Multi-core CPU recommended for optimal parallel performance

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/fastlamp.git
    cd fastlamp
    ```

2.  **Compile the LCM miner:**
    The project requires the `lcm` executable. A `Makefile` is provided in the `lcm53` directory.
    ```bash
    make -C lcm53
    ```
    This will create the `lcm53/lcm` executable.

3.  **Install Python dependencies:**
    This project uses `uv` for dependency management. Install the required packages using:
    ```bash
    uv sync
    ```

## Usage

The main script is `lamp.py`. You can see the full list of options by running `uv run lamp.py --help`.

### Synopsis

```
usage: lamp.py [-h] [-p {fisher,chi,u_test}] [--lcm LCM_PATH] [--max_comb MAX_COMB_SIZE] [-e LOG_FILENAME] [--alternative {greater,less,two.sided}] [-j NUM_JOBS]
               transaction_file value_file alpha
```

### Arguments

-   `transaction_file`: Path to the item-set file. Each row represents a transaction (e.g., a gene), and columns represent items (e.g., transcription factors). A '1' indicates the presence of an item in a transaction, and '0' otherwise.
-   `value_file`: Path to the value file. Each row corresponds to a transaction and contains a numerical or binary value for the statistical test.
-   `alpha`: The significance level (e.g., 0.05).

### Options

-   `-p {fisher,chi,u_test}`, `--pvalue {fisher,chi,u_test}`: The statistical test to use. (default: `fisher`)
-   `--lcm LCM_PATH`: Path to the `lcm` executable. (default: `./lcm53/lcm`)
-   `--max_comb MAX_COMB_SIZE`: The maximum arity (size) of combinations to test. (default: `all`)
-   `--alternative {greater,less,two.sided}`: The alternative hypothesis for the test. (default: `greater`)
-   `-e LOG_FILE`, `--log_file LOG_FILE`: Path to write log output.
-   `-j NUM_JOBS`, `--jobs NUM_JOBS`: Number of parallel processes to use. (default: number of CPU cores)

## Example

You can run LAMP on the provided sample data. This example uses Fisher's exact test to find significant combinations in `sample_item.csv` based on the flags in `sample_expression_over1.csv`.

### Basic Usage
```bash
uv run lamp.py -p fisher sample/sample_item.csv sample/sample_expression_over1.csv 0.05
```

### Parallel Processing (Recommended)
```bash
# Use all available CPU cores (default behavior)
uv run lamp.py -p fisher sample/sample_item.csv sample/sample_expression_over1.csv 0.05

# Use specific number of processes
uv run lamp.py -p fisher -j 4 sample/sample_item.csv sample/sample_expression_over1.csv 0.05
```

### Example Output
```
Read input files ...
Compute the optimal correction factor ... 5
Compute P-values of testable combinations ...
--- Using parallel processing: 12 tests ---
--- Parallel config: 4 processes, chunksize: 3 ---
Output results ...
# LAMP ver. 2.0.3
# item-file: sample/sample_item.csv
# value-file: sample/sample_expression_over1.csv
# significance-level: 0.05
# P-value computing procedure: fisher (greater)
# # of tested elements: 4, # of samples: 15, # of positive samples: 7
# Adjusted significance level: 0.01, Correction factor: 5 (# of target rows >= 5)
# # of significant combinations: 1
Rank    Raw p-value     Adjusted p-value        Combination     Arity   # of target rows        # of positives in the targets
1       0.006993        0.034965        TF1,TF2,TF3     3       5       5
Time (sec.): Computing correction factor 0.003, Enumerating significant combinations 0.001, Total 0.004
```
This output indicates that the combination {TF1, TF2, TF3} was found to be statistically significant, and the parallel processing completed the analysis efficiently.

## Performance Tips

### For Large Datasets
- Use parallel processing (enabled by default) or specify the number of CPU cores with `-j`
- Consider using multiple runs with different random seeds for robustness
- Monitor memory usage for very large datasets

### For Optimal Performance
- Ensure your system has sufficient RAM (recommend 2-4GB per process)
- Use SSD storage for faster I/O operations
- Consider the trade-off between number of processes and memory usage

## Testing

The project includes comprehensive performance benchmarks:

```bash
# Run performance tests
uv run python -m pytest tests/test_lamp.py::TestPerformance -v

# Generate performance comparison graphs
uv run python -c "
from tests.performance_graph_utils import *
from tests.test_lamp import TestPerformance
test = TestPerformance()
results = test.run_performance_analysis('logs')
create_performance_graphs(results, 'logs')
print_performance_summary(results)
"
```

## Original Implementation

This repository is based on the original LAMP implementation by Aika Terada. For comparison with the original version, please refer to the [original repository](https://github.com/a-terada/lamp).

## Dependencies

The project uses modern Python libraries optimized for performance:
- **NumPy**: Efficient numerical operations
- **SciPy**: Statistical functions and distributions  
- **Numba**: JIT compilation for numerical functions
- **Pandas**: Data manipulation and analysis
- **Matplotlib**: Visualization for performance analysis

All dependencies are managed using `uv` and specified in `pyproject.toml`.
