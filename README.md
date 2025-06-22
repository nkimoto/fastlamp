# LAMP: Limitless-Arity Multiple testing Procedure

LAMP is a Python tool for performing rapid multiple hypothesis testing on combinations of features. It identifies statistically significant combinations from a given dataset using Fisher's exact test, the Mann-Whitney U-test, or the Chi-square test, while controlling the Family-Wise Error Rate (FWER).

It is an implementation of the algorithm described in the following paper:

> Minato, S., Uno, T., Tsuda, K., Terada, A., & Sese, J. (2014). A Fast Method of Statistical Assessment for Combinatorial Hypotheses Based on Frequent Itemset Enumeration. In *Machine Learning and Knowledge Discovery in Databases* (pp. 422-436). Springer, Berlin, Heidelberg.

This repository contains a modernized version of the original LAMP software, migrated to Python 3 with updated dependencies and code style.

## Features

- **Multiple Statistical Tests**: Supports Fisher's exact test, Mann-Whitney U-test, and Chi-square test.
- **FWER Control**: Employs a permutation-based method to control the Family-Wise Error Rate, ensuring the reliability of discovered patterns.
- **Efficient Mining**: Utilizes the LCM (Linear time Closed itemset Miner) for efficient discovery of frequent itemsets.

## Requirements

- Python 3.8+
- `uv` for package management
- A C compiler (like `gcc`) to build the `lcm` dependency.

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

The main script is `lamp.py`. You can see the full list of options by running `python3 lamp.py --help`.

### Synopsis

```
usage: lamp.py [-h] [-p {fisher,chi,u_test}] [--lcm LCM_PATH] [--max_comb MAX_COMB_SIZE] [-e LOG_FILENAME] [--alternative {greater,less,two.sided}]
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

## Example

You can run LAMP on the provided sample data. This example uses Fisher's exact test to find significant combinations in `sample_item.csv` based on the flags in `sample_expression_over1.csv`.

```bash
python3 lamp.py -p fisher sample/sample_item.csv sample/sample_expression_over1.csv 0.05
```

### Example Output
```
Read input files ...
Compute the optimal correction factor ... 5
Compute P-values of testable combinations ...
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
Time (sec.): Computing correction factor 0.003, Enumerating significant combinations 0.000, Total 0.003
```
This output indicates that the combination {TF1, TF2, TF3} was found to be statistically significant.
