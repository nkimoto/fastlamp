import argparse
import numpy as np
import pandas as pd
import os


def generate_data(
    num_samples=200,
    num_items=100,
    density=0.2,
    num_significant_items=5,
    value_type="continuous",
):
    """
    Generate benchmark data compatible with both optimized and legacy LAMP.

    Parameters:
    - num_samples: Number of samples (genes for legacy LAMP compatibility)
    - num_items: Number of items (features)
    - density: Proportion of 1s in item matrix
    - num_significant_items: Number of items that should be significant
    - value_type: 'continuous' for u-test, 'binary' for fisher/chi-square tests
    """

    print(f"Generating benchmark data with {num_samples} samples, {num_items} items...")
    print(f"Value type: {value_type}, Significant items: {num_significant_items}")

    np.random.seed(42)  # For reproducibility

    # Generate sample names (genes for legacy LAMP compatibility)
    sample_names = [f"Gene_{i + 1}" for i in range(num_samples)]

    # Generate item names
    item_names = [f"Item_{i + 1}" for i in range(num_items)]

    # Create item matrix (binary: 0 or 1)
    item_matrix = np.random.choice(
        [0, 1], size=(num_samples, num_items), p=[1 - density, density]
    )

    # Create significant pattern for some items
    significant_items = np.random.choice(
        num_items, size=num_significant_items, replace=False
    )

    # Generate values based on item presence for significant items
    if value_type == "binary":
        # For binary values (Fisher/Chi-square tests)
        base_prob = 0.3  # Base probability of being positive
        enhanced_prob = 0.8  # Enhanced probability when significant items are present

        values = np.zeros(num_samples)
        for i in range(num_samples):
            # Check if any significant items are present
            has_significant = np.any(item_matrix[i, significant_items])
            prob = enhanced_prob if has_significant else base_prob
            values[i] = np.random.choice([0, 1], p=[1 - prob, prob])
    else:
        # For continuous values (U-test)
        base_mean = 2.0
        enhancement = 3.0
        noise_std = 1.5

        values = np.zeros(num_samples)
        for i in range(num_samples):
            # Count significant items present
            sig_count = np.sum(item_matrix[i, significant_items])
            mean_value = base_mean + (enhancement * sig_count / len(significant_items))
            values[i] = np.random.normal(mean_value, noise_std)

    # Get parent directory path
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Save item matrix - use #gene instead of #sample for Legacy LAMP compatibility
    item_df = pd.DataFrame(item_matrix, columns=item_names)
    item_df.insert(0, "#gene", sample_names)
    item_df.to_csv(os.path.join(parent_dir, "sample/benchmark_item.csv"), index=False)

    # Save values - use #gene instead of #sample for Legacy LAMP compatibility
    value_df = pd.DataFrame({"#gene": sample_names, "expression": values})
    value_df.to_csv(os.path.join(parent_dir, "sample/benchmark_value.csv"), index=False)

    print("Generated files:")
    print(f"  - sample/benchmark_item.csv: {num_samples} samples x {num_items} items")
    print(
        f"  - sample/benchmark_value.csv: {num_samples} samples with {value_type} values"
    )
    print(f"  - Significant items (indices): {significant_items}")

    return significant_items


def main():
    parser = argparse.ArgumentParser(description="Generate benchmark data for LAMP.")
    parser.add_argument(
        "--samples", type=int, default=500, help="Number of samples (rows)."
    )
    parser.add_argument(
        "--items", type=int, default=200, help="Number of items (columns)."
    )
    parser.add_argument(
        "--density", type=float, default=0.1, help="Density of items in the item file."
    )
    parser.add_argument(
        "--significant",
        type=int,
        default=10,
        help="Number of items to make significant.",
    )
    parser.add_argument(
        "--value-type",
        type=str,
        default="continuous",
        choices=["continuous", "binary"],
        help="Type of data for the value file.",
    )
    args = parser.parse_args()

    generate_data(
        args.samples, args.items, args.density, args.significant, args.value_type
    )


if __name__ == "__main__":
    main()
