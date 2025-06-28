#!/usr/bin/env python3
"""
Performance graph utilities for LAMP - Graph creation only
"""

import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import os


def create_performance_graphs(results, output_dir="logs"):
    """Create performance graphs from test results"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Set up matplotlib style
    plt.style.use(
        "seaborn-v0_8" if "seaborn-v0_8" in plt.style.available else "default"
    )

    # Create a figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(
        "LAMP Performance Analysis: Optimized vs Legacy", fontsize=16, fontweight="bold"
    )

    # Graph 1: Sample size scaling
    if results["sample_scaling"]["sizes"]:
        ax1 = axes[0, 0]
        sizes = results["sample_scaling"]["sizes"]
        optimized_times = results["sample_scaling"]["optimized_times"]
        legacy_times = results["sample_scaling"]["legacy_times"]

        # Plot optimized LAMP
        ax1.plot(
            sizes,
            optimized_times,
            "bo-",
            linewidth=2,
            markersize=6,
            label="Optimized LAMP",
        )

        # Plot legacy LAMP (only if we have valid data)
        if any(t > 0 for t in legacy_times):
            valid_sizes = [s for s, t in zip(sizes, legacy_times) if t > 0]
            valid_legacy_clean = [t for t in legacy_times if t > 0]

            ax1.plot(
                valid_sizes,
                valid_legacy_clean,
                "ro-",
                linewidth=2,
                markersize=6,
                label="Legacy LAMP",
            )

        ax1.set_xlabel("Number of Samples")
        ax1.set_ylabel("Execution Time (seconds)")
        ax1.set_title("Scalability: Sample Size vs Execution Time")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Add speedup annotation
        if any(t > 0 for t in legacy_times) and len(optimized_times) > 0:
            avg_speedup = np.mean(
                [
                    legacy_time / optimized_time
                    for legacy_time, optimized_time in zip(
                        legacy_times, optimized_times
                    )
                    if legacy_time > 0 and optimized_time > 0
                ]
            )
            ax1.text(
                0.05,
                0.95,
                f"Avg Speedup: {avg_speedup:.1f}x",
                transform=ax1.transAxes,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.7),
            )

    # Graph 2: Item count scaling
    if results["item_scaling"]["counts"]:
        ax2 = axes[0, 1]
        counts = results["item_scaling"]["counts"]
        optimized_times = results["item_scaling"]["optimized_times"]
        legacy_times = results["item_scaling"]["legacy_times"]

        # Plot optimized LAMP
        ax2.plot(
            counts,
            optimized_times,
            "go-",
            linewidth=2,
            markersize=6,
            label="Optimized LAMP",
        )

        # Plot legacy LAMP (only if we have valid data)
        if any(t > 0 for t in legacy_times):
            valid_counts = [c for c, t in zip(counts, legacy_times) if t > 0]
            valid_legacy_clean = [t for t in legacy_times if t > 0]

            ax2.plot(
                valid_counts,
                valid_legacy_clean,
                "ro-",
                linewidth=2,
                markersize=6,
                label="Legacy LAMP",
            )

        ax2.set_xlabel("Number of Items")
        ax2.set_ylabel("Execution Time (seconds)")
        ax2.set_title("Scalability: Item Count vs Execution Time")
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        # Add speedup annotation
        if any(t > 0 for t in legacy_times) and len(optimized_times) > 0:
            avg_speedup = np.mean(
                [
                    legacy_time / optimized_time
                    for legacy_time, optimized_time in zip(
                        legacy_times, optimized_times
                    )
                    if legacy_time > 0 and optimized_time > 0
                ]
            )
            ax2.text(
                0.05,
                0.95,
                f"Avg Speedup: {avg_speedup:.1f}x",
                transform=ax2.transAxes,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.7),
            )

    # Graph 3: Density scaling
    if results["density_scaling"]["densities"]:
        ax3 = axes[1, 0]
        densities = results["density_scaling"]["densities"]
        optimized_times = results["density_scaling"]["optimized_times"]
        legacy_times = results["density_scaling"]["legacy_times"]

        # Plot optimized LAMP
        ax3.plot(
            densities,
            optimized_times,
            "mo-",
            linewidth=2,
            markersize=6,
            label="Optimized LAMP",
        )

        # Plot legacy LAMP (only if we have valid data)
        if any(t > 0 for t in legacy_times):
            valid_densities = [d for d, t in zip(densities, legacy_times) if t > 0]
            valid_legacy_clean = [t for t in legacy_times if t > 0]

            ax3.plot(
                valid_densities,
                valid_legacy_clean,
                "ro-",
                linewidth=2,
                markersize=6,
                label="Legacy LAMP",
            )

        ax3.set_xlabel("Data Density")
        ax3.set_ylabel("Execution Time (seconds)")
        ax3.set_title("Scalability: Data Density vs Execution Time")
        ax3.grid(True, alpha=0.3)
        ax3.legend()

        # Add speedup annotation
        if any(t > 0 for t in legacy_times) and len(optimized_times) > 0:
            avg_speedup = np.mean(
                [
                    legacy_time / optimized_time
                    for legacy_time, optimized_time in zip(
                        legacy_times, optimized_times
                    )
                    if legacy_time > 0 and optimized_time > 0
                ]
            )
            ax3.text(
                0.05,
                0.95,
                f"Avg Speedup: {avg_speedup:.1f}x",
                transform=ax3.transAxes,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.7),
            )

    # Graph 4: Speedup comparison (grouped by category)
    ax4 = axes[1, 1]

    # Show speedup ratios grouped by test type
    if (
        results["sample_scaling"]["sizes"]
        and results["item_scaling"]["counts"]
        and results["density_scaling"]["densities"]
    ):
        categories = []
        speedup_values = []
        category_colors = []

        # Sample scaling speedups
        if any(t > 0 for t in results["sample_scaling"]["legacy_times"]):
            sample_speedups = [
                legacy_time / optimized_time
                for legacy_time, optimized_time in zip(
                    results["sample_scaling"]["legacy_times"],
                    results["sample_scaling"]["optimized_times"],
                )
                if legacy_time > 0 and optimized_time > 0
            ]
            sample_sizes = [
                s
                for s, t in zip(
                    results["sample_scaling"]["sizes"],
                    results["sample_scaling"]["legacy_times"],
                )
                if t > 0
            ]

            if sample_speedups:
                for i, (size, speedup) in enumerate(zip(sample_sizes, sample_speedups)):
                    categories.append(f"{size}s")
                    speedup_values.append(speedup)
                    category_colors.append("skyblue")

        # Item scaling speedups
        if any(t > 0 for t in results["item_scaling"]["legacy_times"]):
            item_speedups = [
                legacy_time / optimized_time
                for legacy_time, optimized_time in zip(
                    results["item_scaling"]["legacy_times"],
                    results["item_scaling"]["optimized_times"],
                )
                if legacy_time > 0 and optimized_time > 0
            ]
            item_counts = [
                c
                for c, t in zip(
                    results["item_scaling"]["counts"],
                    results["item_scaling"]["legacy_times"],
                )
                if t > 0
            ]

            if item_speedups:
                for i, (count, speedup) in enumerate(zip(item_counts, item_speedups)):
                    categories.append(f"{count}i")
                    speedup_values.append(speedup)
                    category_colors.append("lightgreen")

        # Density scaling speedups
        if any(t > 0 for t in results["density_scaling"]["legacy_times"]):
            density_speedups = [
                legacy_time / optimized_time
                for legacy_time, optimized_time in zip(
                    results["density_scaling"]["legacy_times"],
                    results["density_scaling"]["optimized_times"],
                )
                if legacy_time > 0 and optimized_time > 0
            ]
            densities = [
                d
                for d, t in zip(
                    results["density_scaling"]["densities"],
                    results["density_scaling"]["legacy_times"],
                )
                if t > 0
            ]

            if density_speedups:
                for i, (density, speedup) in enumerate(
                    zip(densities, density_speedups)
                ):
                    categories.append(f"{density:.1f}d")
                    speedup_values.append(speedup)
                    category_colors.append("lightcoral")

        if speedup_values:
            x_pos = range(len(speedup_values))
            bars = ax4.bar(
                x_pos,
                speedup_values,
                color=category_colors,
                alpha=0.8,
                edgecolor="black",
                linewidth=0.5,
            )

            # Add horizontal line at y=1 (no speedup)
            ax4.axhline(
                y=1,
                color="red",
                linestyle="--",
                alpha=0.7,
                linewidth=2,
                label="No speedup (1.0x)",
            )

            # Customize appearance
            ax4.set_xlabel("Test Configuration", fontsize=10, fontweight="bold")
            ax4.set_ylabel(
                "Speedup\n(Legacy ÷ Optimized)", fontsize=10, fontweight="bold"
            )
            ax4.set_title(
                "Optimized LAMP Speedup vs Legacy LAMP", fontsize=12, fontweight="bold"
            )
            ax4.set_xticks(x_pos)
            ax4.set_xticklabels(categories, rotation=45, ha="right", fontsize=8)
            ax4.grid(True, alpha=0.3, axis="y")

            # Add value labels on top of bars
            for i, (bar, value) in enumerate(zip(bars, speedup_values)):
                height = bar.get_height()
                ax4.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height + max(speedup_values) * 0.01,
                    f"{value:.1f}x",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    fontweight="bold",
                )

            # Add legend for colors
            from matplotlib.patches import Patch

            legend_elements = [
                Patch(facecolor="skyblue", label="Sample Size (s)"),
                Patch(facecolor="lightgreen", label="Item Count (i)"),
                Patch(facecolor="lightcoral", label="Data Density (d)"),
                ax4.get_lines()[0],  # The horizontal line
            ]
            ax4.legend(handles=legend_elements, loc="upper left", fontsize=8)

            # Add summary statistics
            avg_speedup = np.mean(speedup_values)
            max_speedup = max(speedup_values)
            min_speedup = min(speedup_values)

            summary_text = f"Overall Performance:\nAvg: {avg_speedup:.1f}x faster\nBest: {max_speedup:.1f}x faster\nWorst: {min_speedup:.1f}x faster"
            ax4.text(
                0.98,
                0.98,
                summary_text,
                transform=ax4.transAxes,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.9),
                verticalalignment="top",
                horizontalalignment="right",
                fontsize=9,
            )

            # Set y-axis to start from 0 for better visualization
            ax4.set_ylim(bottom=0, top=max(speedup_values) * 1.15)

        else:
            ax4.text(
                0.5,
                0.5,
                "No Legacy Comparison Data Available\n\nLegacy LAMP benchmarks were not run\nor failed to execute properly.",
                ha="center",
                va="center",
                transform=ax4.transAxes,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.5),
                fontsize=12,
            )
            ax4.set_title(
                "Speedup Comparison\n(Legacy Data Not Available)",
                fontsize=14,
                fontweight="bold",
            )
            ax4.set_xlim(0, 1)
            ax4.set_ylim(0, 1)

    # Adjust layout and save
    plt.tight_layout()

    # Save the graph
    graph_path = os.path.join(output_dir, "lamp_performance_analysis.png")
    plt.savefig(graph_path, dpi=300, bbox_inches="tight")
    print(f"\n📊 Performance graph saved: {graph_path}")

    # Save individual graphs as well
    graph_configs = [
        ("sample_scaling", axes[0, 0]),
        ("item_scaling", axes[0, 1]),
        ("density_scaling", axes[1, 0]),
        ("speedup_comparison", axes[1, 1]),
    ]

    for i, (title, ax) in enumerate(graph_configs):
        individual_fig = plt.figure(figsize=(8, 6))
        individual_ax = individual_fig.add_subplot(111)

        # Copy the content from the original axis
        if title == "speedup_comparison":
            # Handle bar chart differently - recreate the improved chart
            categories = []
            speedup_values = []
            category_colors = []

            # Check if we have speedup data
            if (
                results["sample_scaling"]["sizes"]
                and results["item_scaling"]["counts"]
                and results["density_scaling"]["densities"]
            ):
                # Sample scaling speedups
                if any(t > 0 for t in results["sample_scaling"]["legacy_times"]):
                    sample_speedups = [
                        legacy_time / optimized_time
                        for legacy_time, optimized_time in zip(
                            results["sample_scaling"]["legacy_times"],
                            results["sample_scaling"]["optimized_times"],
                        )
                        if legacy_time > 0 and optimized_time > 0
                    ]
                    sample_sizes = [
                        s
                        for s, t in zip(
                            results["sample_scaling"]["sizes"],
                            results["sample_scaling"]["legacy_times"],
                        )
                        if t > 0
                    ]

                    if sample_speedups:
                        for i, (size, speedup) in enumerate(
                            zip(sample_sizes, sample_speedups)
                        ):
                            categories.append(f"{size}s")
                            speedup_values.append(speedup)
                            category_colors.append("skyblue")

                # Item scaling speedups
                if any(t > 0 for t in results["item_scaling"]["legacy_times"]):
                    item_speedups = [
                        legacy_time / optimized_time
                        for legacy_time, optimized_time in zip(
                            results["item_scaling"]["legacy_times"],
                            results["item_scaling"]["optimized_times"],
                        )
                        if legacy_time > 0 and optimized_time > 0
                    ]
                    item_counts = [
                        c
                        for c, t in zip(
                            results["item_scaling"]["counts"],
                            results["item_scaling"]["legacy_times"],
                        )
                        if t > 0
                    ]

                    if item_speedups:
                        for i, (count, speedup) in enumerate(
                            zip(item_counts, item_speedups)
                        ):
                            categories.append(f"{count}i")
                            speedup_values.append(speedup)
                            category_colors.append("lightgreen")

                # Density scaling speedups
                if any(t > 0 for t in results["density_scaling"]["legacy_times"]):
                    density_speedups = [
                        legacy_time / optimized_time
                        for legacy_time, optimized_time in zip(
                            results["density_scaling"]["legacy_times"],
                            results["density_scaling"]["optimized_times"],
                        )
                        if legacy_time > 0 and optimized_time > 0
                    ]
                    densities = [
                        d
                        for d, t in zip(
                            results["density_scaling"]["densities"],
                            results["density_scaling"]["legacy_times"],
                        )
                        if t > 0
                    ]

                    if density_speedups:
                        for i, (density, speedup) in enumerate(
                            zip(densities, density_speedups)
                        ):
                            categories.append(f"{density:.1f}d")
                            speedup_values.append(speedup)
                            category_colors.append("lightcoral")

                if speedup_values:
                    x_pos = range(len(speedup_values))
                    bars = individual_ax.bar(
                        x_pos,
                        speedup_values,
                        color=category_colors,
                        alpha=0.8,
                        edgecolor="black",
                        linewidth=0.5,
                    )

                    # Add horizontal line at y=1 (no speedup)
                    individual_ax.axhline(
                        y=1,
                        color="red",
                        linestyle="--",
                        alpha=0.7,
                        linewidth=2,
                        label="No speedup (1.0x)",
                    )

                    # Customize appearance
                    individual_ax.set_xticks(x_pos)
                    individual_ax.set_xticklabels(
                        categories, rotation=45, ha="right", fontsize=8
                    )
                    individual_ax.grid(True, alpha=0.3, axis="y")

                    # Add value labels on top of bars
                    for i, (bar, value) in enumerate(zip(bars, speedup_values)):
                        height = bar.get_height()
                        individual_ax.text(
                            bar.get_x() + bar.get_width() / 2.0,
                            height + max(speedup_values) * 0.01,
                            f"{value:.1f}x",
                            ha="center",
                            va="bottom",
                            fontsize=8,
                            fontweight="bold",
                        )

                    # Add legend for colors
                    from matplotlib.patches import Patch

                    legend_elements = [
                        Patch(facecolor="skyblue", label="Sample Size (s)"),
                        Patch(facecolor="lightgreen", label="Item Count (i)"),
                        Patch(facecolor="lightcoral", label="Data Density (d)"),
                        individual_ax.get_lines()[0]
                        if individual_ax.get_lines()
                        else None,  # The horizontal line
                    ]
                    individual_ax.legend(
                        handles=[elem for elem in legend_elements if elem is not None],
                        loc="upper left",
                        fontsize=8,
                    )

                    # Add summary statistics
                    avg_speedup = np.mean(speedup_values)
                    max_speedup = max(speedup_values)
                    min_speedup = min(speedup_values)

                    summary_text = f"Overall Performance:\nAvg: {avg_speedup:.1f}x faster\nBest: {max_speedup:.1f}x faster\nWorst: {min_speedup:.1f}x faster"
                    individual_ax.text(
                        0.98,
                        0.98,
                        summary_text,
                        transform=individual_ax.transAxes,
                        bbox=dict(
                            boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.9
                        ),
                        verticalalignment="top",
                        horizontalalignment="right",
                        fontsize=9,
                    )

                    # Set y-axis to start from 0 for better visualization
                    individual_ax.set_ylim(bottom=0, top=max(speedup_values) * 1.15)

                else:
                    individual_ax.text(
                        0.5,
                        0.5,
                        "No Legacy Comparison Data Available\n\nLegacy LAMP benchmarks were not run\nor failed to execute properly.",
                        ha="center",
                        va="center",
                        transform=individual_ax.transAxes,
                        bbox=dict(
                            boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.5
                        ),
                        fontsize=12,
                    )
                    individual_ax.set_xlim(0, 1)
                    individual_ax.set_ylim(0, 1)
        else:
            # Handle line plots
            for line in ax.get_lines():
                individual_ax.plot(
                    line.get_xdata(),
                    line.get_ydata(),
                    color=line.get_color(),
                    marker=line.get_marker(),
                    linestyle=line.get_linestyle(),
                    linewidth=line.get_linewidth(),
                    markersize=line.get_markersize(),
                    label=line.get_label(),
                )

        if title != "speedup_comparison":
            individual_ax.set_xlabel(ax.get_xlabel())
            individual_ax.set_ylabel(ax.get_ylabel())
            individual_ax.set_title(ax.get_title())
            individual_ax.grid(True, alpha=0.3)
            if ax.get_legend():
                individual_ax.legend()
        else:
            # For speedup comparison, set labels manually as they might not be set in the original axis
            individual_ax.set_xlabel(
                "Test Configuration", fontsize=10, fontweight="bold"
            )
            individual_ax.set_ylabel(
                "Speedup\n(Legacy ÷ Optimized)", fontsize=10, fontweight="bold"
            )
            individual_ax.set_title(
                "Optimized LAMP Speedup vs Legacy LAMP", fontsize=12, fontweight="bold"
            )

        individual_path = os.path.join(output_dir, f"lamp_performance_{title}.png")

        # Adjust layout to prevent clipping
        try:
            individual_fig.tight_layout()
        except Exception:
            pass  # Skip if tight_layout fails

        individual_fig.savefig(individual_path, dpi=300, bbox_inches="tight")
        plt.close(individual_fig)
        print(f"📊 Individual graph saved: {individual_path}")

    plt.close(fig)

    return graph_path


def print_performance_summary(results):
    """Print a comprehensive summary of performance test results with statistical analysis"""
    print("\n" + "=" * 80)
    print("COMPREHENSIVE PERFORMANCE ANALYSIS (Optimized vs Legacy LAMP)")
    print("=" * 80)

    # Sample scaling summary
    if results["sample_scaling"]["sizes"]:
        sizes = results["sample_scaling"]["sizes"]
        optimized_times = results["sample_scaling"]["optimized_times"]
        legacy_times = results["sample_scaling"]["legacy_times"]

        print(
            f"\n📊 1. Sample Size Scaling Analysis ({min(sizes)} to {max(sizes)} samples):"
        )
        print(
            f"   📈 Optimized LAMP: {min(optimized_times):.4f}s to {max(optimized_times):.4f}s"
        )

        # Growth rate analysis
        if len(optimized_times) > 1:
            growth_rate = (max(optimized_times) / min(optimized_times)) / (
                max(sizes) / min(sizes)
            )
            print(
                f"   📊 Growth efficiency: {growth_rate:.3f} (lower is better, 1.0 = linear)"
            )

            # Correlation analysis
            correlation = np.corrcoef(sizes, optimized_times)[0, 1]
            print(f"   🔗 Correlation coefficient: {correlation:.3f}")

            # Complexity analysis (approximate)
            if correlation > 0.8:  # Strong linear correlation
                complexity = "O(n) - Linear"
            elif correlation > 0.6:
                complexity = "O(n log n) - Log-linear"
            else:
                complexity = "O(n²) or higher - Polynomial/Exponential"
            print(f"   ⚡ Estimated complexity: {complexity}")

        if any(t > 0 for t in legacy_times):
            valid_legacy = [t for t in legacy_times if t > 0]
            print(
                f"   📉 Legacy LAMP: {min(valid_legacy):.4f}s to {max(valid_legacy):.4f}s"
            )

            # Calculate speedup statistics
            speedups = [
                legacy_time / optimized_time
                for legacy_time, optimized_time in zip(legacy_times, optimized_times)
                if legacy_time > 0 and optimized_time > 0
            ]
            if speedups:
                avg_speedup = np.mean(speedups)
                max_speedup = max(speedups)
                min_speedup = min(speedups)
                std_speedup = np.std(speedups)
                print(
                    f"   🚀 Speedup: Avg={avg_speedup:.2f}x, Max={max_speedup:.2f}x, Min={min_speedup:.2f}x, Std={std_speedup:.2f}"
                )
        else:
            print("   📉 Legacy LAMP: No valid data")

    # Item scaling summary
    if results["item_scaling"]["counts"]:
        counts = results["item_scaling"]["counts"]
        optimized_times = results["item_scaling"]["optimized_times"]
        legacy_times = results["item_scaling"]["legacy_times"]

        print(
            f"\n📊 2. Item Count Scaling Analysis ({min(counts)} to {max(counts)} items):"
        )
        print(
            f"   📈 Optimized LAMP: {min(optimized_times):.4f}s to {max(optimized_times):.4f}s"
        )

        # Growth rate analysis
        if len(optimized_times) > 1:
            growth_rate = (max(optimized_times) / min(optimized_times)) / (
                max(counts) / min(counts)
            )
            print(f"   📊 Growth efficiency: {growth_rate:.3f}")

            # Correlation analysis
            correlation = np.corrcoef(counts, optimized_times)[0, 1]
            print(f"   🔗 Correlation coefficient: {correlation:.3f}")

            # Exponential growth check
            if len(counts) > 3:
                try:
                    log_times = np.log(optimized_times)
                    exp_correlation = np.corrcoef(counts, log_times)[0, 1]
                    if exp_correlation > 0.9:
                        complexity = "O(2^n) - Exponential"
                    elif correlation > 0.9:
                        complexity = "O(n) - Linear"
                    else:
                        complexity = "O(n²) - Polynomial"
                    print(f"   ⚡ Estimated complexity: {complexity}")
                except Exception:
                    pass

        if any(t > 0 for t in legacy_times):
            valid_legacy = [t for t in legacy_times if t > 0]
            print(
                f"   📉 Legacy LAMP: {min(valid_legacy):.4f}s to {max(valid_legacy):.4f}s"
            )

            speedups = [
                legacy_time / optimized_time
                for legacy_time, optimized_time in zip(legacy_times, optimized_times)
                if legacy_time > 0 and optimized_time > 0
            ]
            if speedups:
                avg_speedup = np.mean(speedups)
                max_speedup = max(speedups)
                min_speedup = min(speedups)
                std_speedup = np.std(speedups)
                print(
                    f"   🚀 Speedup: Avg={avg_speedup:.2f}x, Max={max_speedup:.2f}x, Min={min_speedup:.2f}x, Std={std_speedup:.2f}"
                )
        else:
            print("   📉 Legacy LAMP: No valid data")

    # Density scaling summary
    if results["density_scaling"]["densities"]:
        densities = results["density_scaling"]["densities"]
        optimized_times = results["density_scaling"]["optimized_times"]
        legacy_times = results["density_scaling"]["legacy_times"]

        print(
            f"\n📊 3. Density Scaling Analysis ({min(densities):.2f} to {max(densities):.2f}):"
        )
        print(
            f"   📈 Optimized LAMP: {min(optimized_times):.4f}s to {max(optimized_times):.4f}s"
        )

        # Growth rate analysis (density often shows exponential behavior)
        if len(optimized_times) > 1:
            growth_rate = (max(optimized_times) / min(optimized_times)) / (
                max(densities) / min(densities)
            )
            print(f"   📊 Growth efficiency: {growth_rate:.3f}")

            # Power law analysis
            if all(d > 0 for d in densities) and all(t > 0 for t in optimized_times):
                try:
                    log_densities = np.log(densities)
                    log_times = np.log(optimized_times)
                    power_correlation = np.corrcoef(log_densities, log_times)[0, 1]
                    if power_correlation > 0.8:
                        z_power = np.polyfit(log_densities, log_times, 1)
                        power_exponent = z_power[0]
                        print(f"   🔗 Power law correlation: {power_correlation:.3f}")
                        print(
                            f"   ⚡ Power law exponent: {power_exponent:.2f} (time ∝ density^{power_exponent:.2f})"
                        )
                except Exception:
                    pass

        if any(t > 0 for t in legacy_times):
            valid_legacy = [t for t in legacy_times if t > 0]
            print(
                f"   📉 Legacy LAMP: {min(valid_legacy):.4f}s to {max(valid_legacy):.4f}s"
            )

            speedups = [
                legacy_time / optimized_time
                for legacy_time, optimized_time in zip(legacy_times, optimized_times)
                if legacy_time > 0 and optimized_time > 0
            ]
            if speedups:
                avg_speedup = np.mean(speedups)
                max_speedup = max(speedups)
                min_speedup = min(speedups)
                std_speedup = np.std(speedups)
                print(
                    f"   🚀 Speedup: Avg={avg_speedup:.2f}x, Max={max_speedup:.2f}x, Min={min_speedup:.2f}x, Std={std_speedup:.2f}"
                )
        else:
            print("   📉 Legacy LAMP: No valid data")


if __name__ == "__main__":
    print(
        "This module provides graph creation utilities for LAMP performance analysis."
    )
    print(
        "Use the functions create_performance_graphs() and print_performance_summary()"
    )
    print("from your test code to generate performance analysis graphs.")
