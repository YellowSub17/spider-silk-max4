"""
Plot unnormalized raw intensity curves for all protein-vs-blank pairs.
This shows why normalization is critical: without it, instrumental/geometric
differences overwhelm the true protein signal.
"""

import os
import matplotlib.pyplot as plt
import numpy as np

import compare_groups as cg
import src


def plot_unnormalized_pairs(out_dir="figures/01_normalization_overview"):
    """Plot all 6 protein-vs-blank pairs with NO normalization applied."""
    os.makedirs(out_dir, exist_ok=True)

    # Same 6 pairs as in protein_vs_buffer.py
    pairs = [
        (f"{name} vs {cg.BLANK_OF[name]} blank", name, cg.BLANK_OF[name])
        for name in cg.GROUPS
        if cg.SAMPLE_INFO[name]["protein"] and cg.BLANK_OF[name] is not None
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
    axes = axes.flatten()

    for ax, (title, protein_name, blank_name) in zip(axes, pairs):
        # Load raw data WITHOUT any normalization
        protein = cg.load_group(protein_name)
        blank = cg.load_group(blank_name)

        # Get kept frames (quality-filtered by r² streak metric, but NOT normalized)
        keep_p = cg.kept_mask(protein)
        keep_b = cg.kept_mask(blank)

        # Raw intensity curves (no normalization applied)
        qs = protein.qs
        Is_p_raw = protein.Is[keep_p]
        Is_b_raw = blank.Is[keep_b]

        I_p_mean = Is_p_raw.mean(axis=0)
        I_b_mean = Is_b_raw.mean(axis=0)
        std_p = Is_p_raw.std(axis=0)
        std_b = Is_b_raw.std(axis=0)

        # Plot
        ax.plot(qs, I_p_mean, color="tab:red", linewidth=2, label=f"{protein_name} (n={keep_p.sum()})")
        ax.fill_between(qs, I_p_mean - std_p, I_p_mean + std_p, color="tab:red", alpha=0.2, linewidth=0)
        ax.plot(qs, I_b_mean, color="tab:blue", linewidth=2, label=f"{blank_name} (n={keep_b.sum()})")
        ax.fill_between(qs, I_b_mean - std_b, I_b_mean + std_b, color="tab:blue", alpha=0.2, linewidth=0)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(True, alpha=0.3, which="both")

    # Shared labels
    fig.text(0.5, 0.02, "q [1/A]", ha="center", fontsize=12)
    fig.text(0.02, 0.5, "Raw Intensity (unnormalized)", va="center", rotation="vertical", fontsize=12)

    fig.suptitle("Raw (Unnormalized) Scattering Curves for All Protein-vs-Blank Pairs",
                 fontsize=14, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0.03, 0.05, 1, 0.96])

    out_path = os.path.join(out_dir, "sample_blank_curves__raw_unnormalized.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}")

    # Also create a side-by-side comparison: one unnormalized pair + three normalized versions
    _plot_normalization_comparison(pairs[2], out_dir)  # Y2F_h2o_kpi (strongest signal)


def _plot_normalization_comparison(pair_info, out_dir):
    """Show one pair (Y2F_h2o_kpi) in all 4 states: raw + 3 normalizations."""
    title, protein_name, blank_name = pair_info

    protein = cg.load_group(protein_name)
    blank = cg.load_group(blank_name)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharex=True)
    qs = protein.qs

    # --- Panel 1: Raw unnormalized ---
    keep_p = cg.kept_mask(protein)
    keep_b = cg.kept_mask(blank)
    I_p_raw = protein.Is[keep_p].mean(axis=0)
    I_b_raw = blank.Is[keep_b].mean(axis=0)
    std_p_raw = protein.Is[keep_p].std(axis=0)
    std_b_raw = blank.Is[keep_b].std(axis=0)

    axes[0].plot(qs, I_p_raw, color="tab:red", linewidth=2.5, label=protein_name)
    axes[0].fill_between(qs, I_p_raw - std_p_raw, I_p_raw + std_p_raw, color="tab:red", alpha=0.2)
    axes[0].plot(qs, I_b_raw, color="tab:blue", linewidth=2.5, label=blank_name)
    axes[0].fill_between(qs, I_b_raw - std_b_raw, I_b_raw + std_b_raw, color="tab:blue", alpha=0.2)
    axes[0].set_title("RAW (unnormalized)\n[Very different scales!]", fontweight="bold", fontsize=11)
    axes[0].legend(fontsize=9)

    # --- Panels 2-4: Three normalized versions ---
    for i, (norm_name, norm_label) in enumerate([
        ("linear_saxs", "linear_saxs norm"),
        ("neutral", "neutral norm"),
        ("i0", "i0 norm"),
    ]):
        protein_norm = cg.load_group(protein_name)
        blank_norm = cg.load_group(blank_name)

        qs_p, Is_p, _, _ = cg.kept_frames(protein_norm, norm=norm_name)
        qs_b, Is_b, _, _ = cg.kept_frames(blank_norm, norm=norm_name)

        I_p_mean = Is_p.mean(axis=0)
        I_b_mean = Is_b.mean(axis=0)
        std_p = Is_p.std(axis=0)
        std_b = Is_b.std(axis=0)

        axes[i+1].plot(qs_p, I_p_mean, color="tab:red", linewidth=2.5, label=protein_name)
        axes[i+1].fill_between(qs_p, I_p_mean - std_p, I_p_mean + std_p, color="tab:red", alpha=0.2)
        axes[i+1].plot(qs_b, I_b_mean, color="tab:blue", linewidth=2.5, label=blank_name)
        axes[i+1].fill_between(qs_b, I_b_mean - std_b, I_b_mean + std_b, color="tab:blue", alpha=0.2)
        axes[i+1].set_title(f"{norm_label}\n[Curves now comparable]", fontweight="bold", fontsize=11)
        axes[i+1].legend(fontsize=9)

    # All share log scales
    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(src.det.FULL_QRANGE)
        ax.grid(True, alpha=0.3, which="both")

    fig.text(0.5, 0.02, "q [1/A]", ha="center", fontsize=12)
    fig.text(0.02, 0.5, "Intensity", va="center", rotation="vertical", fontsize=12)

    fig.suptitle(f"Why Normalization Matters: {title}", fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0.03, 0.05, 1, 0.96])

    out_path = os.path.join(out_dir, "normalization_comparison__y2f_h2o_kpi.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    plot_unnormalized_pairs()
    plt.show()
