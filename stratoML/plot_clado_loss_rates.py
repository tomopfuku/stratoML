import sys
import os

#mpl_config_dir = os.environ.get("MPLCONFIGDIR")
#if mpl_config_dir is None or not os.access(mpl_config_dir, os.W_OK):
#    os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def usage():
    print(
        "usage: "
        + sys.argv[0]
        + " <mcmc samples csv> <global cladogenesis rate> [output plot]"
    )


def main():
    if len(sys.argv) not in (3, 4):
        usage()
        sys.exit(1)

    chain_file = sys.argv[1]
    try:
        global_clado_rate = float(sys.argv[2])
    except ValueError:
        print(f"global cladogenesis rate must be numeric; received {sys.argv[2]!r}")
        sys.exit(1)

    samples = pd.read_csv(chain_file)
    required_columns = ["param_1", "param_2"]
    missing_columns = [col for col in required_columns if col not in samples.columns]
    if missing_columns:
        print("MCMC chain is missing expected columns: " + ", ".join(missing_columns))
        sys.exit(1)

    loss_rates = samples["param_1"].to_numpy(dtype=np.float64)
    clado_loss_rates = samples["param_2"].to_numpy(dtype=np.float64)
    clado_loss_rates = clado_loss_rates * 10.0 * global_clado_rate

    loss_median = np.median(loss_rates)
    clado_loss_median = np.median(clado_loss_rates)
    total_loss_rates = loss_rates + clado_loss_rates
    valid_totals = total_loss_rates > 0.0
    clado_proportions = clado_loss_rates[valid_totals] / total_loss_rates[valid_totals]
    clado_proportion_median = np.median(clado_proportions)

    fig, ax = plt.subplots(figsize=(9, 5))
    bins = np.histogram_bin_edges(np.concatenate([loss_rates, clado_loss_rates]), bins="auto")
    ax.hist(
        loss_rates,
        bins=bins,
        alpha=0.55,
        density=True,
        label=f"Anagenetic loss rate median={loss_median:.5g}",
        color="steelblue",
    )
    ax.hist(
        clado_loss_rates,
        bins=bins,
        alpha=0.55,
        density=True,
        label=f"Cladogenetic loss rate median={clado_loss_median:.5g}",
        color="darkorange",
    )
    ax.axvline(loss_median, color="steelblue", linestyle=":", linewidth=2)
    ax.axvline(clado_loss_median, color="darkorange", linestyle=":", linewidth=2)
    ax.text(
        0.98,
        0.95,
        f"Median cladogenetic proportion: {clado_proportion_median:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
    )

    ax.set_xlabel("Rate")
    ax.set_ylabel("Density")
    ax.legend()
    fig.tight_layout()

    if len(sys.argv) == 4:
        fig.savefig(sys.argv[3], dpi=300)
        print(f"Saved plot to '{sys.argv[3]}'.")
    else:
        plt.show()


if __name__ == "__main__":
    main()
