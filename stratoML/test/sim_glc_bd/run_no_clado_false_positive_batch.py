#!/usr/bin/env python3

import argparse
import math
import os
import random
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_sim_likelihood_batch as batch


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-sims", type=int, default=100, help="Number of simulated matrices to analyze.")
    parser.add_argument("--out-dir", default="sim_no_clado_false_positive_batch", help="Output directory for matrices, logs, tables, and plots.")
    parser.add_argument("--tree", default="sim.tre", help="Tree file passed to both simulation and likelihood scripts.")
    parser.add_argument("--ranges", default="sim.ranges", help="Stratigraphic ranges file passed to both scripts.")
    parser.add_argument("--num-traits", type=int, default=40, help="Number of characters to simulate.")
    parser.add_argument("--max-states", type=int, default=3, help="Maximum number of states passed to simulate_mfc.py.")
    parser.add_argument("--min-gain-rate", type=float, default=0.02, help="Minimum simulated gain rate.")
    parser.add_argument("--max-gain-rate", type=float, default=0.15, help="Maximum simulated gain rate.")
    parser.add_argument("--strat-model", default="hr97", help="Stratigraphic model flag passed to main_single_tree_like3.py.")
    parser.add_argument("--morph-model", default="mfc2", help="Morphologic model flag passed to main_single_tree_like3.py.")
    parser.add_argument("--force", action="store_true", help="Regenerate matrices and rerun likelihoods even if outputs already exist.")
    parser.add_argument("--stop-on-error", action="store_true", help="Stop at the first failed simulation or likelihood run.")
    return parser.parse_args(argv)


def draw_no_clado_simulation_parameters(args):
    gain_rate = random.uniform(args.min_gain_rate, args.max_gain_rate)
    return {
        "true_gain_rate": gain_rate,
        "true_loss_rate": gain_rate,
        "true_clado_fraction": 0.0,
    }


def is_true(value):
    return value in (True, "True", "true", "1")


def plot_false_positive_rate(rows, plot_stem):
    total = len(rows)
    false_positives = sum(is_true(row.get("clado_lower_aic")) for row in rows)
    false_positive_rate = false_positives / total if total else math.nan
    not_false_positive_rate = 1.0 - false_positive_rate if total else math.nan

    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    ax.bar(
        ["No-cladogenesis simulations"],
        [false_positive_rate],
        label="False positive: cladogenesis lower AIC",
        color="#d62728",
    )
    ax.bar(
        ["No-cladogenesis simulations"],
        [not_false_positive_rate],
        bottom=[false_positive_rate],
        label="No false positive",
        color="#4c78a8",
    )
    ax.set_ylim(0, 1)
    ax.set_ylabel("Proportion of simulations")
    ax.set_title(f"False positive rate = {false_positive_rate:.3f} ({false_positives}/{total})")
    ax.legend(loc="upper right")
    batch.save_plot(fig, plot_stem)
    plt.close(fig)


def plot_false_positive_clado_fraction_histogram(rows, plot_stem):
    false_positive_rows = [
        row for row in rows
        if is_true(row.get("clado_lower_aic"))
    ]
    clado_fractions = [
        batch.as_float(row, "reconstructed_clado_fraction")
        for row in false_positive_rows
    ]
    clado_fractions = [
        value for value in clado_fractions
        if math.isfinite(value)
    ]

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    if clado_fractions:
        ax.hist(clado_fractions, bins=12, color="#d62728", edgecolor="black", alpha=0.8)
    else:
        ax.text(
            0.5,
            0.5,
            "No false positives",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    ax.set_xlabel("Estimated cladogenesis fraction")
    ax.set_ylabel("Count")
    ax.set_title("Estimated cladogenesis fractions among false positives")
    batch.save_plot(fig, plot_stem)
    plt.close(fig)


def make_plots(summary_path, plots_dir):
    rows = batch.load_successful_rows(summary_path)
    if not rows:
        print("No successful likelihood runs to plot.")
        return
    plots_dir.mkdir(parents=True, exist_ok=True)
    plot_false_positive_rate(rows, plots_dir / "false_positive_rate")
    plot_false_positive_clado_fraction_histogram(
        rows,
        plots_dir / "false_positive_clado_fraction_histogram",
    )


def main(argv=None):
    args = parse_args(argv)
    if args.min_gain_rate > args.max_gain_rate:
        raise ValueError("--min-gain-rate cannot exceed --max-gain-rate")

    base_dir = batch.work_dir()
    root_dir = batch.repo_root()
    output_dir = base_dir / args.out_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "simulation_likelihood_summary.csv"
    plots_dir = output_dir / "plots"

    if summary_path.exists() and not args.force:
        print(f"Found existing summary at {summary_path}; recreating plots only.")
        make_plots(summary_path, plots_dir)
        print(f"Wrote plots to {plots_dir}")
        return

    rows = []
    python = sys.executable
    simulate_script = root_dir / "simulate_mfc.py"
    likelihood_script = root_dir / "main_single_tree_like3.py"

    for sim_index in range(1, args.num_sims + 1):
        replicate_dir = output_dir / f"replicate_{sim_index:03d}"
        replicate_dir.mkdir(parents=True, exist_ok=True)
        matrix_path = replicate_dir / "matrix.fa"
        params_path = replicate_dir / "simulated_parameters.json"
        sim_stderr = replicate_dir / "simulate.stderr.log"
        like_stdout = replicate_dir / "likelihood.stdout.log"
        like_stderr = replicate_dir / "likelihood.stderr.log"

        if args.force or not matrix_path.exists():
            sim_params = draw_no_clado_simulation_parameters(args)
        elif params_path.exists():
            sim_params = batch.read_simulation_parameters(params_path)
        else:
            raise FileNotFoundError(
                f"{matrix_path} already exists but {params_path} is missing; "
                "rerun with --force to regenerate matched simulation parameters."
            )

        row = {
            "sim_index": sim_index,
            "status": "failed",
            "matrix_path": str(matrix_path),
            "true_gain_rate": sim_params["true_gain_rate"],
            "true_loss_rate": sim_params["true_loss_rate"],
            "true_clado_fraction": sim_params["true_clado_fraction"],
            "error": "",
        }

        try:
            if args.force or not matrix_path.exists():
                sim_command = [
                    python,
                    str(simulate_script),
                    args.tree,
                    args.ranges,
                    str(args.num_traits),
                    str(args.max_states),
                    str(sim_params["true_clado_fraction"]),
                    "--gain-rate",
                    str(sim_params["true_gain_rate"]),
                    "--loss-rate",
                    str(sim_params["true_loss_rate"]),
                ]
                result = batch.run_command(sim_command, base_dir, matrix_path, sim_stderr)
                if result.returncode != 0:
                    raise RuntimeError(f"simulate_mfc.py failed with exit code {result.returncode}")
                batch.write_simulation_parameters(params_path, sim_params)

            if args.force or not like_stdout.exists():
                like_command = [
                    python,
                    str(likelihood_script),
                    "--trees",
                    args.tree,
                    "--traits",
                    str(matrix_path),
                    "--strat-data",
                    args.ranges,
                    "--strat-model",
                    args.strat_model,
                    "--morph-model",
                    args.morph_model,
                ]
                result = batch.run_command(like_command, base_dir, like_stdout, like_stderr)
                if result.returncode != 0:
                    raise RuntimeError(f"main_single_tree_like3.py failed with exit code {result.returncode}")

            parsed = batch.parse_likelihood_output(like_stdout.read_text())
            row.update(parsed)
            row["status"] = "ok"
            print(f"replicate {sim_index:03d}: ok")
        except Exception as exc:
            row["error"] = str(exc)
            print(f"replicate {sim_index:03d}: failed: {exc}")
            if args.stop_on_error:
                rows.append(row)
                batch.write_summary(rows, summary_path)
                raise

        rows.append(row)
        batch.write_summary(rows, summary_path)

    make_plots(summary_path, plots_dir)
    print(f"Wrote summary to {summary_path}")
    print(f"Wrote plots to {plots_dir}")


if __name__ == "__main__":
    main()
