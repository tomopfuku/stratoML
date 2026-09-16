#!/usr/bin/env python3

import argparse
import csv
import json
import math
import os
import random
import re
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-sims", type=int, default=30, help="Number of simulated matrices to analyze.")
    parser.add_argument("--out-dir", default="sim_likelihood_batch", help="Output directory for matrices, logs, tables, and plots.")
    parser.add_argument("--tree", default="sim.tre", help="Tree file passed to both simulation and likelihood scripts.")
    parser.add_argument("--ranges", default="sim.ranges", help="Stratigraphic ranges file passed to both scripts.")
    parser.add_argument("--num-traits", type=int, default=40, help="Number of characters to simulate.")
    parser.add_argument("--max-states", type=int, default=3, help="Maximum number of states passed to simulate_mfc.py.")
    parser.add_argument("--min-clado-fraction", type=float, default=0.3, help="Minimum simulated cladogenetic substitution fraction.")
    parser.add_argument("--max-clado-fraction", type=float, default=0.7, help="Maximum simulated cladogenetic substitution fraction.")
    parser.add_argument("--min-gain-rate", type=float, default=0.02, help="Minimum simulated gain rate.")
    parser.add_argument("--max-gain-rate", type=float, default=0.15, help="Maximum simulated gain rate.")
    parser.add_argument("--strat-model", default="hr97", help="Stratigraphic model flag passed to main_single_tree_like3.py.")
    parser.add_argument("--morph-model", default="mfc2", help="Morphologic model flag passed to main_single_tree_like3.py.")
    parser.add_argument("--force", action="store_true", help="Regenerate matrices and rerun likelihoods even if outputs already exist.")
    parser.add_argument("--stop-on-error", action="store_true", help="Stop at the first failed simulation or likelihood run.")
    return parser.parse_args(argv)


def repo_root():
    return Path(__file__).resolve().parents[2]


def work_dir():
    return Path(__file__).resolve().parent


def numbers_from_line(line):
    return [float(match.group(0)) for match in FLOAT_RE.finditer(line)]


def parse_likelihood_output(text):
    parsed = {}
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("BDS rates:"):
            nums = numbers_from_line(line)
            if nums:
                parsed["bds_lambda"] = nums[0]
        elif line.startswith("no jump params:"):
            nums = numbers_from_line(line)
            if len(nums) >= 3:
                parsed["reconstructed_gain_rate"] = nums[0]
                parsed["reconstructed_loss_rate"] = nums[1]
                parsed["reconstructed_clado_fraction"] = nums[2] * 10.0
        elif line.startswith("GAIN:"):
            nums = numbers_from_line(line)
            if nums:
                parsed["reconstructed_gain_rate"] = nums[0]
        elif line.startswith("LOSS:"):
            nums = numbers_from_line(line)
            if nums:
                parsed["reconstructed_loss_rate"] = nums[0]
        elif line.startswith("LAMBDA SUB:"):
            nums = numbers_from_line(line)
            if nums:
                parsed["reconstructed_lambda_sub"] = nums[0]
        elif line.startswith("no jump AIC:"):
            nums = numbers_from_line(line)
            if nums:
                parsed["clado_aic"] = nums[0]
        elif line.startswith("no cladogenesis AIC:"):
            nums = numbers_from_line(line)
            if nums:
                parsed["no_clado_aic"] = nums[0]
        elif line.startswith("cladogenesis AIC weight:"):
            nums = numbers_from_line(line)
            if nums:
                parsed["clado_aic_weight"] = nums[0]
        elif line.startswith("no cladogenesis AIC weight:"):
            nums = numbers_from_line(line)
            if nums:
                parsed["no_clado_aic_weight"] = nums[0]

    if "reconstructed_clado_fraction" not in parsed:
        if parsed.get("bds_lambda") and "reconstructed_lambda_sub" in parsed:
            parsed["reconstructed_clado_fraction"] = parsed["reconstructed_lambda_sub"] / parsed["bds_lambda"]

    required = [
        "reconstructed_gain_rate",
        "reconstructed_loss_rate",
        "reconstructed_clado_fraction",
        "clado_aic",
        "no_clado_aic",
    ]
    missing = [key for key in required if key not in parsed]
    if missing:
        raise ValueError("Could not parse likelihood output fields: " + ", ".join(missing))

    parsed["clado_lower_aic"] = parsed["clado_aic"] < parsed["no_clado_aic"]
    parsed["no_clado_lower_aic"] = parsed["no_clado_aic"] < parsed["clado_aic"]
    parsed["aic_tie"] = parsed["clado_aic"] == parsed["no_clado_aic"]
    return parsed


def run_command(command, cwd, stdout_path, stderr_path):
    with open(stdout_path, "w") as stdout_file, open(stderr_path, "w") as stderr_file:
        return subprocess.run(
            command,
            cwd=cwd,
            stdout=stdout_file,
            stderr=stderr_file,
            text=True,
        )


def draw_simulation_parameters(args):
    gain_rate = random.uniform(args.min_gain_rate, args.max_gain_rate)
    clado_fraction = random.uniform(args.min_clado_fraction, args.max_clado_fraction)
    return {
        "true_gain_rate": gain_rate,
        "true_loss_rate": gain_rate,
        "true_clado_fraction": clado_fraction,
    }


def read_simulation_parameters(path):
    with open(path) as infile:
        return json.load(infile)


def write_simulation_parameters(path, params):
    with open(path, "w") as outfile:
        json.dump(params, outfile, indent=2, sort_keys=True)


def write_summary(rows, summary_path):
    fieldnames = [
        "sim_index",
        "status",
        "matrix_path",
        "true_gain_rate",
        "true_loss_rate",
        "true_clado_fraction",
        "reconstructed_gain_rate",
        "reconstructed_loss_rate",
        "reconstructed_clado_fraction",
        "reconstructed_lambda_sub",
        "bds_lambda",
        "clado_aic",
        "no_clado_aic",
        "clado_aic_weight",
        "no_clado_aic_weight",
        "clado_lower_aic",
        "no_clado_lower_aic",
        "aic_tie",
        "error",
    ]
    with open(summary_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_successful_rows(summary_path):
    if not summary_path.exists():
        return []
    rows = []
    with open(summary_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row.get("status") == "ok":
                rows.append(row)
    return rows


def as_float(row, key):
    value = row.get(key, "")
    if value == "" or value is None:
        return math.nan
    return float(value)


def save_plot(fig, plot_stem):
    fig.savefig(plot_stem.with_suffix(".png"), dpi=300)
    fig.savefig(plot_stem.with_suffix(".svg"))


def plot_rate_scatters(rows, plot_stem):
    metrics = [
        ("gain_rate", "Gain rate"),
        ("loss_rate", "Loss rate"),
        ("clado_fraction", "Cladogenesis fraction"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)

    for ax, (metric, label) in zip(axes, metrics):
        true_key = f"true_{metric}"
        recon_key = f"reconstructed_{metric}"
        x = [as_float(row, true_key) for row in rows]
        y = [as_float(row, recon_key) for row in rows]
        ax.scatter(x, y, alpha=0.75, edgecolor="black", linewidth=0.25)
        vals = [val for val in x + y if math.isfinite(val)]
        if vals:
            lo = min(vals)
            hi = max(vals)
            if lo == hi:
                pad = max(abs(lo) * 0.1, 0.01)
                lo -= pad
                hi += pad
            ax.plot([lo, hi], [lo, hi], color="0.35", linestyle="--", linewidth=1)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
        ax.set_title(label)
        ax.set_xlabel("Simulated")
        ax.set_ylabel("Reconstructed")

    fig.suptitle("Simulated vs reconstructed parameters")
    save_plot(fig, plot_stem)
    plt.close(fig)


def plot_aic_bar(rows, plot_stem):
    total = len(rows)
    clado_lower = sum(row.get("clado_lower_aic") in (True, "True", "true", "1") for row in rows)
    not_clado_lower = total - clado_lower
    clado_prop = clado_lower / total if total else 0.0
    not_clado_prop = not_clado_lower / total if total else 0.0

    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    ax.bar(["AIC comparison"], [clado_prop], label="Cladogenesis lower AIC", color="#4c78a8")
    ax.bar(
        ["AIC comparison"],
        [not_clado_prop],
        bottom=[clado_prop],
        label="No cladogenesis lower or tied",
        color="#f58518",
    )
    ax.set_ylim(0, 1)
    ax.set_ylabel("Proportion of simulations")
    ax.legend(loc="upper right")
    ax.set_title("AIC support across simulations")
    save_plot(fig, plot_stem)
    plt.close(fig)


def make_plots(summary_path, plots_dir):
    rows = load_successful_rows(summary_path)
    if not rows:
        print("No successful likelihood runs to plot.")
        return
    plots_dir.mkdir(parents=True, exist_ok=True)
    plot_rate_scatters(rows, plots_dir / "simulated_vs_reconstructed_rates")
    plot_aic_bar(rows, plots_dir / "aic_model_support_proportions")


def main(argv=None):
    args = parse_args(argv)
    if args.min_gain_rate > args.max_gain_rate:
        raise ValueError("--min-gain-rate cannot exceed --max-gain-rate")
    if args.min_clado_fraction > args.max_clado_fraction:
        raise ValueError("--min-clado-fraction cannot exceed --max-clado-fraction")

    base_dir = work_dir()
    root_dir = repo_root()
    output_dir = base_dir / args.out_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "simulation_likelihood_summary.csv"
    plots_dir = output_dir / "plots"

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
            sim_params = draw_simulation_parameters(args)
        elif params_path.exists():
            sim_params = read_simulation_parameters(params_path)
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
                #print(" ".join(sim_command))
                #exit()
                result = run_command(sim_command, base_dir, matrix_path, sim_stderr)
                if result.returncode != 0:
                    raise RuntimeError(f"simulate_mfc.py failed with exit code {result.returncode}")
                write_simulation_parameters(params_path, sim_params)

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
                #print(" ".join(like_command))
                result = run_command(like_command, base_dir, like_stdout, like_stderr)
                if result.returncode != 0:
                    raise RuntimeError(f"main_single_tree_like3.py failed with exit code {result.returncode}")

            parsed = parse_likelihood_output(like_stdout.read_text())
            row.update(parsed)
            row["status"] = "ok"
            print(f"replicate {sim_index:03d}: ok")
        except Exception as exc:
            row["error"] = str(exc)
            print(f"replicate {sim_index:03d}: failed: {exc}")
            if args.stop_on_error:
                rows.append(row)
                write_summary(rows, summary_path)
                raise

        rows.append(row)
        write_summary(rows, summary_path)

    make_plots(summary_path, plots_dir)
    print(f"Wrote summary to {summary_path}")
    print(f"Wrote plots to {plots_dir}")


if __name__ == "__main__":
    main()
