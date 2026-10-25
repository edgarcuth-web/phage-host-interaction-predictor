"""Verify sampling-control numerics and project hosted transfer from real outputs."""
import argparse
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
from phage_cocktail.cohort import parse_dna_fasta, reverse_complement
from phage_cocktail.host_baseline import predict_host_model
from phage_cocktail.host_study import host_features, verify_resources
from phage_cocktail.evo2_host import host_windows
from phage_cocktail.paired_baseline import canonical_fourmers


def verify(repo, source, run, repeat=None):
    report = read_json_object(run / "report.json")
    plan = read_json_object(run / "protocol.json")
    for name, digest in report["output_sha256"].items():
        if sha256_file(run / name) != digest:
            raise ValueError("changed sampling output")
    resources = verify_resources(plan, source)
    metrics = pd.read_csv(run / "per_host_metrics.csv", float_precision="round_trip")
    historical = pd.read_csv(repo / "data/studies/results/genome_host_retrieval.v1/per_host_metrics.csv",
                             float_precision="round_trip")
    retained = metrics[metrics.model.isin(historical.model.unique())].reset_index(drop=True)
    pd.testing.assert_frame_equal(retained, historical, check_exact=True)
    scores = pd.read_csv(run / "predictions.csv.gz", float_precision="round_trip")
    splits = pd.read_csv(run / "host_splits.csv")
    scenarios = splits.scenario.drop_duplicates().tolist()
    metadata = pd.read_csv(resources[("table", "metadata")], sep=";", index_col=0)
    defenses = pd.read_csv(resources[("table", "defenses")], sep=";", index_col=0)
    features, _ = host_features(sorted(splits.bacterium.unique()), metadata, defenses,
                                plan["categorical_features"])
    genome = pd.read_csv(run / "host_genome_features.csv", index_col="bacterium", float_precision="round_trip")
    frames = {"host_logistic": features, "host_4mer_logistic": genome,
              "hybrid_4mer_logistic": features.join(genome)}
    for budget in plan["sampling_control"]["budgets"]:
        loaded = pd.read_csv(run / f"host_window{budget}_features.csv", index_col="bacterium", float_precision="round_trip")
        # The sampler constructs frames from np.stack (C order), whereas CSV
        # reading builds column blocks. Restore the original layout to reproduce
        # scalar BLAS dot products bit for bit; do not change scores or weights.
        frames[f"host_window{budget}_4mer_logistic"] = pd.DataFrame(
            np.array(loaded.to_numpy(), order="C"), index=loaded.index, columns=loaded.columns)
    reloaded = 0
    with threadpool_limits(limits=1):
        for path in sorted(run.glob("model_*.json")):
            parts = path.stem.split("_", 2)
            scenario = scenarios[int(parts[1])]
            method = parts[2] if len(parts) == 3 else "host_logistic"
            target = scores[scores.scenario.eq(scenario) & scores.model.eq(method)]
            model = read_json_object(path)
            np.testing.assert_array_equal(predict_host_model(model, target, frames[method]), target.score)
            reloaded += 1
    # Independent scalar counting on the first three frozen hosts, using recorded
    # coordinates rather than the sampler. No interaction outcomes select these.
    coordinates = pd.read_csv(run / "sampling_windows.csv")
    pilot = read_json_object(repo / "data/studies/results/nvidia_evo2_host_pilot.v1/report.json")
    checks = 0
    for host in sorted(splits.bacterium.unique())[:3]:
        contigs = sorted(min(c, reverse_complement(c)) for c in parse_dna_fasta(
            resources[("host_genome", host)].read_text(encoding="ascii")))
        for budget in plan["sampling_control"]["budgets"]:
            rows = coordinates[coordinates.bacterium.eq(host) & coordinates.window_budget.eq(budget)]
            sequences = []
            for row in rows.itertuples(index=False):
                s = contigs[row.contig_ordinal][row.start:row.start + row.length]
                s = "".join(b if b in "ACGT" else "N" for b in s)
                assert sha256(s.encode()).hexdigest() == row.input_sha256
                sequences.append(s)
            fasta = "".join(f">window_{i}\n{s}\n" for i,s in enumerate(sequences)).encode()
            np.testing.assert_array_equal(canonical_fourmers(fasta), frames[f"host_window{budget}_4mer_logistic"].loc[host])
            if budget == 1:
                parent = next(j for j in pilot["jobs"] if j["input"]["host"] == host and j["input"]["strand"] == "forward")
                assert rows.iloc[0].input_sha256 == parent["input"]["input_sha256"]
            checks += 1
    observed = pd.read_csv(run / "observed_pairs.csv")
    panels = scores[scores.scenario.eq("grouped_host") & scores.partition.eq("test")]
    enumerated = 0
    for (method, host), values in panels.groupby(["model", "bacterium"]):
        labels = observed[observed.bacterium.eq(host)].set_index("phage").infection.to_dict()
        if len(labels) != plan["expected_candidate_phages"]:
            continue
        boundary = values.score.nlargest(2).iloc[-1]
        higher = values.loc[values.score.gt(boundary), "phage"].tolist()
        tied = values.loc[values.score.eq(boundary), "phage"].tolist()
        trials = [any(labels[p] for p in higher + list(chosen))
                  for chosen in combinations(tied, 2-len(higher))]
        expected = sum(trials)/len(trials)
        actual = metrics.loc[metrics.scenario.eq("grouped_host") & metrics.partition.eq("test")
                             & metrics.model.eq(method) & metrics.bacterium.eq(host), "hit@2"].item()
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-15)
        enumerated += 1
    if repeat is not None:
        other = read_json_object(repeat / "report.json")
        if report != other or any(sha256_file(repeat / n) != d for n,d in report["output_sha256"].items()):
            raise ValueError("repeat report or payload differs")
    reference_path = repo / "data/studies/results/nvidia_evo2_window_feasibility.v1.json"
    reference = read_json_object(reference_path)["response"]
    audit = pd.read_csv(run / "sampling_audit.csv")
    planning = []
    for budget, group in audit.groupby("window_budget"):
        calls = int(group.selected_windows.sum()) * 2
        planning.append({"strategy": f"full_windows_{budget}", "window_budget": int(budget), "hosts": len(group), "two_strand_calls": calls,
                         "hosts_with_fewer_windows_than_budget": int(group.selected_windows.lt(budget).sum()),
                         "median_assembly_fraction": float(group.sampled_assembly_fraction.median()),
                         "estimated_json_bytes_no_reuse_credit": calls * reference["response_bytes"],
                         "serial_server_seconds_projection": calls * reference["elapsed_ms"] / 1000,
                         "reference_report_sha256": sha256_file(reference_path),
                         "scope": "Projection from one real 8192bp response, not a measured full run, invoice or schedule"})
    # Refine the old upper-bound projection as well: its original 32-window
    # protocol includes short tails, unlike the full-window composition control.
    legacy_counts, legacy_bases, legacy_fractions = [], [], []
    sizes = audit[audit.window_budget.eq(1)].set_index("bacterium").assembly_bases.to_dict()
    for host in sorted(splits.bacterium.unique()):
        windows = host_windows(resources[("host_genome", host)].read_bytes(), window_bp=8192, max_windows_per_host=32)
        bases = sum(len(s) for _,_,s in windows)
        legacy_counts.append(len(windows)); legacy_bases.append(bases)
        legacy_fractions.append(bases / sizes[host])
    equivalents = 2 * sum(legacy_bases) / 8192
    planning.append({"strategy": "original_32_including_short_tails", "window_budget": 32,
                     "hosts": len(legacy_counts), "two_strand_calls": 2 * sum(legacy_counts),
                     "hosts_with_fewer_windows_than_budget": sum(n < 32 for n in legacy_counts),
                     "median_assembly_fraction": float(np.median(legacy_fractions)),
                     "estimated_json_bytes_no_reuse_credit": int(round(equivalents * reference["response_bytes"])),
                     "serial_server_seconds_projection": equivalents * reference["elapsed_ms"] / 1000,
                     "reference_report_sha256": sha256_file(reference_path),
                     "scope": "Linear input-base projection including short tails; unmeasured service behavior, not an invoice or schedule"})
    return ({"status": "passed", "historical_six_methods_per_host_exact": True,
             "saved_models_reloaded_with_exact_scores": reloaded,
             "independent_scalar_sampled_features": checks, "pilot_one_window_inputs_match": True,
             "independent_primary_top2_tie_panels": enumerated,
             "second_run_all_report_and_output_hashes_equal": repeat is not None,
             "report_sha256": sha256_file(run / "report.json"), "new_nvidia_inference_calls": 0},
            pd.DataFrame(planning))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--repeat", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--budget-out", type=Path, required=True)
    args = parser.parse_args()
    qa, planning = verify(Path(__file__).resolve().parents[1], args.data_root, args.run, args.repeat)
    atomic_write_json(qa, args.out)
    planning.to_csv(args.budget_out, index=False)
    print(qa)


if __name__ == "__main__":
    main()
