from genespeciestreeverdict.benchmark import (
    BenchmarkConfig,
    benchmark_summary,
    load_benchmark_config,
    run_manifest_benchmark,
    run_native_benchmark,
)


def test_native_benchmark_known_truth_states(tmp_path):
    config = BenchmarkConfig(
        seed=17,
        replicates=1,
        taxa=(8,),
        scenarios=(
            "concordant_singlecopy",
            "multicopy_duplication",
            "missing_taxa",
        ),
    )
    outdir = tmp_path / "benchmark"
    records = run_native_benchmark(config, outdir)

    assert len(records) == 3
    by_scenario = {record.scenario: record for record in records}
    assert by_scenario["concordant_singlecopy"].observed_status == "PASS"
    assert by_scenario["multicopy_duplication"].observed_status == "RESOLVE"
    assert by_scenario["missing_taxa"].observed_status == "REVIEW"

    summary = benchmark_summary(records)
    assert summary["resolve_sensitivity_multicopy_truth"] == 1.0
    assert summary["false_resolve_rate_no_duplication_truth"] == 0.0

    for filename in (
        "locus_results.tsv",
        "truth_manifest.tsv",
        "failure_cases.tsv",
        "scenario_metrics.tsv",
        "class_metrics.tsv",
        "confusion_matrix.tsv",
        "performance_summary.tsv",
        "benchmark_report.md",
        "scenario_accuracy.png",
    ):
        assert (outdir / filename).exists()


def test_native_benchmark_is_seed_reproducible(tmp_path):
    config = BenchmarkConfig(
        seed=99,
        replicates=2,
        taxa=(8,),
        scenarios=("branch_length_only", "singlecopy_discordance"),
    )
    first = run_native_benchmark(config, tmp_path / "first")
    second = run_native_benchmark(config, tmp_path / "second")

    first_signature = [
        (record.benchmark_id, record.seed, record.observed_status, record.normalized_rf)
        for record in first
    ]
    second_signature = [
        (record.benchmark_id, record.seed, record.observed_status, record.normalized_rf)
        for record in second
    ]
    assert first_signature == second_signature


def test_truth_manifest_round_trip(tmp_path):
    config = BenchmarkConfig(
        seed=31,
        replicates=1,
        taxa=(8,),
        scenarios=("concordant_singlecopy", "multicopy_duplication"),
    )
    native_dir = tmp_path / "native"
    native = run_native_benchmark(config, native_dir)
    replay = run_manifest_benchmark(native_dir / "truth_manifest.tsv", tmp_path / "replay")

    assert [record.observed_status for record in replay] == [
        record.observed_status for record in native
    ]
    assert all(record.acceptable for record in replay)


def test_yaml_config_loading(tmp_path):
    config_path = tmp_path / "benchmark.yaml"
    config_path.write_text(
        "seed: 123\n"
        "replicates: 2\n"
        "taxa: [8, 12]\n"
        "scenarios: [concordant_singlecopy, multicopy_duplication]\n"
        "missing_fraction: 0.3\n",
        encoding="utf-8",
    )
    config = load_benchmark_config(config_path)

    assert config.seed == 123
    assert config.replicates == 2
    assert config.taxa == (8, 12)
    assert config.scenarios == ("concordant_singlecopy", "multicopy_duplication")
    assert config.missing_fraction == 0.3
