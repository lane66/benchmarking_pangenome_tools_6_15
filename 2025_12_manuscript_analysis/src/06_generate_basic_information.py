from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
PANGENOME_PARSED_DIR = PROJECT_DIR / "results" / "pangenome_parsed"
OUTPUT_DIR = PANGENOME_PARSED_DIR / "06_aggregated_data"
OUTPUT_FILE = OUTPUT_DIR / "basic_information.csv"

RUNS = [
    ("02_all_assembly", "0%"),
    ("03_11per_complete", "11%"),
    ("04_50per_complete", "50%"),
    ("05_89per_complete", "89%"),
    ("01_all_complete", "100%"),
]

METHOD_BY_TOOL = {
    "panaroo": "Gene annotation and clustering",
    "roary": "Gene annotation and clustering",
    "PPanGGOLiN": "Gene annotation and clustering",
    "ggcaller": "Gene annotation and clustering",
    "bifrost": "Unitig calling",
    "cuttlefish": "Unitig calling",
    "pangraph": "Multiple sequence alignments",
}

TOOL_ORDER = [
    "panaroo",
    "roary",
    "PPanGGOLiN",
    "ggcaller",
    "bifrost",
    "cuttlefish",
    "pangraph",
]


def read_run_summary(run_dir, pct_complete):
    run_path = PANGENOME_PARSED_DIR / run_dir
    basic_path = run_path / "Tool_basic_informations.tsv"
    length_path = run_path / "length_summary_stats.csv"

    basic = pd.read_csv(basic_path, sep="\t")
    length = pd.read_csv(length_path)

    summary = basic.merge(
        length,
        left_on="Method",
        right_on="Tool",
        how="left",
        suffixes=("", "_length"),
    )

    summary["Length Median"] = summary["Median"]
    summary["Q3-Q1"] = summary["Q3"] - summary["Q1"]
    summary["Percentage of Complete genome "] = pct_complete
    summary["Method_name"] = summary["Method"].map(METHOD_BY_TOOL).fillna("Unknown")
    summary["Tool_order"] = summary["Method"].map(
        {tool: i for i, tool in enumerate(TOOL_ORDER)}
    )
    summary["Run_order"] = RUNS.index((run_dir, pct_complete))

    return summary[
        [
            "Method",
            "Node Count",
            "Edge Count",
            "Subgraph Count",
            "Length Median",
            "Q3",
            "Q3-Q1",
            "Percentage of Complete genome ",
            "Method_name",
            "Run_order",
            "Tool_order",
        ]
    ].rename(columns={"Method": "Tool", "Method_name": "Method"})


def main():
    all_summaries = pd.concat(
        [read_run_summary(run_dir, pct_complete) for run_dir, pct_complete in RUNS],
        ignore_index=True,
    )

    all_summaries = all_summaries.sort_values(["Run_order", "Tool_order", "Tool"])
    all_summaries = all_summaries.drop(columns=["Run_order", "Tool_order"])

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_summaries.to_csv(OUTPUT_FILE, index=False)

    print(f"Wrote {len(all_summaries)} rows to {OUTPUT_FILE}")
    print(all_summaries["Tool"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
