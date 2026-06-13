import os
import sqlite3

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

from load_data import DB_PATH

CELL_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
SIGNIFICANCE_THRESHOLD = 0.05
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


# ---------------------------------------------------------------------------
# Part 2: Relative frequency of each cell population per sample
# ---------------------------------------------------------------------------

def get_frequency_table() -> pd.DataFrame:
    """Return a long-form table of relative cell population frequencies.

    Columns: sample, total_count, population, count, percentage
    """
    conn = sqlite3.connect(DB_PATH)
    samples_df = pd.read_sql_query(
        f"SELECT sample, {', '.join(CELL_POPULATIONS)} FROM samples",
        conn,
    )
    conn.close()

    samples_df["total_count"] = samples_df[CELL_POPULATIONS].sum(axis=1)

    freq_df = samples_df.melt(
        id_vars=["sample", "total_count"],
        value_vars=CELL_POPULATIONS,
        var_name="population",
        value_name="count",
    )

    freq_df["percentage"] = (freq_df["count"] / freq_df["total_count"] * 100).round(2)

    return freq_df[["sample", "total_count", "population", "count", "percentage"]]


def save_frequency_table() -> None:
    """Persist the frequency table to the database as the cell_frequencies table."""
    freq_df = get_frequency_table()
    conn = sqlite3.connect(DB_PATH)
    freq_df.to_sql("cell_frequencies", conn, if_exists="replace", index=False)
    conn.close()
    print(f"Saved {len(freq_df):,} rows to cell_frequencies table.")


# ---------------------------------------------------------------------------
# Part 3: Statistical analysis — responders vs non-responders
# ---------------------------------------------------------------------------

TIME_POINTS = [0, 7, 14]


def _get_melanoma_miraclib_pbmc(time_point: int = 0) -> pd.DataFrame:
    """Return per-subject frequency table filtered to melanoma / miraclib / PBMC
    at a specific time point, enriched with response metadata.

    One row per subject per population — independent observations for statistics.
    """
    conn = sqlite3.connect(DB_PATH)
    meta_df = pd.read_sql_query(
        f"""
        SELECT s.sample, t.subject, t.response
        FROM samples s
        JOIN treatments t   ON s.treatment_id = t.id
        JOIN subjects sub   ON t.subject = sub.subject
        WHERE sub.condition             = 'melanoma'
          AND t.treatment               = 'miraclib'
          AND s.sample_type             = 'PBMC'
          AND s.time_from_treatment_start = {time_point}
          AND t.response                IN ('yes', 'no')
        """,
        conn,
    )
    conn.close()

    freq_df = get_frequency_table()
    merged = freq_df.merge(meta_df, on="sample")
    return merged


def plot_responder_boxplots(time_point: int = 0, save: bool = True) -> plt.Figure:
    """Boxplot of per-subject cell population frequencies split by responder status
    at the given time point."""
    df = _get_melanoma_miraclib_pbmc(time_point=time_point)

    fig, axes = plt.subplots(1, len(CELL_POPULATIONS), figsize=(18, 6), sharey=False)
    fig.suptitle(
        f"Cell Population Frequencies: Responders vs Non-Responders\n"
        f"(Melanoma · Miraclib · PBMC · Day {time_point})",
        fontsize=13,
        fontweight="bold",
    )

    colors = {"yes": "#4C9BE8", "no": "#E8744C"}

    for ax, pop in zip(axes, CELL_POPULATIONS):
        groups = [
            df.loc[df["population"] == pop, "percentage"][df["response"] == grp]
            for grp in ("yes", "no")
        ]
        bp = ax.boxplot(groups, patch_artist=True, widths=0.5)
        for patch, grp in zip(bp["boxes"], ("yes", "no")):
            patch.set_facecolor(colors[grp])
        ax.set_title(pop, fontsize=10)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Responder", "Non-responder"], fontsize=8)
        ax.set_ylabel("Frequency (%)" if pop == CELL_POPULATIONS[0] else "")

    plt.tight_layout()

    if save:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        path = os.path.join(OUTPUT_DIR, f"boxplot_responders_day{time_point}.png")
        fig.savefig(path, dpi=150)
        print(f"Boxplot saved to: {path}")

    return fig


def get_significance_table(time_point: int = 0) -> pd.DataFrame:
    """Mann-Whitney U test per population at a given time point.

    Unit of observation is one subject (not one sample), ensuring independence.
    """
    df = _get_melanoma_miraclib_pbmc(time_point=time_point)
    results = []

    for pop in CELL_POPULATIONS:
        pop_df = df[df["population"] == pop]
        responders     = pop_df[pop_df["response"] == "yes"]["percentage"]
        non_responders = pop_df[pop_df["response"] == "no"]["percentage"]

        stat, p_value = mannwhitneyu(responders, non_responders, alternative="two-sided")
        results.append({
            "population": pop,
            "n_responders": len(responders),
            "n_non_responders": len(non_responders),
            "mean_responders": round(responders.mean(), 2),
            "mean_non_responders": round(non_responders.mean(), 2),
            "p_value": round(p_value, 4),
            "significant": p_value < SIGNIFICANCE_THRESHOLD,
        })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Part 4: Data subset analysis — melanoma PBMC baseline miraclib samples
# ---------------------------------------------------------------------------

BASE_QUERY = """
    SELECT s.sample, sub.subject, sub.project, sub.sex, t.response
    FROM samples s
    JOIN treatments t   ON s.treatment_id = t.id
    JOIN subjects sub   ON t.subject = sub.subject
    WHERE sub.condition             = 'melanoma'
      AND s.sample_type             = 'PBMC'
      AND s.time_from_treatment_start = 0
      AND t.treatment               = 'miraclib'
"""


def get_baseline_samples() -> pd.DataFrame:
    """All melanoma PBMC baseline samples from miraclib-treated subjects."""
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(BASE_QUERY, conn)
    conn.close()
    return df


def get_samples_per_project() -> pd.DataFrame:
    """Count of baseline samples per project."""
    df = get_baseline_samples()
    return (
        df.groupby("project")
        .agg(n_samples=("sample", "count"))
        .reset_index()
    )


def get_subjects_by_response() -> pd.DataFrame:
    """Count of distinct subjects by responder status."""
    df = get_baseline_samples()
    return (
        df.groupby("response")
        .agg(n_subjects=("subject", "nunique"))
        .reset_index()
    )


def get_subjects_by_sex() -> pd.DataFrame:
    """Count of distinct subjects by sex."""
    df = get_baseline_samples()
    return (
        df.groupby("sex")
        .agg(n_subjects=("subject", "nunique"))
        .reset_index()
    )


if __name__ == "__main__":
    print("=== Part 2: Cell Population Frequencies ===")
    freq_table = get_frequency_table()
    print(freq_table.to_string(index=False))
    save_frequency_table()

    print("\n=== Part 3: Responder vs Non-Responder Significance ===")
    for tp in TIME_POINTS:
        print(f"\n-- Day {tp} --")
        print(get_significance_table(time_point=tp).to_string(index=False))

    for tp in TIME_POINTS:
        plot_responder_boxplots(time_point=tp, save=True)

    print("\n=== Part 4: Melanoma PBMC Baseline Miraclib Subset ===")
    print(f"\nTotal baseline samples: {len(get_baseline_samples())}")

    print("\n-- Samples per project --")
    print(get_samples_per_project().to_string(index=False))

    print("\n-- Subjects by response --")
    print(get_subjects_by_response().to_string(index=False))

    print("\n-- Subjects by sex --")
    print(get_subjects_by_sex().to_string(index=False))
