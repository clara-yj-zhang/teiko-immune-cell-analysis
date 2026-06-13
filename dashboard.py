import sqlite3

import plotly.express as px
import streamlit as st

from analysis import (
    CELL_POPULATIONS,
    TIME_POINTS,
    _get_melanoma_miraclib_pbmc,
    get_baseline_samples,
    get_frequency_table,
    get_samples_per_project,
    get_significance_table,
    get_subjects_by_response,
    get_subjects_by_sex,
)
from load_data import DB_PATH

st.set_page_config(page_title="Immune Cell Analysis", layout="wide")
st.title("Immune Cell Analysis Dashboard")
st.caption("Clinical trial analysis of immune cell populations — Loblaw Bio")


def db_query(sql: str):
    conn = sqlite3.connect(DB_PATH)
    try:
        import pandas as pd
        return pd.read_sql_query(sql, conn)
    finally:
        conn.close()


def get_tables_and_views():
    rows = db_query(
        "SELECT name, type FROM sqlite_master WHERE type IN ('table','view') ORDER BY type, name"
    )
    return {row["name"]: row["type"] for _, row in rows.iterrows()}


# ---------------------------------------------------------------------------
# Sidebar navigation
# ---------------------------------------------------------------------------
st.sidebar.title("Navigation")
section = st.sidebar.radio(
    "Go to",
    [
        "Part 2: Initial Analysis - Data Overview",
        "Part 3: Statistical Analysis",
        "Part 4: Data Subset Analysis",
        "Browse Tables & Views",
        "Custom SQL",
    ],
)

# ---------------------------------------------------------------------------
# Part 2: Frequency Table
# ---------------------------------------------------------------------------
if section == "Part 2: Initial Analysis - Data Overview":
    st.header("Part 2: Initial Analysis - Data Overview")

    freq_df = get_frequency_table()

    col1, col2 = st.columns(2)
    with col1:
        selected_samples = st.multiselect(
            "Filter by sample",
            options=sorted(freq_df["sample"].unique()),
            default=[],
            placeholder="All samples",
        )
    with col2:
        selected_pops = st.multiselect(
            "Filter by population",
            options=CELL_POPULATIONS,
            default=[],
            placeholder="All populations",
        )

    filtered_freq = freq_df.copy()
    if selected_samples:
        filtered_freq = filtered_freq[filtered_freq["sample"].isin(selected_samples)]
    if selected_pops:
        filtered_freq = filtered_freq[filtered_freq["population"].isin(selected_pops)]

    st.dataframe(filtered_freq.sort_values("sample"), use_container_width=True, hide_index=True)

# ---------------------------------------------------------------------------
# Part 3: Responder vs Non-Responder Analysis
# ---------------------------------------------------------------------------
elif section == "Part 3: Statistical Analysis":
    st.header("Part 3: Statistical Analysis")
    st.caption("Melanoma · Miraclib · PBMC · Per subject (one observation per subject per time point)")

    time_point = st.radio(
        "Time point (days from treatment start)",
        options=TIME_POINTS,
        index=2,
        horizontal=True,
        format_func=lambda d: f"Day {d}",
    )

    sig_df = get_significance_table(time_point=time_point)

    st.subheader("Statistical Significance (Mann-Whitney U Test)")
    st.dataframe(
        sig_df.style.apply(
            lambda col: ["background-color: #d4edda" if v else "" for v in col]
            if col.name == "significant"
            else [""] * len(col),
            axis=0,
        ),
        use_container_width=True,
        hide_index=True,
    )

    st.subheader("Boxplots by Cell Population")
    box_col1, box_col2, box_col3 = st.columns(3)
    with box_col1:
        box_view = st.radio("View", ["All populations", "Single population"], horizontal=True, key="box_view")
    with box_col2:
        selected_pop = st.selectbox("Population", CELL_POPULATIONS, disabled=(box_view == "All populations"))
    with box_col3:
        box_time = st.selectbox("Time point (days)", TIME_POINTS, key="box_time")

    all_pop_data = _get_melanoma_miraclib_pbmc(time_point=box_time)

    if box_view == "All populations":
        fig_box = px.box(
            all_pop_data,
            x="response",
            y="percentage",
            color="response",
            facet_col="population",
            color_discrete_map={"yes": "#4C9BE8", "no": "#E8744C"},
            labels={"response": "Response", "percentage": "Frequency (%)"},
            title=f"All Cell Populations — Responders vs Non-Responders (Day {box_time})",
            category_orders={"response": ["yes", "no"], "population": CELL_POPULATIONS},
        )
        fig_box.update_layout(showlegend=False)
        fig_box.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    else:
        pop_data = all_pop_data[all_pop_data["population"] == selected_pop]
        fig_box = px.box(
            pop_data,
            x="response",
            y="percentage",
            color="response",
            color_discrete_map={"yes": "#4C9BE8", "no": "#E8744C"},
            labels={"response": "Response", "percentage": "Frequency (%)"},
            title=f"{selected_pop} — Responders vs Non-Responders (Day {box_time})",
            category_orders={"response": ["yes", "no"]},
        )
        fig_box.update_layout(showlegend=False)

    st.plotly_chart(fig_box, use_container_width=True)

    # Findings summary
    st.subheader("Findings Summary")
    st.caption(f"Based on Mann-Whitney U test at Day {time_point} · α = 0.05")

    sig_pops = sig_df[sig_df["significant"]]
    non_sig_pops = sig_df[~sig_df["significant"]]

    if sig_pops.empty:
        st.info(
            f"At Day {time_point}, **no cell populations** show a statistically significant "
            f"difference in relative frequencies between responders and non-responders "
            f"(all p-values > 0.05). There is insufficient statistical evidence at this "
            f"time point to identify immune cell predictors of miraclib response in melanoma patients."
        )
    else:
        for _, row in sig_pops.iterrows():
            direction = "higher" if row["mean_responders"] > row["mean_non_responders"] else "lower"
            st.success(
                f"**{row['population']}** — significantly {direction} in responders "
                f"({row['mean_responders']}% vs {row['mean_non_responders']}%, "
                f"p = {row['p_value']}, n = {row['n_responders']} responders / "
                f"{row['n_non_responders']} non-responders)"
            )
        if not non_sig_pops.empty:
            non_sig_list = ", ".join(
                f"{r['population']} (p={r['p_value']})"
                for _, r in non_sig_pops.iterrows()
            )
            st.info(f"No significant difference at this time point: {non_sig_list}")


elif section == "Part 4: Data Subset Analysis":
    st.header("Part 4: Data Subset Analysis")

    # --- Interactive filter table ---
    st.subheader("Sample Explorer")
    all_conditions   = sorted(db_query("SELECT DISTINCT condition FROM subjects")["condition"].dropna())
    all_sample_types = sorted(db_query("SELECT DISTINCT sample_type FROM samples")["sample_type"].dropna())
    all_treatments   = sorted(db_query("SELECT DISTINCT treatment FROM treatments")["treatment"].dropna())
    all_times        = sorted(db_query("SELECT DISTINCT time_from_treatment_start FROM samples")["time_from_treatment_start"].dropna().astype(int))

    fc1, fc2, fc3, fc4 = st.columns(4)
    with fc1:
        f_condition = st.multiselect("Condition", all_conditions, default=["melanoma"])
    with fc2:
        f_sample_type = st.multiselect("Sample type", all_sample_types, default=["PBMC"])
    with fc3:
        f_treatment = st.multiselect("Treatment", all_treatments, default=["miraclib"])
    with fc4:
        f_time = st.multiselect("Time (days)", all_times, default=[0])

    clauses = []
    if f_condition:
        vals = ", ".join(f"'{v}'" for v in f_condition)
        clauses.append(f"condition IN ({vals})")
    if f_sample_type:
        vals = ", ".join(f"'{v}'" for v in f_sample_type)
        clauses.append(f"sample_type IN ({vals})")
    if f_treatment:
        vals = ", ".join(f"'{v}'" for v in f_treatment)
        clauses.append(f"treatment IN ({vals})")
    if f_time:
        vals = ", ".join(str(v) for v in f_time)
        clauses.append(f"time_from_treatment_start IN ({vals})")

    where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
    display_cols = "sample, subject, project, condition, sample_type, time_from_treatment_start, treatment, response, sex"
    filter_df = db_query(f"SELECT {display_cols} FROM sample_details {where} ORDER BY sample")
    st.caption(f"{len(filter_df)} samples match the current filters")
    st.dataframe(filter_df, use_container_width=True, hide_index=True)

    st.divider()
    st.subheader("Subset Summary")

    filter_caption = " · ".join(
        ([", ".join(f_condition)] if f_condition else [])
        + ([", ".join(f_sample_type)] if f_sample_type else [])
        + ([", ".join(f_treatment)] if f_treatment else [])
        + ([f"Time {', '.join(str(t) for t in f_time)}"] if f_time else [])
    ) or "All samples"
    st.caption(filter_caption)

    st.metric("Total samples", len(filter_df))

    col_a, col_b, col_c = st.columns(3)

    with col_a:
        st.subheader("Samples per Project")
        proj_df = (
            filter_df.groupby("project")
            .agg(n_samples=("sample", "count"))
            .reset_index()
        )
        st.dataframe(proj_df, use_container_width=True, hide_index=True)
        fig_proj = px.bar(proj_df, x="project", y="n_samples", text="n_samples",
                          color="project", title="Samples per Project")
        fig_proj.update_traces(textposition="outside")
        fig_proj.update_layout(showlegend=False)
        st.plotly_chart(fig_proj, use_container_width=True)

    with col_b:
        st.subheader("Subjects by Response")
        resp_df = (
            filter_df.dropna(subset=["response"])
            .groupby("response")
            .agg(n_subjects=("subject", "nunique"))
            .reset_index()
        )
        st.dataframe(resp_df, use_container_width=True, hide_index=True)
        fig_resp = px.pie(resp_df, names="response", values="n_subjects",
                          color="response",
                          color_discrete_map={"yes": "#4C9BE8", "no": "#E8744C"},
                          title="Subjects by Response")
        st.plotly_chart(fig_resp, use_container_width=True)

    with col_c:
        st.subheader("Subjects by Sex")
        sex_df = (
            filter_df.dropna(subset=["sex"])
            .groupby("sex")
            .agg(n_subjects=("subject", "nunique"))
            .reset_index()
        )
        st.dataframe(sex_df, use_container_width=True, hide_index=True)
        fig_sex = px.pie(sex_df, names="sex", values="n_subjects",
                         color="sex",
                         color_discrete_map={"M": "#6C8EBF", "F": "#D6A0C0"},
                         title="Subjects by Sex")
        st.plotly_chart(fig_sex, use_container_width=True)

# ---------------------------------------------------------------------------
# Explorer: Browse Tables & Views
# ---------------------------------------------------------------------------
elif section == "Browse Tables & Views":
    st.header("Browse Tables & Views")

    objects = get_tables_and_views()
    selected = st.selectbox(
        "Select a table or view",
        list(objects.keys()),
        format_func=lambda n: f"{n}  [{objects[n]}]",
    )

    count_df = db_query(f'SELECT COUNT(*) AS row_count FROM "{selected}"')
    total_rows = count_df["row_count"].iloc[0]

    col1, col2 = st.columns([3, 1])
    with col1:
        st.subheader(selected)
    with col2:
        st.metric("Total rows", f"{total_rows:,}")

    with st.expander("Column info"):
        col_info = db_query(f'PRAGMA table_info("{selected}")')
        st.dataframe(col_info[["name", "type", "notnull", "pk"]], hide_index=True)

    st.markdown("**Filters**")
    col_info_full = db_query(f'PRAGMA table_info("{selected}")')
    columns = col_info_full["name"].tolist()

    filter_col, filter_val = st.columns(2)
    with filter_col:
        filter_column = st.selectbox("Filter column", ["(none)"] + columns)
    with filter_val:
        filter_value = st.text_input("Filter value (exact match)", "")

    limit = st.slider("Rows to display", min_value=10, max_value=1000, value=100, step=10)

    where_clause = ""
    if filter_column != "(none)" and filter_value:
        where_clause = f'WHERE "{filter_column}" = \'{filter_value}\''

    df = db_query(f'SELECT * FROM "{selected}" {where_clause} LIMIT {limit}')
    st.dataframe(df, use_container_width=True, hide_index=True)
    st.caption(f"Showing {len(df):,} of {total_rows:,} rows")

# ---------------------------------------------------------------------------
# Explorer: Custom SQL
# ---------------------------------------------------------------------------
elif section == "Custom SQL":
    st.header("Custom SQL Query")

    with st.expander("Tables & Views reference"):
        objects = get_tables_and_views()
        for name, obj_type in objects.items():
            col_info = db_query(f'PRAGMA table_info("{name}")')
            cols = ", ".join(col_info["name"].tolist())
            st.markdown(f"**{name}** `[{obj_type}]`  \n`{cols}`")

    # Interactive filter builder
    with st.expander("Interactive Filter Builder", expanded=False):
        st.markdown("Build a `WHERE` clause from `sample_details` without writing SQL.")
        fc1, fc2, fc3, fc4 = st.columns(4)
        with fc1:
            f_condition = st.multiselect("condition", ["melanoma", "carcinoma", "healthy"])
        with fc2:
            f_treatment = st.multiselect("treatment", ["miraclib", "phauximab", "none"])
        with fc3:
            f_response = st.multiselect("response", ["yes", "no"])
        with fc4:
            f_sex = st.multiselect("sex", ["M", "F"])
        fc5, fc6 = st.columns(2)
        with fc5:
            f_sample_type = st.multiselect("sample_type", ["PBMC", "WB"])
        with fc6:
            f_time = st.multiselect("time_from_treatment_start", [0, 7, 14])
        f_limit = st.number_input("LIMIT", min_value=1, max_value=10000, value=50)

        clauses = []
        if f_condition:    clauses.append(f"condition IN ({', '.join(repr(v) for v in f_condition)})")
        if f_treatment:    clauses.append(f"treatment IN ({', '.join(repr(v) for v in f_treatment)})")
        if f_response:     clauses.append(f"response IN ({', '.join(repr(v) for v in f_response)})")
        if f_sex:          clauses.append(f"sex IN ({', '.join(repr(v) for v in f_sex)})")
        if f_sample_type:  clauses.append(f"sample_type IN ({', '.join(repr(v) for v in f_sample_type)})")
        if f_time:         clauses.append(f"time_from_treatment_start IN ({', '.join(str(v) for v in f_time)})")

        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        generated_sql = f"SELECT *\nFROM sample_details\n{where}\nLIMIT {f_limit}"
        st.code(generated_sql, language="sql")
        if st.button("Use this query", key="use_generated"):
            st.session_state["prefilled_sql"] = generated_sql

    default_sql = st.session_state.get("prefilled_sql", """\
SELECT ROUND(AVG(s.b_cell), 2) AS avg_b_cell
FROM samples s
JOIN treatments t   ON s.treatment_id = t.id
JOIN subjects sub   ON t.subject = sub.subject
WHERE sub.condition                 = 'melanoma'
  AND sub.sex                       = 'M'
  AND t.response                    = 'yes'
  AND s.time_from_treatment_start   = 0""")

    sql = st.text_area("SQL", value=default_sql, height=180)

    if st.button("Run", type="primary"):
        try:
            result = db_query(sql)
            st.success(f"{len(result):,} rows returned")
            st.dataframe(result, use_container_width=True, hide_index=True)
            csv = result.to_csv(index=False).encode()
            st.download_button("Download CSV", csv, "query_result.csv", "text/csv")
        except Exception as e:
            st.error(f"Query error: {e}")


