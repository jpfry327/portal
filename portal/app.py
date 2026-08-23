"""Chordoma Multi-omic Atlas — companion portal prototype.

Run:  streamlit run portal/app.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

import portal_data as D

st.set_page_config(page_title="Chordoma Multi-omic Atlas", page_icon="🧬", layout="wide")

SUBTYPE_ORDER = ["C1", "C2", "C3", "C4"]


def subtype_legend():
    cols = st.columns(4)
    blurbs = {
        "C1": "Hypermetabolic — loss of notochordal differentiation, nucleotide biosynthesis. Poor prognosis.",
        "C2": "Neural proliferative — MYCN-driven, highly proliferative, neural/GABAergic programs. Poor prognosis.",
        "C3": "Immunogenic — immune/inflammatory infiltration. Intermediate outcomes.",
        "C4": "LYST-activated — hypomethylated enhancer landscape, LYST super-enhancer. Best prognosis.",
    }
    for c, (s, blurb) in zip(cols, blurbs.items()):
        with c:
            st.markdown(
                f"<div style='border-left: 5px solid {D.SUBTYPE_COLORS[s]}; padding-left: 8px'>"
                f"<b>{D.SUBTYPE_NAMES[s]}</b><br><small>{blurb}</small></div>",
                unsafe_allow_html=True,
            )


def box_strip(df: pd.DataFrame, value_col: str, title: str, ylab: str) -> go.Figure:
    """Manuscript-style violin: viridis fill, black outline, inner grey box."""
    df = df[df.subtype.isin(SUBTYPE_ORDER)]
    fig = go.Figure()
    for s in SUBTYPE_ORDER:
        g = df[df.subtype == s]
        if g.empty:
            continue
        fig.add_trace(go.Violin(
            x=g["subtype"], y=g[value_col], name=s,
            fillcolor=D.SUBTYPE_COLORS[s], opacity=0.95,
            line=dict(color="black", width=1.5),
            box=dict(visible=True, width=0.22, fillcolor="#DCDCDC",
                     line=dict(color="black", width=1.2)),
            meanline_visible=False, points=False, spanmode="hard",
            hoverinfo="y+name",
        ))
    p_txt = kruskal_p_text(df, value_col)
    if p_txt:
        fig.add_annotation(x=0.03, y=0.98, xref="paper", yref="paper",
                           text=p_txt, showarrow=False, font=dict(size=15))
    fig.update_layout(
        showlegend=False, yaxis_title=ylab, xaxis_title="", height=420,
        title=title, margin=dict(t=48, b=24), template="simple_white",
        xaxis=dict(categoryorder="array", categoryarray=SUBTYPE_ORDER),
    )
    return fig


def kruskal_p_text(df: pd.DataFrame, value_col: str) -> str:
    groups = [g[value_col].dropna().values for _, g in df.groupby("subtype") if len(g) >= 3]
    if len(groups) < 2:
        return ""
    _, p = stats.kruskal(*groups)
    return "p < 0.0001" if p < 1e-4 else f"p = {p:.3g}"


def kruskal_note(df: pd.DataFrame, value_col: str) -> str:
    groups = [g[value_col].dropna().values for _, g in df.groupby("subtype") if len(g) >= 3]
    if len(groups) < 2:
        return ""
    h, p = stats.kruskal(*groups)
    return f"Kruskal–Wallis across subtypes: H = {h:.1f}, p = {p:.2e}"


# ---------------------------------------------------------------- pages
def page_overview():
    st.title("Chordoma Multi-omic Atlas")
    st.markdown(
        "Companion portal to *An integrative multi-omic classification of chordoma*: "
        "four molecular subgroups (C1–C4) defined by integrated DNA methylation + RNA "
        "expression in 135 discovery chordomas and reproduced in three validation cohorts. "
        "Query gene expression, proteomics, and methylation by subtype and cohort."
    )
    subtype_legend()
    st.divider()

    s = D.samples()
    c1, c2 = st.columns([1, 1])
    with c1:
        st.subheader("Samples by cohort and subtype")
        counts = (
            s.groupby(["cohort", "subtype"]).size().reset_index(name="n")
        )
        counts["cohort"] = counts["cohort"].map(D.COHORT_LABELS)
        fig = px.bar(
            counts, x="cohort", y="n", color="subtype", barmode="stack",
            category_orders={"subtype": SUBTYPE_ORDER},
            color_discrete_map=D.SUBTYPE_COLORS,
        )
        fig.update_layout(height=400, xaxis_title="", yaxis_title="samples",
                          legend_title="subtype")
        st.plotly_chart(fig, use_container_width=True)
        st.caption(
            "Discovery subtypes are assigned by multi-omic clustering; validation "
            "cohorts are classifier predictions from RNA (and methylation where available)."
        )
    with c2:
        st.subheader("Tumor microenvironment (discovery)")
        props = D.cell_proportions()
        cell_cols = [c for c in props.columns if c not in ("sample_id", "subtype")]
        long = props.melt(id_vars=["sample_id", "subtype"], value_vars=cell_cols,
                          var_name="cell_type", value_name="fraction")
        fig = px.box(
            long, x="cell_type", y="fraction", color="subtype",
            category_orders={"subtype": SUBTYPE_ORDER},
            color_discrete_map=D.SUBTYPE_COLORS,
        )
        fig.update_layout(height=400, xaxis_title="", yaxis_title="BayesPrism fraction")
        st.plotly_chart(fig, use_container_width=True)
        st.caption("Bulk deconvolution (BayesPrism) cell-type fractions per subtype.")

    st.subheader("Sample table")
    tbl = s.copy()
    tbl["cohort"] = tbl["cohort"].map(D.COHORT_LABELS)
    st.dataframe(tbl, use_container_width=True, height=300)
    st.download_button("Download sample table (CSV)", tbl.to_csv(index=False),
                       "chordoma_atlas_samples.csv", "text/csv")


def page_expression():
    st.title("Gene expression by subtype")
    genes = D.gene_list()
    default = genes.index("LYST") if "LYST" in genes else 0
    gene = st.selectbox("Gene", genes, index=default,
                        help="Symbol as used in each cohort's expression matrix")
    cohorts = st.multiselect(
        "Cohorts", D.RNA_COHORTS, default=D.RNA_COHORTS,
        format_func=lambda c: D.COHORT_LABELS[c],
    )
    if not cohorts:
        st.info("Pick at least one cohort.")
        return

    tabs = st.tabs([D.COHORT_LABELS[c] for c in cohorts])
    rows = []
    for tab, cohort in zip(tabs, cohorts):
        with tab:
            df = D.expression(gene, cohort)
            if df.empty:
                st.warning(f"{gene} not present in the {cohort} expression matrix.")
                continue
            unit = "log2 RPKM" if cohort == "discovery" else "normalized log expression"
            fig = box_strip(df, "value", f"{gene} — {D.COHORT_LABELS[cohort]}", unit)
            st.plotly_chart(fig, use_container_width=True)
            note = kruskal_note(df, "value")
            if note:
                st.caption(note)
            rows.append(df.assign(gene=gene))

    if rows:
        allc = pd.concat(rows, ignore_index=True)
        st.download_button(
            f"Download {gene} values (all selected cohorts, CSV)",
            allc.to_csv(index=False), f"{gene}_expression.csv", "text/csv",
        )

    prot = D.protein(gene)
    if not prot.empty:
        st.divider()
        st.subheader(f"{gene} protein abundance (Zhang proteogenomics)")
        c1, c2 = st.columns([1, 1])
        with c1:
            fig = box_strip(prot, "value", f"{gene} protein by subtype", "corrected log abundance")
            st.plotly_chart(fig, use_container_width=True)
            note = kruskal_note(prot, "value")
            if note:
                st.caption(note)
        with c2:
            pairs = D.rna_protein_pairs(gene)
            if len(pairs) >= 5:
                r, p = stats.spearmanr(pairs["rna"], pairs["prot"])
                fig = px.scatter(
                    pairs, x="rna", y="prot", color="subtype",
                    color_discrete_map=D.SUBTYPE_COLORS, hover_data=["sample_id"],
                    title=f"RNA vs protein — Spearman ρ = {r:.2f}, p = {p:.1e} (n={len(pairs)})",
                    trendline="ols",
                )
                fig.update_layout(height=420, xaxis_title=f"{gene} RNA",
                                  yaxis_title=f"{gene} protein")
                st.plotly_chart(fig, use_container_width=True)


def page_methylation():
    st.title("DNA methylation by subtype")
    st.markdown(
        "Query normalized, batch-corrected **M-values** for any of the 756,448 array "
        "probes in the discovery cohort (n=161 arrays incl. replicates)."
    )
    probe = st.text_input("Probe ID (Illumina cg identifier)", value="cg26928153")
    if not probe:
        return
    df = D.methylation(probe.strip())
    if df.empty:
        st.warning(f"Probe '{probe}' not found. Probes look like cg26928153.")
        return
    c1, c2 = st.columns([2, 1])
    with c1:
        fig = box_strip(df, "value", f"{probe} methylation by subtype", "M-value")
        st.plotly_chart(fig, use_container_width=True)
        note = kruskal_note(df, "value")
        if note:
            st.caption(note)
    with c2:
        st.metric("Samples with data", len(df))
        st.metric("Median M-value", f"{df.value.median():.2f}")
        beta = 2 ** df.value / (2 ** df.value + 1)
        st.metric("Median beta (approx.)", f"{beta.median():.2f}")
    st.download_button(f"Download {probe} values (CSV)", df.to_csv(index=False),
                       f"{probe}_methylation.csv", "text/csv")


def page_correlation():
    st.title("Gene–gene correlation")
    genes = D.gene_list()
    c1, c2, c3 = st.columns(3)
    with c1:
        gene_a = st.selectbox("Gene A", genes, index=genes.index("LYST") if "LYST" in genes else 0)
    with c2:
        gene_b = st.selectbox("Gene B", genes, index=genes.index("TBXT") if "TBXT" in genes else 1)
    with c3:
        cohort = st.selectbox("Cohort", D.RNA_COHORTS, format_func=lambda c: D.COHORT_LABELS[c])
    method = st.radio("Correlation", ["Spearman", "Pearson"], horizontal=True)

    pairs = D.pair_expression(gene_a, gene_b, cohort)
    if len(pairs) < 5:
        st.warning("Fewer than 5 samples with both genes measured in this cohort.")
        return
    fn = stats.spearmanr if method == "Spearman" else stats.pearsonr
    r, p = fn(pairs["x"], pairs["y"])
    sym = "ρ" if method == "Spearman" else "r"
    fig = px.scatter(
        pairs, x="x", y="y", color="subtype", color_discrete_map=D.SUBTYPE_COLORS,
        category_orders={"subtype": SUBTYPE_ORDER}, hover_data=["sample_id"],
        title=f"{gene_a} vs {gene_b} — {method} {sym} = {r:.2f}, p = {p:.1e} (n={len(pairs)})",
        trendline="ols",
    )
    fig.update_layout(height=520, xaxis_title=gene_a, yaxis_title=gene_b)
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Within-subtype correlations")
    rows = []
    for s, g in pairs.groupby("subtype"):
        if len(g) >= 5:
            r_s, p_s = fn(g["x"], g["y"])
            rows.append({"subtype": s, "n": len(g), sym: round(r_s, 3), "p": f"{p_s:.2e}"})
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True)
    st.download_button("Download paired values (CSV)", pairs.to_csv(index=False),
                       f"{gene_a}_{gene_b}_{cohort}.csv", "text/csv")


SURV_ENDPOINTS = {  # first entry = default (RFS/PFS have more events)
    "discovery": {
        "Recurrence-free survival": ("time_rfs", "status_rfs"),
        "Overall survival": ("time_os", "status_os"),
    },
    "zhang": {
        "Progression-free survival": ("time_pfs", "status_pfs"),
        "Overall survival": ("time_os", "status_os"),
    },
}


def km_figure(df, tcol, scol, groups, color_map, title):
    from lifelines import KaplanMeierFitter
    kmf = KaplanMeierFitter()
    fig = go.Figure()
    for grp in groups:
        g = df[df.group == grp]
        if len(g) < 3:
            continue
        kmf.fit(g[tcol] / 365.25, g[scol], label=grp)
        sf = kmf.survival_function_.reset_index()
        fig.add_trace(go.Scatter(
            x=sf["timeline"], y=sf[grp], mode="lines", name=f"{grp} (n={len(g)})",
            line=dict(shape="hv", width=2.5, color=color_map.get(grp)),
            legendgroup=grp,
        ))
        # ggsurvplot-style "+" marks at censoring times
        cens = g.loc[g[scol] == 0, tcol] / 365.25
        if len(cens):
            cs = kmf.survival_function_at_times(cens).values
            fig.add_trace(go.Scatter(
                x=cens, y=cs, mode="markers", showlegend=False, legendgroup=grp,
                marker=dict(symbol="cross-thin-open", size=7,
                            line=dict(width=1.5, color=color_map.get(grp)),
                            color=color_map.get(grp)),
                hoverinfo="skip",
            ))
    fig.update_layout(
        height=500, xaxis_title="Years", yaxis_title="Survival probability",
        yaxis_range=[0, 1.02], title=title, template="simple_white",
        legend=dict(x=0.99, xanchor="right", y=0.99, yanchor="top"),
    )
    return fig


def page_survival():
    from lifelines.statistics import multivariate_logrank_test
    import survival_stats as SS

    st.title("Survival analysis")
    c1, c2 = st.columns(2)
    with c1:
        cohort = st.selectbox("Cohort", ["discovery", "zhang"],
                              format_func=lambda c: D.COHORT_LABELS[c])
    with c2:
        endpoint = st.selectbox("Endpoint", list(SURV_ENDPOINTS[cohort]))
    tcol, scol = SURV_ENDPOINTS[cohort][endpoint]

    clin = D.clinical(cohort).dropna(subset=[tcol, scol])
    note_extra = ""
    if cohort == "zhang" and tcol == "time_os":
        # Recurrent samples inherit the matched primary's time-to-death (time
        # zero = diagnosis, not this resection) -> OS restricted to primaries.
        clin = clin[clin["tumor_sample_type"] == 1]
        note_extra = (" OS restricted to primary tumors: recurrent samples share the "
                      "matched primary's death time, so time-from-resection is only "
                      "defined for primaries (see paper methods).")
    st.caption(f"{len(clin)} samples with follow-up, {int(clin[scol].sum())} events.{note_extra}")

    mode = st.radio("Stratify by", ["Molecular subtype", "Gene expression"], horizontal=True)

    if mode == "Molecular subtype":
        df = clin.copy()
        df["group"] = df["subtype"]
        order = [s for s in SUBTYPE_ORDER if s in set(df.group)]
        res = multivariate_logrank_test(df[tcol], df["group"], df[scol])
        fig = km_figure(df, tcol, scol, order, D.SUBTYPE_COLORS,
                        f"{endpoint} by subtype — log-rank p = {res.p_value:.2e}")
        st.plotly_chart(fig, use_container_width=True)
        return

    # ---- per-gene: DISPLAY optimal-cutpoint KM, REPORT Cox LR p + HR ----
    genes = D.gene_list()
    gcol, mcol = st.columns([2, 1])
    with gcol:
        gene = st.selectbox("Gene", genes, index=genes.index("LYST") if "LYST" in genes else 0)
    with mcol:
        molecule = st.radio("Molecule", ["RNA", "Protein"], horizontal=True,
                            disabled=cohort != "zhang",
                            help="Protein abundance is available for the Zhang cohort only")
    if cohort != "zhang":
        molecule = "RNA"
    if molecule == "Protein":
        expr = D.protein(gene)[["sample_id", "value"]]
    else:
        expr = D.expression(gene, cohort)[["sample_id", "value"]]
    df = clin.merge(expr, on="sample_id", how="inner").dropna(subset=["value"])
    if len(df) < 20:
        st.warning(f"{gene}: fewer than 20 samples with expression + follow-up in this cohort.")
        return

    cox = SS.cox_continuous(df, tcol, scol)
    cut = SS.optimal_cutpoint(df, tcol, scol)
    if not cox or not cut:
        st.warning("Not enough events for a stable model for this gene/endpoint.")
        return

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Cox HR (per SD)", f"{cox['hr']:.2f}",
              help="Hazard ratio per 1 SD increase in expression, continuous Cox PH model")
    m2.metric("95% CI", f"{cox['hr_lo']:.2f}–{cox['hr_hi']:.2f}")
    m3.metric("Likelihood-ratio p", f"{cox['lr_p']:.2e}")
    m4.metric("n (events)", f"{cox['n']} ({cox['events']})")

    lab = f"{gene} {molecule}"
    df["group"] = np.where(df["value"] > cut["cutpoint"], f"{lab} high", f"{lab} low")
    color_map = {f"{lab} high": "#B2182B", f"{lab} low": "#2166AC"}
    fig = km_figure(
        df, tcol, scol, [f"{lab} high", f"{lab} low"], color_map,
        f"{lab} — {endpoint}, optimal cutpoint {cut['cutpoint']:.2f} "
        f"({cut['n_high']} high / {cut['n_low']} low)",
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "KM at the optimal cutpoint (maximally selected log-rank statistic, scanned "
        "between the 10th–90th expression percentiles, ≥10 per arm) is for display; "
        f"the reported statistics are from Cox PH regression on continuous expression "
        f"(LR p = {cox['lr_p']:.2e}, HR = {cox['hr']:.2f}). The cutpoint log-rank p "
        f"({cut['logrank_p']:.2e}) is optimistically biased by cutpoint selection and "
        "should not be quoted."
    )
    st.download_button("Download per-sample data (CSV)", df.to_csv(index=False),
                       f"{gene}_{cohort}_{tcol}_survival.csv", "text/csv")


def page_celllines():
    st.title("Chordoma cell lines & PDXs")
    st.markdown(
        "RNA-seq TPM across chordoma cell lines and patient-derived xenografts "
        "(CCL Tracker resource). Single gene: median bar ± SD with replicate points. "
        "Multiple genes: per-sample distribution with per-gene dots."
    )
    genes_avail = D.cellline_genes()
    sel = st.multiselect("Gene(s)", genes_avail,
                         default=["TBXT"] if "TBXT" in genes_avail else genes_avail[:1])
    include_pdx = st.toggle("Include PDXs", value=False)
    if not sel:
        return

    frames = [D.cellline_tpm(g).assign(gene=g) for g in sel]
    df = pd.concat(frames, ignore_index=True)
    if not include_pdx:
        df = df[df.sample_type == "Cell Line"]
    if df.empty:
        st.warning("No data for the selected gene(s).")
        return
    df["log_tpm"] = np.log10(df.tpm + 1)

    if len(sel) == 1:
        gene = sel[0]
        stats_df = (df.groupby("sample_id")
                    .agg(med=("tpm", "median"), sd=("tpm", "std"), n=("tpm", "size"))
                    .fillna({"sd": 0}).reset_index()
                    .sort_values("med"))
        stats_df["y"] = np.log10(stats_df.med + 1)
        stats_df["lo"] = np.log10(np.maximum(stats_df.med - stats_df.sd, 0) + 1)
        stats_df["hi"] = np.log10(stats_df.med + stats_df.sd + 1)
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=stats_df.sample_id, y=stats_df.y,
            error_y=dict(type="data", symmetric=False,
                         array=stats_df.hi - stats_df.y, arrayminus=stats_df.y - stats_df.lo,
                         color="#333", thickness=1.2, width=4),
            marker=dict(color="#31688E", opacity=0.75), showlegend=False,
        ))
        fig.add_trace(go.Scatter(
            x=df.sample_id, y=df.log_tpm, mode="markers",
            marker=dict(color="#1a1a2e", size=6, opacity=0.65),
            text=df.replicate, showlegend=False,
        ))
        fig.update_layout(
            template="simple_white", height=520,
            title=f"{gene} — chordoma cell lines" + (" + PDXs" if include_pdx else ""),
            xaxis=dict(categoryorder="array", categoryarray=stats_df.sample_id.tolist(),
                       tickangle=-45),
            yaxis_title="log10(TPM + 1)", margin=dict(b=110),
        )
    else:
        pts = (df.groupby(["gene", "sample_id"])["log_tpm"].mean().reset_index())
        order = (pts.groupby("sample_id")["log_tpm"].min().sort_values().index.tolist())
        fig = go.Figure()
        fig.add_trace(go.Box(
            x=pts.sample_id, y=pts.log_tpm, boxpoints=False,
            fillcolor="rgba(0,0,0,0.06)", line=dict(color="#666", width=1.2),
            showlegend=False,
        ))
        for g in sel:
            sub = pts[pts.gene == g]
            if sub.empty:
                continue
            fig.add_trace(go.Scatter(x=sub.sample_id, y=sub.log_tpm, mode="markers",
                                     name=g, marker=dict(size=8, opacity=0.85)))
        fig.update_layout(
            template="simple_white", height=520,
            title=f"{len(sel)} genes — chordoma cell lines" + (" + PDXs" if include_pdx else ""),
            xaxis=dict(categoryorder="array", categoryarray=order, tickangle=-45),
            yaxis_title="mean log10(TPM + 1)", margin=dict(b=110),
        )
    st.plotly_chart(fig, use_container_width=True)
    st.download_button("Download values (CSV)", df.to_csv(index=False),
                       f"celllines_{'_'.join(sel)[:60]}.csv", "text/csv")


PAGES = {
    "Overview": page_overview,
    "Gene expression": page_expression,
    "Methylation": page_methylation,
    "Correlation": page_correlation,
    "Survival": page_survival,
    "Cell lines": page_celllines,
}

with st.sidebar:
    st.markdown("## 🧬 Chordoma Atlas")
    choice = st.radio("Navigate", list(PAGES), label_visibility="collapsed")
    st.divider()
    st.markdown(
        "<small>Companion resource to the chordoma multi-omic subtyping manuscript. "
        "Data: 135 discovery + 3 validation cohorts (RNA, protein, methylation). "
        "Please cite the paper when using this resource.</small>",
        unsafe_allow_html=True,
    )

PAGES[choice]()
