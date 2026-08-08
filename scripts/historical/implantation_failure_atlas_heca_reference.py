from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path("[local path omitted]")
HECA_PATH = Path("/Volumes/Extreme SSD/reference_atlas/HECA/h5ad/endometriumAtlasV2_cells_with_counts.h5ad")
OUTDIR = ROOT / "analysis" / "05_implantation_failure_atlas" / "heca_reference"
TABDIR = OUTDIR / "tables"


def map_broad_celltype(celltype: str, lineage: str) -> str:
    celltype = str(celltype)
    lineage = str(lineage)
    if lineage == "Endothelial":
        return "Endothelial"
    if lineage == "Immune":
        if "Myeloid" in celltype:
            return "Myeloid"
        return "Lymphoid"
    if lineage == "Mesenchymal":
        if celltype.startswith("dStromal"):
            return "Decidual_Stroma"
        return "Stroma"
    if lineage == "Epithelial":
        if celltype in {"Luminal", "SOX9_functionalis_I", "SOX9_functionalis_II"}:
            return "Luminal_Epi"
        if celltype in {"Glandular", "Glandular_secretory", "preGlandular", "Ciliated"}:
            return "Glandular_Epi"
        return "Epithelial_Other"
    return "Other"


def simplify_stage(stage: str) -> str:
    stage = str(stage)
    if "Secretory Early-Mid" in stage:
        return "Secretory_EarlyMid"
    if "Secretory Early" in stage:
        return "Secretory_Early"
    if "Secretory Mid" in stage:
        return "Secretory_Mid"
    if "Secretory Late" in stage:
        return "Secretory_Late"
    if stage == "Secretory":
        return "Secretory_Mid"
    if "Proliferative" in stage:
        return "Proliferative"
    if "Menstrual" in stage:
        return "Menstrual"
    if "Hormones" in stage:
        return "Hormones"
    return "Other"


def normalize_selected(x, n_counts):
    scaled = x / np.maximum(n_counts[:, None], 1.0) * 1e4
    return np.log1p(scaled)


def main() -> None:
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg_cache")
    TABDIR.mkdir(parents=True, exist_ok=True)

    receptivity_sig = pd.read_csv(ROOT / "analysis/03_rif/round5_maximal/tables/receptivity_signature_genes_round5.csv")
    rif_sig = pd.read_csv(ROOT / "analysis/03_rif/round5_maximal/tables/rif_directional_signature_genes_round5.csv")
    lesion_deg = pd.read_csv(ROOT / "analysis/02_endometriosis/round4_deep/tables/Figure_ENDO_2A_Ectopic_vs_Eutopic_deg.csv")
    adeno_cons = pd.read_csv(ROOT / "analysis/04_adenomyosis/round4_deep/tables/consensus_adenomyosis_signature.csv")

    marker_sets = {
        "Luminal_Epi": ["EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "TACSTD2"],
        "Glandular_Epi": ["EPCAM", "KRT8", "KRT18", "PAEP", "CXCL14", "SCGB2A1"],
        "Stroma": ["DCN", "COL1A1", "COL3A1", "LUM", "COL6A1", "CFD"],
        "Decidual_Stroma": ["IGFBP1", "PRL", "LEFTY2", "FOXO1", "HAND2", "PAEP"],
        "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ESAM"],
        "Lymphoid": ["NKG7", "KLRD1", "CD3D", "TRBC1", "IL7R", "CCL5"],
        "Myeloid": ["LST1", "C1QC", "TYROBP", "AIF1", "FCER1G", "CTSB"],
    }
    selected_genes = sorted(
        {
            *sum(marker_sets.values(), []),
            *receptivity_sig["late_genes"].dropna().tolist(),
            *receptivity_sig["early_genes"].dropna().tolist(),
            *rif_sig["rif_up"].dropna().tolist(),
            *rif_sig["rif_down"].dropna().tolist(),
            *adeno_cons["consensus_up"].dropna().tolist(),
            *adeno_cons["consensus_down"].dropna().tolist(),
            *lesion_deg.loc[(lesion_deg["adj.P.Val"] < 0.05) & (lesion_deg["logFC"] > 0), "gene"].head(30).tolist(),
            *lesion_deg.loc[(lesion_deg["adj.P.Val"] < 0.05) & (lesion_deg["logFC"] < 0), "gene"].head(30).tolist(),
            "IGFBP1", "PRL", "LEFTY2", "FOXO1", "WNT4", "HAND2", "IL15", "SPP1", "GPX3",
            "PGR", "ESR1", "GREB1", "IHH", "HOXA10", "HOXA11", "KLF9", "NR2F2",
        }
    )

    adata = sc.read_h5ad(HECA_PATH, backed="r")
    obs = adata.obs[["celltype", "lineage", "Stage", "Binary Stage", "dataset", "sample", "n_counts"]].copy()
    obs["broad_celltype"] = [map_broad_celltype(ct, ln) for ct, ln in zip(obs["celltype"], obs["lineage"])]
    obs["stage_simple"] = [simplify_stage(x) for x in obs["Stage"]]

    var_index = pd.Index(adata.var_names.astype(str))
    keep_genes = [g for g in selected_genes if g in set(var_index)]
    gene_idx = var_index.get_indexer(keep_genes)
    x = adata[:, gene_idx].X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    expr = normalize_selected(x, obs["n_counts"].to_numpy(dtype=float))
    expr_df = pd.DataFrame(expr, columns=keep_genes)

    ref_long = pd.concat([obs.reset_index(drop=True), expr_df], axis=1)
    ref_long = ref_long[ref_long["broad_celltype"].isin(
        ["Luminal_Epi", "Glandular_Epi", "Stroma", "Decidual_Stroma", "Endothelial", "Lymphoid", "Myeloid"]
    )]

    celltype_summary = (
        ref_long.groupby(["broad_celltype", "stage_simple"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["broad_celltype", "stage_simple"])
    )
    celltype_summary.to_csv(TABDIR / "heca_broad_celltype_stage_counts.csv", index=False)

    centroid = (
        ref_long.groupby("broad_celltype", observed=True)[keep_genes]
        .mean()
        .T.reset_index()
        .rename(columns={"index": "gene"})
    )
    centroid.to_csv(TABDIR / "heca_broad_centroids_logexpr.csv", index=False)

    stage_centroid = (
        ref_long[ref_long["stage_simple"].isin(["Proliferative", "Secretory_Early", "Secretory_EarlyMid", "Secretory_Mid", "Secretory_Late"])]
        .groupby(["broad_celltype", "stage_simple"], observed=True)[keep_genes]
        .mean()
        .reset_index()
    )
    stage_centroid.to_csv(TABDIR / "heca_broad_stage_centroids_logexpr.csv", index=False)

    meta_summary = (
        ref_long.groupby(["broad_celltype", "Binary Stage"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["broad_celltype", "Binary Stage"])
    )
    meta_summary.to_csv(TABDIR / "heca_broad_binary_stage_counts.csv", index=False)

    gene_table = pd.DataFrame({"gene": keep_genes})
    gene_table.to_csv(TABDIR / "heca_selected_genes.csv", index=False)


if __name__ == "__main__":
    main()