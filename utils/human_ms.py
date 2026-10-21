"""Human orthologs of the mtDSB OL gene sets, for scoring human MS single-nucleus data (analysis/16).

Mouse symbols map to human by upper-casing, with a few explicit exceptions (Serpina3n -> SERPINA3, MHC class I,
Tmem173 -> STING1 / TMEM173). Sets:

* ``MTDSB_CORE``  the mtDSB-specific OL markers (ATF4 / ISR + p53 / p21, analysis/10 and /13)
* ``ISR``, ``P53``  the prespecified pathway sets of utils/gene_sets.MECHANISM_SETS
* ``DRIVERS_UP`` / ``DRIVERS_DOWN``  the data-derived mtDSB OL driver genes (results/ol_gene_drivers)
* ``IFNG``        the interferon-gamma programme that marks EAE OLs (analysis/10), as a comparator
* ``IMMUNE_OL``   the immune-visible OL core shared by mtDSB and EAE
* ``MYELIN``      prespecified myelin genes
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from utils.gene_sets import GENE_SETS, MECHANISM_SETS
from utils.ol_response import ISR_GENES

SPECIAL = {"Serpina3n": ["SERPINA3"], "H2-D1": ["HLA-A", "HLA-B", "HLA-C"], "H2-K1": ["HLA-A", "HLA-B", "HLA-C"],
           "Tmem173": ["STING1", "TMEM173"], "Gbp7": ["GBP1", "GBP2", "GBP5"], "Irgm1": ["IRGM"], "Cd74": ["CD74"]}


def to_human(genes) -> list[str]:
    out = []
    for g in genes:
        for h in SPECIAL.get(g, [g.upper()]):
            if h not in out:
                out.append(h)
    return out


def gene_sets(repo: Path) -> dict[str, list[str]]:
    drivers = pd.read_csv(repo / "results" / "ol_gene_drivers" / "tables" / "driver_genes.csv", index_col=0)
    return {
        "MTDSB_CORE": to_human(["Trib3", "Gdf15", "Atf5", "Ddit3", "Cdkn1a", "Eif4ebp1", "Aldh1l2", "Sesn2", "Bbc3"]),
        "ISR": to_human(sorted(set(ISR_GENES) | set(MECHANISM_SETS["ATF4 / ISR"]))),
        "P53": to_human(MECHANISM_SETS["p53 targets"]),
        "DRIVERS_UP": to_human(drivers.index[drivers.direction == "up"]),
        "DRIVERS_DOWN": to_human(drivers.index[drivers.direction == "down"]),
        "IFNG": to_human(["Irf1", "Psmb8", "Psmb9", "Stat1", "Tap1", "Gbp7", "Irgm1", "Cxcl10", "Ciita"]),
        "IMMUNE_OL": to_human(["Serpina3n", "Cd74", "H2-D1", "B2m", "Lgals1", "Hmox1", "C4b"]),
        "MYELIN": to_human(GENE_SETS["MYELIN"]),
    }
