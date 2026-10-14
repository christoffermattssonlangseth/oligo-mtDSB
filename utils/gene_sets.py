"""Prespecified gene sets (from the literature, restricted to genes on the Xenium 5k panel).

Fixed before use in the neighbourhood and staging analyses so that they are not selected from these data.
Genes absent from the panel are listed in comments.

ISR           ATF4 / mitochondrial integrated stress response targets (as utils.ol_response.ISR_GENES)
MYELIN        mature-OL / myelin genes (absent: Plp1, Mobp, Tf, Gjc2, Plekhb1)
DAM           disease-associated microglia, Keren-Shaul et al. 2017 (absent: Clec7a, Spp1, Apoe, Ctsb, Ctsd, Lgals3, Tyrobp)
HOMEOSTATIC   homeostatic microglia
ASTRO_REACTIVE reactive astrocytes, Zamanian 2012 / Liddelow 2017 (absent: Vim, Lcn2, Steap4, Cp, Aspg, Hspb1, Srgn)
ENDO_ACTIVATION endothelial activation / leukocyte adhesion (absent: Lcn2, Lrg1, Il6)
MHC1          MHC class I antigen presentation (absent: H2-Q7)
IFN           type-I interferon response (absent: Rsad2, Bst2, Ifitm3)
PROLIF        cell cycle / proliferation (absent: Ccnb2, Mcm6, Ube2c, Tpx2, Hmgb2)
OPC           oligodendrocyte precursor cells (absent: Lhfpl3; Gpr17 counted under NEWLY_FORMED_OL)
NEWLY_FORMED_OL committed / newly formed OLs (absent: Bcas1, Tmem2, Nfasc, Plp1; Mag left out, also mature)
AXON_INJURY   regeneration-associated genes induced in injured neurons (absent: Sprr1a, Atf3, Hspb1, Npy,
              Klf6, Csrp3; neurofilaments left out because their direction after injury is ambiguous)
"""

from utils.ol_response import ISR_GENES

GENE_SETS = {
    "ISR": list(ISR_GENES),
    "MYELIN": ["Mbp", "Mog", "Mag", "Cldn11", "Cnp", "Mal", "Opalin", "Ugt8a", "Aspa", "Ermn", "Myrf", "Enpp6", "Car2"],
    "DAM": ["Cst7", "Lpl", "Itgax", "Axl", "Cd74", "Ctsl", "Igf1", "Cd9", "Trem2", "Cybb", "H2-Ab1", "H2-Aa", "Gpnmb", "Ank", "Csf1"],
    "HOMEOSTATIC": ["P2ry12", "Tmem119", "Cx3cr1", "Sall1", "Csf1r", "Siglech", "Hexb", "Fcrls", "Gpr34", "Olfml3"],
    "ASTRO_REACTIVE": ["Gfap", "Serpina3n", "C4b", "Cd44", "Osmr", "Timp1", "Socs3", "Cxcl10", "S100a10", "C3", "H2-D1", "Gbp2", "Psmb8"],
    "ENDO_ACTIVATION": ["Vcam1", "Icam1", "Sele", "Selp", "Ackr1", "Cxcl1", "Cxcl12", "Plvap", "Tnfaip3"],
    "MHC1": ["H2-D1", "H2-K1", "B2m", "Tap1", "Psmb8", "Psmb9", "Nlrc5"],
    "IFN": ["Isg15", "Ifit1", "Ifit3", "Irf7", "Stat1", "Oasl2", "Usp18"],
    "PROLIF": ["Mki67", "Top2a", "Cdk1", "Ccnb1", "Pcna", "Mcm2", "Cenpf", "Birc5", "Aurkb"],
    "OPC": ["Pdgfra", "Cspg4", "Ptprz1", "Olig1", "Sox6", "Vcan", "Tnr", "Cspg5"],
    "NEWLY_FORMED_OL": ["Enpp6", "Tcf7l2", "Gpr17", "Fyn", "Bmp4", "Itpr2", "Prom1"],
    "AXON_INJURY": ["Gap43", "Jun", "Sox11", "Stmn2", "Tubb3", "Gadd45a", "Ecel1", "Sema6a", "Gal", "Ngfr",
                    "Lgals1", "Smad1", "Cd44"],
}

IMMATURE_OL = ["Immature Oligodendrocytes I", "Immature Oligodendrocytes II", "Immature Oligodendrocytes III"]
NEURONS = ["Neurons", "Excitatory Neurons", "GABAergic Neurons", "Medium Spiny Neurons"]

# which neighbour readouts are tested in which cell type (neighbourhood analysis)
NEIGHBOUR_READOUTS = {
    "Microglia": ["DAM", "HOMEOSTATIC", "MHC1", "IFN"],
    "Astrocytes": ["ASTRO_REACTIVE", "MHC1", "IFN"],
    "Endothelial Cells": ["ENDO_ACTIVATION", "MHC1", "IFN"],
}

# repair (OL lineage) and neuronal readouts, tested with the same neighbourhood framework
REPAIR_READOUTS = {ct: ["PROLIF", "OPC", "NEWLY_FORMED_OL"] for ct in IMMATURE_OL}
NEURON_READOUTS = {ct: ["AXON_INJURY", "ISR"] for ct in NEURONS}
