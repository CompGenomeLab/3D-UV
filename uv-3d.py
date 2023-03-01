import pandas as pd
import bioframe as bf
import polars as pl
import numpy as np


def rpkm(name: str, count: str, length: str) -> pl.Expr:
    return (
        pl.col(f"{count}")
        / ((pl.col(f"{length}")) / 1000 * (pl.sum(f"{count}") / 1000000))
    ).alias(f"rpkm_{name}")


genes = pl.read_csv(
    "ensembl_hg38_genes.tsv",
    sep="\t",
    has_header=True,
)
pcoding_genes = genes.filter(pl.col("Gene type") == "protein_coding")
pcoding_genes = pcoding_genes.select(
    [
        pl.col("Chromosome/scaffold name").alias("chrom"),
        pl.col("Gene start (bp)").alias("start"),
        pl.col("Gene end (bp)").alias("end"),
    ]
).sort("chrom")
pcoding_genes = pcoding_genes.select(
    [(pl.lit("chr") + pl.col("chrom")).alias("chrom"), pl.col("start"), pl.col("end")]
)

pcoding_genes_merged: pl.DataFrame = pl.from_pandas(bf.merge(pcoding_genes.to_pandas()))


HELA_coords = pl.read_csv(
    "HELA_hg38_coords_sorted.bed",
    sep="\t",
    has_header=False,
    new_columns=["chrom", "start", "end", "x", "y", "z"],
).with_columns(
    [
        (pl.col("end") - pl.col("start")).alias("length"),
        (np.sqrt(pl.col("x") ** 2 + pl.col("y") ** 2 + pl.col("z") ** 2)).alias("r"),
    ]
)

HELA = HELA_coords.with_column(
    pl.when(pl.col("r") < 1)
    .then("0-1")
    .when(pl.col("r") < 2)
    .then("1-2")
    .when(pl.col("r") < 3)
    .then("2-3")
    .when(pl.col("r") < 4)
    .then("3-4")
    .when(pl.col("r") < 5)
    .then("3-4")
    .alias("label")
)

# UV damage
HELA_damage_64 = pl.read_csv(
    "ds_64.bed",
    sep="\t",
    has_header=False,
    columns=[0, 1, 2],
    new_columns=["chrom", "start", "end"],
)
HELA_damage_cpd = pl.read_csv(
    "ds_cpd.bed",
    sep="\t",
    has_header=False,
    columns=[0, 1, 2],
    new_columns=["chrom", "start", "end"],
)

# Repair
HELA_repair_64 = pl.read_csv(
    "/xr_64.bed",
    sep="\t",
    has_header=False,
    columns=[0, 1, 2],
    new_columns=["chrom", "start", "end"],
)
HELA_repair_cpd = pl.read_csv(
    "xr_cpd.bed",
    sep="\t",
    has_header=False,
    columns=[0, 1, 2],
    new_columns=["chrom", "start", "end"],
)

## HeLa
# Input simulations
# UV Damage
HELA_input_sim_damage_64 = pl.read_csv(
    "ds_64_sim_input.bed",
    sep="\t",
    has_header=False,
    columns=[0, 1, 2],
    new_columns=["chrom", "start", "end"],
)
HELA_input_sim_damage_cpd = pl.read_csv(
    "ds_cpd_sim_input.bed",
    sep="\t",
    has_header=False,
    columns=[0, 1, 2],
    new_columns=["chrom", "start", "end"],
)

# Repair
HELA_input_sim_repair_64 = pl.read_csv(
    "xr_64_sim_input.bed",
    sep="\t",
    has_header=False,
    columns=[0, 1, 2],
    new_columns=["chrom", "start", "end"],
)
HELA_input_sim_repair_cpd = pl.read_csv(
    "xr_cpd_sim_input.bed",
    sep="\t",
    has_header=False,
    columns=[0, 1, 2],
    new_columns=["chrom", "start", "end"],
)

overlap: pl.DataFrame = pl.from_pandas(
    bf.overlap(HELA.to_pandas(), pcoding_genes_merged.to_pandas(), return_overlap=True)
).select(
    [
        pl.all().exclude("n_intervals_"),
        (pl.col("overlap_end") - pl.col("overlap_start")).alias("genic_length"),
    ]
)

# groups
group = [
    ("damage_64", HELA_damage_64),
    ("damage_cpd", HELA_damage_cpd),
    ("repair_64", HELA_repair_64),
    ("repair_cpd", HELA_repair_cpd),
    ("damage_64_sim", HELA_input_sim_damage_64),
    ("damage_cpd_sim", HELA_input_sim_damage_cpd),
    ("repair_64_sim", HELA_input_sim_repair_64),
    ("repair_cpd_sim", HELA_input_sim_repair_cpd),
]

# get reads overlapping with genes
for g_name, g_data in group:
    print(g_name)

    overlap: pl.DataFrame = pl.from_pandas(
        bf.count_overlaps(
            overlap.to_pandas(),
            g_data.to_pandas(),
            cols1=["chrom", "overlap_start", "overlap_end"],
        )
    ).select(
        [pl.all().exclude("count"), pl.col("count").alias(f"{g_name}_genic_count")]
    )

# group genes in each bead
overlap = overlap.groupby(
    ["chrom", "start", "end", "x", "y", "z", "r", "label", "length"]
).agg(
    [
        pl.sum("genic_length").alias("genic_length"),
        pl.sum("damage_64_genic_count").alias("damage_64_genic_count"),
        pl.sum("damage_cpd_genic_count").alias("damage_cpd_genic_count"),
        pl.sum("repair_64_genic_count").alias("repair_64_genic_count"),
        pl.sum("repair_cpd_genic_count").alias("repair_cpd_genic_count"),
        pl.sum("damage_64_sim_genic_count").alias("damage_64_sim_genic_count"),
        pl.sum("damage_cpd_sim_genic_count").alias("damage_cpd_sim_genic_count"),
        pl.sum("repair_64_sim_genic_count").alias("repair_64_sim_genic_count"),
        pl.sum("repair_cpd_sim_genic_count").alias("repair_cpd_sim_genic_count"),
    ]
)

# reads overlapping with beads
for g_name, g_data in group:
    print(g_name)

    overlap: pl.DataFrame = pl.from_pandas(
        bf.count_overlaps(overlap.to_pandas(), g_data.to_pandas())
    ).select([pl.all().exclude("count"), pl.col("count").alias(f"{g_name}_count")])

# reads in intergenic regions
overlap = overlap.with_columns(
    [
        (pl.col("length") - pl.col("genic_length")).alias("intergenic_length"),
        (pl.col("damage_64_count") - pl.col("damage_64_genic_count")).alias(
            "damage_64_intergenic_count"
        ),
        (pl.col("damage_cpd_count") - pl.col("damage_cpd_genic_count")).alias(
            "damage_cpd_intergenic_count"
        ),
        (pl.col("repair_64_count") - pl.col("repair_64_genic_count")).alias(
            "repair_64_intergenic_count"
        ),
        (pl.col("repair_cpd_count") - pl.col("repair_cpd_genic_count")).alias(
            "repair_cpd_intergenic_count"
        ),
        (pl.col("damage_64_sim_count") - pl.col("damage_64_sim_genic_count")).alias(
            "damage_64_sim_intergenic_count"
        ),
        (pl.col("damage_cpd_sim_count") - pl.col("damage_cpd_sim_genic_count")).alias(
            "damage_cpd_sim_intergenic_count"
        ),
        (pl.col("repair_64_sim_count") - pl.col("repair_64_sim_genic_count")).alias(
            "repair_64_sim_intergenic_count"
        ),
        (pl.col("repair_cpd_sim_count") - pl.col("repair_cpd_sim_genic_count")).alias(
            "repair_cpd_sim_intergenic_count"
        ),
    ]
)

# get rpkm values for each condition
overlap = overlap.with_columns(
    [
        rpkm("damage_64", "damage_64_count", "length"),
        rpkm("damage_64_genic", "damage_64_genic_count", "genic_length"),
        rpkm("damage_64_intergenic", "damage_64_intergenic_count", "intergenic_length"),
        rpkm("damage_cpd", "damage_cpd_count", "length"),
        rpkm("damage_cpd_genic", "damage_cpd_genic_count", "genic_length"),
        rpkm(
            "damage_cpd_intergenic", "damage_cpd_intergenic_count", "intergenic_length"
        ),
        rpkm("repair_64", "repair_64_count", "length"),
        rpkm("repair_64_genic", "repair_64_genic_count", "genic_length"),
        rpkm("repair_64_intergenic", "repair_64_intergenic_count", "intergenic_length"),
        rpkm("repair_cpd", "repair_cpd_count", "length"),
        rpkm("repair_cpd_genic", "repair_cpd_genic_count", "genic_length"),
        rpkm(
            "repair_cpd_intergenic", "repair_cpd_intergenic_count", "intergenic_length"
        ),
        rpkm("damage_64_sim", "damage_64_sim_count", "length"),
        rpkm("damage_64_sim_genic", "damage_64_sim_genic_count", "genic_length"),
        rpkm(
            "damage_64_sim_intergenic",
            "damage_64_sim_intergenic_count",
            "intergenic_length",
        ),
        rpkm("damage_cpd_sim", "damage_cpd_sim_count", "length"),
        rpkm("damage_cpd_sim_genic", "damage_cpd_sim_genic_count", "genic_length"),
        rpkm(
            "damage_cpd_sim_intergenic",
            "damage_cpd_sim_intergenic_count",
            "intergenic_length",
        ),
        rpkm("repair_64_sim", "repair_64_sim_count", "length"),
        rpkm("repair_64_sim_genic", "repair_64_sim_genic_count", "genic_length"),
        rpkm(
            "repair_64_sim_intergenic",
            "repair_64_sim_intergenic_count",
            "intergenic_length",
        ),
        rpkm("repair_cpd_sim", "repair_cpd_sim_count", "length"),
        rpkm("repair_cpd_sim_genic", "repair_cpd_sim_genic_count", "genic_length"),
        rpkm(
            "repair_cpd_sim_intergenic",
            "repair_cpd_sim_intergenic_count",
            "intergenic_length",
        ),
    ]
)

# damage normalized with simulation
A = ["", "_genic", "_intergenic"]
B = ["64", "cpd"]

for a in A:
    for b in B:
        o = o.with_column(
            (pl.col(f"rpkm_damage_{b}{a}") / pl.col(f"rpkm_damage_{b}_sim{a}")).alias(
                f"DS_{b}{a}"
            )
        )


# repair normalized with simulation
A = ["", "_genic", "_intergenic"]
B = ["64", "cpd"]

for a in A:
    for b in B:
        o = o.with_column(
            (pl.col(f"rpkm_repair_{b}{a}") / pl.col(f"rpkm_repair_{b}_sim{a}")).alias(
                f"RS_{b}{a}"
            )
        )

# double normalized repair
A = ["", "_genic", "_intergenic"]
B = ["64", "cpd"]

for a in A:
    for b in B:
        o = o.with_column(
            (pl.col(f"RS_{b}{a}") / pl.col(f"DS_{b}{a}")).alias(f"RR_{b}{a}")
        )
