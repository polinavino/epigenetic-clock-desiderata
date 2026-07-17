# Gene identifier note — GSE314970 rat retina

Source: GEO series GSE314970 supplementary file
`GSE314970_rat_raw_counts_wo_muscle.txt.gz` (first column "GeneName"),
subset to the 85 retina samples.

## Matrix row identifiers (`ensembl_gene_id` column in the counts files)
The row key is the study's "GeneName" field. It is a mixed identifier:
- 26,135 / 30,454 rows (85.8%) are **gene symbols / NCBI gene names**
  (e.g. `Rho`, `Vom2r3`, and `LOC…` NCBI identifiers for uncharacterized genes).
  These are used directly as the gene symbol.
- 4,319 / 30,454 rows (14.2%) are **bare Ensembl Rat gene IDs** (`ENSRNOG…`),
  used as a fallback where the annotation had no symbol. Of these, 218 could be
  resolved to a symbol via Ensembl BioMart (mRatBN7.2); the remainder are
  novel/predicted genes with no current symbol.

Net: 86.5% of rows carry a usable gene symbol.

## Companion file
`GSE314970_gene_symbol_map.csv` — one row per matrix gene:
- `gene_id`      : exact row key used in the counts matrix
- `gene_symbol`  : gene symbol (empty for the ~13.5% unresolved Ensembl IDs)
- `id_type`      : `symbol` (row key already a symbol/name) or `ensembl` (ENSRNOG id)

Symbol resolution for `ENSRNOG` ids used Ensembl BioMart
(`rnorvegicus_gene_ensembl`, attributes `ensembl_gene_id`, `external_gene_name`).

## Values
Raw gene-level counts (RSEM-style; non-integer values occur from
multi-mapped-read fractional assignment). Not normalized. Per-sample library
sizes range ~16.3M–28.6M. 6,890 genes are all-zero across the retina subset.
