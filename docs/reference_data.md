# Controlled reference data

The pipeline treats the routine influenza databases as run inputs, not as mutable implementation details. Large or season-specific references must be supplied explicitly and should be released as a single tested reference bundle.

## Required bundle contents

A routine bundle normally contains:

- `human_HA.fasta` and `human_NA.fasta`;
- `H5_genotype_database.fasta` for workflows that genotype;
- `reassortment_database.fasta`;
- mammalian and inhibition mutation workbooks;
- the `sequence_references/` hierarchy consumed by mutation calling;
- the `nextclade_datasets/` hierarchy described below;
- a reference table used by the wrappers to validate seasonal reference selection.

Keep a bundle-level version, release date, source/accession metadata, and SHA-256 manifest beside these files. Update the bundle in a staging directory, validate it, and switch routine runs to it only after the complete bundle is ready. Do not update files in place while a run is active.

## Nextclade dataset layout

`--nextclade_dataset` points to a directory of local Nextclade datasets. Every dataset directory must contain `pathogen.json`. Directory names follow `<subtype>_<segment>`, for example:

```text
nextclade_datasets/
├── H1N1_HA/
├── H1N1_NA/
├── H3N2_HA/
├── H3N2_NA/
├── VIC_HA/
└── VIC_NA/
```

`B_VIC_<segment>` and `BVIC_<segment>` are accepted alternatives for Victoria-lineage influenza B. H5 datasets use the called H5 subtype, such as `H5N1_HA`. The established H5 path runs Nextclade on HA; other avian segments use the translation workflow.

The pipeline stops if a required local dataset is missing. It does not silently download a moving "latest" Nextclade dataset.

## Runtime records

Each workflow hashes the reference files staged for that run and publishes `pipeline_info/reference_manifest.tsv`. The Nextflow parameter dump records the source paths. Archive both files with the final report.

Seasonal subclade rule archives are currently retrieved from immutable Git commit URLs and recorded in `versions.yml`; their reference accessions are versioned. For a completely network-isolated deployment, mirror those immutable URLs internally or vendor the verified rule bundle and configure the process through an institutional network policy.
