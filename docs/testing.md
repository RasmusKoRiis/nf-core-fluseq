# Testing and acceptance

The repository uses three levels of verification:

1. Python unit tests exercise report transformations, mutation-reference mapping, input contracts, and helper programs.
2. Nextflow and nf-test smoke tests parse the workflows and run small components against committed synthetic fixtures.
3. A private end-to-end acceptance run exercises the human FASTQ workflow with controlled references and representative operational data.

The committed fixtures contain invented records only. Real sequence data, operational sample sheets, controlled reference sets, work directories, and acceptance results must remain outside the Git repository.

## Local checks

Run the checks used by CI:

```bash
python3 -m pip install --requirement requirements-test.txt
NXF_SYNTAX_PARSER=v1 nextflow config . -profile docker,server
NXF_SYNTAX_PARSER=v1 nextflow run . --help
python3 -m pytest -p no:cacheprovider
nf-test test tests/modules/local/surveillance_summary.nf.test modules/local/cat_fastq/tests/main.nf.test
```

The Python suite also validates all local `meta.yml` files and their input/output declarations. These checks use committed synthetic fixtures and a vendored metadata schema; they do not require private data or local reference database paths. Passing them does not establish complete module coverage or biological acceptance. See [nfcore_practices.md](nfcore_practices.md) for the remaining gaps.

The pipeline currently requires the v1 syntax parser when it is run with Nextflow 26.04 or newer.

## Private human FASTQ acceptance run

This is an optional local acceptance procedure, not a prerequisite for contributing documentation or running the public component tests. FASTQ files contain sequencing reads; they do not supply the separate HA/NA references, protein reference sequences, mutation workbooks, and Nextclade datasets needed for a full analysis. Reference paths are needed only on the computer executing that analysis. They can stay in its private parameter file and do not need to be shared or committed.

Create a private YAML file outside the repository containing the controlled reference parameters. It must define at least `ha_database`, `na_database`, `inhibition_mutation_db`, `reassortment_database`, `sequence_references`, and `nextclade_dataset`. Do not put sequence data, sample identifiers, credentials, or local reference paths in a tracked file.

Then run:

```bash
export FLUSEQ_ACCEPTANCE_DATA=/private/path/to/run-data
export FLUSEQ_ACCEPTANCE_PARAMS=/private/path/to/references.yml
export FLUSEQ_ACCEPTANCE_OUTDIR=/private/path/to/acceptance-results
scripts/run_private_acceptance.sh
```

The script requires exactly one CSV file near the private data root and exactly one `fastq_pass` directory beneath it. It refuses to use repository paths for private data, parameters, work files, or results. A successful Nextflow exit and a final CSV in `reporthuman/` are both required.

Set `FLUSEQ_ACCEPTANCE_PROFILE` to override the default `docker` profile. Set `FLUSEQ_ACCEPTANCE_WORKDIR` to retain a known work directory for `-resume`; otherwise the script creates one below `/tmp` and prints its location through Nextflow.
