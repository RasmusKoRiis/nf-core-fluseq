process REPORTHUMAN {
    label 'process_single'



    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    path(subtype)
    path(coverage)
    path(mutation_human)
    path(mutation_inhibtion)
    path(lookup)
    path(nextclade_summary_ha)
    path(nextclade_sample)
    path(mutation_vaccine)
    val runid
    val release_version
    path(filtered_fasta)
    path(irma_depth)
    val seq_instrument
    path(samplesheet)
    path(reassortment_report)
    path(subclade_nomenclature_report)
    
    output:

    path("${runid}.csv"), emit: report
    path("${runid}.fasta"), emit: filtered_fasta
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when


    script:
    """ 
    set -euo pipefail

    # Generate date
    current_date=\$(date '+%Y-%m-%d')


    # Convert CSV to TSV with a real CSV parser so quoted commas are preserved.
    python - ${samplesheet} samplesheet.tsv <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as source, open(sys.argv[2], "w", newline="") as target:
    csv.writer(target, delimiter="\t", lineterminator="\n").writerows(csv.reader(source))
PY

    report.py samplesheet.tsv

    #Add constant parameters to the report
    # Add RunID column
    awk -v runid=${runid} -v OFS=',' '{ if (NR == 1) { print  \$0, "RunID" } else { print  \$0, runid } }' merged_report.csv > ${runid}_temp1.csv

    # Add Instrument ID column
    awk -v seq_instrument=${seq_instrument} -v OFS=',' '{ if (NR == 1) { print  \$0, "Instrument ID" } else { print  \$0, seq_instrument } }' ${runid}_temp1.csv > ${runid}_temp2.csv

    # Add Date column
    awk -v date="\$current_date" -v OFS=',' '{ if (NR == 1) { print \$0, "Date" } else { print \$0, date } }' ${runid}_temp2.csv > ${runid}_temp3.csv

    # Add Release Version column
    awk -v version="${release_version}" -v OFS=',' '{ if (NR == 1) { print \$0, "Release Version" } else { print \$0, version } }' ${runid}_temp3.csv > ${runid}_temp4.csv

    report_QC_calculation.py ${runid}_temp4.csv -o ${runid}.csv

    #Merge all filtered fasta files to one
    cat ${filtered_fasta} > ${runid}.fasta

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    
    """

}
