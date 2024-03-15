
process SUBTYPEFINDER {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'

    conda "bioconda::blast=2.15.0"
    container 'docker.io/ncbi/blast:latest'

    input:
    tuple val(meta), path(ha_fasta), path(na_fasta),path(pa_fasta), path(pb1_fasta),path(pb2_fasta), path(np_fasta),path(ns_fasta), path(m_fasta)
    path(ha_database)
    path(na_database)

    output:
    tuple val(meta), path("*.txt"), emit: subtype
    tuple val(meta), path("*.csv"), emit: subtype_file
    tuple val(meta), path("*.tsv"), emit: subtype_file_test
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    blastn -query $ha_fasta -subject $ha_database -outfmt 6   > "ha_${prefix}.tsv"
    subtype_ha=\$(awk -F '\t' '{split(\$2,a,\"_\"); print a[2]}' "ha_${prefix}.tsv" | head -n 1)
    
    blastn -query $na_fasta -subject $na_database -outfmt 6   > "na_${prefix}.tsv"
    subtype_na=\$(awk -F '\t' '{split(\$2,a,\"_\"); print a[2]}' "na_${prefix}.tsv" | head -n 1)

    subtype=""$meta.id",\$subtype_ha\$subtype_na"
    subtype_txt="\$subtype_ha\$subtype_na"

    echo \$subtype > ""$meta.id"_subtype.csv"
    echo \$subtype_txt > ""$meta.id"_subtype.txt"

    blastn -query $ha_fasta -db nt -out ha.txt -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_hsps 3 -max_target_seqs 5 -remote
    blastn -query $na_fasta -db nt -out na.txt -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_hsps 3 -max_target_seqs 5 -remote
    blastn -query $pb1_fasta -db nt -out pb1.txt -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_hsps 3 -max_target_seqs 5 -remote
    blastn -query $pb2_fasta -db nt -out pb2.txt -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_hsps 3 -max_target_seqs 5 -remote
    blastn -query $pa_fasta -db nt -out pa.txt -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_hsps 3 -max_target_seqs 5 -remote
    blastn -query $ns_fasta -db nt -out ns.txt -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_hsps 3 -max_target_seqs 5 -remote
    blastn -query $np_fasta -db nt -out np.txt -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_hsps 3 -max_target_seqs 5 -remote
    blastn -query $m_fasta -db nt -out m.txt -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_hsps 3 -max_target_seqs 5 -remote

    # Read the first line of the results file
    ha=\$(head -n 1 ha.txt)
    echo \$ha
    
    # Use grep with Perl-compatible regular expressions to find patterns like HxNy
    if [[ \$ha =~ (H[0-9]+N[0-9]+) ]]; then
        subtype_ha_re=\${BASH_REMATCH[1]}
        echo \$subtype_ha_re
    else
        echo "Subtype not found."
    fi

    na=\$(head -n 1 na.txt)
    print(\$na)
    if [[ \$na =~ (H[0-9]+N[0-9]+) ]]; then
        subtype_na_re=\${BASH_REMATCH[1]}
    else
        echo "Subtype not found."
    fi

    pa=\$(head -n 1 pa.txt)
    if [[ \$pa =~ (H[0-9]+N[0-9]+) ]]; then
        subtype_pa_re=\${BASH_REMATCH[1]}
    else
        echo "Subtype not found."
    fi

    pb1=\$(head -n 1 pb1.txt)
    if [[ \$pb1 =~ (H[0-9]+N[0-9]+) ]]; then
        subtype_pb1_re=\${BASH_REMATCH[1]}
    else
        echo "Subtype not found."
    fi

    pb2=\$(head -n 1 pb2.txt)
    if [[ \$pb1 =~ (H[0-9]+N[0-9]+) ]]; then
        subtype_pb1_re=\${BASH_REMATCH[1]}
    else
        echo "Subtype not found."
    fi

    np=\$(head -n 1 np.txt)
    if [[ \$np =~ (H[0-9]+N[0-9]+) ]]; then
        subtype_np_re=\${BASH_REMATCH[1]}
    else
        echo "Subtype not found."
    fi

    ns=\$(head -n 1 ns.txt)
    if [[ \$ns =~ (H[0-9]+N[0-9]+) ]]; then
        subtype_ns_re=\${BASH_REMATCH[1]}
    else
        echo "Subtype not found."
    fi

    m=\$(head -n 1 m.txt)
    if [[ \$m =~ (H[0-9]+N[0-9]+) ]]; then
        subtype_m_re=\${BASH_REMATCH[1]}
    else
        echo "Subtype not found."
    fi

    subtype_re=""$meta.id",HA_\$subtype_ha_re,NA_\$subtype_na_re,PB2_\$subtype_pb2_re,PB1_\$subtype_pb1_re,PA_\$subtype_pa_re,NP_\$subtype_np_re,NS_\$subtype_ns_re,M_\$subtype_m_re,"
    echo \$subtype_re > ""$meta.id"_subtype_re.csv"


    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

}
