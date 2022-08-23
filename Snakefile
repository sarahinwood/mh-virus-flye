#!/usr/bin/env python3

###############
# CONTAINERS ##
###############

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
blast_container = 'docker://ncbi/blast:2.12.0'
flye_container = 'docker://staphb/flye:2.9'
hmmer_container = 'docker://biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'
medaka_container = 'docker://ontresearch/medaka:v1.7.0'
nextpolish_container = 'docker://pvstodghill/nextpolish:1.4.1__2022-07-29'
prodigal_container = 'docker://biocontainers/prodigal:v1-2.6.3-4-deb_cv1'
viralFlye_container = 'library://sinwood/viralflye/viralflye_0.2:0.0.1'

#########
# RULES #
#########

rule target:
    input:
    		# initial MhV contigs blast
    	'output/03_blast_MhV/blastn.outfmt6',
    		# stats
    	'output/05_nextpolish/bb_stats.out',
    		# annotation
    	'data/07_prodigal_blastp/prodigal_blastp.outfmt6',
    	'output/08_hmmer/prodigalPFAM_domtblout.out',
    		# hrs self blast
    	'output/09_hrs_self_blast/blastn.outfmt6'

#####################################
##  blast against self for repeats ##
#####################################

rule hrs_self_blast:
    input:
        viral_contigs = 'output/05_nextpolish/genome.nextpolish.fasta',
        db = 'output/09_hrs_self_blast/mhv.nhr'
    output:
        blast_res = 'output/09_hrs_self_blast/blastn.outfmt6'
    params:
        db = 'output/09_hrs_self_blast/mhv'
    threads:
        20
    log:
        'output/logs/09_hrs_self_blast.log'
    singularity:
    	blast_container
    shell:
        'blastn '
        '-query {input.viral_contigs} '
        '-db {params.db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blast_res} '
        '2>{log}'    

rule make_MhV_blast_db:
    input:
        'output/05_nextpolish/genome.nextpolish.fasta'
    output:
        'output/09_hrs_self_blast/mhv.nhr'
    params:
        db_name = 'mhv',
        db_dir = 'output/09_hrs_self_blast/mhv'
    log:
        'output/logs/09_make_blast_db.log'
    singularity:
    	blast_container
    shell:
        'makeblastdb '
        '-in {input} '
        '-dbtype nucl '
        '-blastdb_version 4 '
        '-title {params.db_name} '
        '-out {params.db_dir} '
        '-parse_seqids '
        '2> {log}'

##############################
## predict genes & annotate ##
##############################

## pfam all genes
##hmmscan
rule hmmscan:
    input:
        peptides_viral_contigs = 'output/06_prodigal/protein_translations.faa',
        db = 'bin/db/trinotate_dbs/Pfam-A.hmm'
    output:
        'output/08_hmmer/prodigalPFAM_domtblout.out'
    threads:
        40
    log:
        'output/logs/08_hmmscan.log'
    singularity:
        hmmer_container
    shell:
        'hmmscan '
        '--cpu {threads} '
        '--domtblout {output} '
        '{input.db} '
        '{input.peptides_viral_contigs} '
        '> {log}'

rule prodigal_blastp:
    input:
        prodigal = 'output/06_prodigal/protein_translations.faa'
    output:
        blastp_res = 'data/07_prodigal_blastp/prodigal_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        40
    singularity:
        blast_container
    log:
        'output/logs/07_prodigal_blastp.log'
    shell:
        'blastp '
        '-query {input.prodigal} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 5e-04 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'

rule prodigal:
    input:
        'output/05_nextpolish/genome.nextpolish.fasta'
    output:
        protein_translations = 'output/06_prodigal/protein_translations.faa',
        nucleotide_seq = 'output/06_prodigal/nucleotide_seq.fasta',
        gene_predictions = 'output/06_prodigal/gene_predictions.gff'
    log:
        'output/logs/06_prodigal.log'
    threads:
        1
    singularity:
        prodigal_container
    shell:
        'prodigal '
        '-i {input} '
        '-a {output.protein_translations} '
        '-d {output.nucleotide_seq} '
        '-f gff '
        '-p meta '
        '-o {output.gene_predictions} '
        '2> {log} '

#####################
## polish assembly ##
#####################

rule bb_stats:
    input:
        assembly = 'output/05_nextpolish/genome.nextpolish.fasta'
    output:
        stats = 'output/05_nextpolish/bb_stats.out'
    log:
        'output/logs/bbstats.log'
    singularity:
        bbduk_container
    shell:
        'stats.sh '
        'in={input.assembly} '
        'out={output.stats} '
        '2> {log}'

# nextpolish - short read polishing
rule nextpolish:
    input:
        assembly = 'output/04_medaka/consensus.fasta',
        config = 'data/nextpolish/hac.cfg'
    output:
        'output/05_nextpolish/genome.nextpolish.fasta'
    log:
        'output/logs/05_nextpolish.log'
    threads:
        40
    singularity:
        nextpolish_container
    shell:
        'nextPolish '
        '{input.config} 2> {log}'

# medaka - long read polishing
rule medaka:
    input:
        assembly = 'output/02_viralFlye/circulars_viralFlye.fasta',
        reads = 'data/reads/nanopore-hac-non-hyp.fq.gz'
    output:
        'output/04_medaka/consensus.fasta'
    params:
        wd = 'output/04_medaka'
    log:
        'output/logs/04_medaka.log'
    threads:
        40
    singularity:
        medaka_container
    shell:
        'medaka_consensus '
        '-i {input.reads} '
        '-d {input.assembly} '
        '-o {params.wd} '
        '-t {threads} '
        '-m r941_min_hac_g507 ' # closest minion 9.4.1 flowcell model to guppy v6.0.0

################################
## assemble & check for virus ##
################################

rule blast_against_mhv:
    input:
        fasta = 'output/02_viralFlye/circulars_viralFlye.fasta',
        db = 'output/03_blast_MhV/initial_MhV.nhr'
    output:
        blast_res = 'output/03_blast_MhV/blastn.outfmt6'
    params:
        db = 'output/03_blast_MhV/initial_MhV'
    threads:
        20
    log:
        'output/logs/03_blast_against_MhV.log'
    singularity:
    	blast_container
    shell:
        'blastn '
        '-query {input.fasta} '
        '-db {params.db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blast_res} '
        '2>{log}' 

rule make_initial_MhV_blast_db:
    input:
        'data/initial_MhV_prodigal/nucleotide_seq.fasta'
    output:
        'output/03_blast_MhV/initial_MhV.nhr'
    params:
        db_name = 'initial_MhV',
        db_dir = 'output/03_blast_MhV/initial_MhV'
    threads:
        20
    log:
        'output/logs/make_blast_db.log'
    singularity:
    	blast_container
    shell:
        'makeblastdb '
        '-in {input} '
        '-dbtype nucl '
        '-blastdb_version 4 '
        '-title {params.db_name} '
        '-out {params.db_dir} '
        '-parse_seqids '
        '2> {log}'

rule viralFlye:
    input:
        flye = 'output/01_flye/assembly.fasta',
        nanopore = 'data/reads/nanopore-hac-non-hyp.fq.gz',
    output:
        'output/02_viralFlye/circulars_viralFlye.fasta'
    params:
        flye_dir = 'output/01_flye',
        wd = 'output/02_viralFlye',
        hmm = '/viralFlye/Pfam-A.hmm.gz'
    log:
        'output/logs/02_viralFlye.log'
    threads:
        40
    singularity:
        viralFlye_container
    shell:
        'source activate viralFlye ; ' # ; allows first command to run and then second to run after with viralFlye conda env activated
        '/viralFlye/viralFlye.py '
        '--dir {params.flye_dir} '
        '--hmm {params.hmm} '
        '--reads {input.nanopore} '
        '--outdir {params.wd} '
        '--threads {threads} '
        '2>{log}'

rule flye_metagenome:
    input:
        nanopore = 'data/reads/nanopore-hac-non-hyp.fq.gz'
    output:
        'output/01_flye/assembly.fasta'
    params:
        wd = 'output/01_flye'
    log:
        'output/logs/01_flye_metagenome.log'
    threads:
        40
    singularity:
        flye_container
    shell:
        'flye '
        '--meta '
        '--nano-raw {input.nanopore} '
        '--genome-size 200k '
        '-o {params.wd} '
        '--threads {threads} '
        '2> {log}'