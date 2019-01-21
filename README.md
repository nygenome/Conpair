# Conpair
Conpair: concordance and contamination estimator for tumor–normal pairs

Conpair is a fast and robust method dedicated for human tumor-normal studies to perform concordance verification (= samples coming from the same individual), as well as cross-individual contamination level estimation in whole-genome and whole-exome sequencing experiments. Importantly, our method of estimating contamination in the tumor samples is not affected by copy number changes and is able to detect contamination levels as low as 0.1%.

* Version: 0.2
* Author: Ewa A Bergmann
* Contact: ewa.a.bergmann@gmail.com

**Required input files:** two bam files (tumor, normal)

**Required software:** GATK 2.3 or later, python 2.7 or higher, scipy, numpy, java

**Required data:** Human genome file (GRCh37 or GRCh38) 

GRCh37:

The fasta file can be downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz  
In order to be able use the fasta file as a reference 2 additional files are required:
`human_g1k_v37.dict`, `human_g1k_v37.fa.fai`  
To create these files please follow: http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference

GRCh38:

The fasta file can be downloaded from: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
In order to be able use the fasta file as a reference 2 additional files are required:
`hg38.dict`, `hg38.fa.fai`
To create these files please follow: http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference


# Manual

**Dependencies:**

* python 2.7 or higher :     [www.python.org](https://www.python.org/)
* GATK 3 or higher :         [www.broadinstitute.org/gatk/download](http://www.broadinstitute.org/gatk/download/)
* java :                     [http://java.com](http://java.com/en/download/)

**Installation:**   

Use pip to install the package:
```
pip install .
```
If you are getting permission issues, you may want to install into a [virtual environment](https://docs.python.org/3/tutorial/venv.html).

To set up GATK, you should either set the environment variable:
```
export GATK_JAR=/your/path/to/GenomeAnalysisTK.jar
```

Or install gatk4 through [bioconda](https://bioconda.github.io/):
```
conda install -c bioconda gatk
```

**Default reference genome:**

To avoid specifying the reference file every time you run Conpair, you can put/link the following files in the specified directory:
`/your/path/to/CONPAIR/genomes/human_g1k_v37.fa`  
`/your/path/to/CONPAIR/genomes/human_g1k_v37.fa.fai`  
`/your/path/to/CONPAIR/genomes/human_g1k_v37.dict`
<br/>

**Most common usage and additional options:**   
To generate pileups (GATK required):
```

run_gatk_pileup_for_sample.py -B TUMOR_bam -O TUMOR_pileup
run_gatk_pileup_for_sample.py -B NORMAL_bam -O NORMAL_pileup


Optional:
--help                              show help message and exit
--reference REFERENCE               reference genome in the fasta format, two additional files (.fai, .dict) located in the same directory as the fasta file are required. You may choose to avoid specifying the reference by following the steps in the "default reference genome" section above.
--markers MARKERS                   the set of preselected genomic positions in the BED format. Default: ${CONPAIR_DIR}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed
--conpair_dir CONPAIR_DIR           path to ${CONPAIR_DIR}
--gatk GATK                         path to GATK JAR [$GATK by default]
--java JAVA                         path to JAVA [java by default]
--temp_dir_java TEMP_DIR_JAVA       java temporary directory to set -Djava.io.tmpdir
--xmx_java  XMX_JAVA                Xmx java memory setting [default: 12g]

```
Verifying concordance between two samples (tumor and normal):
```  

${CONPAIR_DIR}/scripts/verify_concordance.py -T TUMOR_pileup -N NORMAL_pileup


Optional:
--help                              show help message and exit
--outfile OUTFILE                   write output to OUTFILE
--normal_homozygous_markers_only    use only normal homozygous positions to calculate concordance between TUMOR and NORMAL 
--min_cov MIN_COV                   require min of MIN_COV in both TUMOR and NORMAL to use the marker
--min_mapping_quality MIN_MAP_QUAL  do not use reads with mapping qual below MIN_MAP_QUAL [default: 10]
--min_base_quality  MIN_BASE_QUAL   do not use reads with base qual below MIN_BASE_QUAL of a specified position [default: 20]
--markers MARKERS                   the set of preselected genomic positions in the TXT format. Default: ${CONPAIR_DIR}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt

```  
Estimating contamination level in both the tumor and the normal:
```

${CONPAIR_DIR}/scripts/estimate_tumor_normal_contamination.py -T TUMOR_pileup -N NORMAL_pileup


Optional:
--help                              show help message and exit
--outfile OUTFILE                   write output to OUTFILE
--min_mapping_quality MIN_MAP_QUAL  do not use reads with mapping qual below MIN_MAP_QUAL [default: 10] 
--markers MARKERS                   the set of preselected genomic positions in the TXT format. Default: ${CONPAIR_DIR}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt
--conpair_dir CONPAIR_DIR           path to ${CONPAIR_DIR}
--grid  GRID                        grid interval [default: 0.01]

```  

# Output files
**Pileup**  
An example of a pileup file (10 first lines) can be viewed here: ([`pileup.txt`](https://github.com/nygenome/Conpair/blob/master/data/example/pileup/NA12878_normal40x.gatk.pileup.10lines.txt)).

**Concordance**  
An example of a concordance file can be viewed here: ([`concordance.txt`](https://github.com/nygenome/Conpair/blob/master/data/example/concordance/NA12878_tumor80x--NA12878_normal40x.concordance.txt)). 

**Contamination**  
An example of a concordance file can be viewed here: ([`contamination.txt`](https://github.com/nygenome/Conpair/blob/master/data/example/contamination/NA12878_tumor80x--NA12878_normal40x.contamination.txt)). 


# Interpretation  
**Concordance**  
To eliminate the effect of copy number variation on the concordance levels, we recommend using the -H flag. 
If two samples are concordant the expected concordance level should be close to 99-100%.  
For discordant samples concordance level should be close to 40%.  
You can observe slighly lower concordance (80-99%) in presence of contamination and/or copy number changes (if the -H option wasn't used) in at least one of the samples.   
<br/>
**Contamination**   
Even a very low contamination level (such as 0.5%) in the tumor sample will have a severe effect on calling somatic mutations, resulting in decreased specificity. Cross-individual contamination in the normal sample usually has a milder effect on somatic calling.
