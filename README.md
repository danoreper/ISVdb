# ISVdb v1.2
Code to generate the whole genome data stored in the Inbred Strain Variant Database (v1.2). 
In particular, this code imputes genotypes and diplotypes for CC strains, and writes this information out to zipped csv files (along with the founder strains as well). It may potentially become a CRAN package eventually but will reside on Github for now. It is more likely that we will eventually create a non-CRAN R package. 

# Requirements
* OSX or linux
* VCFTools>0.1.13
* R>=3.2.3
* Python>=2.7.12
* R CRAN packages: data.table, reshape2, ggplot2, igraph, stringr
* R Bioconductor packages: Biostrings, IRanges, biomaRt, GenomicRanges, rtracklayer
* Python packages: PyYAML, ipdb

# Download
Git or SVN checkout using the web URL: https://github.com/danoreper/ISVdb.git

# Install
Enter the following at the command line, from ISVDB_LOCATION (the cloned, local ISVDB repository):
ISVDB\_LOCATION\$ bash install.sh 

For now, this install script only downloads data files (that are too large for github), and places them in their expected location relative to source code. Note that the time to complete installation may be lengthy, depending on network speed. VCFTools, R, Python, and Python and R libraries need to be self-installed by the user.

# To create all csv files representing genotypes and diplotyes:
1. Navigate to ISVDB_LOCATION/src as the working directory
2. Enter the following at command line: 
ISVDB\_LOCATION/src\$ R CMD BATCH ./genomerep/variantdb2/dbdriver.R

2a. Running locally will--- realistically --- require hundreds of hours. If you have access to UNC research computing, computation can finish in a few hours on the Killdevil computing cluster. To use Killdevil, enter the following command while logged intto killdevil:

ISVDB\_LOCATION/src$ R CMD BATCH '--args ../config/defaultCluster.yaml' ./genomerep/variantdb/dbdriver.R

2b. Another option to reduce runtime is to run locally, but restrict the size of the database that is created. To do so, edit ISVDB\_LOCATION/config/defaultParams.yaml. In particular, edit the following properties:

* **variantdb:  var\_limit: .na** The max number of variants of each type (indel, snp) that the vcf parsing will bother storing. Meant for testing. .na means build them all.

 
* **variantdb: varcc\_limit: .na** The max number of cc lines that variant db will bother storing. Meant for testing. .na means build all.



* **variantdb:  chr\_range: .na** The chromosomes that will be built. .na means build them all




# Directory structure
* ISVDB\_LOCATION/src: The source code that generates genotype and diplotype files
* ISVDB\_LOCATION/config: configuration files for running locally or on killdevil.
* ISVDB\_LOCATION/data: Data files for generating genotype and diplotype files. Largely empty before instll.sh is run.
* ISVDB\_LOCATION/output/isvdb: Generated genotype and diplotype files, for both exon and whole genome
* ISVDB\_LOCATION/output/isvdb/exon1410: exon information
* ISVDB\_LOCATION/output/isvdb/full1410: whole genome information

Restricting description to exons,
* ISVDB\_LOCATION/output/isvdb/exon1410/genotype: exon genotypes
* ISVDB\_LOCATION/output/isvdb/exon1410/diplotype: exon diplotypes
* ISVDB\_LOCATION/output/isvdb/exon1410/genotype\_sampling: genotype sampling files (useful for simulating crosses)
* ISVDB\_LOCATION/output/isvdb/exon1410/diplotype\_sampling: diplotype sampling files (useful for simulating crosses)

* ISVDB\_LOCATION/v1.1: version 1.1 of this project, which used an SQL representation, and only contained exons+/- 100bp rather than the whole genome.

# Output file organization
Output files are organized in folders according to strain and then chromosome. For example, 
ISVDB\_LOCATION/output/isvdb/exon1410/diplotype contains a folder corresponding to the diplotypes for StrainName1, with a separate file per chromosome.

---StrainName1

------1.txt.tar.gz

------2.txt.tar.gz

------3.txt.tar.gz

------...

------MT.txt.tar.gz

------X.txt.tar.gz

------Y.txt.tar.gz

# Diplotype output format
Each *.txt.tar.gz file in the diplotypes folder is a compressed csv file of diplotypes for a particular strain and chromosome, with 2 header rows:
1) A row specifying the strain and chromosome of the format: strain:strainname,chr:chrname 
2) The fields stored in the csv file. 

The fields consist of the following:
* variant_id: a positive integer ID specifiying the variant. Unique to every variant; chromosomes don't share variant_id. However, this ID is shared across strains, and is consistent between the diplotype dump and the genotype dump files.
* pos: positive integer specifying variant position in bp along chromosome.
* founder\_1: one of the founder haplotypes at this variant for the file-specified strain and chromosome. This is part of an unphased diplotype.
* founder\_2: the second founder haplotypes at this variant for the file-specified strain and chromosome. This is part of an unphased diplotype.
* prob: the probability on [0,1] that the unphased diplotype of the variant at this position is (founder\_1, founder\_2) 
* gene_name: the name of a gene enclosing the variant.


Note that there may be multiple records per variant if the diplotype of the variant is uncertain; in such a case there is one row per non-zero diplotype probability, and the probabilities add approximately to 1. There also may be multiple rows per variant if a variant is enclosed by more than one gene.


# Genotype output format
Each *.txt.tar.gz file in the genotype folder is a compressed csv file of genotypes for a particular strain and chromosome, with 2 header rows:
1) A row specifying the strain and chromosome of the format: strain:strainname,chr:chrname 
2) The fields stored in the csv file. 
The fields consist of the following:
* variant_id: a positive integer ID specifiying the variant. Unique to every variant; chromosomes don't share variant_id. However, this ID is shared across strains, and is consistent between the diplotype dump and the genotype dump files.
* pos: positive integer specifying variant position in bp along chromosome.
* allele\_1: one of the alleles at this variant for the file-specified strain. This is part of an unphased genotype.
* allele\_2: the second alelle at this variant for the file-specified strain. This is part of an unphased diplotype.
* prob: the probability that the unphased diplotype of the variant at this position is (allele_1, allele_2) 
* is\_max: whether this is the max likelihood genotype for this variant in this strain. Redundant with prob, but available for convenience.
* gene\_name: the name of a gene enclosing the variant.
Note that there may multiple records per variant if the diplotype of the variant is uncertain; in such a case there is one row per non-zero diplotype probability, and the probabilities add approximately to 1. There also may be multiple rows per variant if a variant is enclosed by more than one gene.
* transcript\_name: the name of a transcript enclosing the variant.
* consequence\_1: the function consequence of allele_1 on transcript transcript_name, with respect to the B6 reference allele. If the allele_2 is the reference allele, the consequence is "reference". 
* consequence\_2: the function consequence of allele_2 on transcipt ranscript_name, with respect to the B6 reference allele. If allele_2 is the reference allele, the consequence is "reference". 

Note that there may multiple records per variant if the genotype of the variant is uncertain; in such a case there is one row per non-zero genotype probability, and the probabilities add approximately to 1. There also may be multiple rows per variant if a variant is enclosed by more than one gene and/or more than one transcript, as the consequence changes depending on the transcript.

