# genome-signal-tracks

This repo provides tools for creating 

1. Mappability file generation
1. Signal track production

## Mappability file generation

Create mappability tracks for a reference genome, at specific kmer lengths. Mappability is based on unique, exact matches. It is inspired by the mappability pipeline created by Anshul Kundaje for [align2rawsignal](https://code.google.com/p/align2rawsignal/)

### Requirements

1. [ensembl-hive](https://github.com/Ensembl/ensembl-hive) tested with version 2.2
1. A database for use with hive (developed with MySQL 5.5.27, but should work with anything hive supports)
1. A cluster workload mangement system supported by hive. This pipeline was developed with Platform LSF 8.0.1.
1. Perl (tested with 5.12.5)
1. [ bowtie](http://bowtie-bio.sourceforge.net/index.shtml) tested with version 1.1.1 The pipeline is not compatible with bowtie2
1. [samtools](http://www.htslib.org/) tested with version 1.2
1. [bedtools](http://bedtools.readthedocs.org/en/latest/) tested with version 2.22.0
1. [bedGraphToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig) from the [Kent src tree](https://github.com/ENCODE-DCC/kentUtils)
1. Additional perl module dependencies - please see cpanfile in the root of the repository. If you have cpanm, you can install them with `cpanm --installdeps .` 
1.  A reference genome, split into fasta files. One chromosome/contig per file works well, as the each file is loaded into memory by default. The eHive pipeline is tested to work with chromsomes up to the size of human chr1.
1. A bowtie index file for the reference genome
1. A chrom sizes file for use with bedGraphToBigWig - 2 column, tab-separated  text file, listing the reference sequence names and lengths. 

### Installation

1. Get the code: `git clone git@github.com:FAANG/genome-signal-tracks.git`
1. Add  the library to your PERL5LIB:  `export PERL5LIB=$PERL5LIB:genome-signal-tracks/lib`
1. Get all software and data listed in Requirements

### Usage

1. First, initialise the hive pipeline and tell it where to find the software it depends upon:
```
ensembl-hive/scripts/init_pipeline.pl Bio::GenomeSignalTracks::PipeConfig::MappabilityConf \
	-hive_driver mysql \
	-host db_host \
	-port db_port \
	-user db_user \
	-password db_password \
	-bedtools /path/to/bedtools \
	-samtools /path/to/samtools \
	-bowtie /path/to/bowtie \
	-bedGraphToBigWig /path/to/bedGraphToBigWig \
	-fasta_suffix fa
```
The paths to bedtools, samtool and bowtie are required. The fasta suffix is used when searching for fasta files and defaults to fa (i.e. it will look for files matching *.fa). 

Running init_pipeline.pl successfully will provide you with a database URL to use in subsequent steps.

2. Secondly, tell hive what to work on:
```
ensembl-hive/scripts/seed_pipeline.pl \
	-url mysql://username:password@host:port/dbname \
	-logic_name start \
	-input_id "{fasta_dir => '/path/to/fasta/', output_dir => '/path/to/output', kmer_sizes => '35,42,90..100', index_dir => '/path/to/bowtie_index', index_name => 'name_of_index', chrom_list => '/path/to/chrom_list'}"
```
 * __fasta\_dir__ should be a directory containing fasta files  (matching the fasta suffix supplied in step 1).
 * __output\_dir__ is the directory where output files will be created. Working directories for each kmer length will also be creted here
 * __kmer\_size__ a list of kmer lengths to use. These should match the read lengths you expect to work with. This can include lists of values or ranges, so '35,45..50', would cause the pipeline to run for kmer lengths 35,45,45,47,48,49,50.
 * __index\_dir__ should be the directory containing the bowtie index
 * __index\_name__ should be the base name of the bowtie index, located in index\_dir
 * __chrom\_list___ should be a tab separated file listing the sequences in the reference genome, and their length
 
3. Run hive. The controller script (beekeeper) may be running for many hours, and should be run under [gnu screen](http://www.gnu.org/software/screen/) or similar.
```
screen
beekeeper.pl -url mysql://username:password@host:port/dbname -loop
```
4. Once the pipeline has finished, you should have bam, bedGraph and bigBed files (one per kmer size) in the __output\_dir__. Output will be named based on the index_name and kmer size, e.g. Sscrofa102.k42.bg

You can repeat steps 2 and 3 for as many reference genomes and kmer lengths as you require. The pipeline cleans up intermediate files files as it goes along.

## Signal track production

TODO

### Requirements

1. An alignment in BAM format. It must be coordinate sorted. 
1. A mappability track for the reference genome and read length of interest
1. [WiggleTools](https://github.com/Ensembl/Wiggletools), tested with version 1.
