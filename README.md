# mLinker
### Tools for Analyzing Long and Linked Read Sequencing

Table of Contents
------------

  * [Installation](#gh-md-toc)
  * [Description](#description)
  * [Commands](#commands)
  * [Input/Output](#input/output)

Installation
------------

From the mlinker directory run `./build.sh` or install manually. This code requires bamtools, htslib, c++11, and zlib libraries.

  * htslib: https://github.com/samtools/htslib
  * bamtools: https://github.com/pezmaster31/bamtools

Installing htslib locally:
```bash
git clone https://github.com/samtools/htslib
cd htslib; autoheader; autoconf
./configure --prefix=/path/to/mlinker/packages/htslib/
make; make install
cd ..
```

Installing bamtools locally:
```bash
git clone git://github.com/pezmaster31/bamtools.git
cd bamtools; mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/mlinker/packages/bamtools/ ..
make; make install
cd ../..
```

Installing python dependencies and scripts:
```bash
cd <root of mlinker>
virtualenv -p python3 env
source env/bin/activate
./setup.py [develop|install]
```

Make an empty output directory if not present:
```bash
mkdir ./output
```

Description
-----------
mLinker is a suite of C++ tools useful for interpreting long and linked read sequencing of cancer genomes.  The most significant information long read sequencing provides is the local haplotype of a sample.  In cancer cell lines where Aneuplpoidy, Loss of Heterozygosity, and Structural Variation is common, haplotypes can provide a better resolution of the samples karyotype, and clarify the cancer cells genomic evolution.

This program was used to determine Whole Chromosome Haplotypes as detailed here: [biorxiv](https://www.biorxiv.org/content/10.1101/629337v1).  mLinker currently supports 10X, Pacbio, Oxford Nanopore, and HiC sequencing technologies.  A diagram of the germline haplotype phasing workflow is shown below.

<p align="center">
<img src="https://github.com/rwtourdot/mlinker/blob/main/new_linker_flowchart.png" width=900/>
</p>

Associated plotting scripts found in `./R_plots` can be used to determine block energy cutoffs, assess haplotype accuracy, and plot haplotype specific copy number.

<p align="center">
<img src="https://github.com/rwtourdot/mlinker/blob/main/plot_output.png" width=600/>
</p>

Commands
--------

**List Commands**
```bash
./mlinker -h
```
<!--
#### Flags
  * -i: /path/to/input/bamfile
  * -v: /path/to/input/vcffile
  * -u: /path/to/input/svfile
  * -l: /path/to/haplotype_solution_file
  * -t: aneuploid (tumor) coverage file
  * -m: diploid (normal) coverage file
  * -e: sequencing technology (tenx,pacbio,nanopore,illumina)
  * -c: chromosome (chr1,chr2,...,chrX)
  * -s: (optional) second bam file /path/to/second/bamfile
  * -n: (optional) id string for output files
  * -b: (optional) binsize (default is 10kb - 10000)
-->


**Extract Phasing Information from Long and Linked Reads**

This command extracts all long and linked read phasing information from an aligned bamfile given a corresponding vcf file containing heterozygous sites.  The long read technology flag (tenx,pacbio,nanopore,hic) should correspond to the bamfile chosen.

```bash
./mlinker extract -i ./input.bam -v het_sites.vcf -e tenx -c chr21 -n sample_name
```
  * Output is a graph file: graph_variant_{}.dat

The output of this command is a graph_variant file which lists all of the unique hashes associated with each het site base.  This file can be concatenated with other graph_variant files, to combine all phasing information from multiple technologies.


**Solve For Local Haplotype Blocks from Long Read Links**

After extracting all phasing information into graph_variant files and combining or trimming certain hashes the samples haplotype can be solved for by:

```bash
./mlinker solve -i ./graph_variant_{}_chr21.dat -c chr21 -n sample_name
```
  * Output is the haplotype solution file: hap_solution_{}.dat

The haplotype file contains a Block Switch Energy column which can be used to define blocks.  The lower (more negative) the block energy the more likely a heterozygous site is phased correctly. More details on defining blocks can be found in the paper.


**Generate A Whole Chromosome Haplotype Scaffold**

This command takes a haplotype solution file and combines it with hic phasing information in a corresponding graph_variant_hic file to generate a full chromosome haplotype scaffold.  An energy cutoff which defines haplotype blocks is specified by the -e flag.

```bash
./mlinker scaffold -i ./hap_solution_{}_chr21.dat -g ./graph_variant_hic_chr21.dat -e -700 -c chr21 -n sample_name
```
  * Output is a haplotype scaffold file: hap_full_scaffold_{}.dat

The hap_full_scaffold file contains less heterozygous sites than the haplotype solution file but is accurate over the full length of the chromosome.


**Recover and Phase Variants to a Haplotype Scaffold**

This command takes a haplotype scaffold file and phases all variants in a graph file to that scaffold.  This is useful since some small haplotype blocks and their variants are lost when generating a haplotype scaffold.  Though this command was built to recover variants, it can also be used to phase somatic variants to the germline haplotype - given a corresponding graph file.

```bash
./mlinker recover -i ./hap_full_scaffold_{}_chr21.dat -g ./graph_variant_{}_chr21.dat -e tech -c chr21 -n sample_name
```
  * Output is a recovered haplotype file: hap_recovered_{}.dat

The recovered haplotype file is similar to the scaffold file but its additional columns explicitly count the number of links connecting the reference/variant base to haplotype A/B.


**Extract Heterozygous Site Coverage**

This command takes a vcf and bam file and extracts the read coverage of each base. Base counts at each heterozygous site must pass a base quality and map quality cutoff.

```bash
./mlinker coverage -i ./input.bam -v het_sites.vcf -e illumina -c chr21 -n sample_name
```
  * Output is heterozygous site coverage file: het_coverage.dat

 More information on this file is described in the I/O section below.

<!--
#### Phase Germline Haplotypes from Long and Linked Reads

This commmand takes in a vcf file and a long or linked read bam file to compute phased haplotype blocks.  The vcf file should contain all germline heterozygous sites and most likely originates from a paired normal sample.  The bam file or files should be obtained with a long read technology and could originate for tumor, normal, or tumor+normal.

```
./mlinker phase -i ./input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
./mlinker phase -i ./input.bam -s ./second_input.bam -v het_sites.vcf -e pacbio -c chr4 -n august15
```
  * Output is haplotype solution file: haplotype_solution.dat

The output of this command is a file which contains the minimum energy solution to the germline haplotype.  More information on this file is described in the I/O section below.
-->

<!--
#### Extract Heterozygous Site Coverage

This command takes a vcf and bam file and extracts the read coverage of each allele. In order to count a base at a heterozygous site both the base quality and read map quality must pass a cutoff.

```
./mlinker coverage -i ./input.bam -v het_sites.vcf -e illumina -c chr4 -n august15
```
  * Output is heterozygous site coverage file: coverage.dat

The output of this command is a file which contains the read depth coverage of each heterozygous site for both reference and variant bases.  More information on this file is described in the I/O section below.

#### Phase Aneuploid Samples based on Copy Number

Haplotypes can be phased further based on tumor copy number.  Copy number phasing works better in tumor samples where aneuploidy and loss of heterozygocity is prevalent.

```
./mlinker cn_phase -l haplotype_solution.dat -m ./het_coverage_aneuploid.dat -n august15
```
  * Output is a copy number phased haplotype file: cn_haplotype.dat

#### Phase Structural Variants (10X/Nanopore)

Once haplotypes are found, associated Structural Variants can be phased with a 10X or Nanopore tumor sample.  Structual Variant's can be called with a SV caller and converted to the svfile input file format described in the I/O section below.

```
./mlinker sv_phase -i ./input.bam -v het_sites.vcf -e tenx -u ./svfile.dat -n august15
```
  * Output is a phased sv file: sv_phased.dat

#### Create Linked Read Matrix

A local alignment map can be extracted from any long or linked read technology. This map does not contain any allelic information but can more clearly show translocations and inversions.

```
./mlinker matrix -i ./input.bam -e tenx -c chr4 -b 10000 -n august15
```
  * Output is a matrix file: matrix.dat

#### Filter Het Sites by Coverage and Coverage Fraction

In order refine variant calls it can be useful to use allele fraction from a normal sample to extract a more confident subset of variants.  This command takes in a coverage file and filters het sites based on allelic fraction.  More specifically, this command checks if the fraction of variant bases is between 10 - 90 percent and passes it to a filtered output file.  This command can also filter tumor coverage data based on the allele fraction in a corresponding normal sample.

```
./mlinker filter -m normal_coverage.dat -n august15
./mlinker filter -t tumor_coverage.dat -m normal_coverage.dat -n august15
```
  * Output is a filtered coverage file: filtered_coverage.dat

#### Bin 10X Sample by BX tag

10X sequencing coverage data can be refined by looking at unique bx tag density.  This command takes a 10X bam file and bins the genome by coverage of unique bx tag.

```
./mlinker bx_bin -i ./input.bam -e tenx -c chr4 -b 10000 -n august15
```
-->

Input/Output
--------

### Required Files

**.bam file**

The .bam file can originate from Long Read or Linked Read sequencing. The alignment method will depend on the type of sequencing, but a consistent genome reference should be used between samples.

**.vcf file**

The .vcf file can be obtained with GATK (https://software.broadinstitute.org/gatk/) or another variant caller and should contain all heterozygous sites.  Indels will not be considered in mlinkers phasing methods, and SNP's called in centromere and variable regions of the genome should be removed.  .vcf files should have the same reference as the long or linked read .bam file

<!--
**SV file**

There are several Structural Variant callers that work with paired end sequencing. The output of their calls may differ, mlinker takes structural variant input in the form below.

```
RAindex	chr1	pos1	str1	chr2	pos2	str2	TotalCount
20  chr1	1651206	-1	chr1	1717357	1	3
```
-->

### Generated Files

**graph_variant file**

A file which extracts all phasing information in a .bam file.  These files can be concatenated given the same sample and chromosome.
```
variant_id  chrom position  tech_hash readname
46692825_C_C	chr21	46692825	GATAACCGTACTGCTA-1	HMVT3CCXX:6:2107:9790:22739
```

**hap_solution file**

A file that contain all haplotype blocks solved for from a graph_variant file.  The **haplotype** column defines each het sites local haplotype in a block defined by **blockE**.
```
index position  reference_id  alternate_id  haplotype ref_count alt_count spinE blockE  default_block range
0	13000241	chr21_13000241_A_A	chr21_13000241_A_C	1	23	30	-2083.39	-2083.39	0	27362
```

**scaffold file**

A file that contains the final haplotype accurate over the whole chromosome.  The **scaffold_hap** column is a spin vector which defines the germline haplotype of the sample.
```
chrom arm position  reference_id  alternate_id  scaffold_hap  block_hap block spinE blockE
chr21	q	13207111	13207111_T_T	13207111_T_C	-1	block:	1	18	tenx_energy(spin/block):	-168	-168
```

**het_coverage file**

A file that contains base counts at each heterozygous site.
```
variant_id  position  variant:reference Indel Deletion  Gbase Cbase Abase Tbase Total
190204289_T_C	190204289	T:C	I|0	D|0	G|0	C|44	A|0	T|30	74
```
