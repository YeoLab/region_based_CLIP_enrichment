# region_based_CLIP_enrichment
Pipeline for region-based enrichment. This workflow is defined using
CWL v1.0 spec.

### Prerequisites:
- For Yeolab: module already installed
- For everyone else:
    - samtools 1.3+
    - perl 5.22+
        - Statistics::Basic
        - Statistics::Distributions
        - Statistics::R
    - CWL 1.0
    - R (any version that works with Statistics::R)

- You can try running:
    - ```source create_environment``` (installs Perl 5.22)
    - ```source run_perlbrew_perl5.10.1.sh``` (attempts to use perlbrew to install Perl 5.10)

### How to install:
- For Yeolab: module already installed
    - (On TSCC) ```module load eclipregionnormalize``` (doesn't work as of 1/12, C library issues)
    - Please use the 'create environment strategy above, this should still work'
- For everyone else:
    - Please ensure that the following are in your $PATH
    (or ```source add_paths.sh``` in the top directory.
    ```create_environments.sh``` should already do this for you):
        - bin/
        - bin/perl/
        - cwl/

### How to run:
- Fill out the requisite fields in the yaml file (please see: example/wf_region_based_enrichment.204_01.yaml as an example).

### Notes:
- Make sure Perl and its modules are properly installed! You may run into C problems, resulting in negative p-values!
- I suggest running the workflow in its entirety for one sample (if batch processing) to ensure no silent error messages are thrown

#### Description of requisite fields:

##### Essential Input Files
This is the resulting BAM file that is the merged product of multiple
PCR-duplicate-removed IP-associated barcoded files (read 2 only).
```
clipBamFile:
  class: File
  path: inputs/204_01_RBFOX2.merged.r2.bam
```
This is the equivalent size-matched input BAM file that is the merged product of one or more
PCR-duplicate-removed INPUT-associated barcoded files (read 2 only). In the current protocol,
these input files are not assigned any inline barcodes and remain unassigned after
the demultiplexing step.
```
inputBamFile:
  class: File
  path: inputs/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam
```

##### Default Parameters - shouldn't change except between species
By default, this is a gencode-formatted GTF file for the species of interest.
```
gencodeGTFFile:
  class: File
  path: inputs/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
```
This is a parsed UCSC-formatted file parsed from the above GTF file.
```
gencodeTableBrowserFile:
  class: File
  path: inputs/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
```

##### Output files
These files contain read counts for each of the following regions (broad features):
- CDS
- 3utr
- 5utr
- 3utr
- 5utr|3utr (regions that overlap both 5' and 3' UTR)
- intron
- noncoding_exon
- noncoding_intron
- intergenic (regions that do not fall into any category)
The combinedOutput will contain a tabbed concatenation of clipBroadFeatures and inputBroadFeatures:
```
clipBroadFeatures: 204_01_RBFOX2.clip.broadfeaturecounts.txt
inputBroadFeatures: 204_01_RBFOX2.input.broadfeaturecounts.txt
combinedOutput: 204_01_RBFOX2.combined.ReadsByLoc
```
This file contains the log2 fold changes for each region of each gene. NaNs that are
present in these files indicate ALL of the following:
- that less than 10 reads were present in either the CLIP or INPUT sample.
- that less than 10 reads were present in either the CLIP or the EXPECTED INPUT
- that less than 10 reads were present in either the INPUT or the EXPECTED CLIP

where:
- EXPECTED INPUT = clip reads * (input mapped read number / clip mapped read number),
- EXPECTED CLIP = input reads * (clip mapped read number / input mapped read number)

l2fcWithPvalEnr files contain l2fc values followed by a -log10 pvalue via fisher or chisquared test, delimited by '|'
- (ie. 1.1614|8.82646 indicates a log2 fold change of 1.1614 and a -log10p value of 8.82646).
```
l2fc: 204_01_RBFOX2.l2fc.txt
l2fcWithPvalEnr: 204_01_RBFOX2.l2fc_significant_regioncalls.txt
```

##### Misc output files (may be useful as inputs to another workflow in the future)
Specify these output files that will contain the mapped read number for your BAM files.
```
clipMappedReadNum: 204_01_RBFOX2.clip.mapped_read_num
inputMappedReadNum: 204_01_RBFOX2.input.mapped_read_num
```