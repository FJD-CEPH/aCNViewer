#aCNViewer
comprehensive genome-wide visualization of absolute copy number and copy neutral variations

**Contact:** Victor Renault / Alexandre How-Kit (aCNViewer@cephb.fr)

##Table of contents
- [Dependencies](#dependencies)
- [Overview](#overview)
- [Usage](#usage)
- [Tutorial](#tutorial):
  + Quick start:
    * [Affymetrix example](#affymetrix)
    * [ASCAT CNVs example](#testaffy2)
    * [Sequenza CNVs example](#testsequenzacnvs)
  + [Affymetrix](#affymetrix)
  + [Illumina](#illumina)
  + [NGS](#ngs)
  + [Processing CNV files](#processing-cnv-file)
  + Other plots:
    * [Plot heatmaps](#plotheatmaps)
    * [Plot dendrograms](#plotdendrograms)
    * [plot all](#plotall)
  


***

##Dependencies:

* APT ([Affymetrix Power Tools](http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx#1_2)) if you plan to process Affymetrix SNP arrays

* a recent version of R with [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) installed for generating the different graphical outputs:
  + [ASCAT](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) (will be automatically installed if not already installed) if you are analyzing raw SNP array data
  + [Sequenza](https://cran.r-project.org/web/packages/sequenza/index.html) if you are analyzing paired (tumor / normal) bams

* Python with version &ge; 2.5


##Overview:
![Overview of CNViewer:](/img/aCNViewer.png?raw=true "Overview of aCNViewer")  


##Usage:



***


##Tutorial

###Processing raw data

####Affymetrix

Let's respectively call `aCNViewer_TEST_DATA` and `BIN_DIR` the location where respectively [aCNViewer_TEST_DATA.tar.gz](https://drive.google.com/file/d/0B9ZcXWVM-9y1SDktTTBjVVd1ZVk/view?usp=sharing) and [APT archive]((http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx#1_2)) have been uncompressed into.

#####Requirements:

* Test data set [aCNViewer_TEST_DATA.tar.gz](https://drive.google.com/file/d/0B9ZcXWVM-9y1SDktTTBjVVd1ZVk/view?usp=sharing)

* Download [Affymetrix Power Tools](http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx#1_2) from Affymetrix website and uncompress it into `BIN_DIR`



#####TestAffy: generate a quantitative stacked histogram from CEL files with a window size of 2Mbp
`python aCNViewer.py -f CEL_DIR -c CHR_SIZE_FILE -t OUTPUT_DIR --histogram 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR --gcFile ASCAT_GC_FILE --platform AFFY_PLATFORM -l AFFY_LIB_DIR --gw6Dir GW6_DIR`<br>
where:
* `CEL_DIR` is the folder containing ".cel" ou ".cel.gz" files
* `CHR_SIZE_FILE`: a tab-delimited file with 2 columns respectively chromosome name and chromosome length
* `WINDOW_SIZE`: segment size in bp. Please note that alternatively, `-p PERCENTAGE` can be used instead of `-w WINDOW_SIZE` in order to set the segment size in percentage of chromosome length where `PERCENTAGE` is a floating number between 0 and 100
* `CENTROMERE_FILE`: : file giving the centromere bounds. Can be generated using `curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/BUILD/database/cytoBand.txt.gz" | gunzip -c | grep acen > centro_build.txt`
* `ASCAT_GC_FILE`: GC content file necessary for ASCAT GC correction when analyzing SNP array data. Please check [ASCAT website](https) for available GC content files
* `AFFY_PLATFORM`: name of ASCAT supported Affymetrix platform with a GC content file available ("Affy250k_sty", "Affy250k_nsp", "Affy500k" or "AffySNP6"). Please refer to [ASCAT website](//www.crick.ac.uk/peter-van-loo/software/ASCAT) for more details
* `AFFY_LIB_DIR`: Affymetrix library file downloadable from [Affymetrix website](http://www.affymetrix.com/support/technical/byproduct.affx?cat=dnaarrays)
* `GW6_DIR` refers to the folder where [gw6.tar.gz](http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz) has been uncompressed into. This archive contains different programs and files necessary to process Affymetrix SNP array

Here is an example:

`python aCNViewer.py -f aCNViewer_TEST_DATA/snpArrays250k_sty/ -c aCNViewer_TEST_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --histogram 1 -C aCNViewer_TEST_DATA/snpArrays250k_sty/centro.txt -w 2000000 -b BIN_DIR --gcFile aCNViewer_TEST_DATA/snpArrays250k_sty/GC_Affy250k.txt --platform Affy250k_sty -l aCNViewer_TEST_DATA/snpArrays250k_sty/LibFiles/ --gw6Dir aCNViewer_TEST_DATA/snpArrays250k_sty/gw6/`

If ASCAT is not installed and if you want to install it into a custom folder, please add the following option to the previous command line: `--rLibDir RLIB`

Compare `OUTPUT_DIR2/lrr_baf.segments_merged_hist_2000000nt.png` with `aCNViewer_TEST_DATA/snpArrays250k_sty/expectedResults/test1/lrr_baf.segments_merged_hist_2000000nt.png`


####Illumina

#####TestIlluminaWithoutNormalization

`python aCNViewer.py -f REPORT_FILES -c CHR_SIZE_FILE --histogram 1 -m 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR [--sampleList SAMPLE_TO_PROCESS_FILE] -n 0 --probeFile PROBE_POS_FILE --platform ILLUMINA_PLATFORM -g ASCAT_GC_FILE`<br>
where:
  * `REPORT_FILES` is the list of Illumina final report files to process specified either as a comma-separated string with all the report files to process or as a directory containing these files. Each Illumina final report file should contain at least the following columns:
    - `SNP Name`
    - `Sample ID`
    - `Log R Ratio`
    - `B Allele Freq`
  * `CHR_SIZE_FILE`: a tab-delimited file with 2 columns respectively chromosome name and chromosome length
  * `CENTROMERE_FILE`: file giving the centromere bounds. Can be generated using `curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/BUILD/database/cytoBand.txt.gz" | gunzip -c | grep acen > centro_build.txt`
  * `WINDOW_SIZE`: segment size in bp. Please note that alternatively, `-p PERCENTAGE` can be used instead of `-w WINDOW_SIZE` in order to set the segment size in percentage of chromosome length where `PERCENTAGE` is a floating number between 0 and 100
  * `PROBE_POS_FILE`: file listing the probes used on the SNP array with their genomic position. The file is tab-delimited with the following columns:
    - `Name`
    - `Chr`
    - `MapInfo`
  * `ILLUMINA_PLATFORM`: name of ASCAT supported Illumina platform with a GC content file available ("Illumina660k" or "HumanOmniExpress12"). Please refer to [ASCAT website](//www.crick.ac.uk/peter-van-loo/software/ASCAT) for more details
  * `ASCAT_GC_FILE`: GC content file necessary for ASCAT GC correction when analyzing SNP array data. Please check [ASCAT website](https) for available GC content files
  * `SAMPLE_TO_PROCESS_FILE`: optional, used to specify list of samples to process in one of the following formats:
    - a comma-separated string listing all the samples to process
    - the name of text file with one line per sample to process
    - the name of a Python dump file with the extension ".pyDump"

#####TestIlluminaWithTQN

Same command as above with `-n 1` instead of `-n 0`.


####NGS

#####ExampleSequenza: 

######CreateSequenzaCNVs

Here are the different options available, with the most relevant option first, to generate Sequenza CNVs when analyzing paired bam files:

1. run Sequenza without intermediary mpileup file creation and by chromosomes (this maximizes space and time):<br>
`python aCNViewer.py -P sequenza -r REF_FILE -b BIN_DIR -d BAM_DIR [--pattern BAM_FILE_PATTERN] -o TARGET_DIR --sampleFile SAMPLE_FILE --createMpileUp 0 -n NB_THREADS --byChr 1 -M MEMORY [> LOG_FILE]`<br>
where:
  * `BAM_DIR` refers to the location of the bam files
  * `BAM_FILE_PATTERN` is optional and is used when there are several bams associated to one sample. The default value is ".bam"
  * `SAMPLE_FILE` is used to list all the pairs of samples to process. It is a tab-delimited file with the following columns in the order below:
    - `idvdName`
    - `sampleName`
    - `seqFile`
    - `patientType` with value "N" for normal or "T" for tumor
  * `NB_THREADS` indicates the number of threads that will be used to create chromosomal seqz files for each sample pair
  * `MEMORY` in GB to run Sequenza. The default value is 8 (GB) ans should work for most WES analysis
  * `LOG_FILE` is optional and will keep a record of all the submitted jobs if specified


2. run Sequenza with intermediary mpileup file creation and by chromosomes (use this option only if Sequenza is freezing on some chromosomes):
the command is the same as above with `createMpileUp` set to `1`:<br>
`python aCNViewer.py -P sequenza -r REF_FILE -b BIN_DIR -d BAM_DIR [--pattern BAM_FILE_PATTERN] -o TARGET_DIR --sampleFile SAMPLE_FILE --createMpileUp 1 -n NB_THREADS --byChr 1 -M MEMORY [> LOG_FILE]`


3. run Sequenza without intermediary mpileup file creation on the whole bams (not recommended as it is slower than the previous 2 steps):
the command is the same as in 1 with `byChr` set to 0 or removed from the command and without `-n`:<br>
`python aCNViewer.py -P sequenza -r REF_FILE -b BIN_DIR -d BAM_DIR [--pattern BAM_FILE_PATTERN] -o TARGET_DIR --sampleFile SAMPLE_FILE --createMpileUp 0 [--byChr 0] -M MEMORY [> LOG_FILE]`


4. run Sequenza with intermediary mpileup file creation on the whole bams (only if step 1, 2 and 3 failed):
the command is the same as above with `createMpileUp` set to `1`:<br>
`python aCNViewer.py -P sequenza -r REF_FILE -b BIN_DIR -d BAM_DIR [--pattern BAM_FILE_PATTERN] -o TARGET_DIR --sampleFile SAMPLE_FILE --createMpileUp 1 [--byChr 0] -M MEMORY [> LOG_FILE]`

Once Sequenza CNVs have been generated sucessfully, you can proceed to the [graph generation](#testsequenzacnvs).


###Processing CNV file

Both examples below require to download [aCNViewer_TEST_DATA.tar.gz]().

####TestAffy2: generate quantitative stacked histogram from ASCAT segment files with a window size of 2Mbp

`python aCNViewer.py -f ASCAT_SEGMENT_FILE -c CHR_SIZE_FILE -t OUTPUT_DIR --histogram 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR`<br>
where:
* `ASCAT_SEGMENT_FILE`: ASCAT segment file (`ascat.output$segments` obtained by running `ascat.runAscat`) with the following columns:
  + `sample`
  + `chr`
  + `startpos`
  + `endpos`
  + `nMajor`
  + `nMinor`
* <a id="chrSize"></a>`CHR_SIZE_FILE`: a tab-delimited file with 2 columns respectively chromosome name and chromosome length
* `CENTROMERE_FILE`: file giving the centromere bounds. Can be generated using `curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/BUILD/database/cytoBand.txt.gz" | gunzip -c | grep acen > centro_build.txt`
* `WINDOW_SIZE`: segment size in bp. Please note that alternatively, `-p PERCENTAGE` can be used instead of `-w WINDOW_SIZE` in order to set the segment size in percentage of chromosome length where `PERCENTAGE` is a floating number between 0 and 100

An example can be found below:

`python aCNViewer.py -f aCNViewer_TEST_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_TEST_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --histogram 1 -C aCNViewer_TEST_DATA/snpArrays250k_sty/centro.txt -w 2000000 -b BIN_DIR`

Compare `OUTPUT_DIR/GSE9845_lrr_baf.segments_merged_hist_2000000nt.png` with `aCNViewer_TEST_DATA/snpArrays250k_sty/expectedResults/test2/GSE9845_lrr_baf.segments_merged_hist_2000000nt.png`


####TestSequenzaCNVs: generate quantitative stacked from Sequenza results with a window size of 2Mbp
`python aCNViewer.py -f SEQUENZA_RES_DIR --fileType Sequenza -c CHR_SIZE_FILE -t TARGET_DIR --histogram 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR`<br>
where:
* `SEQUENZA_RES_DIR` is the folder containing Sequenza results (`*_segments.txt`)
* <a href="#chrSize">`CHR_SIZE_FILE`</a>

An example can be found below:

`python aCNViewer.py -f aCNViewer_TEST_DATA/wes/ --fileType Sequenza -c aCNViewer_TEST_DATA/wes/hg19.chrom.sizes -t TARGET_DIR --histogram 1 -C aCNViewer_TEST_DATA/wes/centro_hg19.txt -w 2000000 -b BIN_DIR`

Compare generated histogram with `aCNViewer_TEST_DATA/wes/expectedResults/ascat_merged_hist_2000000nt.png`


###Other plots

All the different grapgs below can be generated using the 

####PlotHeatmaps




####PlotDendrograms

####PlotAll
