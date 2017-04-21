# aCNViewer
comprehensive genome-wide visualization of absolute copy number and copy neutral variations

**Contact:** Victor Renault / Alexandre How-Kit (aCNViewer@cephb.fr)

## Table of contents
- [Installation](#installation)
- [Overview](#overview)
- [Tutorial](#tutorial):
    * Processing SNP array data:
      + [Affymetrix](#affymetrix)
      + [Illumina](#illumina)
    * [Processing NGS data](#ngs)
    * [Processing CNV file](#processing-cnv-file)
- [Output files](#outputfiles)
  

***

## Installation

### Docker installation

The easiest way to install aCNViewer is to install the [Docker application](https://hub.docker.com/r/fjdceph/acnviewer/) (supports multi-threading but not computer clusters which are better suited for processing NGS bams):
`docker pull fjdceph/acnviewer`



### Installation from source

aCNViewer can also be installed from its source by:
1. downloading aCNViewer's data (includes test data sets and some of the third-party softwares listed in the [dependencies](#dependencies) section [APT and tQN])
1. installing the dependencies listed [below](#dependencies).
1. downloading the github source code from this page


#### Dependencies:

* APT ([Affymetrix Power Tools](http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx#1_2)) if you plan to process raw Affymetrix SNP arrays (to uncompress into `BIN_DIR`)

* a recent version of R with [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) installed for generating the different graphical outputs:
  + [ASCAT](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) (will be automatically installed if not already installed) if you are analyzing raw SNP array data
  + [Sequenza](https://cran.r-project.org/web/packages/sequenza/index.html) if you are analyzing paired (tumor / normal) bams
  + [plotrix](https://cran.r-project.org/web/packages/plotrix/index.html) for plotting dendrograms (will be automatically installed if not already installed)
  + [gplots](https://cran.r-project.org/web/packages/gplots/index.html)
  + [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

* [tQN](http://cbbp.thep.lu.se/~markus/software/tQN/tQN-1.1.2.zip) if you plan to process raw Illumina SNP arrays (to uncompress into `BIN_DIR`) and run tQN normalisation. If the cluster file for the Illumina SNP array you plan to analyze is not in the tQN lib folder, you can download additional cluster files from [here](http://cbbp.thep.lu.se/~markus/software/tQN/)

* Python with version &ge; 2.5


## Overview:
![Overview of aCNViewer:](/img/aCNViewer.png?raw=true "Overview of aCNViewer")  



***


## Tutorial

### Requirements:

Let's call `aCNViewer_DATA` the location where the test data set [aCNViewer_DATA.tar.gz](https://www.cng.fr/genodata/pub/LIVER/aCNViewer_DATA.tar.gz) has been uncompressed into and `BIN_DIR` the folder containing all third-party softwares located in `aCNViewer_DATA/bin`. Download the test data set [aCNViewer_DATA.tar.gz](https://www.cng.fr/genodata/pub/LIVER/aCNViewer_DATA.tar.gz).


### Processing SNP array data

#### Affymetrix

##### TestAffy: generate a quantitative stacked histogram from CEL files with a window size of 2Mbp

<p style="background-color:green;">Docker command:</p>

`python aCNViewer/code/aCNViewer.py -f aCNViewer_DATA/snpArrays250k_sty/ -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --histogram 1 -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 -b BIN_DIR --gcFile aCNViewer_DATA/snpArrays250k_sty/GC_Affy250k.txt --platform Affy250k_sty -l aCNViewer_DATA/snpArrays250k_sty/LibFiles/ --gw6Dir aCNViewer_DATA/snpArrays250k_sty/gw6/`

If ASCAT is not installed and if you want to install it into a custom folder, please add the following option to the previous command line: `--rLibDir RLIB`

The histogram `OUTPUT_DIR/lrr_baf.segments_merged_hist_2000000nt.png` can also be found in `aCNViewer_DATA/snpArrays250k_sty/expectedResults/test1/lrr_baf.segments_merged_hist_2000000nt.png`.

Please note that while generating the histogram, the following files are created as well:
* `*_10pc_ploidy.txt` is a matrix of segments of 10% chromosomal length for all samples. The last column indicates the calculated ploidy which corresponds to the most frequent ploidy
* `*_samples.txt` contains, for each segment, the percentage of samples belonging to each CNV key value and the list of these samples


==**Here is the full command:**==

`python aCNViewer.py -f CEL_DIR -c CHR_SIZE_FILE -t OUTPUT_DIR --histogram 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR --gcFile ASCAT_GC_FILE --platform AFFY_PLATFORM -l AFFY_LIB_DIR --gw6Dir GW6_DIR [--lohToPlot LOH_TO_PLOT] [--rColorFile RCOLOR_FILE]`<br>
where:
* `CEL_DIR` is the folder containing ".cel" ou ".cel.gz" files
* <a id="chrSize"></a>`CHR_SIZE_FILE`: a tab-delimited file with 2 columns respectively chromosome name and chromosome length
* <a id="windowSize"></a>`WINDOW_SIZE`: segment size in bp. Please note that alternatively, `-p PERCENTAGE` can be used instead of `-w WINDOW_SIZE` in order to set the segment size in percentage of chromosome length where `PERCENTAGE` is a floating number between 0 and 100
* <a id="centromereFile"></a>`CENTROMERE_FILE`: file giving the centromere bounds. Can be generated using `curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/BUILD/database/cytoBand.txt.gz" | gunzip -c | grep acen > centro_build.txt`
* `ASCAT_GC_FILE`: GC content file necessary for ASCAT GC correction when analyzing SNP array data. Please check [ASCAT website](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) for available GC content files
* `AFFY_PLATFORM`: name of ASCAT supported Affymetrix platform with a GC content file available ("Affy250k_sty", "Affy250k_nsp", "Affy500k" or "AffySNP6"). Please refer to [ASCAT website](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) for more details
* `AFFY_LIB_DIR`: Affymetrix library file downloadable from [Affymetrix website](http://www.affymetrix.com/support/technical/byproduct.affx?cat=dnaarrays)
* `GW6_DIR` refers to the folder where [gw6.tar.gz](http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz) has been uncompressed into. This archive contains different programs and files necessary to process Affymetrix SNP array
* <a id="lohToPlot"></a>`LOH_TO_PLOT`: histogram option for LOH plotting. Values should be one of "cn-LOH" for plotting cn-LOH only, "LOH" for LOH only or "both" for cn-LOH and LOH. The default value is "cn-LOH"
* <a id="rColorFile"></a> `RCOLOR_FILE`: colors in histograms (section "[histogram]". If defined, should contain exactly 10 colors [one per line] corresponding to CNV values in the following order: "&le; -4", "-3", "-2", "-1", "1", "2", "3", "4", "5", "&ge; 6"), dendrograms (section "[group]". If defined, should contain at least the same number of colors than the number of distinct values for the phenotypic / clinical feature of interest) and heatmaps (sections "[chr]" [if defined, should contain 22 colors corresponding to chromosomes 1 to 22], "[group]" and "[heatmap]" [if defined, should contain 10 colors [one per line] corresponding to CNV values in the following order: "0", "1", "2", "3", "4", "5", "6", "7", "8", "&ge; 9"]) can be redefined in that file. An example can be found [here](/img/rColor.txt).


#### Illumina

##### TestIlluminaWithoutNormalization

`python aCNViewer.py -f REPORT_FILES -c CHR_SIZE_FILE --histogram 1 -m 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR [--sampleList SAMPLE_TO_PROCESS_FILE] -n 0 --probeFile PROBE_POS_FILE --platform ILLUMINA_PLATFORM -g ASCAT_GC_FILE [--lohToPlot LOH_TO_PLOT] [--rColorFile RCOLOR_FILE]`<br>
where:
  * `REPORT_FILES` is the list of Illumina final report files to process specified either as a comma-separated string with all the report files to process or as a directory containing these files. Each Illumina final report file should contain at least the following columns:
    - `SNP Name`
    - `Sample ID`
    - `Log R Ratio`
    - `B Allele Freq`
  * <a href="#chrSize">`CHR_SIZE_FILE`</a>
  * <a href="#centromereFile">`CENTROMERE_FILE`</a>
  * <a href="#windowSize">`WINDOW_SIZE`</a>
  * `PROBE_POS_FILE`: file listing the probes used on the SNP array with their genomic position. The file is tab-delimited with the following columns:
    - `Name`
    - `Chr`
    - `MapInfo`
  * `ILLUMINA_PLATFORM`: name of ASCAT supported Illumina platform with a GC content file available ("Illumina660k" or "HumanOmniExpress12"). Please refer to [ASCAT website](//www.crick.ac.uk/peter-van-loo/software/ASCAT) for more details
  * `ASCAT_GC_FILE`: GC content file necessary for ASCAT GC correction when analyzing SNP array data. Please check [ASCAT website](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) for available GC content files
  * `SAMPLE_TO_PROCESS_FILE`: optional, used to specify list of samples to process in one of the following formats:
    - a comma-separated string listing all the samples to process
    - the name of text file with one line per sample to process
    - the name of a Python dump file with the extension ".pyDump"
  * <a href="#lohToPlot">`LOH_TO_PLOT`</a>
  * <a href="#rColorFile">`RCOLOR_FILE`</a>

##### TestIlluminaWithTQN

Same command as above with `-n 1` instead of `-n 0`.


#### NGS

##### ExampleSequenza: 

###### CreateSequenzaCNVs

This step requires an access to cluster (supported clusters are SGE, SLURM, MOAB and LSF. Tests have been sucessfully made on SGE and SLURM environments). Here are the different options available, with the most relevant option first, to generate Sequenza CNVs when analyzing paired bam files:

1. run Sequenza without intermediary mpileup file creation and by chromosomes (this maximizes space and time):<br>
`python aCNViewer.py -P sequenza -r REF_FILE -b BIN_DIR -D BAM_DIR [--pattern BAM_FILE_PATTERN] -t TARGET_DIR --sampleFile SAMPLE_FILE --createMpileUp 0 -n NB_THREADS --byChr 1 [-M MEMORY] [> LOG_FILE]`<br>
where:
  * `BAM_DIR` refers to the location of the bam files
  * `BAM_FILE_PATTERN` is optional and is used when there are several bams associated to one sample. The default value is ".bam"
  * `SAMPLE_FILE` is used to list all the pairs of samples to process. It is a tab-delimited file with the following columns in the order below:
    - `idvdName`
    - `sampleName`
    - `seqFile`
    - `patientType` with value "N" for normal or "T" for tumor
  * `NB_THREADS` indicates the number of threads that will be used to create chromosomal seqz files for each sample pair
  * `MEMORY`: optional argument specifying allocated memory in GB to run Sequenza. The default value is 8 (GB) and should work for most WES analysis
  * `LOG_FILE` is optional and will keep a record of all the submitted jobs if specified


2. run Sequenza with intermediary mpileup file creation and by chromosomes (use this option only if Sequenza is freezing on some chromosomes):
the command is the same as above with `createMpileUp` set to `1`:<br>
`python aCNViewer.py -P sequenza -r REF_FILE -b BIN_DIR -D BAM_DIR [--pattern BAM_FILE_PATTERN] -t TARGET_DIR --sampleFile SAMPLE_FILE --createMpileUp 1 -n NB_THREADS --byChr 1 [-M MEMORY] [> LOG_FILE]`


3. run Sequenza without intermediary mpileup file creation on the whole bams (not recommended as it is slower than the previous 2 steps):
the command is the same as in 1 with `byChr` set to 0 or removed from the command and without `-n`:<br>
`python aCNViewer.py -P sequenza -r REF_FILE -b BIN_DIR -D BAM_DIR [--pattern BAM_FILE_PATTERN] -t TARGET_DIR --sampleFile SAMPLE_FILE --createMpileUp 0 [--byChr 0] [-M MEMORY] [> LOG_FILE]`


4. run Sequenza with intermediary mpileup file creation on the whole bams (only if step 1, 2 and 3 failed):
the command is the same as above with `createMpileUp` set to `1`:<br>
`python aCNViewer.py -P sequenza -r REF_FILE -b BIN_DIR -D BAM_DIR [--pattern BAM_FILE_PATTERN] -t TARGET_DIR --sampleFile SAMPLE_FILE --createMpileUp 1 [--byChr 0] [-M MEMORY] [> LOG_FILE]`

Once Sequenza CNVs have been generated sucessfully, you can proceed to the [graph generation](#testsequenzacnvs).


### Processing CNV file

At the moment, ASCAT segment file and Sequenza results can be used as an input to aCNViewer. It is possible however to feed aCNViewer with CNV results from any other softwares as explained in the section [below] (#OtherCNVformats).


Both examples below require to download [aCNViewer_DATA.tar.gz](http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz).

#### TestAffy2
Generate quantitative stacked histogram from ASCAT segment files with a window size of 2Mbp:<br>

[aCNViewer_DATA.tar.gz](http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz) is required to run this example.

`python aCNViewer.py -f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --histogram 1 -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 -b BIN_DIR`

The generated histogram `OUTPUT_DIR/GSE9845_lrr_baf.segments_merged_hist_2000000nt.png` can be found in `aCNViewer_DATA/snpArrays250k_sty/expectedResults/test2/GSE9845_lrr_baf.segments_merged_hist_2000000nt.png`
![quantitative stacked histogram example:](/img/GSE9845_lrr_baf.segments_merged_hist_2000000nt.png?raw=true "Quantitative stacked histogram example")


==**Here is the full command:**==

`python aCNViewer.py -f ASCAT_SEGMENT_FILE -c CHR_SIZE_FILE -t OUTPUT_DIR --histogram 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR [--lohToPlot LOH_TO_PLOT] [--rColorFile RCOLOR_FILE]`<br>
where:
* <a id="ascatSegmentFile"></a>`ASCAT_SEGMENT_FILE`: ASCAT segment file (`ascat.output$segments` obtained by running `ascat.runAscat`) with the following columns:
  + `sample`
  + `chr`
  + `startpos`
  + `endpos`
  + `nMajor`
  + `nMinor`
* <a href="#chrSize">`CHR_SIZE_FILE`</a>
* <a href="#centromereFile">`CENTROMERE_FILE`</a>
* <a href="#windowSize">`WINDOW_SIZE`</a>
* <a href="#lohToPlot">`LOH_TO_PLOT`</a>
* <a href="#rColorFile">`RCOLOR_FILE`</a>


#### TestSequenzaCNVs
Generate quantitative stacked histogram from Sequenza results with a window size of 2Mbp:<br>

[aCNViewer_DATA.tar.gz](http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz) is required to run this example.

`python aCNViewer.py -f aCNViewer_DATA/wes/ --fileType Sequenza -c aCNViewer_DATA/wes/hg19.chrom.sizes -t TARGET_DIR --histogram 1 -C aCNViewer_DATA/wes/centro_hg19.txt -w 2000000 -b BIN_DIR`

The histogram `TARGET_DIR/ascat_merged_hist_2000000nt.png` can be found in `aCNViewer_DATA/wes/expectedResults/ascat_merged_hist_2000000nt.png`.


==**Here is the full command:**==

`python aCNViewer.py -f SEQUENZA_RES_DIR --fileType Sequenza -c CHR_SIZE_FILE -t TARGET_DIR --histogram 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR [--lohToPlot LOH_TO_PLOT] [--rColorFile RCOLOR_FILE]`<br>
where:
* `SEQUENZA_RES_DIR` is the folder containing Sequenza results (`*_segments.txt`)
* <a href="#chrSize">`CHR_SIZE_FILE`</a>
* <a href="#centromereFile">`CENTROMERE_FILE`</a>
* <a href="#windowSize">`WINDOW_SIZE`</a>
* <a href="#lohToPlot">`LOH_TO_PLOT`</a>
* <a href="#rColorFile">`RCOLOR_FILE`</a>


#### OtherCNVformats

CNV results from any software can be processed by aCNViewer if formatted in the ASCAT segment file format i.e. a tab-delimited file with the following columns: 
* `sample`
* `chr`
* `startpos`
* `endpos`
* `nMajor`
* `nMinor`

If there is only a global CNV value `v` (and this no allele-specific CNV value), `nMajor` and `nMinor` should take the value of `v / 2`.



### Other plots
All the different examples below **require a sample file** with at least one phenotypic / clinical information in order separate samples according to a given phenotypic / clinical group. All the previous commands can be adjusted to plot the chosen graph below by appending to each command the sample file and the chosen phenotypic / clinical column name. All the examples below will start from ASCAT segment file but could start from any previously described input. [aCNViewer_DATA.tar.gz](http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz) is required to run all the following examples.

#### PlotHeatmaps

<u>Heatmap with relative copy number values:</u>
`python aCNViewer.py -f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --heatmap 1 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b BIN_DIR --chrLegendPos 0,.55 --groupLegendPos .9,1.05 --useRelativeCopyNbForClustering 1`

<u>Heatmap with regions ordered by genomic positions (only clustering on samples):</u>
`python aCNViewer.py -f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --heatmap 1 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b BIN_DIR --chrLegendPos 0,.55 --groupLegendPos .9,1.05 --useRelativeCopyNbForClustering 1 --keepGenomicPosForHistogram 1`

<u>Heatmap with copy number values:</u>
`python aCNViewer.py -f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --heatmap 1 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b BIN_DIR --chrLegendPos 0,.55 --groupLegendPos .9,1.05`

See [Heatmap example](/img/matrix_None2000000nt_heatmap_BCLC_stage_None.pdf)


==**Full command:**==

`python aCNViewer.py -f ASCAT_SEGMENT_FILE -c CHR_SIZE_FILE -t OUTPUT_DIR --heatmap 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR --sampleFile SAMPLE_FILE -G PHENOTYPIC_COLUMN_NAME [--labRow LAB_ROW] [--labCol LAB_COL] [--cexCol CEX_COL] [--cexRow CEX_ROW] [--height HEIGHT] [--width WIDTH] [--margins MARGINS] [--hclust HCLUST] [--groupLegendPos GROUP_LEGEND_POS] [--chrLegendPos CHR_LEGEND_POS] [--rColorFile RCOLOR_FILE] [--useRelativeCopyNbForClustering USE_RELATIVE_COPY_NB_FOR_CLUSTERING]`<br>
where:
* <a href="#ascatSegmentFile">`ASCAT_SEGMENT_FILE`</a>
* <a href="#chrSize">`CHR_SIZE_FILE`</a>
* <a href="#centromereFile">`CENTROMERE_FILE`
* <a href="#windowSize">`WINDOW_SIZE`
* <a id="sampleFile"></a>`SAMPLE_FILE`: a tab-delimited file that should contain at least a column `Sample` with the name of each sample and another column with the phenotypic / clinical feature. This file can contain a `sampleAlias` which will be used as the official sample id if provided.
* <a id="phenotypicColumnName"></a>`PHENOTYPIC_COLUMN_NAME` refers to the name of the column of the phenotypic / clinical feature of interest in `SAMPLE_FILE`. If you omit this parameter, one plot per feature defined in `SAMPLE_FILE` will be generated
* `LAB_ROW` is an optional parameter telling whether heatmap's row names (chromosomal regions) should be shown. The default value is `0`
* `LAB_COL` is an optional parameter telling whether heatmap's column names (sample names) should be shown. The default value is `1`
* `CEX_COL` is an optional parameter setting `cexCol` for heatmaps. The default value is `0.7`. See R heatmap.2 documentation for more details
* `CEX_ROW` is an optional parameter setting `cexRow` for heatmaps. The default value is `0.45`. See R heatmap.2 documentation for more details
* `HEIGHT` is an optional parameter setting `height` for heatmaps. The default value is `12`. for heatmaps.
* `WIDTH` is an optional parameter setting `width` for heatmaps. The default value is `10`. See R heatmap.2 documentation for more details
* `MARGINS` is an optional parameter setting `margins` as a comma-separated string for heatmaps. The default value is `5,5`. See R heatmap.2 documentation for more details
* <a id="hclust"></a>`HCLUST` is an optional parameter setting `hclust` method for heatmaps / dendrograms. See R heatmap.2 documentation for more details
* `GROUP_LEGEND_POS` is an optional parameter setting the phenotypic / clinical feature legend's position within the heatmap. The default value is `topright` and can be changed to coordinates (for example `0.1,0.5`) or in R specified logical location (`top`, `bottom`, `left`, `right`, etc)
* `CHR_LEGEND_POS` is an optional parameter setting the chromosome legend's position within the heatmap. The default value is `bottomleft` and can be changed to coordinates (for example `0.1,0.5`) or in R specified logical location (`top`, `bottom`, `left`, `right`, etc)
* <a href="#rColorFile">`RCOLOR_FILE`</a>
* `USE_RELATIVE_COPY_NB_FOR_CLUSTERING` is an optional parameter specifying whether the CNV matrix used for the heatmap should be relative copy number values or not. The default value is `0`


#### PlotDendrograms

`python aCNViewer.py -f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --dendrogram 1 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b BIN_DIR -u 1`

![dendrogram example:](/img/matrix_None2000000nt_dendro_BCLC_stage.png?raw=true "Dendrogram example")


==**Full command:**==

`python aCNViewer.py -f ASCAT_SEGMENT_FILE -c CHR_SIZE_FILE -t OUTPUT_DIR --dendrogram 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR --sampleFile SAMPLE_FILE -G PHENOTYPIC_COLUMN_NAME [-u USE_SHAPE] [--hclust HCLUST] [--rColorFile RCOLOR_FILE]`<br>
where:
* <a href="#ascatSegmentFile">`ASCAT_SEGMENT_FILE`</a>
* <a href="#chrSize">`CHR_SIZE_FILE`</a>
* <a href="#centromereFile">`CENTROMERE_FILE`
* <a href="#windowSize">`WINDOW_SIZE`
* <a href="#sampleFile">`SAMPLE_FILE`</a>
* <a href="#phenotypicColumnName">`PHENOTYPIC_COLUMN_NAME`</a>
* `USE_SHAPE` is optional and tells whether colored shaped leaves should be used instead of sample names. The default value is `0`
* <a href="#hclust">`HCLUST`</a>
* <a href="#rColorFile">`RCOLOR_FILE`</a>


#### PlotAll
If you want to plot all available graphical representations (a quantitative stacked histogram, a heatmap and a dendrogram):<br>
`python aCNViewer.py -f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --dendrogram 1 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b BIN_DIR -u 1 --plotAll 1`


==**Full command:**==

`python aCNViewer.py -f ASCAT_SEGMENT_FILE -c CHR_SIZE_FILE -t OUTPUT_DIR --plotAll 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b BIN_DIR --sampleFile SAMPLE_FILE -G PHENOTYPIC_COLUMN_NAME [OPTIONS]`<br>
where:
* <a href="#ascatSegmentFile">`ASCAT_SEGMENT_FILE`</a>
* <a href="#chrSize">`CHR_SIZE_FILE`</a>
* <a href="#centromereFile">`CENTROMERE_FILE`
* <a href="#windowSize">`WINDOW_SIZE`
* <a href="#sampleFile">`SAMPLE_FILE`</a>
* <a href="#phenotypicColumnName">`PHENOTYPIC_COLUMN_NAME`</a>
* `OPTIONS` refers to all the options defined for heatmaps / dendrograms and quantitative stacked histograms


### OutputFiles

#### ASCAT

When processing raw SNP array data with aCNViewer, ASCAT is used to calculate CNV profiles. These results are saved into a folder named `ASCAT` in the user selected target directory with the following files:
* `*.segments.txt`: file containing ASCAT predicted CNV segments
* `*.ascatInfo.txt`: file containing the following ASCAT values for all the samples: `aberrantcellfraction`, `goodnessOfFit`, `psi` and `ploidy`
* `*.png`: the various ASCAT graphical outputs:

| File | Description |
| --- | --- |
| .ASCATprofile.png | genome-wide representation of ASCAT CNVs |
| .ASPCF.png | results of segmentation using Allele-Specific Piecewise Constant Fitting |
| .rawprofile.png | genome-wide representation of raw ASCAT CNVs |
| .sunrise.png | sunrise plot showing the optimal solution of tumor ploidy and percentage of aberrant tumor |
| .tumour.png | representation of LogR and BAF values |
| tumorSep*.png	 |  |
| .ascatInfo.txt | ASCAT values of aberrantcellfraction, goodnessOfFit, psi and ploidy for all samples |
| .segments.txt | list of all CNVs with the copy number for each allele |


#### HistogramOutputs

When generating histograms, 3 text files with the suffix `_samples.txt` will be created along:
* one with all the genomic segments
* one with only the LOH events (file with suffix `_loh_samples.txt`)
* one with only the cn-LOH (file with suffix `_cnLoh_samples.txt`)

Each file is in the same format with the following columns:
* `CNV key`: the relative copy number value compared to the tumor ploidy
* `chrName`
* `start`
* `segmentLength`
* `percentage`: the percentage of samples with the relative ploidy value in `CNV key` for the segment (`chrName`, `start`, `start + segmentLength`)
* `samples`: the list of the samples falling in the above category


#### Dendrograms and heatmaps

2 folders (`relCopyNb` and `rawCopyNb`) will be created and will respectively contain graphs generated from relative copy number values and raw copy number values.
