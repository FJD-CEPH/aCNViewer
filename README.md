# aCNViewer
comprehensive genome-wide visualization of absolute copy number and copy neutral variations

**Contact:** Victor Renault / Alexandre How-Kit (aCNViewer@cephb.fr)

## Table of contents
- [Installation](#installation)
- [Overview](#overview)
- [Tutorial](#tutorial):
    * [Glossary](#glossary)
    * Processing SNP array data:
      + [Affymetrix](#affymetrix)
        - [Test using Ascat results](#testaffyascat) (**start here if you want an overview of all the plotting options and different usage scenarios**)
        - [Test using Affymetrix Cel files](#testaffycel)
      + [Illumina](#illumina)
    * [Processing NGS data](#ngs)
      + [from paired (tumor / normal) bams](#testsequenzaraw)
      + [from Sequenza segment files](#testsequenzacnvs)
    * [Processing CNV file](#processing-cnv-files)
      + [from ASCAT](#testaffyascat)
      + [from Sequenza](#testsequenzacnvs)
      + [from other tools](#othercnvformats)
- [Output files](#outputfiles)
    * [Dendrograms and heatmaps](#dendrograms-and-heatmaps)
  

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
  + [Sequenza](https://cran.r-project.org/web/packages/sequenza/index.html) and [samtools](https://sourceforge.net/projects/samtools/files/latest/download?source=files) if you are analyzing paired (tumor / normal) bams
  + [plotrix](https://cran.r-project.org/web/packages/plotrix/index.html) for plotting dendrograms (will be automatically installed if not already installed)
  + [gplots](https://cran.r-project.org/web/packages/gplots/index.html)
  + [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

* [tQN](http://cbbp.thep.lu.se/~markus/software/tQN/tQN-1.1.2.zip) if you plan to process raw Illumina SNP arrays (to uncompress into <a href="#binDir">`BIN_DIR`</a>) and run tQN normalisation. If the cluster file for the Illumina SNP array you plan to analyze is not in the tQN lib folder, you can download additional cluster files from [here](http://cbbp.thep.lu.se/~markus/software/tQN/)

* [GISTIC](http://portals.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=216&p=t) if you want to have an advanced statistical way to prioritize regions of interest. Create a folder named `GISTIC_VERSION` in `BIN_DIR` and uncompress the GISTIC archive into it. Follow the instructions listed in `INSTALL.txt` at the root of the GISTIC folder in order to install MATLAB Component Runtime required by GISTIC and set the associated environment variables (`LD_LIBRARY_PATH` and `XAPPLRESDIR`).

* Python with version &ge; 2.7


## Overview:
![Overview of aCNViewer:](/img/aCNViewer.png?raw=true "Overview of aCNViewer")  



***


## Tutorial

For all the examples below, there is a folder named `expectedResults` in the folder associated with the `-f` option.


### Glossary:

Let's call:
- `aCNViewer_DATA` the location where the test data set [aCNViewer_DATA.tar.gz](https://www.cng.fr/genodata/pub/LIVER/aCNViewer_DATA.tar.gz) has been uncompressed into
- <a id="binDir">`BIN_DIR`</a> the folder containing all third-party softwares located in `aCNViewer_DATA/bin`.
- <a id="dockerOrPython">`DOCKER_OR_PYTHON`</a> refers to the fact that `docker run fjdceph/acnviewer` or `python aCNViewer/code/aCNViewer.py` can be used as a prefix to run aCNViewer depending on the chosen installation method.


### Requirements:

Download the test data set [aCNViewer_DATA.tar.gz](https://www.cng.fr/genodata/pub/LIVER/aCNViewer_DATA.tar.gz).


### Processing SNP array data

#### Affymetrix

##### TestAffyAscat
Generate all available plots from ASCAT segment files with a window size of 2Mbp:<br>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt`

The generated histogram `OUTPUT_DIR/GSE9845_lrr_baf.segments_merged_hist_2000000nt.png` can be found in `aCNViewer_DATA/snpArrays250k_sty/expectedResults/test2/GSE9845_lrr_baf.segments_merged_hist_2000000nt.png`
![quantitative stacked histogram example:](/img/GSE9845_lrr_baf.segments_merged_hist_2000000nt.png?raw=true "Quantitative stacked histogram example")

Here are other typical plots you may be interested in:

<u>Quantitative histogram with GISTIC results:</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 -b aCNViewer_DATA/bin --runGISTIC 1 --refBuild hg18`


<u>Heatmap with relative copy number values:</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --PLOT_ALL 0 --heatmap 1 --dendrogram 0 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b aCNViewer_DATA/bin --chrLegendPos 0,.55 --groupLegendPos .9,1.05 --useRelativeCopyNbForClustering 1`


<u>Heatmap with regions ordered by genomic positions (only clustering on samples):</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --PLOT_ALL 0 --heatmap 1 --dendrogram 0 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b aCNViewer_DATA/bin --chrLegendPos 0,.55 --groupLegendPos .9,1.05 --useRelativeCopyNbForClustering 1 --keepGenomicPosForHistogram 1`


[u]Heatmap with copy number values:[/u]

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --PLOT_ALL 0 --heatmap 1 --dendrogram 0 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b aCNViewer_DATA/bin --chrLegendPos 0,.55 --groupLegendPos .9,1.05`

See [Heatmap example](/img/matrix_None2000000nt_heatmap_BCLC_stage_None.pdf)


<u>Heatmap with copy number values:</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR --dendrogram 1 -G "BCLC stage" -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt -b aCNViewer_DATA/bin -u 1`

![dendrogram example:](/img/matrix_None2000000nt_dendro_BCLC_stage.png?raw=true "Dendrogram example")



==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f ASCAT_SEGMENT_FILE -c CHR_SIZE_FILE -t OUTPUT_DIR -C CENTROMERE_FILE -w WINDOW_SIZE -b` <a href="#binDir">`BIN_DIR`</a> `[--lohToPlot LOH_TO_PLOT] [--histogram HISTOGRAM] [--sampleFile SAMPLE_FILE -G PHENOTYPIC_COLUMN_NAME --rColorFile RCOLOR_FILE --plotAll PLOT_ALL --outputFormat OUTPUT_FORMAT] [--heatmap HEATMAP --labRow LAB_ROW --labCol LAB_COL --cexCol CEX_COL --cexRow CEX_ROW --height HEIGHT --width WIDTH --margins MARGINS --hclust HCLUST --groupLegendPos GROUP_LEGEND_POS --chrLegendPos CHR_LEGEND_POS --useRelativeCopyNbForClustering USE_RELATIVE_COPY_NB_FOR_CLUSTERING --keepGenomicPosForHistogram KEEP_GENOMIC_POS] [--dendrogram DENDROGRAM --useShape USE_SHAPE] [--runGISTIC RUN_GISTIC --refBuild REF_BUILD --geneGistic GENE_GISTIC --smallMem SMALL_MEM --broad BROAD --brLen BR_LEN --conf CONF --armPeel ARM_PEEL --saveGene SAVE_GENE --gcm GCM]`<br>
where:
* <a id="ascatSegmentFile"></a>`ASCAT_SEGMENT_FILE`: ASCAT segment file (`ascat.output$segments` obtained by running `ascat.runAscat`) with the following columns:
  + `sample`
  + `chr`
  + `startpos`
  + `endpos`
  + `nMajor`
  + `nMinor`
* <a id="chrSize"></a>`CHR_SIZE_FILE`: a tab-delimited file with 2 columns respectively chromosome name and chromosome length
* <a id="centromereFile"></a>`CENTROMERE_FILE`: file giving the centromere bounds. Can be generated using `curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/BUILD/database/cytoBand.txt.gz" | gunzip -c | grep acen > centro_build.txt`
* <a id="windowSize"></a>`WINDOW_SIZE`: segment size in bp. Please note that alternatively, `-p PERCENTAGE` can be used instead of `-w WINDOW_SIZE` in order to set the segment size in percentage of chromosome length where `PERCENTAGE` is a floating number between 0 and 100
* <a id="lohToPlot"></a>`LOH_TO_PLOT`: histogram option for LOH plotting. Values should be one of "cn-LOH" for plotting cn-LOH only, "LOH" for LOH only, "both" for cn-LOH and LOH or "none" to disable this feature. The default value is "cn-LOH".


<a id="generalPlotOptions"></a>**The following options are general plotting options:**
* <a id="sampleFile"></a>`SAMPLE_FILE`: a tab-delimited file that should contain a column named `Sample` with the name of each sample and at least another column with the phenotypic / clinical feature. This file can contain a `sampleAlias` which will be used as the official sample id if provided.
* <a id="phenotypicColumnName"></a>`PHENOTYPIC_COLUMN_NAME` is optional and refers to the name of the column of the phenotypic / clinical feature of interest in `SAMPLE_FILE`. If you omit this parameter, one plot per feature defined in `SAMPLE_FILE` will be generated.
* <a id="rColorFile"></a> `RCOLOR_FILE`: colors in histograms (section "[histogram]". If defined, should contain exactly 10 colors [one per line] corresponding to CNV values in the following order: "&le; -4", "-3", "-2", "-1", "1", "2", "3", "4", "5", "&ge; 6"), dendrograms (section "[group]". If defined, should contain at least the same number of colors than the number of distinct values for the phenotypic / clinical feature of interest) and heatmaps (sections "[chr]" [if defined, should contain 22 colors corresponding to chromosomes 1 to 22], "[group]" and "[heatmap]" [if defined, should contain 10 colors [one per line] corresponding to CNV values in the following order: "0", "1", "2", "3", "4", "5", "6", "7", "8", "&ge; 9"]) can be redefined in that file. An example can be found [here](/img/rColor.txt).
* `HISTOGRAM`: specify whether an histogram should be generated. The default value is 1.
* `PLOT_ALL`: specify whether all available plots should be generated. The default value is 1.
* `OUTPUT_FORMAT` 


<a id="gisticOptions"></a>**The following options are GISTIC options** (more details can be found [here](ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm)):
* `RUN_GISTIC`: specify whether to run GISTIC in order to have a statistical way to prioritize regions of interest. The default value is `0`.
* `REF_BUILD`: the genome build used to generate the CNV segments (allowed values are `hg16`, `hg17`, `hg18` and `hg19`)
* `GENE_GISTIC`: tell whether gene GISTIC algorithm should be used to calculate the significance of deletions at a gene level instead of a marker level. The default value is `1`.
* `SMALL_MEM`: tell GISTIC whether to use memory compression at the cost of a longer runtime. The default value is `0`.
* `BROAD`: tell GISTIC to run the broad-level analysis as well. The default value is `1`.
* `BR_LEN`: set GISTIC's `broad_len_cutoff`. The default value is `0.5`.
* `CONF`: set the confidence level used to calculate the region containing a driver. The default value is `0.9`.
* `ARM_PEEL`: set GISTIC's `arm_peeloff`. The default value is `1`.
* `SAVE_GENE`: tell GISTIC whether to save gene tables. The default value is `1`.
* `GCM`: set GISTIC's `gene_collapse_method`. The default value is `extreme`.


<a id="heatmapDendroOptions"></a>**The following options are mainly specific to heatmaps while a few are related to dendrograms:**

* `HEATMAP` is an optional parameter used only if `PLOT_ALL` is set to `0` to tell whether to plot heatmaps or not. The default value is `1`
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
* `USE_RELATIVE_COPY_NB_FOR_CLUSTERING` is an optional parameter specifying whether the CNV matrix used for the heatmap should be relative copy number values or not. The default value is `0`. If `PLOT_ALL` is `1` then plots for both values of `USE_RELATIVE_COPY_NB_FOR_CLUSTERING` will be generated.
* `KEEP_GENOMIC_POS` is optional and will keep the segmented genome in its original position if set to `1` and not cluster segments according to sample CNV patterns (the default value is `0`).
* `DENDROGRAM` is an optional dendrogram parameter used only if `PLOT_ALL` is set to `0` to tell whether to plot dendrograms or not. The default value is `1`
* `USE_SHAPE` is an optional dendrogram parameter and if set to `1` (default value) will replace sample labels with colored shapes in the leaves of the dendrogram(s).


##### TestAffyCel

**Generate a quantitative stacked histogram from CEL files (subset of [data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9845) of hepatocellular carcinomas with hepatitis C virus etiology used in Chiang et al. Cancer Res, 2008) with a window size of 2Mbp:**

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/ -c aCNViewer_DATA/snpArrays250k_sty/hg18.chrom.sizes -t OUTPUT_DIR -C aCNViewer_DATA/snpArrays250k_sty/centro.txt -w 2000000 -b aCNViewer_DATA/bin/ --platform Affy250k_sty -l aCNViewer_DATA/snpArrays250k_sty/LibFiles/ --gw6Dir aCNViewer_DATA/snpArrays250k_sty/gw6/`

If ASCAT is not installed (i.e you are not using the [docker](https://hub.docker.com/r/fjdceph/acnviewer/) application) and if you want to install it into a custom R library folder, please add the following option to the previous command line: `--rLibDir RLIB`.

The histogram `OUTPUT_DIR/lrr_baf.segments_merged_hist_2000000nt.png` can also be found in `aCNViewer_DATA/snpArrays250k_sty/expectedResults/test1/lrr_baf.segments_merged_hist_2000000nt.png`.


==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f CEL_DIR -c CHR_SIZE_FILE -t OUTPUT_DIR -C CENTROMERE_FILE -w WINDOW_SIZE -b` <a href="#binDir">`BIN_DIR`</a> `--platform AFFY_PLATFORM -l AFFY_LIB_DIR --gw6Dir GW6_DIR [--gcFile ASCAT_GC_FILE] [--lohToPlot LOH_TO_PLOT]` [[`GENERAL_PLOT_OPTIONS`](#generalPlotOptions)] [`[GISTIC_OPTIONS]`](#gisticOptions) [`[HEATMAP_DENDRO_OPTIONS]`](#heatmapDendroOptions)<br>
where:
* `CEL_DIR` is the folder containing ".cel" ou ".cel.gz" files
* <a href="#chrSize">`CHR_SIZE_FILE`</a>
* <a href="#windowSize">`WINDOW_SIZE`</a>
* <a href="#centromereFile">`CENTROMERE_FILE`</a>
* `AFFY_PLATFORM`: name of ASCAT supported Affymetrix platform with a GC content file available ("Affy250k_sty", "Affy250k_nsp", "Affy500k" or "AffySNP6"). Please refer to [ASCAT website](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) for more details
* `AFFY_LIB_DIR`: Affymetrix library file downloadable from [Affymetrix website](http://www.affymetrix.com/support/technical/byproduct.affx?cat=dnaarrays)
* `GW6_DIR` refers to the folder where [gw6.tar.gz](http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz) has been uncompressed into. This archive contains different programs and files necessary to process Affymetrix SNP array
* <a id="ascatGcFile">`ASCAT_GC_FILE`</a>: GC content file necessary for ASCAT GC correction when analyzing SNP array data. This parameter is optional as its value will be automatically deduced from the value of `AFFY_PLATFORM`. Please check [ASCAT website](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) for available GC content files
* <a href="#lohToPlot">`LOH_TO_PLOT`</a>


#### Illumina

##### TestIllu660k

Generate a quantitative stacked histogram from [raw Illumina data](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE47357&format=file&file=GSE47357%5FMatrix%5Fsignal%5F660w%2Etxt%2Egz) from non-Hodgkin lymphoma patients used in Yang F *et al*. PLoS One 2014 with a window size of 2Mbp:

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrayIllu660k/GSE47357_Matrix_signal_660w.txt.gz -c aCNViewer_DATA/wes/hg19.chrom.sizes -C aCNViewer_DATA/wes/centro_hg19.txt -w 2000000 -b aCNViewer_DATA/bin/ --probeFile aCNViewer_DATA/snpArrayIllu660k/Human660W-Quad_v1_H_SNPlist.txt --platform Illumina660k -t ILL_OUT --beadchip "human660w-quad"`

==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f ILLU_FILES -c CHR_SIZE_FILE -C CENTROMERE_FILE -w WINDOW_SIZE -b` <a href="#binDir">`BIN_DIR`</a> `[--sampleList SAMPLE_TO_PROCESS_FILE] --probeFile PROBE_POS_FILE --platform ILLUMINA_PLATFORM [--beadchip BEADCHIP] [-g ASCAT_GC_FILE] [--lohToPlot LOH_TO_PLOT]` [[`GENERAL_PLOT_OPTIONS`](#generalPlotOptions)] [`[GISTIC_OPTIONS]`](#gisticOptions) [`[HEATMAP_DENDRO_OPTIONS]`](#heatmapDendroOptions)<br>
where:
  * `ILLU_FILES` can either be the list of Illumina final report files to process specified either as a comma-separated string with all the report files to process or as a directory containing these files. Each Illumina final report file should contain at least the following columns:
    - `SNP Name`
    - `Sample ID`
    - `Log R Ratio`
    - `B Allele Freq`
  
    Alternatively, it can be the raw Illumina files with at least the following columns:
    - `ID`
    - `SAMPLE1.X`
    - `SAMPLE1.Y`
    - ...
    - `SAMPLEn.X`
    - `SAMPLEn.Y`
    
  * <a href="#chrSize">`CHR_SIZE_FILE`</a>
  * <a href="#centromereFile">`CENTROMERE_FILE`</a>
  * <a href="#windowSize">`WINDOW_SIZE`</a>
  * `PROBE_POS_FILE`: file listing the probes used on the SNP array with their genomic position. The file is tab-delimited with the following columns:
    - `Name`
    - `Chr`
    - `MapInfo`
  * `ILLUMINA_PLATFORM`: name of ASCAT supported Illumina platform with a GC content file available ("Illumina660k" or "HumanOmniExpress12"). Please refer to [ASCAT website](//www.crick.ac.uk/peter-van-loo/software/ASCAT) for more details
  * <a href="#ascatGcFile">`ASCAT_GC_FILE`</a>
  * `SAMPLE_TO_PROCESS_FILE`: optional, used to specify list of samples to process in one of the following formats:
    - a comma-separated string listing all the samples to process
    - the name of text file with one line per sample to process
    - the name of a Python dump file with the extension ".pyDump"
  * <a href="#lohToPlot">`LOH_TO_PLOT`</a>


#### NGS

Sequenza is used to process NGS **paired (tumor / normal) bams** and produce CNV segments. These segments are then used by aCNViewer to produce the different available outputs. This step is best executed on a computer cluster (supported clusters are <a id="supportedClusters">SGE, SLURM, MOAB and LSF</a>. Tests have been successfully made on SGE and SLURM clusters) but will work on a single machine as well (although it will be much slower).

##### testSequenzaRaw

Generate a quantitative histogram from paired (tumor / normal) bams (if you run this example on a <a href="#supportedClusters">supported computer cluster</a> which do not support array jobs then you will need to run the command below a second time after the first set of jobs are finished in order to generate the different plots you are interested in):

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/wes/bams/ -c aCNViewer_DATA/wes/hg19.chrom.sizes -t WES_OUT -C aCNViewer_DATA/wes/centro_hg19.txt -w 2000000 -b aCNViewer_DATA/bin/ --fileType Sequenza --samplePairFile aCNViewer_DATA/wes/bams/sampleFile.txt -r aCNViewer_DATA/wes/Homo_sapiens_assembly19.fasta`

==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f BAM_DIR -c CHR_SIZE_FILE -t OUTPUT_DIR -C CENTROMERE_FILE -w WINDOW_SIZE -b` <a href="#binDir">`BIN_DIR`</a> `--fileType Sequenza --samplePairFile SAMPLE_PAIR_FILE -r REF_FILE [--byChr 1] [-n NB_THREADS] [--createMpileUp CREATE_MPILEUP] [--pattern BAM_FILE_PATTERN] [-M MEMORY]` [[`GENERAL_PLOT_OPTIONS`](#generalPlotOptions)] [`[GISTIC_OPTIONS]`](#gisticOptions) [`[HEATMAP_DENDRO_OPTIONS]`](#heatmapDendroOptions)<br>
where:
  * `BAM_DIR` is the folder containing the paired bam files
  * `BAM_FILE_PATTERN` is an optional parameter which default value is `.bam`
  * <a href="#chrSize">`CHR_SIZE_FILE`</a>
  * <a href="#centromereFile">`CENTROMERE_FILE`</a>
  * <a href="#windowSize">`WINDOW_SIZE`</a>
  * `SAMPLE_PAIR_FILE` is a tab-delimited file with the following three column names:
    - `idvdName`
    - `sampleName`
    - `type` which should either be `T` for tumoral samples or `N` for normal samples
  * `REF_FILE` is the reference file in fasta format used to generate the bam files
  * `BY_CHR` is an optional parameter to indicate whether Sequenza should create seqz (Sequenza intermediate file) files by chromosome or not (the default value is `1`)
  * `NB_CPUS` is an optional parameter specifying the number of cores which will be used for each sample pair to create chromosomal seqz files if `BY_CHR` has been set to `1`. If aCNViewer is ran on a <a href="#supportedClusters">supported computer cluster</a> master node, jobs will be submitted to the cluster. Otherwise, multi-threading will be used run Sequenza.
  * `CREATE_MPILEUP` is an optional parameter telling Sequenza whether to create intermediate mpileup files when generating results. The default value is `1` and it is recommended not to change its value as Sequenza may freeze in some cases when set to `0`.
  * `MEMORY`: optional argument specifying allocated memory in GB to run Sequenza when using a computer cluster. The default value is 8 (GB) and should work for most WES analysis


#### TestSequenzaCNVs
Generate quantitative stacked histogram from Sequenza results with a window size of 2Mbp:<br>

[aCNViewer_DATA.tar.gz](http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz) is required to run this example.

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/wes/ --fileType Sequenza -c aCNViewer_DATA/wes/hg19.chrom.sizes -t TARGET_DIR --histogram 1 -C aCNViewer_DATA/wes/centro_hg19.txt -w 2000000 -b aCNViewer_DATA/bin`

The histogram `TARGET_DIR/ascat_merged_hist_2000000nt.png` can be found in `aCNViewer_DATA/wes/expectedResults/ascat_merged_hist_2000000nt.png`.


==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f SEQUENZA_RES_DIR --fileType Sequenza -c CHR_SIZE_FILE -t TARGET_DIR --histogram 1 -C CENTROMERE_FILE -w WINDOW_SIZE -b` <a href="#binDir">`BIN_DIR`</a> `[--lohToPlot LOH_TO_PLOT] [--rColorFile RCOLOR_FILE]`<br>
where:
* `SEQUENZA_RES_DIR` is the folder containing Sequenza results (`*_segments.txt`)
* <a href="#chrSize">`CHR_SIZE_FILE`</a>
* <a href="#centromereFile">`CENTROMERE_FILE`</a>
* <a href="#windowSize">`WINDOW_SIZE`</a>
* <a href="#lohToPlot">`LOH_TO_PLOT`</a>
* <a href="#rColorFile">`RCOLOR_FILE`</a>
  
  
### Processing CNV files

At the moment, ASCAT segment file and Sequenza results can be used as an input to aCNViewer. It is possible however to feed aCNViewer with CNV results from any other softwares as explained in the section [below] (#OtherCNVformats).


Both examples below require to download [aCNViewer_DATA.tar.gz](http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz).


#### [Ascat file generated from Affymetrix SNP arrays](#testaffyascat)

#### [Sequenza segments](#testsequenzacnvs)

#### OtherCNVformats

CNV results from any software can be processed by aCNViewer if formatted in the ASCAT segment file format i.e. a tab-delimited file with the following columns: 
* `sample`
* `chr`
* `startpos`
* `endpos`
* `nMajor`
* `nMinor`

If there is only a global CNV value `v` (and this no allele-specific CNV value), `nMajor` and `nMinor` can take any value as long as `nMajor + nMinor = v`. When plotting the quantitative histogram, add option `--lohToPlot none` to disable LOH plotting.


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


The following files are created as well:
* `*_10pc_ploidy.txt` is a matrix of segments of 10% chromosomal length for all samples. The last column indicates the calculated ploidy which corresponds to the most frequent ploidy

* `*.R` is the R script used to create the various graphical representations. You can modify and re-run this script if you want to further customize your graphical outputs and aCNViewer do not propose the customizations you are looking for.


##### Dendrograms and heatmaps

2 folders (`relCopyNb` and `rawCopyNb`) will be created and will respectively contain graphs generated from relative copy number values and raw copy number values.
