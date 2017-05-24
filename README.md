# aCNViewer
comprehensive genome-wide visualization of absolute copy number and copy neutral variations

**Contact:** Victor Renault / Alexandre How-Kit (aCNViewer@cephb.fr)

aCNViewer (Absolute CNV Viewer) is a tool which allows the visualization of absolute CNVs and cn-LOH across a group of cancer samples. aCNViewer proposes three graphical representations : dendrograms, bidimensional heatmaps allowing the visualization of chromosomal regions sharing similar abnormality patterns and quantitative stacked histograms facilitating the identification of recurrent absolute CNVs and cn-LOH. aCNViewer include a complete pipeline allowing the processing of raw data from SNP array and whole exome/genome sequencing experiments using respectively ASCAT and Sequenza algorithms to generate absolute CNV and cn-LOH data used for the graphical outputs.


## Table of contents
- [Installation](#installation)
- [Overview](#overview)
- [Tutorial](#tutorial):
    * [Glossary](#glossary)
    * Processing SNP array data:
      + [Affymetrix](#affymetrix)
        - [Test using Ascat results](#testaffyascat) (**start here if you want an overview of all the plotting options and different usage scenarios**):
          * [all plots](#allPlots)
          * [use custom colors for plots](#customColors)
          * [GISTIC example](#gisticExample)
          * [heatmap using relative copy number values](#heatmapRel)
          * [heatmap with CNVs ordered by genomic position](heatmapGenPos)
          * [heatmap using copy number values](#heatmapCNV)
          * [dendrogram example](#dendroExample)
          * [custom graphical output formats](#outputFormatExamples)
        - [Test using Affymetrix Cel files](#testaffycel)
      + [Illumina](#illumina)
    * [Processing NGS data](#ngs)
      + [from paired (tumor / normal) bams](#testsequenzaraw)
      + [from Sequenza segment files](#testsequenzacnvs)
    * [Processing CNV file](#processing-cnv-files)
      + [from ASCAT](#testaffyascat)
      + [from PennCNV](#penncnv)
      + [from Sequenza](#testsequenzacnvs)
      + [from other tools](#othercnvformats)
- [Output files](#outputfiles)
    * [ASCAT](#ascat)
    * [GISTIC](#gistic)
    * [Sequenza](#sequenza)
    * [Histogram](HistogramOutputs)
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

* a recent version of R (version &ge; 3.2) with [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) installed for generating the different graphical outputs:
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

The results of all the examples below can be found in `aCNViewer_DATA/allTests` in their respective target folder. All examples of this tutorial can be run at once using: <a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-P testAll -t TARGET_DIR [--fastTest 0 --smallMem 0 --runGistic 1]`. 

If `--fastTest` is set to `1`, only tests which run in a *reasonable* amount of time will be run (all tests except [Illumina SNP array](#illumina), [paired bams with Sequenza](#testsequenzaraw), [GISTIC](#gistic) and [Affymetrix SNP arrays from CEL files](#testaffycel)). If `--runGistic` is `1`, GISTIC will be tested and if `--smallMem` is set `1`, GISTIC will run in small memory mode and will only require about 5GB of RAM vs 50GB of RAM at the expense of a longer running time.


### Glossary:

Let's call:
- `aCNViewer_DATA` the location where the test data set [aCNViewer_DATA.tar.gz](https://www.cng.fr/genodata/pub/LIVER/aCNViewer_DATA.tar.gz) has been uncompressed into
- <a id="binDir">`BIN_DIR`</a> the folder containing all third-party softwares located in `aCNViewer_DATA/bin`.
- <a id="dockerOrPython">`DOCKER_OR_PYTHON`</a> refers to the fact that `docker run fjdceph/acnviewer` or `python aCNViewer/code/aCNViewer.py` can be used as a prefix to run aCNViewer depending on the chosen installation method.


### Requirements:

Download the test data set [aCNViewer_DATA.tar.gz (~5GB and ~20GB uncompressed)](https://www.cng.fr/genodata/pub/LIVER/aCNViewer_DATA.tar.gz). In terms of computing resources: if you plan to:
- run Sequenza on paired bam files, an access to a computer cluster is highly recommended as even though aCNViewer will be able to process your data in multi-threading mode, it may take quite a long time depending on the number of sample pairs to analyze
- run GISTIC in order to have a robust statistical way to prioritize recurrent regions of interest, a machine with at least 50GB of RAM is necessary with `--smallMem 0` and 5GB with `--smallMem 1` (this option will make GISTIC run substantially longer)



### Processing SNP array data

#### Affymetrix

##### TestAffyAscat
<a id="allPlots"></a>Generate all available plots from ASCAT segment files using base resolution for the quantitative histograms and using a window size of 2Mbp for the other plots:<br>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt`

![quantitative stacked histogram example:](/img/GSE9845_lrr_baf.segments_hg18_cov_hist.png?raw=true "Quantitative stacked histogram example")

![Histogram of heterozygous / homozygous CNVs:](/img/GSE9845_lrr_baf.segments_hg18_cov_hist_hetHom.png?raw=true "Histogram of heterozygous / homozygous CNVs")

Here are other typical plots you may be interested in:

<a id="customColors"></a><u>Customize colors:</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_RCOLOR --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt --rColorFile aCNViewer_DATA/rColor.txt`


<a id="gisticExample"></a><u>Quantitative histogram with GISTIC results:</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_GISTIC --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --runGISTIC 1`


<a id="heatmapRel"></a><u>Heatmap of relative copy number values only for the clinical feature `BCLC stage` with the chromosome legend position set at `0,.55` i.e. at the left-most of the graph and at 55% on the y axis and the group legend position set at `.9,1.05` (basically at the top right corner):</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_HEATMAP1 --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt --plotAll 0 --heatmap 1 --dendrogram 0 -G "BCLC stage" --chrLegendPos 0,.55 --groupLegendPos .9,1.05 --useRelativeCopyNbForClustering 1`

![Heatmap of relative copy number values using the clinical feature `BCLC stage`:](/img/matrix_relCopyNb_2000000nt_heatmap_BCLC_stage_ward_10.png?raw=true "Heatmap of relative copy number values using the clinical feature `BCLC stage`")



<a id="heatmapGenPos"></a><u>Heatmap with regions ordered by genomic positions (only clustering on samples):</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_HEATMAP_GENPOS --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt --plotAll 0 --heatmap 1 --dendrogram 0 -G "BCLC stage" --chrLegendPos 0,.55 --groupLegendPos .9,1.05 --useRelativeCopyNbForClustering 1 --keepGenomicPosForHistogram 1`

![Heatmap of relative copy number values with regions ordered by genomic positions using the clinical feature `BCLC stage`:](/img/matrix_relCopyNb_2000000nt_heatmap_BCLC_stage_ward_11.png?raw=true "Heatmap of relative copy number values with regions ordered by genomic positions using the clinical feature `BCLC stage`")


<a id="heatmapCNV"></a><u>Heatmap with copy number values:</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_HEATMAP2 --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt --plotAll 0 --heatmap 1 --dendrogram 0 -G "BCLC stage" --chrLegendPos 0,.55 --groupLegendPos .9,1.05`

![Heatmap of copy number values using the clinical feature `BCLC stage`:](/img/matrix_rawCopyNb_2000000nt_heatmap_BCLC_stage_ward_00.png?raw=true "Heatmap of copy number values using the clinical feature `BCLC stage`")


<a id="dendroExample"></a><u>Dendrogram with copy number values:</u>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_DENDRO --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt --plotAll 0 --heatmap 0 --dendrogram 1 -G "BCLC stage" -u 1`

![Dendrogram of copy number values using the clinical feature `BCLC stage`:](/img/matrix_rawCopyNb_2000000nt_dendro_BCLC_stage_rawCopyNb.png?raw=true "Dendrogram of copy number values using the clinical feature `BCLC stage`")


<a id="outputFormatExamples"></a><u>Customize output formats:</u>

- all outputs set to pdf: <a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_PDF --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt --outputFormat pdf`

- all output set to jpg: <a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_PDF --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt --outputFormat jpg`

- heatmaps set to `bmp`, histograms to `tiff` and dendrograms to `pdf` with the R plot parameters `width=10,height=8`: `-f aCNViewer_DATA/snpArrays250k_sty/GSE9845_lrr_baf.segments.txt -t TEST_AFFY_OTHER_OUT --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --sampleFile aCNViewer_DATA/snpArrays250k_sty/GSE9845_clinical_info2.txt --outputFormat "heat:bmp;hist:tiff;dend:pdf(width=10,height=8)"`


==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f ASCAT_SEGMENT_FILE --refBuild REF_BUILD -b` <a href="#binDir">`BIN_DIR`</a> `[--histogram HISTOGRAM --lohToPlot LOH_TO_PLOT --useFullResolutionForHist USE_FULL_RESOLUTION_FOR_HIST] [-c CHR_SIZE_FILE -t OUTPUT_DIR -C CENTROMERE_FILE -w WINDOW_SIZE --sampleFile SAMPLE_FILE -G PHENOTYPIC_COLUMN_NAME --rColorFile RCOLOR_FILE --plotAll PLOT_ALL --outputFormat OUTPUT_FORMAT --ploidyFile PLOIDY_FILE] [--heatmap HEATMAP --labRow LAB_ROW --labCol LAB_COL --cexCol CEX_COL --cexRow CEX_ROW --height HEIGHT --width WIDTH --margins MARGINS --hclust HCLUST --groupLegendPos GROUP_LEGEND_POS --chrLegendPos CHR_LEGEND_POS --useRelativeCopyNbForClustering USE_RELATIVE_COPY_NB_FOR_CLUSTERING --keepGenomicPosForHistogram KEEP_GENOMIC_POS] [--dendrogram DENDROGRAM --useShape USE_SHAPE] [--runGISTIC RUN_GISTIC --geneGistic GENE_GISTIC --smallMem SMALL_MEM --broad BROAD --brLen BR_LEN --conf CONF --armPeel ARM_PEEL --saveGene SAVE_GENE --gcm GCM]`<br>
where:
* <a id="ascatSegmentFile"></a>`ASCAT_SEGMENT_FILE`: ASCAT segment file (`ascat.output$segments` obtained by running `ascat.runAscat`) with the following columns:
  + `sample`
  + `chr`
  + `startpos`
  + `endpos`
  + `nMajor`
  + `nMinor`
* `REF_BUILD`: the genome build used to generate the CNV segments (`hg18` and `hg19` are currently supported. If you want to add another build `BUILD`, please add a folder in `BUILD` in `aCNViewer_DATA/genomes` containing at least a tab-delimited file named `BUILD.chrom.sizes` with each chromosome name and length and a tab-delimited file named `BUILD.centro.txt` with the centromere positions by chr [this file can be generated using `curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/BUILD/database/cytoBand.txt.gz" | gunzip -c | grep acen > centro_build.txt`])

<a id="generalPlotOptions"></a>**The following options are general plotting options:**
* <a id="chrSize"></a>`CHR_SIZE_FILE`: a tab-delimited file with 2 columns respectively chromosome name and chromosome length. When `REF_BUILD` is set, `CHR_SIZE_FILE` is automatically set to `aCNViewer_DATA/genomes/REF_BUILD.chrom.sizes`.
* <a id="centromereFile"></a>`CENTROMERE_FILE`: file giving the centromere bounds. Can be generated using `curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/BUILD/database/cytoBand.txt.gz" | gunzip -c | grep acen > centro_build.txt`. When `REF_BUILD` is set, `CENTROMERE_FILE` is automatically set to `aCNViewer_DATA/genomes/REF_BUILD.centro.txt`.
* <a id="windowSize"></a>`WINDOW_SIZE`: segment size in bp. Please note that alternatively, `-p PERCENTAGE` can be used instead of `-w WINDOW_SIZE` in order to set the segment size in percentage of chromosome length where `PERCENTAGE` is a floating number between 0 and 100. If `WINDOW_SIZE` and `PERCENTAGE` are null then `WINDOW_SIZE` is set to 2Mb by default.
* <a id="sampleFile"></a>`SAMPLE_FILE`: a tab-delimited file that should contain a column named `Sample` with the name of each sample and at least another column with the phenotypic / clinical feature. This file can contain a `sampleAlias` which will be used as the official sample id if provided.
* <a id="phenotypicColumnName"></a>`PHENOTYPIC_COLUMN_NAME` is optional and refers to the name of the column of the phenotypic / clinical feature of interest in `SAMPLE_FILE`. If you omit this parameter, one plot per feature defined in `SAMPLE_FILE` will be generated.
* <a id="rColorFile"></a> `RCOLOR_FILE`: colors in histograms (section "[histogram]". If defined, should contain exactly 10 colors [one per line] corresponding to CNV values in the following order: "&le; -4", "-3", "-2", "-1", "1", "2", "3", "4", "5", "&ge; 6"), dendrograms (section "[group]". If defined, should contain at least the same number of colors than the number of distinct values for the phenotypic / clinical feature of interest) and heatmaps (sections "[chr]" [if defined, should contain 22 colors corresponding to chromosomes 1 to 22], "[group]" and "[heatmap]" [if defined, should contain 10 colors [one per line] corresponding to CNV values in the following order: "0", "1", "2", "3", "4", "5", "6", "7", "8", "&ge; 9"]) can be redefined in that file. An example can be found [here](/img/rColor.txt).
* `PLOT_ALL`: specify whether all available plots should be generated. The default value is `1`.
* `OUTPUT_FORMAT`: allow to customize output formats for the different types of available plots (histograms, heatmaps and dendrograms). Examples of use can be found [above](#outputFormatExamples). The default value is `hist:png(width=4000,height=1800,res=300);hetHom:png(width=4000,height=1800,res=300);dend:png(width=4000,height=2200,res=300);heat:pdf(width=10,height=12)`.
* `PLOIDY_FILE`: custom ploidy values for each sample. Can either be a tab-delimited file with at least 2 columns: "sample" and "ploidy" or an integer which will set the same ploidy to all samples. By default, the ploidy is calculated using the CNV file segmented in fragments of 10% of chromosomal length and its value will be the most represented CNV value for each sample.

<a id="histogramOptions"></a>**The following options are histogram specific**:
* `HISTOGRAM`: specify whether an histogram should be generated. The default value is `0` but its value is overriden to `1` when option `--plotAll 1` is set.
* <a id="lohToPlot"></a>`LOH_TO_PLOT`: histogram option for LOH plotting. Values should be one of "cn-LOH" for plotting cn-LOH only, "LOH" for LOH only, "both" for cn-LOH and LOH or "none" to disable this feature. The default value is "cn-LOH".
* <a id="useFullRes"></a>`USE_FULL_RESOLUTION_FOR_HIST`: tell whether to plot histogram using full resolution i.e. CNVs are not segmented according to a user-defined length. The default value is `1`. If `0`, the resolution of the plot will be given by either <a href="#windowSize">`WINDOW_SIZE`</a> or <a href="#windowSize">`PERCENTAGE`</a>.


<a id="gisticOptions"></a>**The following options are GISTIC options** (more details can be found [here](ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm)):
* `RUN_GISTIC`: specify whether to run GISTIC in order to have a statistical way to prioritize regions of interest. The default value is `0`.
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
* `GROUP_LEGEND_POS` is an optional parameter setting the phenotypic / clinical feature legend's position within the heatmap. The default value is `topright` and can be changed to coordinates (for example `0.1,0.5` which will put the legend at 10% of the total width of the graph on the x axis and 50% of the total height of the graph on the y axis i.e. in the middle of the y axis) or in R specified logical location (`top`, `bottom`, `left`, `right`, etc)
* `CHR_LEGEND_POS` is an optional parameter setting the chromosome legend's position within the heatmap. The default value is `bottomleft` and can be changed to coordinates (for example `0.1,0.5` which will put the legend at 10% of the total width of the graph on the x axis and 50% of the total height of the graph on the y axis i.e. in the middle of the y axis) or in R specified logical location (`top`, `bottom`, `left`, `right`, etc)
* <a href="#rColorFile">`RCOLOR_FILE`</a>
* `USE_RELATIVE_COPY_NB_FOR_CLUSTERING` is an optional parameter specifying whether the CNV matrix used for the heatmap should be relative copy number values or not. The default value is `0`. If `PLOT_ALL` is `1` then plots for both values of `USE_RELATIVE_COPY_NB_FOR_CLUSTERING` will be generated.
* `KEEP_GENOMIC_POS` is optional and will keep the segmented genome in its original position if set to `1` and not cluster segments according to sample CNV patterns (the default value is `0`).
* `DENDROGRAM` is an optional dendrogram parameter used only if `PLOT_ALL` is set to `0` to tell whether to plot dendrograms or not. The default value is `1`
* `USE_SHAPE` is an optional dendrogram parameter and if set to `1` (default value) will replace sample labels with colored shapes in the leaves of the dendrogram(s).


##### TestAffyCel

**Generate a quantitative stacked histogram from CEL files (subset of [data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9845) of hepatocellular carcinomas with hepatitis C virus etiology used in Chiang et al. Cancer Res, 2008) with a window size of 2Mbp:**

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrays250k_sty/ -t TEST_AFFY_CEL --refBuild hg18 -w 2000000 -b aCNViewer_DATA/bin --platform Affy250k_sty -l aCNViewer_DATA/snpArrays250k_sty/LibFiles/`

If ASCAT is not installed (i.e you are not using the [docker](https://hub.docker.com/r/fjdceph/acnviewer/) application) and if you want to install it into a custom R library folder, please add the following option to the previous command line: `--rLibDir RLIB`.

The histogram `OUTPUT_DIR/lrr_baf.segments_merged_hist_2000000nt.png` can also be found in `aCNViewer_DATA/snpArrays250k_sty/expectedResults/test1/lrr_baf.segments_merged_hist_2000000nt.png`.


==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f CEL_DIR --refBuild` [`REF_BUILD`](#refbuild) `-t OUTPUT_DIR -b` <a href="#binDir">`BIN_DIR`</a> `--platform AFFY_PLATFORM -l AFFY_LIB_DIR [--gw6Dir GW6_DIR] [--gcFile ASCAT_GC_FILE]` [[`GENERAL_PLOT_OPTIONS`](#generalPlotOptions)] [[`HISTOGRAM_OPTIONS`]](#histogramOptions) [`[GISTIC_OPTIONS]`](#gisticOptions) [`[HEATMAP_DENDRO_OPTIONS]`](#heatmapDendroOptions)<br>
where:
* `CEL_DIR` is the folder containing ".cel" ou ".cel.gz" files
* `AFFY_PLATFORM`: name of ASCAT supported Affymetrix platform with a GC content file available ("Affy250k_sty", "Affy250k_nsp", "Affy500k" or "AffySNP6"). Please refer to [ASCAT website](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) for more details
* `AFFY_LIB_DIR`: Affymetrix library file downloadable from [Affymetrix website](http://www.affymetrix.com/support/technical/byproduct.affx?cat=dnaarrays)
* `GW6_DIR` is optional and refers to the folder where [gw6.tar.gz](http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz) has been uncompressed into. This archive contains different programs and files necessary to process Affymetrix SNP array and has been uncompressed into `aCNViewer_DATA/bin/PennCNV/gw6/` (default value).
* <a id="ascatGcFile">`ASCAT_GC_FILE`</a>: GC content file necessary for ASCAT GC correction when analyzing SNP array data. This parameter is optional as its value will be automatically deduced from the value of `AFFY_PLATFORM`. Please check [ASCAT website](https://www.crick.ac.uk/peter-van-loo/software/ASCAT) for available GC content files. It is also possible to [create custom GC file](https://github.com/Crick-CancerGenomics/ascat/tree/master/gcProcessing).


#### Illumina

##### TestIllu660k

Generate a quantitative stacked histogram from [raw Illumina data](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE47357&format=file&file=GSE47357%5FMatrix%5Fsignal%5F660w%2Etxt%2Egz) from non-Hodgkin lymphoma patients used in Yang F *et al*. PLoS One 2014 with a window size of 2Mbp:

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/snpArrayIllu660k/GSE47357_Matrix_signal_660w.txt.gz -t TEST_ILLU --refBuild hg19 -w 2000000 -b aCNViewer_DATA/bin --probeFile aCNViewer_DATA/snpArrayIllu660k/Human660W-Quad_v1_H_SNPlist.txt --platform Illumina660k --beadchip "human660w-quad"`

==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f ILLU_FILES --refBuild` [`REF_BUILD`](#refbuild) `-b` <a href="#binDir">`BIN_DIR`</a> `[--sampleList SAMPLE_TO_PROCESS_FILE] --probeFile PROBE_POS_FILE --platform ILLUMINA_PLATFORM [--beadchip BEADCHIP] [-g ASCAT_GC_FILE]` [[`GENERAL_PLOT_OPTIONS`](#generalPlotOptions)] [[`HISTOGRAM_OPTIONS`]](#histogramOptions) [`[GISTIC_OPTIONS]`](#gisticOptions) [`[HEATMAP_DENDRO_OPTIONS]`](#heatmapDendroOptions)<br>
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

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/wes/bams/ -t TEST_WES_RAW --refBuild hg19 -w 2000000 -b aCNViewer_DATA/bin --fileType Sequenza --samplePairFile aCNViewer_DATA/wes/bams/sampleFile.txt`

==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f BAM_DIR -t OUTPUT_DIR --refBuild` [`REF_BUILD`](#refbuild) `-b` <a href="#binDir">`BIN_DIR`</a> `--fileType Sequenza --samplePairFile SAMPLE_PAIR_FILE [-r REF_FILE] [--byChr 1] [-n NB_THREADS] [--createMpileUp CREATE_MPILEUP] [--pattern BAM_FILE_PATTERN] [-M MEMORY]` [[`GENERAL_PLOT_OPTIONS`](#generalPlotOptions)] [[`HISTOGRAM_OPTIONS`]](#histogramOptions) [`[GISTIC_OPTIONS]`](#gisticOptions) [`[HEATMAP_DENDRO_OPTIONS]`](#heatmapDendroOptions)<br>
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
  * `REF_FILE` is the reference file in fasta format used to generate the bam files. When `REF_BUILD` is set, `REF_FILE` is automatically set to the fasta file present in `aCNViewer_DATA/genomes/REF_BUILD`.
  * `BY_CHR` is an optional parameter to indicate whether Sequenza should create seqz (Sequenza intermediate file) files by chromosome or not (the default value is `1`)
  * `NB_THREADS` is an optional parameter specifying the number of cores which will be used for each sample pair to create chromosomal seqz files if `BY_CHR` has been set to `1`. If aCNViewer is ran on a <a href="#supportedClusters">supported computer cluster</a> master node, jobs will be submitted to the cluster. Otherwise, multi-threading will be used run Sequenza.
  * `CREATE_MPILEUP` is an optional parameter telling Sequenza whether to create intermediate mpileup files when generating results. The default value is `1` and it is recommended not to change its value as Sequenza may freeze in some cases when set to `0`.
  * `MEMORY`: optional argument specifying allocated memory in GB to run Sequenza when using a computer cluster. The default value is 8 (GB) and should work for most WES analysis


#### TestSequenzaCNVs
Generate quantitative stacked histogram from Sequenza results with a window size of 2Mbp:<br>

[aCNViewer_DATA.tar.gz](http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz) is required to run this example.

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/wes/ -t TEST_WES_SEQUENZA --refBuild hg19 -w 2000000 -b aCNViewer_DATA/bin --fileType Sequenza`

The histogram `TARGET_DIR/ascat_merged_hist_2000000nt.png` can be found in `aCNViewer_DATA/wes/expectedResults/ascat_merged_hist_2000000nt.png`.


==**Here is the full command:**==

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f SEQUENZA_RES_DIR --fileType Sequenza -t TARGET_DIR --refBuild` [`REF_BUILD`](#refbuild) `-b` <a href="#binDir">`BIN_DIR`</a> [[`GENERAL_PLOT_OPTIONS`](#generalPlotOptions)] [[`HISTOGRAM_OPTIONS`]](#histogramOptions) [`[GISTIC_OPTIONS]`](#gisticOptions) [`[HEATMAP_DENDRO_OPTIONS]`](#heatmapDendroOptions)<br>
where:
* `SEQUENZA_RES_DIR` is the folder containing Sequenza results (`*_segments.txt`)
  
  
### Processing CNV files

At the moment, ASCAT segment file, PennCNV and Sequenza results can be used as an input to aCNViewer. It is possible however to feed aCNViewer with CNV results from any other softwares as explained in the section [below] (#OtherCNVformats).


Both examples below require to download [aCNViewer_DATA.tar.gz](http://www.cephb.fr/tools/aCNViewer/aCNViewer_DATA.tar.gz).


#### [Ascat file generated from Affymetrix SNP arrays](#testaffyascat)

#### PennCNV

Generate quantitative stacked histogram from PennCNV results (79 samples from [Hapmap3](ftp://ftp.ncbi.nlm.nih.gov/hapmap/raw_data/hapmap3_affy6.0/Broad_hapmap3_r2_Affy6_cels_excluded.tgz)):<br>

<a href="#dockerOrPython">`DOCKER_OR_PYTHON`</a> `-f aCNViewer_DATA/pennCNV/hapmap3.rawcnv -t TEST_PENN_CNV --refBuild hg18 -b aCNViewer_DATA/bin --lohToPlot none`

#### [Sequenza segments](#testsequenzacnvs)



#### OtherCNVformats

CNV results from any software can be processed by aCNViewer if formatted in the ASCAT segment file format i.e. a tab-delimited file with the following columns: 
* `sample`
* `chr`
* `startpos`
* `endpos`
* `nMajor`
* `nMinor`

The result file should be sorted according to the following ordered column names: `sample`, `chr`, `startpos`, `endpos` and chromosome names in the `chr` column should not contain the prefix `chr` so `chr1` should appear as `1`. If there is only a global CNV value `v` (and this no allele-specific CNV value), `nMajor` and `nMinor` can take any value as long as `nMajor + nMinor = v`. When plotting the quantitative histogram, add option `--lohToPlot none` to disable LOH plotting.


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
| tumorSep*.png	 | plot of BAF values |
| .ascatInfo.txt | ASCAT values of aberrantcellfraction, goodnessOfFit, psi and ploidy for all samples |
| .segments.txt | list of all CNVs with the copy number for each allele |


#### GISTIC

For the full list of GISTIC output files, please refer to the section `Output Files` of the following [website](ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm). Here are the main output files of interest:

| File | Description |
| --- | --- |
| broad_significance_results.txt | The list of broad events with related q-values and frequencies |
| all_lesions.conf_*.txt | the list of all focal events along with their level of significance |
| amp_* | list of all focal amplification events |
| del_* | list of all focal deletion events |


#### Sequenza

The Sequenza results of each sample pair are stored in a folder named `TUMOR_NORMAL_sequenza` in the `sequenza` folder and contains the following files:

| File | Description |
| --- | --- |
| `*_segments.txt` | predicted CNVs |
| `*_CP_contours.pdf`, `*_confints_CP.txt` & `*_model_fit.pdf` | inferred cellularity and ploidy |
| `*_alternative_fit.pdf` & `*_alternative_solutions.txt` | alternative inferred cellularities and ploidies |
| `*_chromosome_view.pdf` | chromosome view with mutations, BAF, depth ratio and segments |
| `*_genome_view.pdf` | genome view of all the CNVs |
| `*_mutations.txt` | list of detected mutations |
| `*_CN_bars.pdf` | frequency of all the copy number values |

For more information about Sequenza output files, please refer to its [user guide](https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.pdf).


#### HistogramOutputs

When generating histograms, 3 text files with the suffix `_samples.txt` will be created along:
* one with all the genomic segments
* one with only the LOH events (file with suffix `_loh_samples.txt`)
* one with only the cn-LOH (file with suffix `_cnLoh_samples.txt`)

Each file is in the same format with the following columns:
* `CNV key`: the relative copy number value compared to the tumor ploidy
* `chrName`
* `start`: the middle of the segment so the real start is `start - segmentLength / 2`
* `segmentLength`: length of the current segment
* `percentage`: the percentage of samples with the relative ploidy value in `CNV key` for the segment (`chrName`, [`start - segmentLength / 2`, `start + segmentLength / 2`])
* `samples`: the list of the samples falling in the above category


The following files are created as well:
* `*_10pc_ploidy.txt` is a matrix of segments of 10% chromosomal length for all samples. The last column indicates the calculated ploidy which corresponds to the most frequent ploidy

* `*.R` are R scripts used to create the various graphical representations. You can modify and re-run these scripts if you want to further customize your graphical outputs and if aCNViewer do not propose the customizations you are looking for.


##### Dendrograms and heatmaps

2 folders (`relCopyNb` and `rawCopyNb`) will be created and will respectively contain graphs generated from relative copy number values and raw copy number values.


