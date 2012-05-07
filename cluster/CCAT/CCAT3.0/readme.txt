User's guide for CCAT (version 3.0)
Author: Xu Han
Contact: hanxu@gis.a-star.edu.sg; sungk@gis.a-star.edu.sg

1. Installing CCAT 

   CCAT was writen in standard C lauguage. This package is the Unix version of CCAT. To compile CCAT in your Unix system, please go to the directory <CCAT path>/src and run:
./make
This command will generate the executable of CCAT under <CCAT path>/bin.

2. Starting with CCAT

The command line for running CCAT is:

<CCAT path>/bin/CCAT  <ChIP library read file name>  <control library read file name> <chromosome length file name>  <config file name>  <project name>

<CCAT path>: the path of the directory with CCAT executable;

<ChIP library read file name>  : The uniquely mapped ChIP-seq reads in ChIP library, organized in .bed format;

<control library read file name>  : The uniquely mapped ChIP-seq reads in control library, organized in .bed format;

<chromosome length file name>  : Descriptions of chromosomes;

<config file name>  : The configuration file;

<project name> : The name of project, assigned by the user.

To save the intermediate information in the process, one may direct the standard output to a log file.

Example:

Go to the /example directory, type:

../bin/CCAT ES_CTCF_chr19.bed ES_GFP_chr19.bed genome_length_mm8_chr19.txt config_TF.txt ES_CTCF_chr19 > ES_CTCF_chr19.log

This command line will generate 4 files:

ES_CTCF_chr19.significant.peak: The peaks that are identified to be significant;

ES_CTCF_chr19.significant.region: The regions that contain significant peaks. Some regions may contain multiple peaks and the most significant peak will be reported;

ES_CTCF_chr19.top100000: The top 100000 peaks. Sometimes you need this file to take a look at the peaks below the cut-off threshold. The number of output peaks is configurable in the config file;

ES_CTCF_chr19.log: the log-file containing the intermediate information.

3. File format

a)	ChIP-seq read file (.bed)

<chromosome> <start of read> <end of read> <not used, fill with any character> <not used, fill with any character> <strand of read>
Example:
chr1 12345 12381 0 0 +
chr1 12400 12436 0 0 -

	(Note: users are required to convert the output of mapping software into the bed format before running CCAT)

b)	chromosome length file

<chromosome name> <length>
      Example:
      chr1    197069962
      chr2    181976762

c)	configuration file

<name of configuration> <value>

Please see section 4 for the details of configurations.

d)	output file of CCAT 

Both the .peak and .region files take the same format

<chromosome> <position of the peak> <start of region> <end of region> <read
counts in ChIP library> <read counts in control library> <fold-change score>
<local FDR>

In the .top100000.peak file, column 2,3,4 take the same value for the peaks below cut-off threshold.


4. Configurations of CCAT

There are 9 parameters in the configuration file of CCAT. 

fragmentSize: the length of DNA fragment, which can be determined either by gel-shift or computationally. Default = 200 bps; 

slidingWinSize: the size of sliding window. The value of this parameter
depends on how sparse the reads are distributed in the signal regions. For the application of transcription factor binding, the default setting is 300 bps; For the application of histone modifications, the default setting is 1000 bps.

movingStep: the step of window sliding. A small value will result in high resolution but heavy computation. 10 bps for TF application and 50 bps for histone modifications.

isStrandSensitiveMode: If set to be 1, the peaks will be determined to be the transition from sense strand to anti-sense strand; if set to be 0, the peaks will be determined to be the local maximum of read-enrichment profile

minCount: the minimum number of read counts at the peak. This parameter is used in the filtering step prior to peak-finding, which reduce the computational time of processing.

outputNum: the number of the peaks reported in .top<n>.peak file, <n> is its value and is set to be 100000 by default.

randomSeed: the random seed. The results will be identical with the same random seeds.

minScore: the minimum score. In CCAT3.0, we use maximum-likelihood estimate of fold change. default = 5.0.

bootstrapPass: the number of passes in the bootstrapping process. default = 50.
 
For the convenience of users, we provide two pre-compiled configuration files in the CCAT directory. config_TF.txt for the applications of TF binding, and config_histone.txt for the applications of histone modification.



      

