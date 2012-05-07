/**************************************************************************************************
ccat.h: header file, user defined data structures in the software of CCAT
Author: Xu Han (hanxu@gis.a-star.edu.sg)
Date: 15/07/2008
***************************************************************************************************/

#ifndef CCAT_H
#define CCAT_H
#endif

#include "math_api.h"

#define MAX_CHROM_NAME_LEN	100
#define MAX_CHROM_NUM		100

//bin structure. L1: positive library with specific antibody; L2: negative control library with no antibody (WCE) or non-specific antibody (eg. GFP)
typedef struct
{
	int l1Counts;   //L1 counts of tags 
	int l2Counts;   //L2 counts of tags 
}BIN_STRUCT;

typedef struct
{
	int chromIndex;  //Index of chromosome
	int start;       //start of the peak region
	int end;         //end
	int peak;        //peak location
	int l1Count;     //read counts in L1
	int l2Count;     //read counts in L2
	int reSampledL1Count; //resampled read counts in L1
	int reSampledL2Count; //resampled read counts in L2
	double foldChange;  //fold change
	double qValue;      //q value
	int isSignificant;  //1 if the peak is signficant, otherwise 0
}PEAK_STRUCT;

//information of chromosome
typedef struct
{
	char chromName[MAX_CHROM_NAME_LEN];        //name
	int chromIndex;                            //index
	int chromSize;                             //size
	int l1PosTagNum;                           //number of L1 positive (sense) reads
	int *l1PosTags;                            //location of L1 positive reads
	int l1NegTagNum;                           //number of L1 negative (anti-sense) reads
	int *l1NegTags;                            //location of L1 negative reads
	int l2PosTagNum;                           //number of L2 positive (sense) reads
	int *l2PosTags;                            //location of L2 positive reads
	int l2NegTagNum;                           //number of L2 negative (anti-sense) reads
	int *l2NegTags;                            //location of L2 negative reads
	int l1PeakNum;                               //number of candidate peaks
	PEAK_STRUCT *l1Peaks;                        //candidate peaks
	int l2PeakNum;
	PEAK_STRUCT *l2Peaks;
}CHROM_INFO;

//configuration
typedef struct
{
	//experimental setting	
	int fragmentSize;                          
	int slidingWinSize;
	int movingStep;
	int minCount;
	int isStrandSensitiveMode;
	int outputNum;
	int randomSeed;
	double minScore;
	int bootstrapPass;
	double smoothingFactor;
} CONFIG_STRUCT;
