/**************************************************************************************************
CCAT.c: main file for CCAT analysis
Author: Xu Han (hanxu@gis.a-star.edu.sg)
Date: 15/11/2008
***************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "ccat.h"

#define MAX_COUNT 10000

//ReadConfigFile: read the configuration file into config. Return the number of parameters read.
extern int ReadConfigFile(char *fileName, CONFIG_STRUCT *config); 

//GetChromInfo: Read chromosome name and length, return the number of chromosomes read
extern int GetChromInfo(char *fileName, CHROM_INFO *chroms);

//GetTags: Read the tags from bed file of L1 and L2 to the structure of chroms. Return the total number of tags read.
extern int GetTags(char *l1FileName, char *l2FileName, CHROM_INFO *chroms, int chromNum);

//ComputeNoiseRate: compute the noise rate
extern double ComputeNoiseRate(CHROM_INFO *chroms, int chromNum, double *l1Ratio, double *l2Ratio, CONFIG_STRUCT config);

//SignificanceAnalysis: perform significance analysis, will be called by the main routine of CCAT
extern int SignificanceAnalysis(CHROM_INFO *chroms, int chromNum, double l1Ratio, double l2Ratio, int maxL1Count, int maxL2Count,CONFIG_STRUCT config);

//PeakFinding: peak processing, will be called by the main routine of CCAT
extern int PeakFinding(CHROM_INFO *chroms, int chromNum, double l1Ratio,double l2Ratio, int *maxL1Count, int *maxL2Count, CONFIG_STRUCT config);

//OutputFiles: output the results
extern int OutputFiles(CHROM_INFO *chroms, int chromNum, char *projectName, CONFIG_STRUCT config);

//PreProcessing: sort and get unique tags
int PreProcessing(CHROM_INFO *chroms, int chromNum);

//Free the memory in CHROM_INFO structures
void FreeChromMem(CHROM_INFO *chroms, int chromNum);

//sort tag list and get unique tags
int SortAndGetUniqueTags(int *tags, int tagNum);

int main(int argc, char* argv[])
{
	char lib1FileName[1000], lib2FileName[1000], chromFileName[1000],configFileName[2000], projectName[1000];
	CONFIG_STRUCT config;
	CHROM_INFO chroms[MAX_CHROM_NUM];
	int chromNum;
	double noiseRate, l1Ratio, l2Ratio;
	int maxL1Count, maxL2Count;
	int i;

	if (argc!=6)
	{
		printf("Usage: <library 1 tag file name> <library 2 tag file name> <chromosome length file name> <config file name> <project name>\n");
		return 0;
	}

	strcpy(lib1FileName, argv[1]);
	strcpy(lib2FileName, argv[2]);
	strcpy(chromFileName, argv[3]);
	strcpy(configFileName, argv[4]);
	strcpy(projectName, argv[5]);

	chromNum = GetChromInfo(chromFileName, chroms);

	if (chromNum<=0)
	{
		return 0;
	}                                                                                                                                                                                                                             

	printf("chromosome length information read. chromNum = %d!\n", chromNum);

	for (i=0;i<chromNum;i++)
	{
		chroms[i].l1PosTagNum = 0;
		chroms[i].l2PosTagNum = 0;
		chroms[i].l1NegTagNum = 0;
		chroms[i].l2NegTagNum = 0;
		chroms[i].l1PeakNum = 0;
		chroms[i].l2PeakNum = 0;
	}

	if (ReadConfigFile(configFileName, &config)<0)
	{
		printf("config file not complet. Default setting used.\n");
	}

	printf("config file read.\n");

	printf("fragmentSize = %d\n", config.fragmentSize);
	printf("isStrandSensitiveMode = %d\n", config.isStrandSensitiveMode);
	printf("slidingWinSize = %d\n", config.slidingWinSize);
	printf("movingStep = %d\n", config.movingStep);
	printf("outputNum = %d\n", config.outputNum);
	printf("minCount = %d\n", config.minCount);
	printf("minScore = %f\n", config.minScore);
	printf("bootstrapPass = %d\n", config.bootstrapPass);
	printf("randSeed = %d\n", config.randomSeed);

	printf("reading tag files......\n");
	
	if (!GetTags(lib1FileName, lib2FileName, chroms, chromNum))
	{
		printf("tag file error!\n");
		return 0;
	}

	printf("tag file read.\n");

	printf("pre-processing......\n");

	if (PreProcessing(chroms, chromNum)<0)
	{
		printf("Preprocessing failed!\n CCAT aborted!\n");
		return -1;
	}

	printf("pre-processing finished.\n");

	printf("estimating noise rate......\n");

	srand(config.randomSeed);

	noiseRate = ComputeNoiseRate(chroms,chromNum, &l1Ratio,&l2Ratio, config);

	printf("noise rate = %f\n", noiseRate);

	printf("peak-finding......\n");

	if (PeakFinding(chroms,chromNum,l1Ratio,l2Ratio, &maxL1Count, &maxL2Count, config)<=0)
	{
		return -1;
	}

	printf("peak-finding finished.\n");

	printf("Significance Analysis......\n");

	if (SignificanceAnalysis(chroms, chromNum, l1Ratio, l2Ratio, maxL1Count, maxL2Count, config)<0)
	{
		printf("DPS paramter esitmation failed!\n");
		return -1;
	}

	printf("Significance analysis finished.\n");

	printf("saving results.......\n");

	OutputFiles(chroms, chromNum, projectName, config);

	printf("results saved.\n");

	FreeChromMem(chroms,chromNum);

	printf("CCAT process completed!\n");
	return 0;
}

//PreProcessing: sort and get unique tags
int PreProcessing(CHROM_INFO *chroms, int chromNum)
{
	int i;

	for (i=0;i<chromNum;i++)
	{
		if ((chroms[i].l1PosTagNum = SortAndGetUniqueTags(chroms[i].l1PosTags,chroms[i].l1PosTagNum))<0)
		{
			return -1;
		}

		if ((chroms[i].l2PosTagNum = SortAndGetUniqueTags(chroms[i].l2PosTags,chroms[i].l2PosTagNum))<0)
		{
			return -1;
		}

		if ((chroms[i].l1NegTagNum = SortAndGetUniqueTags(chroms[i].l1NegTags,chroms[i].l1NegTagNum))<0)
		{
			return -1;
		}

		if ((chroms[i].l2NegTagNum = SortAndGetUniqueTags(chroms[i].l2NegTags,chroms[i].l2NegTagNum))<0)
		{
			return -1;
		}
	}

	return 1;
}
//Free the memory in CHROM_INFO structures
void FreeChromMem(CHROM_INFO *chroms, int chromNum)
{
	int i;

	for (i=0;i<chromNum;i++)
	{
		if (chroms[i].l1PosTagNum)
		{
			free(chroms[i].l1PosTags);
		}

		if (chroms[i].l2PosTagNum)
		{
			free(chroms[i].l2PosTags);
		}

		if (chroms[i].l1NegTagNum)
		{
			free(chroms[i].l1NegTags);
		}

		if (chroms[i].l2NegTagNum)
		{
			free(chroms[i].l2NegTags);
		}

		if (chroms[i].l1PeakNum)
		{
			free(chroms[i].l1Peaks);
		}

		if (chroms[i].l2PeakNum)
		{
			free(chroms[i].l2Peaks);
		}
	}
}

//sort tag list and get unique tags
int SortAndGetUniqueTags(int *tags, int tagNum)
{
	int i;
	int *tmpA;
	int count;

	if (tagNum<=0)
	{
		return tagNum;
	}
	
	tmpA =(int *)malloc(tagNum*sizeof(int));

	if (!tmpA)
	{
		printf("not enough memory!\n");
		return -1;
	}

	memcpy(tmpA, tags, tagNum*sizeof(int));

	QuicksortI(tmpA, 0, tagNum-1);

	count = 0;

	for (i=0;i<tagNum;i++)
	{
		if (i>0)
		{
			if (tmpA[i]==tmpA[i-1])
			{
				continue;
			}
		}

		tags[count] = tmpA[i];
		count++;
	}

	free(tmpA);

	return count;
}
