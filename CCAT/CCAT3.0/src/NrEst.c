/**************************************************************************************************
NrEst.c: source file, noise rate estimation
Author: Xu Han (hanxu@gis.a-star.edu.sg)
Date: 15/11/2008
***************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include "ccat.h"

#define MIN_COUNT 10000
#define BIN_SIZE 1000
#define MAX_ITERATION 20

//ComputeNoiseRate: compute the noise rate
double ComputeNoiseRate(CHROM_INFO *chroms, int chromNum, double *l1Ratio, double *l2Ratio, CONFIG_STRUCT config);

//ComputeNoiseRateInOneIteration: compute the noise rate by one iteration
int ComputeNoiseRateInOneIteration(CHROM_INFO *chroms, int chromNum, double *noiseRate, CONFIG_STRUCT config);

//CountFragmentsInOneChrom: sample portion of the tags by ratio, and add fragment counts to the bins
void CountFragmentsInOneChrom(CHROM_INFO *chrom, BIN_STRUCT *bins, int binNum, double l1Ratio, double l2Ratio, CONFIG_STRUCT config);

//ComputeSampleRatio: compute the ratio of reSampling based on the noiseRate
void ComputeSampleRatio(double noiseRate, int l1TagNum, int l2TagNum, double *l1Ratio, double *l2Ratio);

//Allocate memory for the bins
int AllocBinMem(CHROM_INFO *chroms, int chromNum, BIN_STRUCT **bins, int *binNum);

//Free memory of bins
void FreeBinMem(BIN_STRUCT **bins, int *binNum, int chromNum);

//Free the memory in CHROM_INFO structures
void FreeChromMem(CHROM_INFO *chroms, int chromNum);

//ComputeNoiseRate: compute the noise rate
double ComputeNoiseRate(CHROM_INFO *chroms, int chromNum, double *l1Ratio, double *l2Ratio, CONFIG_STRUCT config)
{
	int i;
	double nr, prevNr;
	double sum;
	int count, flag;
	int l1TagNum = 0, l2TagNum = 0;

	nr = 1.0;
	prevNr = 1.0;

	sum = 0;
	count = 0;
	flag = 0;

	for (i=0;i<MAX_ITERATION;i++)
	{
		if (ComputeNoiseRateInOneIteration(chroms, chromNum,&nr,config)<0)
		{
			printf("noise rate estimation failed. noise rate set to be 1.0\n");
			nr = 1.0;
			break;
		}

		if (nr>1.0)
		{
			printf("noise rate estimation failed. noise rate set to be 1.0\n");
			nr = 1.0;
			break;
		}

		if (nr>prevNr)
		{
			flag = 1;
		}

		prevNr = nr;

		if (flag)
		{
			sum+=nr;
			count++;
		}

		printf("iteration %d: nr=%f\n", i, nr); 
	}

	if (count>0)
	{
		nr = sum/count;
	}

	for (i=0;i<chromNum;i++)
	{
		l1TagNum += chroms[i].l1PosTagNum+chroms[i].l1NegTagNum;
		l2TagNum += chroms[i].l2PosTagNum+chroms[i].l2NegTagNum;
	}

	ComputeSampleRatio(nr, l1TagNum, l2TagNum, l1Ratio,l2Ratio);

	return nr;
}

//ComputeNoiseRate: compute the noise rate
int ComputeNoiseRateInOneIteration(CHROM_INFO *chroms, int chromNum, double *noiseRate, CONFIG_STRUCT config)
{
	int i,j;
	int l1TagNum=0, l2TagNum = 0;
	double l1Ratio, l2Ratio;
	int thresh, tmpL1Count, tmpL2Count, sum1, sum2;
	int *hist;
	BIN_STRUCT **bins;
	int *binNum;

	for (i=0;i<chromNum;i++)
	{
		l1TagNum += chroms[i].l1PosTagNum+chroms[i].l1NegTagNum;
		l2TagNum += chroms[i].l2PosTagNum+chroms[i].l2NegTagNum;
	}

	ComputeSampleRatio(*noiseRate, l1TagNum, l2TagNum, &l1Ratio, &l2Ratio);
	
	hist = (int *)malloc(l2TagNum*sizeof(int));

	if (!hist)
	{
		printf("not enough memory!\n");
		return -1;
	}

	bins = (BIN_STRUCT **)malloc(chromNum*sizeof(BIN_STRUCT *));
	binNum = (int *)malloc(chromNum*sizeof(int));

	if (AllocBinMem(chroms, chromNum, bins, binNum)<0)
	{
		return -1;
	}

	memset(hist, 0, l2TagNum*sizeof(int));

	for (i=0;i<chromNum;i++)
	{
		CountFragmentsInOneChrom(chroms+i, bins[i], binNum[i], l1Ratio, l2Ratio, config);

		for (j=0;j<binNum[i];j++)
		{
			tmpL1Count = bins[i][j].l1Counts;
			tmpL2Count = bins[i][j].l2Counts;

			hist[tmpL1Count<tmpL2Count?tmpL2Count:tmpL1Count]+= tmpL2Count;
		}
	}
	
	sum1 = 0;

	for (i=0;i<l2TagNum;i++)
	{
		sum1 +=hist[i];

		if (((double)sum1>MIN_COUNT)&&((double)sum1>l2TagNum*l2Ratio*0.5))
		{
			break;
		}
		
	}

	free(hist);

	if (i==l2TagNum)
	{
		printf("not enough reads for noise rate estimation.\n");
		return -1;
	}

	thresh = i;

	sum1 = 0;
	sum2 = 0;

	for (i=0;i<chromNum;i++)
	{
		for (j=0;j<binNum[i];j++)
		{
			tmpL1Count = bins[i][j].l1Counts;
			tmpL2Count = bins[i][j].l2Counts;

			if ((tmpL1Count<tmpL2Count?tmpL2Count:tmpL1Count)<=thresh)
			{
				sum1 += tmpL1Count;
				sum2 += tmpL2Count;
			}
		}
	}

	(*noiseRate)*=(double)sum1/sum2;

	FreeBinMem(bins,binNum, chromNum);

	return 1;
}

//CountFragmentsInOneChrom: sample portion of the tags by ratio, and add fragment counts to the bins
void CountFragmentsInOneChrom(CHROM_INFO *chrom, BIN_STRUCT *bins, int binNum, double l1Ratio, double l2Ratio, CONFIG_STRUCT config)
{
	int j;
	int tmpIndex;

	//count fragments for L1
	for (j=0;j<binNum;j++)
	{
		bins[j].l1Counts = 0;
		bins[j].l2Counts = 0;
	}

	for (j=0;j<chrom->l1PosTagNum;j++)
	{
		if (rand()>RAND_MAX*l1Ratio)
		{
			continue;
		}

		tmpIndex = (chrom->l1PosTags[j]+config.fragmentSize/2)/BIN_SIZE;
		
		tmpIndex = tmpIndex>=binNum?binNum-1:tmpIndex;
		tmpIndex = tmpIndex<0?0:tmpIndex;

		bins[tmpIndex].l1Counts++;
	}

	for (j=0;j<chrom->l1NegTagNum;j++)
	{
		if (rand()>RAND_MAX*l1Ratio)
		{
			continue;
		}

		tmpIndex = (chrom->l1NegTags[j]-config.fragmentSize/2)/BIN_SIZE;
		
		tmpIndex = tmpIndex>=binNum?binNum-1:tmpIndex;
		tmpIndex = tmpIndex<0?0:tmpIndex;

		bins[tmpIndex].l1Counts++;
	}

	for (j=0;j<chrom->l2PosTagNum;j++)
	{
		if (rand()>RAND_MAX*l2Ratio)
		{
			continue;
		}
		
		tmpIndex = (chrom->l2PosTags[j]+config.fragmentSize/2)/BIN_SIZE;
		
		tmpIndex = tmpIndex>=binNum?binNum-1:tmpIndex;
		tmpIndex = tmpIndex<0?0:tmpIndex;

		bins[tmpIndex].l2Counts++;
	}

	for (j=0;j<chrom->l2NegTagNum;j++)
	{
		if (rand()>RAND_MAX*l2Ratio)
		{
			continue;
		}

		tmpIndex = (chrom->l2NegTags[j]-config.fragmentSize/2)/BIN_SIZE;
		
		tmpIndex = tmpIndex>=binNum?binNum-1:tmpIndex;
		tmpIndex = tmpIndex<0?0:tmpIndex;

		bins[tmpIndex].l2Counts++;
	}
}

//ComputeSampleRatio: compute the ratio of reSampling based on the noiseRate
void ComputeSampleRatio(double noiseRate, int l1TagNum, int l2TagNum, double *l1Ratio, double *l2Ratio)
{
	if (l1TagNum*noiseRate > l2TagNum)
	{
		*l1Ratio = (double)(l2TagNum)/(l1TagNum*noiseRate);
		*l2Ratio = 1.0;
	}
	else
	{
		*l2Ratio = (double)(l1TagNum*noiseRate)/l2TagNum;
		*l1Ratio = 1.0;
	}
}


//Allocate memory for the bins
int AllocBinMem(CHROM_INFO *chroms, int chromNum, BIN_STRUCT **bins, int *binNum)
{
	int i;

	for (i=0;i<chromNum;i++)
	{
		if (chroms[i].chromSize>0)
		{
			binNum[i] = chroms[i].chromSize/BIN_SIZE+1;

			bins[i] = (BIN_STRUCT *)malloc(binNum[i]*sizeof(BIN_STRUCT));

			if (!bins[i])
			{
				printf("not enough memory!\n");
				return -1;
			}
		}
	}

	return 1;
}

//Free memory of bins
void FreeBinMem(BIN_STRUCT **bins, int *binNum, int chromNum)
{
	int i;

	for (i=0;i<chromNum;i++)
	{
		if (binNum[i]>0)
		{
			free(bins[i]);
		}
	}

	free(binNum);
	free(bins);
}
