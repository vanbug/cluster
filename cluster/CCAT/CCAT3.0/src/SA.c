/**************************************************************************************************
SA.c: source file, computing parameters and determine cut-off threshold of fold change, post-processing
Author: Xu Han (hanxu@gis.a-star.edu.sg)
Date: 15/11/2008
***************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "ccat.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define Q_VALUE_STEP_NUM 1000

//SignificanceAnalysis: perform significance analysis, will be called by the main routine of CCAT
int SignificanceAnalysis(CHROM_INFO *chroms, int chromNum, double l1Ratio, double l2Ratio, int maxL1Count, int maxL2Count,CONFIG_STRUCT config);

//ComputeThreshold: compute the local FDR, return the FDR at the cut-off threshold;
double ComputeLocalFDR(PEAK_STRUCT *l1Peaks, int l1PeakNum, PEAK_STRUCT *l2Peaks, int l2PeakNum, double *q, double *value, CONFIG_STRUCT config);

//PostProcessing: post-processing of the peaks: boostrapping for fold-change calculation, region identification, peak refinement
int PostProcessing(CHROM_INFO *chroms, int chromNum, double l1Ratio, double l2Ratio, int maxL1Count, int maxL2Count, CONFIG_STRUCT config);

//ProcessOneChrom: process one chromosome
int ProcessOneChrom(CHROM_INFO *chrom, double l1Ratio, double l2Ratio, CONFIG_STRUCT config);

//GenRegionProfile: generate profile of significant regions. 
int GenRegionProfile(int *l1Profile, int *l2Profile, int *region, int len, double l1Ratio, double l2Ratio, CONFIG_STRUCT config);

//BootstrapFoldChange: compute the fold-change by bootstrapping
double BootstrapFoldChange(int l1Count, int l2Count, double l1Ratio, double l2Ratio, int bootstrapPass);

//Compute smoothing parameter based on negative binomial distribution
double ComputeSmoothingParameter(PEAK_STRUCT *peaks, int peakNum);

//Compute the fold-change
double ComputeFoldChange(int count1, int count2);

double *lookUpTable;
int *flag;
int row, column;
double q[Q_VALUE_STEP_NUM], value[Q_VALUE_STEP_NUM];
double smoothingFactor;
int binCount, tagCount;

//SignificanceAnalysis: perform significance analysis, will be called by the main routine of CCAT
int SignificanceAnalysis(CHROM_INFO *chroms, int chromNum, double l1Ratio, double l2Ratio, int maxL1Count, int maxL2Count,CONFIG_STRUCT config)
{
	PEAK_STRUCT *l1Peaks, *l2Peaks;
	int l1PeakNum, l2PeakNum;
	int i,j;
	double threshFDR;
	int sum = 0;

	//Get peaks

	l1PeakNum = 0;
	l2PeakNum = 0;

	for (i=0;i<chromNum;i++)
	{
		l1PeakNum += chroms[i].l1PeakNum;
		l2PeakNum += chroms[i].l2PeakNum;
	}

	if ((l1PeakNum<=0)||(l2PeakNum<=0))
	{
		printf("no peak detected, CCAT aborted!\n");
		return -1;
	}

	l1Peaks = (PEAK_STRUCT *)malloc(l1PeakNum*sizeof(PEAK_STRUCT));
	l2Peaks = (PEAK_STRUCT *)malloc(l2PeakNum*sizeof(PEAK_STRUCT));

	l1PeakNum = 0;
	l2PeakNum = 0;

	//re-sampling

	for (i=0;i<chromNum;i++)
	{
		for (j=0;j<chroms[i].l1PeakNum;j++)
		{
			memcpy(l1Peaks+l1PeakNum, chroms[i].l1Peaks+j, sizeof(PEAK_STRUCT));
			l1PeakNum++;
		}

		for (j=0;j<chroms[i].l2PeakNum;j++)
		{
			memcpy(l2Peaks+l2PeakNum, chroms[i].l2Peaks+j, sizeof(PEAK_STRUCT));
			l2PeakNum++;
		}
	}
	
	printf("estimating paramters......\n");

	smoothingFactor = ComputeSmoothingParameter(l1Peaks, l1PeakNum);

	if (smoothingFactor<0)
	{
		printf("smoothingFactor estimation error, set to be 1.\n");
		smoothingFactor = 1.0;
	}
	else
	{
		printf("smoothing factor = %f\n", smoothingFactor);
	}

	printf("estimating local FDR.....\n");

	threshFDR = ComputeLocalFDR(l1Peaks, l1PeakNum, l2Peaks, l2PeakNum,q, value, config);

	printf("FDR at the threshold %f is %f. ", config.minScore, threshFDR);

	if (threshFDR>0.2)
	{
		printf("You may want to adjust the threshold by inspecting the output files.\n");
	}
	else
	{
		printf("\n");
	}

	free(l1Peaks);
	free(l2Peaks);

	printf("post-processing......\n");

	PostProcessing(chroms, chromNum, l1Ratio, l2Ratio, maxL1Count, maxL2Count,config);

	printf("post-processing finished.\n");

	return 1;
}

//ComputeThreshold: compute the local FDR, return the FDR at the cut-off threshold;
double ComputeLocalFDR(PEAK_STRUCT *l1Peaks, int l1PeakNum, PEAK_STRUCT *l2Peaks, int l2PeakNum, double *q, double *value, CONFIG_STRUCT config)
{
	int posCount, negCount;
	int i, j, tmpIndex;
	double score;
	double *p1, *p2;
	double threshFDR;

	p1 = (double *)malloc(l1PeakNum*sizeof(double));
	p2 = (double *)malloc(l2PeakNum*sizeof(double));

	for (i=0;i<l1PeakNum;i++)
	{
		p1[i] = ComputeFoldChange(l1Peaks[i].reSampledL1Count, l1Peaks[i].reSampledL2Count);
	}

	for (i=0;i<l2PeakNum;i++)
	{
		p2[i] = ComputeFoldChange(l2Peaks[i].reSampledL1Count, l2Peaks[i].reSampledL2Count);
	}
	
	QuicksortF(p1, 0, l1PeakNum-1);
	QuicksortF(p2, 0, l2PeakNum-1);

	for (j=0;j<Q_VALUE_STEP_NUM; j++)
	{
		q[j] = (double)(j+1)/Q_VALUE_STEP_NUM;

		score = 0;
		tmpIndex = -1;

		for (i=l1PeakNum-1;i>=0;i--)
		{
			if (i<l1PeakNum-1)
			{
				if (p1[i]-p1[i+1]<0.00000001)
				{
					continue;
				}
			}

			posCount = i+1;
			negCount = bTreeSearchingF(p1[i]-0.00000001, p2, 0, l2PeakNum-1)+1;

			if (posCount*q[j]-negCount>score)
			{
				score = posCount*q[j]-negCount;
				tmpIndex = i;
			}
		}

		if (tmpIndex == -1)
		{
			value[j] = DBL_MAX;
		}
		else
		{
			value[j] = p1[tmpIndex];
		}
	}

	threshFDR = 0.0;

	for (i=1;i<Q_VALUE_STEP_NUM;i++)
	{
		if (value[i]>value[i-1])
		{
			value[i] = value[i-1];
		}

		if (value[i]>config.minScore)
		{
			threshFDR = (double)i/Q_VALUE_STEP_NUM;
		}
	}

	free(p1);
	free(p2);

	return threshFDR;
}

//PostProcessing: post-processing of the peaks: boostrapping for fold-change calculation, region identification, peak refinement
int PostProcessing(CHROM_INFO *chroms, int chromNum, double l1Ratio, double l2Ratio, int maxL1Count, int maxL2Count, CONFIG_STRUCT config)
{
	int i;

	row = maxL1Count+1;
	column = maxL2Count+1;

	lookUpTable = (double *)malloc(row*column*sizeof(double));
	flag = (int *)malloc(row*column*sizeof(int));

	if ((!lookUpTable)||(!flag))
	{
		printf("not enough memory!\n");
		return -1;
	}

	memset(flag, 0, maxL1Count*maxL2Count*sizeof(int));

	for (i=0;i<chromNum;i++)
	{
		printf("%s......", chroms[i].chromName);

		if (ProcessOneChrom(chroms+i, l1Ratio, l2Ratio, config)<0)
		{
			free(lookUpTable);
			free(flag);
			return -1;
		}

		printf("finished.\n");
	}

	free(lookUpTable);
	free(flag);

	return 1;
}


//ProcessOneChrom: process one chromosome
int ProcessOneChrom(CHROM_INFO *chrom, double l1Ratio, double l2Ratio, CONFIG_STRUCT config)
{
	int i;
	int *region, *profile1,*profile2;
	int profileLen = (chrom->chromSize)/config.movingStep+1;
	int tmpIndex, tmpStart, tmpEnd, lastStart, lastEnd;

	//allocate memory for the profile
	region = (int *)malloc(profileLen*sizeof(int));
	profile1 = (int *)malloc(profileLen*sizeof(int));
	profile2 = (int *)malloc(profileLen*sizeof(int));

	if ((!region)||(!profile1)||(!profile2))
	{
		printf("not enough memory!\n");
		return -1;
	}

	memset(region, 0, profileLen*sizeof(int));
	memset(profile1, 0, profileLen*sizeof(int));
	memset(profile2, 0, profileLen*sizeof(int));

	//generate profile, assign reads from both ChIP and control libraries to the profile.

	for (i=0;i<chrom->l1PosTagNum;i++)
	{
		if (chrom->l1PosTags[i]>=chrom->chromSize)
		{
			break;
		}

		tmpIndex = (chrom->l1PosTags[i]+config.fragmentSize/2)/config.movingStep;

		tmpIndex = tmpIndex<0?0:tmpIndex;
		tmpIndex = tmpIndex>=profileLen?profileLen-1:tmpIndex;

		profile1[tmpIndex]++;
	}

	for (i=0;i<chrom->l1NegTagNum;i++)
	{
		if (chrom->l1NegTags[i]>=chrom->chromSize)
		{
			break;
		}

		tmpIndex = (chrom->l1NegTags[i]-config.fragmentSize/2)/config.movingStep;

		tmpIndex = tmpIndex<0?0:tmpIndex;
		tmpIndex = tmpIndex>=profileLen?profileLen-1:tmpIndex;

		profile1[tmpIndex]++;
	}

	for (i=0;i<chrom->l2PosTagNum;i++)
	{
		if (chrom->l2PosTags[i]>=chrom->chromSize)
		{
			break;
		}

		tmpIndex = (chrom->l2PosTags[i]+config.fragmentSize/2)/config.movingStep;

		tmpIndex = tmpIndex<0?0:tmpIndex;
		tmpIndex = tmpIndex>=profileLen?profileLen-1:tmpIndex;

		profile2[tmpIndex]++;
	}

	for (i=0;i<chrom->l2NegTagNum;i++)
	{
		if (chrom->l2NegTags[i]>=chrom->chromSize)
		{
			break;
		}

		tmpIndex = (chrom->l2NegTags[i]-config.fragmentSize/2)/config.movingStep;

		tmpIndex = tmpIndex<0?0:tmpIndex;
		tmpIndex = tmpIndex>=profileLen?profileLen-1:tmpIndex;

		profile2[tmpIndex]++;
	}

	SmoothProfile(profile1, profileLen, config.slidingWinSize/config.movingStep/2);

	SmoothProfile(profile2, profileLen, config.slidingWinSize/config.movingStep/2);

	GenRegionProfile(profile1,profile2,region,profileLen,l1Ratio,l2Ratio,config);

	lastStart = -1;
	lastEnd = -1;

	for (i=0;i<chrom->l1PeakNum;i++)
	{
		chrom->l1Peaks[i].foldChange = lookUpTable[chrom->l1Peaks[i].l1Count*column+chrom->l1Peaks[i].l2Count];

		tmpIndex = bTreeSearchingF(chrom->l1Peaks[i].foldChange-0.00000001, value,0, Q_VALUE_STEP_NUM-1);

		chrom->l1Peaks[i].qValue = q[tmpIndex];

		if (chrom->l1Peaks[i].foldChange>config.minScore)
		{
			chrom->l1Peaks[i].isSignificant = 1;
		}
		else
		{
			chrom->l1Peaks[i].isSignificant = 0;
		}

		if (chrom->l1Peaks[i].isSignificant)
		{
			if (chrom->l1Peaks[i].peak<=lastEnd)
			{
				tmpStart = lastStart;
				tmpEnd = lastEnd;
			}
			else
			{
				for (tmpStart = chrom->l1Peaks[i].peak-1; tmpStart>=0; tmpStart--)
				{
					if (!region[tmpStart])
					{
						break;
					}
				}

				tmpStart++;

				for (tmpEnd = chrom->l1Peaks[i].peak+1; tmpEnd<profileLen;tmpEnd++)
				{
					if (!region[tmpEnd])
					{
						break;
					}
				}

				tmpEnd--;

				lastStart = tmpStart;
				lastEnd = tmpEnd;
			}

			chrom->l1Peaks[i].start = tmpStart*config.movingStep;
			chrom->l1Peaks[i].end = (tmpEnd+1)*config.movingStep;
		}
		else
		{
			chrom->l1Peaks[i].start = chrom->l1Peaks[i].peak*config.movingStep+config.movingStep/2;
			chrom->l1Peaks[i].end = chrom->l1Peaks[i].peak*config.movingStep+config.movingStep/2;
		}

		chrom->l1Peaks[i].peak = chrom->l1Peaks[i].peak*config.movingStep+config.movingStep/2;
	}

	free(profile1);
	free(profile2);
	free(region);

	return 1;
}


//GenRegionProfile: generate profile of significant regions. 
int GenRegionProfile(int *l1Profile, int *l2Profile, int *region, int len, double l1Ratio, double l2Ratio, CONFIG_STRUCT config)
{
	int i,j;
	double tmpF;
	int halfWinSize = config.slidingWinSize/config.movingStep/2;
	int *tmpRegion;

	tmpRegion = (int *)malloc(len*sizeof(int));

	if (!tmpRegion)
	{
		return -1;
	}

	for (i=0;i<len;i++)
	{
		if ((l1Profile[i]>=row)||(l2Profile[i]>=column))
		{
			printf("somethings wrong.\n");
			return -1;
		}

		if (flag[l1Profile[i]*column+l2Profile[i]])
		{
			tmpRegion[i] = lookUpTable[l1Profile[i]*column+l2Profile[i]]>config.minScore?1:0;
			continue;
		}

		tmpF = BootstrapFoldChange(l1Profile[i],l2Profile[i], l1Ratio,l2Ratio,config.bootstrapPass);

		if (tmpF>config.minScore)
		{
			tmpRegion[i] = 1;
		}
		else
		{
			tmpRegion[i] = 0;
		}

		lookUpTable[l1Profile[i]*column+l2Profile[i]] = tmpF;
		flag[l1Profile[i]*column+l2Profile[i]] = 1;
	}

	memcpy(region, tmpRegion, len*sizeof(int));

	for (i=1;i<len-1;i++)
	{
		if ((tmpRegion[i]==1)&&(tmpRegion[i-1]==0))
		{
			for (j=i-1;j>=i-halfWinSize;j--)
			{
				if ((tmpRegion[j])||(j<0))
				{
					break;
				}

				region[j] = 1;
			}
		}

		if ((tmpRegion[i]==1)&&(tmpRegion[i+1]==0))
		{
			for (j=i+1;j<=i+halfWinSize;j++)
			{
				if ((tmpRegion[j])||(j>=len))
				{
					break;
				}

				region[j] = 1;
			}
		}
	}

	free(tmpRegion);

	return 1;
} 

//BootstrapFoldChange: compute the fold-change by bootstrapping
double BootstrapFoldChange(int l1Count, int l2Count, double l1Ratio, double l2Ratio, int bootstrapPass)
{
	int count1, count2;
	double *tmpA, foldChange;
	int i,j;

	tmpA = (double *)malloc(bootstrapPass*sizeof(double));

	for (i=0;i<bootstrapPass;i++)
	{
		count1 = 0;

		for (j=0;j<l1Count;j++)
		{
			if (rand()>RAND_MAX*l1Ratio)
			{
				continue;
			}

			count1++;
		}

		count2 = 0;

		for (j=0;j<l2Count;j++)
		{
			if (rand()>RAND_MAX*l2Ratio)
			{
				continue;
			}

			count2++;
		}

		if (count1+count2==0)
		{
			tmpA[i] = 0;
		}
		else
		{
			tmpA[i] = ComputeFoldChange(count1, count2);
		}
	}

	foldChange = 0.0;

	for (i=0;i<bootstrapPass;i++)
	{
		foldChange += tmpA[i];
	}

	free(tmpA);

	return foldChange/bootstrapPass;
}

//Compute smoothing parameter based on negative binomial distribution
double ComputeSmoothingParameter(PEAK_STRUCT *peaks, int peakNum)
{
	int i;
	int hist[1000];
	double m;

	memset(hist, 0, 1000*sizeof(int));

	for (i=0;i<peakNum;i++)
	{
		if (peaks[i].reSampledL2Count<1000)
		{
			hist[peaks[i].reSampledL2Count]++;
			binCount ++;
			tagCount += peaks[i].reSampledL2Count;
		}
	}

	FitNegBinomDist(hist, 1000, &m, &(smoothingFactor), 0.1, 1000.0);

	return smoothingFactor;
}

//Compute the fold-change
double ComputeFoldChange(int count1, int count2)
{
	return (double)count1/(count2+smoothingFactor)*(tagCount+smoothingFactor*binCount)/tagCount;
}