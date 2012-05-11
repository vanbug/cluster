/**************************************************************************************************
PeakFinding.c: source file, search for peaks
Author: Xu Han (hanxu@gis.a-star.edu.sg)
Date: 15/11/2008
***************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "ccat.h"

#define TRANSITION_SIZE 30

//PeakFinding: peak processing, will be called by the main routine of CCAT
int PeakFinding(CHROM_INFO *chroms, int chromNum, double l1Ratio,double l2Ratio, int *maxL1Count, int *maxL2Count, CONFIG_STRUCT config);

//AssignTagInChrom: add tag counts to the peaks
void AssignTagInChrom(CHROM_INFO *chrom, CONFIG_STRUCT config);

//GetPeaksInOneChrom1: strand-insensitive mode: get peak location from tags for one chromosome
int GetPeaksInOneChrom1(CHROM_INFO *chrom, double l1Ratio, double l2Ratio, int *maxL1Count, int *maxL2Count, CONFIG_STRUCT config);

//GetPeaksInOneChrom2: strand-sensitive mode: get peak location from tags for one chromosome
int GetPeaksInOneChrom2(CHROM_INFO *chrom, double l1Ratio, double l2Ratio, int *maxL1Count, int *maxL2Count, CONFIG_STRUCT config);

//GetLocalMaxima: find the local maxima in the profile
int GetLocalMaxima(int *profile, int profileSize, PEAK_STRUCT *peaks, int minDist, int minCount);

//GetTransition: get the transition from sense reads to anti-sense reads 
int GetTransition(int *profile, int profileSize, PEAK_STRUCT *peaks, int minDist);

//PeakFinding: peak processing, will be called by the main routine of CCAT
int PeakFinding(CHROM_INFO *chroms, int chromNum, double l1Ratio,double l2Ratio, int *maxL1Count, int *maxL2Count, CONFIG_STRUCT config)
{
	int i;
	int peakNum, tmpI;

	*maxL1Count = 0;
	*maxL2Count = 0;
	peakNum = 0;

	for (i=0;i<chromNum;i++)
	{
		printf("%s......", chroms[i].chromName);

		if (config.isStrandSensitiveMode)
		{
			tmpI = GetPeaksInOneChrom2(chroms+i,l1Ratio,l2Ratio,maxL1Count, maxL2Count, config);
		}
		else 
		{
			tmpI = GetPeaksInOneChrom1(chroms+i,l1Ratio,l2Ratio,maxL1Count, maxL2Count, config);
		}

		if (tmpI<0)
		{
			return -1;
		}

		printf("%d candidate peaks.\n", tmpI);

		peakNum += tmpI;
	}

	return peakNum;
}

//GetPeaksInOneChrom1: strand-insensitive mode: get peak location from tags for one chromosome
int GetPeaksInOneChrom1(CHROM_INFO *chrom, double l1Ratio, double l2Ratio, int *maxL1Count, int *maxL2Count, CONFIG_STRUCT config)
{
	int i;
	int *profile1,*profile2, *rsProfile1, *rsProfile2;
	int profileLen = (chrom->chromSize)/config.movingStep+1;
	int tmpIndex;
	PEAK_STRUCT *tmpPeaks;
	int tmpPeakNum;
	int minDist; 

	//allocate memory for the profile
	profile1 = (int *)malloc(profileLen*sizeof(int));
	profile2 = (int *)malloc(profileLen*sizeof(int));
	rsProfile1 = (int *)malloc(profileLen*sizeof(int));
	rsProfile2 = (int *)malloc(profileLen*sizeof(int));

	if ((!profile1)||(!profile2)||(!rsProfile1)||(!rsProfile2))
	{
		printf("not enough memory!\n");
		return -1;
	}

	memset(profile1, 0, profileLen*sizeof(int));
	memset(profile2, 0, profileLen*sizeof(int));
	memset(rsProfile1, 0, profileLen*sizeof(int));
	memset(rsProfile2, 0, profileLen*sizeof(int));

	//generate profile, assign reads from both ChIP and control libraries to the profile. Reads are re-sampled based on noise rate

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

		if (rand()>RAND_MAX*l1Ratio)
		{
			continue;
		}

		rsProfile1[tmpIndex]++;
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

		if (rand()>RAND_MAX*l1Ratio)
		{
			continue;
		}

		rsProfile1[tmpIndex]++;
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

		if (rand()>RAND_MAX*l2Ratio)
		{
			continue;
		}

		rsProfile2[tmpIndex]++;
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

		if (rand()>RAND_MAX*l2Ratio)
		{
			continue;
		}

		rsProfile2[tmpIndex]++;
	}

	SmoothProfile(profile1, profileLen, config.slidingWinSize/config.movingStep/2);

	SmoothProfile(profile2, profileLen, config.slidingWinSize/config.movingStep/2);

	SmoothProfile(rsProfile1, profileLen, config.slidingWinSize/config.movingStep/2);

	SmoothProfile(rsProfile2, profileLen, config.slidingWinSize/config.movingStep/2);

	for (i=0;i<profileLen;i++)
	{
		if (profile1[i]>*maxL1Count)
		{
			*maxL1Count = profile1[i];
		}

		if (profile2[i]>*maxL2Count)
		{
			*maxL2Count = profile2[i];
		}
	}

	//find peaks from profile

	tmpPeaks = (PEAK_STRUCT *)malloc(((chrom->chromSize)/config.movingStep+1)*sizeof(PEAK_STRUCT));

	if (!tmpPeaks)
	{
		printf("not enough memory!\n");
		return -1;
	}

	minDist = config.slidingWinSize/config.movingStep+1;

	tmpPeakNum = GetLocalMaxima(rsProfile1,profileLen,tmpPeaks,minDist,config.minCount);

	chrom->l1PeakNum = 0;

	for (i=0;i<tmpPeakNum;i++)
	{
		tmpPeaks[i].l1Count = profile1[tmpPeaks[i].peak];
		tmpPeaks[i].l2Count = profile2[tmpPeaks[i].peak];
		tmpPeaks[i].reSampledL1Count = rsProfile1[tmpPeaks[i].peak];
		tmpPeaks[i].reSampledL2Count = rsProfile2[tmpPeaks[i].peak];

		if ((tmpPeaks[i].reSampledL1Count>=config.minCount))
		{
			(chrom->l1PeakNum) ++;
		}
	}

	if (chrom->l1PeakNum<=0)
	{
		free(profile1);
		free(profile2);
		free(rsProfile1);
		free(rsProfile2);
		free(tmpPeaks);
		return 0;
	}

	chrom->l1Peaks = (PEAK_STRUCT *)malloc((chrom->l1PeakNum)*sizeof(PEAK_STRUCT));

	if (!(chrom->l1Peaks))
	{
		printf("not enough memory!\n");
		return -1;
	}

	chrom->l1PeakNum = 0;

	for (i=0;i<tmpPeakNum;i++)
	{
		if (tmpPeaks[i].reSampledL1Count>=config.minCount)
		{
			memcpy(chrom->l1Peaks+chrom->l1PeakNum, tmpPeaks+i,sizeof(PEAK_STRUCT));
			(chrom->l1PeakNum) ++;
		}
	}

	tmpPeakNum = GetLocalMaxima(rsProfile2,profileLen,tmpPeaks,minDist,config.minCount);

	chrom->l2PeakNum = 0;

	for (i=0;i<tmpPeakNum;i++)
	{
		tmpPeaks[i].l1Count = profile2[tmpPeaks[i].peak];
		tmpPeaks[i].l2Count = profile1[tmpPeaks[i].peak];
		tmpPeaks[i].reSampledL1Count = rsProfile2[tmpPeaks[i].peak];
		tmpPeaks[i].reSampledL2Count = rsProfile1[tmpPeaks[i].peak];

		if ((tmpPeaks[i].reSampledL1Count>=config.minCount))
		{
			(chrom->l2PeakNum) ++;
		}
	}

	if (chrom->l2PeakNum<=0)
	{
		free(profile1);
		free(profile2);
		free(rsProfile1);
		free(rsProfile2);
		free(tmpPeaks);
		return 0;
	}

	chrom->l2Peaks = (PEAK_STRUCT *)malloc((chrom->l2PeakNum)*sizeof(PEAK_STRUCT));

	if (!(chrom->l2Peaks))
	{
		printf("not enough memory!\n");
		return -1;
	}

	chrom->l2PeakNum = 0;

	for (i=0;i<tmpPeakNum;i++)
	{
		if (tmpPeaks[i].reSampledL1Count>=config.minCount)
		{
			memcpy(chrom->l2Peaks+chrom->l2PeakNum, tmpPeaks+i,sizeof(PEAK_STRUCT));
			(chrom->l2PeakNum) ++;
		}
	}

	free(profile1);
	free(profile2);
	free(rsProfile1);
	free(rsProfile2);
	free(tmpPeaks);

	return chrom->l1PeakNum;
}

//GetPeaksInOneChrom2: strand-sensitive mode: get peak location from tags for one chromosome
int GetPeaksInOneChrom2(CHROM_INFO *chrom, double l1Ratio, double l2Ratio, int *maxL1Count, int *maxL2Count, CONFIG_STRUCT config)
{
	int i;
	int *profile1,*profile2, *rsProfile1, *rsProfile2, *diffProfile1, *diffProfile2;
	int profileLen = (chrom->chromSize)/config.movingStep+1;
	int tmpIndex;
	PEAK_STRUCT *tmpPeaks;
	int tmpPeakNum;
	int minDist; 

	//allocate memory for the profile
	profile1 = (int *)malloc(profileLen*sizeof(int));
	profile2 = (int *)malloc(profileLen*sizeof(int));
	rsProfile1 = (int *)malloc(profileLen*sizeof(int));
	rsProfile2 = (int *)malloc(profileLen*sizeof(int));
	diffProfile1 = (int *)malloc(profileLen*sizeof(int));
	diffProfile2 = (int *)malloc(profileLen*sizeof(int));

	if ((!profile1)||(!profile2)||(!rsProfile1)||(!rsProfile2)||(!diffProfile1)||(!diffProfile2))
	{
		printf("not enough memory!\n");
		return -1;
	}

	memset(profile1, 0, profileLen*sizeof(int));
	memset(profile2, 0, profileLen*sizeof(int));
	memset(rsProfile1, 0, profileLen*sizeof(int));
	memset(rsProfile2, 0, profileLen*sizeof(int));
	memset(diffProfile1, 0, profileLen*sizeof(int));
	memset(diffProfile2, 0, profileLen*sizeof(int));

	//generate profile, assign reads from both ChIP and control libraries to the profile. Reads are re-sampled based on noise rate

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

		if (rand()>RAND_MAX*l1Ratio)
		{
			continue;
		}

		rsProfile1[tmpIndex]++;

		tmpIndex = (chrom->l1PosTags[i])/config.movingStep;

		tmpIndex = tmpIndex<0?0:tmpIndex;
		tmpIndex = tmpIndex>=profileLen?profileLen-1:tmpIndex;

		diffProfile1[tmpIndex]++;
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

		if (rand()>RAND_MAX*l1Ratio)
		{
			continue;
		}

		rsProfile1[tmpIndex]++;

		tmpIndex = (chrom->l1NegTags[i])/config.movingStep;

		tmpIndex = tmpIndex<0?0:tmpIndex;
		tmpIndex = tmpIndex>=profileLen?profileLen-1:tmpIndex;
		diffProfile1[tmpIndex]--;
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

		if (rand()>RAND_MAX*l2Ratio)
		{
			continue;
		}

		rsProfile2[tmpIndex]++;
		
		tmpIndex = (chrom->l2PosTags[i])/config.movingStep;

		tmpIndex = tmpIndex<0?0:tmpIndex;
		tmpIndex = tmpIndex>=profileLen?profileLen-1:tmpIndex;
		diffProfile2[tmpIndex]++;
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

		if (rand()>RAND_MAX*l2Ratio)
		{
			continue;
		}

		rsProfile2[tmpIndex]++;

		tmpIndex = (chrom->l2NegTags[i])/config.movingStep;

		tmpIndex = tmpIndex<0?0:tmpIndex;
		tmpIndex = tmpIndex>=profileLen?profileLen-1:tmpIndex;
		diffProfile2[tmpIndex]--;
	}

	SmoothProfile(profile1, profileLen, config.slidingWinSize/config.movingStep/2);

	SmoothProfile(profile2, profileLen, config.slidingWinSize/config.movingStep/2);

	SmoothProfile(rsProfile1, profileLen, config.slidingWinSize/config.movingStep/2);

	SmoothProfile(rsProfile2, profileLen, config.slidingWinSize/config.movingStep/2);

	SmoothProfile(diffProfile1, profileLen, TRANSITION_SIZE/config.movingStep/2);

	SmoothProfile(diffProfile2, profileLen, TRANSITION_SIZE/config.movingStep/2);

	for (i=0;i<profileLen;i++)
	{
		if (profile1[i]>*maxL1Count)
		{
			*maxL1Count = profile1[i];
		}

		if (profile2[i]>*maxL2Count)
		{
			*maxL2Count = profile2[i];
		}
	}

	//find peaks from profile

	tmpPeaks = (PEAK_STRUCT *)malloc(((chrom->chromSize)/config.movingStep+1)*sizeof(PEAK_STRUCT));

	if (!tmpPeaks)
	{
		printf("not enough memory!\n");
		return -1;
	}

	minDist = config.slidingWinSize/config.movingStep+1;

	tmpPeakNum = GetTransition(diffProfile1,profileLen,tmpPeaks,minDist);

	chrom->l1PeakNum = 0;

	
	for (i=0;i<tmpPeakNum;i++)
	{
		tmpPeaks[i].l1Count = profile1[tmpPeaks[i].peak];
		tmpPeaks[i].l2Count = profile2[tmpPeaks[i].peak];
		tmpPeaks[i].reSampledL1Count = rsProfile1[tmpPeaks[i].peak];
		tmpPeaks[i].reSampledL2Count = rsProfile2[tmpPeaks[i].peak];

		if ((tmpPeaks[i].reSampledL1Count>=config.minCount))
		{
			(chrom->l1PeakNum) ++;
		}
	}

	if (chrom->l1PeakNum<=0)
	{
		free(profile1);
		free(profile2);
		free(rsProfile1);
		free(rsProfile2);
		free(diffProfile1);
		free(diffProfile2);
		free(tmpPeaks);
		return 0;
	}

	chrom->l1Peaks = (PEAK_STRUCT *)malloc((chrom->l1PeakNum)*sizeof(PEAK_STRUCT));

	if (!(chrom->l1Peaks))
	{
		free(profile1);
		free(profile2);
		free(rsProfile1);
		free(rsProfile2);
		free(diffProfile1);
		free(diffProfile2);
		printf("not enough memory!\n");
		return -1;
	}

	chrom->l1PeakNum = 0;

	for (i=0;i<tmpPeakNum;i++)
	{
		if (tmpPeaks[i].reSampledL1Count>=config.minCount)
		{
			memcpy(chrom->l1Peaks+chrom->l1PeakNum, tmpPeaks+i,sizeof(PEAK_STRUCT));
			(chrom->l1PeakNum) ++;
		}
	}

	tmpPeakNum = GetTransition(diffProfile2,profileLen,tmpPeaks,minDist);

	chrom->l2PeakNum = 0;

	for (i=0;i<tmpPeakNum;i++)
	{
		tmpPeaks[i].l1Count = profile2[tmpPeaks[i].peak];
		tmpPeaks[i].l2Count = profile1[tmpPeaks[i].peak];
		tmpPeaks[i].reSampledL1Count = rsProfile2[tmpPeaks[i].peak];
		tmpPeaks[i].reSampledL2Count = rsProfile1[tmpPeaks[i].peak];

		if ((tmpPeaks[i].reSampledL1Count>=config.minCount))
		{
			(chrom->l2PeakNum) ++;
		}
	}

	if (chrom->l2PeakNum<=0)
	{
		free(profile1);
		free(profile2);
		free(rsProfile1);
		free(rsProfile2);
		free(diffProfile1);
		free(diffProfile2);
		free(tmpPeaks);
		return 0;
	}

	chrom->l2Peaks = (PEAK_STRUCT *)malloc((chrom->l2PeakNum)*sizeof(PEAK_STRUCT));

	if (!(chrom->l2Peaks))
	{
		free(profile1);
		free(profile2);
		free(rsProfile1);
		free(rsProfile2);
		free(diffProfile1);
		free(diffProfile2);
		printf("not enough memory!\n");
		return -1;
	}

	chrom->l2PeakNum = 0;

	for (i=0;i<tmpPeakNum;i++)
	{
		if (tmpPeaks[i].reSampledL1Count>=config.minCount)
		{
			memcpy(chrom->l2Peaks+chrom->l2PeakNum, tmpPeaks+i,sizeof(PEAK_STRUCT));
			(chrom->l2PeakNum) ++;
		}
	}

	free(profile1);
	free(profile2);
	free(rsProfile1);
	free(rsProfile2);
	free(diffProfile1);
	free(diffProfile2);
	free(tmpPeaks);

	return chrom->l1PeakNum;
}

//GetLocalMaxima: find the local maxima in the profile
int GetLocalMaxima(int *profile, int profileSize, PEAK_STRUCT *peaks, int minDist, int minCount)
{
	int tmpStart, tmpEnd, tmpStart1,tmpEnd1;
	int i,j;
	int peakNum;
	int flag;

	peakNum = 0;

	tmpStart = -1;
	tmpEnd = -1;

	for (i=1;i<profileSize;i++)
	{
		if (profile[i]>=minCount)
		{
			if (profile[i]>profile[i-1])
			{
				tmpStart = i;
				tmpEnd = i;
			}

			if (profile[i]==profile[i-1])
			{
				tmpEnd = i;
			}

			if (profile[i]<profile[i-1])
			{
				if ((tmpStart==-1)||(tmpEnd==-1))
				{
					continue;
				}
				
				peaks[peakNum].peak = (tmpStart+tmpEnd)/2;

				tmpStart1 = tmpStart-minDist<0?0:tmpStart-minDist;
				tmpEnd1 = tmpEnd+minDist>profileSize-1?profileSize-1:tmpEnd+minDist;

				flag = 1;
				for (j=tmpStart1;j<tmpStart;j++)
				{
					if (profile[j]>=profile[peaks[peakNum].peak])
					{
						flag = 0;
						break;
					}
				}

				if (!flag)
				{
					tmpStart = -1;
					tmpEnd = -1;
					continue;
				}

				for (j=tmpEnd1;j>tmpEnd;j--)
				{
					if (profile[j]>profile[peaks[peakNum].peak])
					{
						flag = 0;
						break;
					}
				}

				if (!flag)
				{
					tmpStart = -1;
					tmpEnd = -1;
					continue;
				}

				peakNum++;

				tmpStart = -1;
				tmpEnd = -1;
			}
		}
		else
		{
			if ((tmpStart==-1)||(tmpEnd==-1))
			{
				continue;
			}
				
			peaks[peakNum].peak = (tmpStart+tmpEnd)/2;

			tmpStart1 = tmpStart-minDist<0?0:tmpStart-minDist;
			tmpEnd1 = tmpEnd+minDist>profileSize-1?profileSize-1:tmpEnd+minDist;

			flag = 1;
			for (j=tmpStart1;j<tmpStart;j++)
			{
				if (profile[j]>=profile[peaks[peakNum].peak])
				{
					flag = 0;
					break;
				}
			}

			if (!flag)
			{
				tmpStart = -1;
				tmpEnd = -1;
				continue;
			}

			for (j=tmpEnd1;j>tmpEnd;j--)
			{
				if (profile[j]>profile[peaks[peakNum].peak])
				{
					flag = 0;
					break;
				}
			}

			if (!flag)
			{
				tmpStart = -1;
				tmpEnd = -1;
				continue;
			}

			peakNum++;

			tmpStart = -1;
			tmpEnd = -1;
		}
	}
	
	return peakNum;
}

//GetTransition: get the transition from sense reads to anti-sense reads 
int GetTransition(int *profile, int profileSize, PEAK_STRUCT *peaks, int minDist)
{
	int tmpStart, tmpEnd;
	int i,j,k;
	int peakNum, tmpPeakNum;
	int *isFiltered;
	int sum;
	PEAK_STRUCT *tmpPeak;

	peakNum = 0;

	for (i=1;i<profileSize;i++)
	{
		if (profile[i]>0)
		{
			tmpStart = i;

			for (j=i+1;j<profileSize;j++)
			{
				if (profile[j]>0)
				{
					tmpStart = j;
				}

				if (profile[j]<0)
				{
					break;
				}
			}

			if (j-tmpStart<minDist)
			{
				if ((tmpStart+j)%2)
				{
					if (profile[tmpStart]>-profile[j])
					{
						peaks[peakNum].peak = (tmpStart+j)/2;
					}
					else
					{
						peaks[peakNum].peak = (tmpStart+j)/2+1;
					}
				}
				else
				{
					peaks[peakNum].peak = (tmpStart+j)/2;
				}

				if ((peakNum>0)&&(peaks[peakNum].peak-peaks[peakNum-1].peak<minDist))
				{
					sum = 0;

					for (k=peaks[peakNum-1].peak+1;k<peaks[peakNum].peak;k++)
					{
						sum += profile[k];
					}

					if (sum>0)
					{
						memcpy(peaks+peakNum-1,peaks+peakNum, sizeof(PEAK_STRUCT));
					}
				}
				else
				{
					peakNum++;
				}
			}

			i=j;
		}
	}

	if (peakNum==0)
	{
		return 0;
	}

	isFiltered = (int *)malloc(peakNum*sizeof(int));

	if (!isFiltered)
	{
		printf("not enough memory!\n");
		return -1;
	}

	for (i=0;i<peakNum;i++)
	{
		isFiltered[i] = 0;

		tmpStart = peaks[i].peak-minDist-1>=0?peaks[i].peak-minDist-1:0;
		
		tmpEnd = peaks[i].peak+minDist+1<profileSize?peaks[i].peak+minDist+1:profileSize-1;

		sum = 0;

		for (j=tmpStart;j<peaks[i].peak;j++)
		{
			sum += profile[j];
		}

		if (sum<=0)
		{
			isFiltered[i] = 1;
		}

		sum = 0;

		for (j=peaks[i].peak+1;j<=tmpEnd;j++)
		{
			sum += profile[j];
		}

		if (sum>=0)
		{
			isFiltered[i] = 1;
		}
	}

	tmpPeak = (PEAK_STRUCT *)malloc(peakNum*sizeof(PEAK_STRUCT));

	if (!tmpPeak)
	{
		printf("not enough memory!\n");
		return -1;
	}
	
	memcpy(tmpPeak, peaks, peakNum*sizeof(PEAK_STRUCT));

	tmpPeakNum = peakNum;
	peakNum = 0;

	for (i=0;i<tmpPeakNum;i++)
	{
		if (!isFiltered[i])
		{
			memcpy(peaks+peakNum, tmpPeak+i,sizeof(PEAK_STRUCT));
			peakNum++;
		}
	}

	free(tmpPeak);
	free(isFiltered);

	return peakNum;
}