/**************************************************************************************************
file_io.c: source file, read input data from file or save output data to files
Author: Xu Han (hanxu@gis.a-star.edu.sg)
Date: 25/11/2008
***************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "ccat.h"

//ReadConfigFile: read the configuration file into config. Return the number of parameters read.
int ReadConfigFile(char *fileName, CONFIG_STRUCT *config); 

//GetChromInfo: Read chromosome name and length, return the number of chromosomes read
int GetChromInfo(char *fileName, CHROM_INFO *chroms);

//GetTags: Read the tags from bed file of L1 and L2 to the structure of chroms. Return the total number of tags read.
int GetTags(char *l1FileName, char *l2FileName, CHROM_INFO *chroms, int chromNum);

//OutputFiles: output the results
int OutputFiles(CHROM_INFO *chroms, int chromNum, char *projectName, CONFIG_STRUCT config);

//ChromToIndex: Convert the string of chromosome name to the chromosome index, return the index
int ChromToIndex(char *chromName, CHROM_INFO *chrom, int chromNum);

//IndexToChrom: convert the chromosome index to string of chromosome name. 
char *IndexToChrom(int index, char *chromName, int len, CHROM_INFO *chroms, int chromNum);

//Quicksort peaks by foldChange
void QuicksortPeaks(PEAK_STRUCT *a, int lo, int hi);

//ReadConfigFile: read the configuration file into config. Return the number of parameters read.
int ReadConfigFile(char *fileName, CONFIG_STRUCT *config)
{
	FILE *fh;
	char tmpStr[1000];
	int count = 0;

	//default setting
	config->fragmentSize = 200;
	config->slidingWinSize = 300;
	config->movingStep = 10;
	config->isStrandSensitiveMode = 1;
	config->minCount = 3;
	config->outputNum = 100000;
	config->randomSeed = 123456;
	config->minScore = 5.0;
	config->bootstrapPass = 50;

	//read configurations

	fh = (FILE *)fopen(fileName, "r");

	if (!fh)
	{
		return 0;
	}
	
	fscanf(fh, "%s", tmpStr);

	while (!feof(fh))
	{
		//experimental settings
		if (!strcmp(tmpStr, "fragmentSize"))
		{
			fscanf(fh, "%s", tmpStr);
			config->fragmentSize = atoi(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}

		else if (!strcmp(tmpStr, "isStrandSensitiveMode"))
		{
			fscanf(fh, "%s", tmpStr);
			config->isStrandSensitiveMode = atoi(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}

		else if (!strcmp(tmpStr, "movingStep"))
		{
			fscanf(fh, "%s", tmpStr);
			config->movingStep = atoi(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}

		else if (!strcmp(tmpStr, "slidingWinSize"))
		{
			fscanf(fh, "%s", tmpStr);
			config->slidingWinSize = atoi(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}

		else if (!strcmp(tmpStr, "minCount"))
		{
			fscanf(fh, "%s", tmpStr);
			config->minCount = atoi(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}

		else if (!strcmp(tmpStr, "outputNum"))
		{
			fscanf(fh, "%s", tmpStr);
			config->outputNum = atoi(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}
		
		else if (!strcmp(tmpStr, "randomSeed"))
		{
			fscanf(fh, "%s", tmpStr);
			config->randomSeed = atoi(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}

		else if (!strcmp(tmpStr, "minScore"))
		{
			fscanf(fh, "%s", tmpStr);
			config->minScore = atof(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}

		else if (!strcmp(tmpStr, "bootstrapPass"))
		{
			fscanf(fh, "%s", tmpStr);
			config->bootstrapPass = atoi(tmpStr);
			fscanf(fh, "%s", tmpStr);
			count++;
			continue;
		}

		else
		{
			fscanf(fh, "%s", tmpStr);
			fscanf(fh, "%s", tmpStr);
		}
	}

	//validate settings

	if (config->fragmentSize<0)
	{
		printf("config:fragmentSize not valid! Default setting chosen.\n");
		config->fragmentSize = 200;
	}

	if (config->slidingWinSize<1)
	{
		printf("config:slidingWinSize not valid! Default setting chosen.\n");
		config->slidingWinSize = 300;
	}

	if (config->movingStep<1)
	{
		printf("config:movingStep not valid! Default setting chosen.\n");
		config->movingStep = 10;
	}

	if (config->outputNum<0)
	{
		printf("config:outputNum not valid! Default setting chosen.\n");
		config->outputNum = 100000;
	}

	if (config->minCount<=1)
	{
		printf("config:minCount not valid! Default setting chosen.\n");
		config->minCount = 2;
	}

	if (config->minScore<0)
	{
		printf("config:minScore not valid! Default setting chosen.\n");
		config->minScore = 5.0;
	}

	if (config->bootstrapPass<1)
	{
		printf("config:bootstrapPass not valid! Default setting chosen.\n");
		config->bootstrapPass = 50;
	}

	fclose(fh);
	return count;
}

//GetChromInfo: Read chromosome name and length, return the number of chromosomes read
int GetChromInfo(char *fileName, CHROM_INFO *chroms)
{
	FILE *fh;
	char tmpStr[1000];
	int i;
	int chromNum;

	fh = (FILE *)fopen(fileName, "r");

	if (!fh)
	{
		return -1;
	}

	chromNum = 0;

	fscanf(fh, "%s", tmpStr);

	while (!feof(fh))
	{
		fscanf(fh, "%s", tmpStr);
		fscanf(fh, "%s", tmpStr);
		chromNum++;
	}

	fclose(fh);

	memset(chroms,0, chromNum*sizeof(CHROM_INFO));

	fh = (FILE *)fopen(fileName, "r");

	if (!fh)
	{
		return -1;
	}

	for (i=0;i<chromNum;i++)
	{
		fscanf(fh, "%s", chroms[i].chromName);
		fscanf(fh, "%s", tmpStr);
		chroms[i].chromIndex = i;
		chroms[i].chromSize = atoi(tmpStr);
	}

	fclose(fh);
	return chromNum;
}

//GetTags: Read the tags from bed file of L1 and L2 to the structure of chroms. Return the total number of tags read.
int GetTags(char *l1FileName, char *l2FileName, CHROM_INFO *chroms, int chromNum)
{
	FILE *fh;
	char tmpChrom[1000],tmpStrand[10],tmpStr[1000];
	int tmpIndex, tmpStart, tmpEnd;
	int i;
	int l1TagNum,l2TagNum;

	fh = fopen(l1FileName, "r");

	if (!fh)
	{
		printf("Cannot open %s or file format error!\n", l1FileName);
		return 0;
	}

	//Get the total number of L1 tags for each chrom.

	fscanf(fh, "%s", tmpChrom); //chrom name

	for (i=0;i<chromNum;i++)
	{
		chroms[i].l2PosTagNum=0;
		chroms[i].l2NegTagNum=0;
	}

	l1TagNum = 0;

	while (!feof(fh))
	{
		fscanf(fh, "%d",&tmpStart); //start
		fscanf(fh, "%d",&tmpEnd);   //end
		fscanf(fh, "%s",tmpStr);    //unused
		fscanf(fh, "%s", tmpStr);   //unused
		fscanf(fh, "%s", tmpStrand);//strand

		tmpIndex = ChromToIndex(tmpChrom, chroms, chromNum);

		if (tmpIndex<0)
		{
			fscanf(fh, "%s", tmpChrom);
			continue;
		}

		if (tmpStrand[0]=='+')
		{
			l1TagNum++;
			chroms[tmpIndex].l1PosTagNum++;
		}
		else if (tmpStrand[0]=='-')
		{
			l1TagNum++;
			chroms[tmpIndex].l1NegTagNum++;
		}
		else
		{
			printf("%s format error!\n", l1FileName);
			fclose(fh);
			return -1;
		}
		fscanf(fh, "%s", tmpChrom); //chrom name
	}

	fclose(fh);

	//Get the total number of tags in L2 for each chrom

	fh = fopen(l2FileName, "r");

	if (!fh)
	{
		printf("Cannot open %s or file format error!\n", l2FileName);
		return 0;
	}

	fscanf(fh, "%s", tmpChrom);

	for (i=0;i<chromNum;i++)
	{
		chroms[i].l2PosTagNum=0;
		chroms[i].l2NegTagNum=0;
	}

	l2TagNum = 0;

	while (!feof(fh))
	{
		fscanf(fh, "%d",&tmpStart);
		fscanf(fh, "%d",&tmpEnd);
		fscanf(fh, "%s",tmpStr);
		fscanf(fh, "%s", tmpStr);
		fscanf(fh, "%s", tmpStrand);

		tmpIndex = ChromToIndex(tmpChrom, chroms, chromNum);

		if (tmpIndex<0)
		{
			fscanf(fh, "%s", tmpChrom);
			continue;
		}

		if (tmpStrand[0]=='+')
		{
			l2TagNum++;
			chroms[tmpIndex].l2PosTagNum++;
		}
		else if (tmpStrand[0]=='-')
		{
			l2TagNum++;
			chroms[tmpIndex].l2NegTagNum++;
		}
		else
		{
			printf("%s format error!\n", l2FileName);
			fclose(fh);
			return -1;
		}
		fscanf(fh, "%s", tmpChrom); //chrom name
	}

	fclose(fh);

	//Allocate memory

	for (i=0;i<chromNum;i++)
	{
		if (chroms[i].l1PosTagNum>0)
		{
			chroms[i].l1PosTags = (int *)malloc(chroms[i].l1PosTagNum*sizeof(int));
		}

		if (chroms[i].l1NegTagNum>0)
		{
			chroms[i].l1NegTags = (int *)malloc(chroms[i].l1NegTagNum*sizeof(int));
		}

		if (chroms[i].l2PosTagNum>0)
		{
			chroms[i].l2PosTags = (int *)malloc(chroms[i].l2PosTagNum*sizeof(int));
		}

		if (chroms[i].l2NegTagNum>0)
		{
			chroms[i].l2NegTags = (int *)malloc(chroms[i].l2NegTagNum*sizeof(int));
		}
	}

	//Read L1 tag information to the structure of CHROM_INFO

	fh = fopen(l1FileName, "r");

	for (i=0;i<chromNum;i++)
	{
		chroms[i].l1PosTagNum=0;
		chroms[i].l1NegTagNum=0;
	}

	fscanf(fh, "%s", tmpChrom);

	while (!feof(fh))
	{
		fscanf(fh, "%d",&tmpStart);
		fscanf(fh, "%d",&tmpEnd);
		fscanf(fh, "%s",tmpStr);
		fscanf(fh, "%s", tmpStr);
		fscanf(fh, "%s", tmpStrand);

		tmpIndex = ChromToIndex(tmpChrom, chroms, chromNum);

		if (tmpIndex<0)
		{
			fscanf(fh, "%s", tmpChrom);
			continue;
		}

		if (tmpStrand[0]=='+')
		{
			chroms[tmpIndex].l1PosTags[chroms[tmpIndex].l1PosTagNum] = tmpStart;
			chroms[tmpIndex].l1PosTagNum++;
		}
		else if (tmpStrand[0]=='-')
		{
			chroms[tmpIndex].l1NegTags[chroms[tmpIndex].l1NegTagNum] = tmpEnd;
			chroms[tmpIndex].l1NegTagNum++;
		}
		else
		{
			printf("%s format error!\n", l1FileName);
			fclose(fh);
			return -1;
		}

		fscanf(fh, "%s", tmpChrom); //chrom name
	}

	fclose(fh);

	//Read L2 tag information to the structure of CHROM_INFO

	fh = fopen(l2FileName, "r");

	for (i=0;i<chromNum;i++)
	{
		chroms[i].l2PosTagNum=0;
		chroms[i].l2NegTagNum=0;
	}

	fscanf(fh, "%s", tmpChrom);

	while (!feof(fh))
	{
		fscanf(fh, "%d",&tmpStart);
		fscanf(fh, "%d",&tmpEnd);
		fscanf(fh, "%s",tmpStr);
		fscanf(fh, "%s", tmpStr);
		fscanf(fh, "%s", tmpStrand);

		tmpIndex = ChromToIndex(tmpChrom, chroms, chromNum);

		if (tmpIndex<0)
		{
			fscanf(fh, "%s", tmpChrom);
			continue;
		}

		if (tmpStrand[0]=='+')
		{
			chroms[tmpIndex].l2PosTags[chroms[tmpIndex].l2PosTagNum] = tmpStart;
			chroms[tmpIndex].l2PosTagNum++;
		}
		else if (tmpStrand[0]=='-')
		{
			chroms[tmpIndex].l2NegTags[chroms[tmpIndex].l2NegTagNum] = tmpEnd;
			chroms[tmpIndex].l2NegTagNum++;
		}
		else
		{
			printf("%s format error!\n", l2FileName);
			fclose(fh);
			return -1;
		}

		fscanf(fh, "%s", tmpChrom); //chrom name
	}

	fclose(fh);

	printf("%d tags in L1, %d tags in L2.\n", l1TagNum, l2TagNum);

	return l1TagNum+l2TagNum;
}

//OutputFiles: output the results
int OutputFiles(CHROM_INFO *chroms, int chromNum, char *projectName, CONFIG_STRUCT config)
{
	PEAK_STRUCT *peaks;
	int peakNum;
	int i,j;
	FILE *fh;
	char fileName[1000];
	int startChrom;
	int sumL1, sumL2;

	peakNum = 0;
	sumL1 = 0;
	sumL2 = 0;

	for (i=0;i<chromNum;i++)
	{
		peakNum+=chroms[i].l1PeakNum;
		sumL1 += chroms[i].l1PosTagNum+chroms[i].l1NegTagNum;
		sumL2 += chroms[i].l2PosTagNum+chroms[i].l2NegTagNum;
	}

	if (peakNum<=0)
	{
		return 0;
	}

	peaks = (PEAK_STRUCT *)malloc(peakNum*sizeof(PEAK_STRUCT));

	if (!peaks)
	{
		printf("not enough memory!\n");
		return -1;
	}

	peakNum = 0;

	for (i=0;i<chromNum;i++)
	{
		for (j=0;j<chroms[i].l1PeakNum;j++)
		{
			memcpy(peaks+peakNum, chroms[i].l1Peaks+j, sizeof(PEAK_STRUCT));
			peaks[peakNum].chromIndex = i;
			peakNum++;
		}
	}

	QuicksortPeaks(peaks, 0, peakNum-1);

	//output significant peaks
	sprintf(fileName, "%s.significant.peak", projectName);

	fh = fopen(fileName, "w");

	if (!fh)
	{
		free(peaks);
		printf("%s cannot open!\n", fileName);
		return -1;
	}

	for (i=0;i<peakNum;i++)
	{
		if (peaks[i].isSignificant)
		{
			fprintf(fh, "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",
				chroms[peaks[i].chromIndex].chromName,
				peaks[i].peak,
				peaks[i].start,
				peaks[i].end,
				peaks[i].l1Count,
				peaks[i].l2Count,
				peaks[i].foldChange,
				peaks[i].qValue);
		}
	}

	fclose(fh);

	//output top peaks
	sprintf(fileName, "%s.top%d.peak", projectName, config.outputNum);

	fh = fopen(fileName, "w");

	if (!fh)
	{
		free(peaks);
		printf("%s cannot open!\n", fileName);
		return -1;
	}

	for (i=0;i<peakNum;i++)
	{
		if (i<config.outputNum)
		{
			fprintf(fh, "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",
				chroms[peaks[i].chromIndex].chromName,
				peaks[i].peak,
				peaks[i].start,
				peaks[i].end,
				peaks[i].l1Count,
				peaks[i].l2Count,
				peaks[i].foldChange,
				peaks[i].qValue);
		}
	}

	fclose(fh);

	//output regions

	for (startChrom=0;startChrom<chromNum;startChrom++)
	{
		if (chroms[startChrom].l1PeakNum>0)
		{
			break;
		}
	}

	memcpy(peaks,chroms[startChrom].l1Peaks, sizeof(PEAK_STRUCT));
	peaks[0].chromIndex = startChrom;
	peakNum = 1;

	for (i=startChrom;i<chromNum;i++)
	{
		for (j=0;j<chroms[i].l1PeakNum;j++)
		{
			if ((i==startChrom)&&(j==0))
			{
				continue;
			}

			if ((i==peaks[peakNum-1].chromIndex)&&(chroms[i].l1Peaks[j].start==peaks[peakNum-1].start))
			{
				if (chroms[i].l1Peaks[j].foldChange>peaks[peakNum-1].foldChange)
				{
					memcpy(peaks+peakNum-1,chroms[i].l1Peaks+j, sizeof(PEAK_STRUCT));
					peaks[peakNum-1].chromIndex = i;
				}
			}
			else
			{
				memcpy(peaks+peakNum,chroms[i].l1Peaks+j, sizeof(PEAK_STRUCT));
				peaks[peakNum].chromIndex = i;
				peakNum++;
			}
		}
	}

	QuicksortPeaks(peaks, 0, peakNum-1);

	sprintf(fileName, "%s.significant.region", projectName);

	fh = fopen(fileName, "w");

	if (!fh)
	{
		free(peaks);
		printf("%s cannot open!\n", fileName);
		return -1;
	}

	for (i=0;i<peakNum;i++)
	{
		if (peaks[i].isSignificant)
		{
			fprintf(fh, "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",
				chroms[peaks[i].chromIndex].chromName,
				peaks[i].peak,
				peaks[i].start,
				peaks[i].end,
				peaks[i].l1Count,
				peaks[i].l2Count,
				peaks[i].foldChange,
				peaks[i].qValue);
		}
	}

	fclose(fh);

	free(peaks);

	return 1;
}

//ChromToIndex: Convert the string of chromosome name to the chromosome index, return the index
int ChromToIndex(char *chromName, CHROM_INFO *chroms, int chromNum)
{
	int i;

	for (i=0;i<chromNum;i++)
	{
		if (!strcmp(chromName, chroms[i].chromName))
		{
			break;
		}
	}

	if (i<chromNum)
	{
		return i;
	}
	else
	{
		return -1;
	}
}

//IndexToChrom: convert the chromosome index to string of chromosome name. 
char *IndexToChrom(int index, char *chromName, int len, CHROM_INFO *chroms, int chromNum)
{
	if ((index<0)||(index>=chromNum))
	{
		return 0;
	}

	if ((int)(strlen(chroms[index].chromName))>=len)
	{
		return 0;
	}

	strcpy(chromName, chroms[index].chromName);

	return chromName;
}

//Quicksort peaks by foldChange
void QuicksortPeaks(PEAK_STRUCT *a, int lo, int hi)
{
	int i=lo, j=hi;
	double x=a[(lo+hi)/2].foldChange;
	PEAK_STRUCT h;

	if (hi<lo)
	{
		return;
	}

    //  partition
    while (i<=j)
    {    
		while ((a[i].foldChange>x)&&(i<=j))
		{
			i++;
		}
		while ((a[j].foldChange<x)&&(i<=j))
		{
			j--;
		}
        if (i<=j)
        {
			memcpy(&h, &(a[i]), sizeof(PEAK_STRUCT));
			memcpy(&(a[i]), &(a[j]), sizeof(PEAK_STRUCT));
			memcpy(&(a[j]), &h, sizeof(PEAK_STRUCT));
            i++; j--;
        }
    } 

    //  recursion
    if (lo<j) QuicksortPeaks(a, lo, j);
    if (i<hi) QuicksortPeaks(a, i, hi);
}

