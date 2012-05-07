/**************************************************************************************************
math_api.c: source file: the APIs for the mathematical calculation in CCAT
Author: Xu Han (hanxu@gis.a-star.edu.sg)
Date: 15/07/2008
***************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "math_api.h"

//Overlap: Caculate the overlap of two regions. Return 0 if no overlap
int Overlap(int start1, int end1, int start2, int end2)
{
	int start,end;

	start = start1>start2?start1:start2;
	end = end1>end2?end2:end1;

	if (end<start)
	{
		return 0;
	}
	else
	{
		return end-start+1;
	}
}

//GetMaxPos: Get the maximum position in a profile.
int GetMaxPos(double *profile, int len)
{
	int i, start,end,index;
	double maxValue;

	index = 0;
	maxValue = profile[index];

	for (i=0;i<len;i++)
	{
		if (profile[i]>maxValue)
		{
			index = i;
			maxValue = profile[i];
		}
	}

	for (start = index-1;start>=0;start--)
	{
		if (profile[start]<maxValue-0.0000000001)
		{
			break;
		}
	}

	for (end = index+1; end<len;end++)
	{
		if (profile[end]<maxValue-0.0000000001)
		{
			break;
		}
	}

	return (start+end)/2;
}

//SmoothProfile: Smooth a profile by summing over a sliding window
void SmoothProfile(int *array, int len, int halfWinSize)
{
	int i;
	int sum;
	int *tmp;

	if (len<=halfWinSize)
	{
		return;
	}

	tmp = (int *)malloc(len*sizeof(int));

	sum = 0;
	for (i=0;i<halfWinSize*2+1;i++)
	{
		sum += array[i];
	}

	for (i=0;i<=halfWinSize;i++)
	{
		tmp[i] = sum;
	}

	for (i=halfWinSize+1;i<len-halfWinSize;i++)
	{
		tmp[i] = tmp[i-1]-array[i-halfWinSize-1]+array[i+halfWinSize];
	}

	for (i=len-halfWinSize;i<len;i++)
	{
		tmp[i] = tmp[i-1];
	}

	memcpy(array,tmp,len*sizeof(int));

	free(tmp);
}

//Quicksort an array in real values, in descending order
void QuicksortF(double *a, int lo, int hi)
{
	int i=lo, j=hi;
	double x=a[(lo+hi)/2];
	double h;

	if (hi<lo)
	{
		return;
	}

    //  partition
    while (i<=j)
    {    
		while ((a[i]>x)&&(i<=j))
		{
			i++;
		}
		while ((a[j]<x)&&(i<=j))
		{
			j--;
		}
        if (i<=j)
        {
			h = a[i];
			a[i] = a[j];
			a[j] = h;
            i++; j--;
        }
    } 

    //  recursion
    if (lo<j) QuicksortF(a, lo, j);
    if (i<hi) QuicksortF(a, i, hi);
}

//Quicksort an array in integer values, in ascending order
void QuicksortI(int *a, int lo, int hi)
{
	int i=lo, j=hi;
	int x=a[(lo+hi)/2];
	int h;

	if (hi<lo)
	{
		return;
	}

    //  partition
    while (i<=j)
    {    
		while ((a[i]<x)&&(i<=j))
		{
			i++;
		}
		while ((a[j]>x)&&(i<=j))
		{
			j--;
		}
        if (i<=j)
        {
			h = a[i];
			a[i] = a[j];
			a[j] = h;
            i++; j--;
        }
    } 

    //  recursion
    if (lo<j) QuicksortI(a, lo, j);
    if (i<hi) QuicksortI(a, i, hi);
}

//BTreeSearchingF: Searching value in array, which was organized in ascending order previously
int  bTreeSearchingF(double value, double *a, int lo, int hi)
{
	if (value>=a[lo])
	{
		return lo;
	}

	if (value<=a[hi])
	{
		return hi;
	}

	if (hi-lo<=1)
	{
		return lo;
	}

	if (value<=a[(lo+hi)/2])
	{
		return bTreeSearchingF(value, a, (lo+hi)/2, hi);
	}
	else
	{
		return bTreeSearchingF(value, a,  lo, (lo+hi)/2);
	}
}

//BTreeSearchingI: Searching value in array, which was organized in ascending order previously
int  bTreeSearchingI(int value, int *a, int lo, int hi)
{
	if (value<=a[lo])
	{
		return lo;
	}

	if (value>=a[hi])
	{
		return hi;
	}

	if (hi-lo<=1)
	{
		return lo;
	}

	if (value>=a[(lo+hi)/2])
	{
		return bTreeSearchingI(value, a, (lo+hi)/2, hi);
	}
	else
	{
		return bTreeSearchingI(value, a,  lo, (lo+hi)/2);
	}
}

//ComputeMean: compute mean of the array
double ComputeMedian(double *array, int len)
{
	QuicksortF(array, 0, len-1);

	if (len%2)
	{
		return array[len/2];
	}
	else
	{
		return (array[len/2-1]+array[len/2])/2;
	}
}

//fit a negative binomial distribution
int FitNegBinomDist(int *hist, int histLen, double *m, double *k, double k_min, double k_max)
{
	int i,j;
	double sum;
	double *histNorm;
	double tmp_k_min, tmp_k_max, tmpK, step;
	double error,minError;
	double bestFit, prevBestFit;

	if ((k_max<=k_min)||(k_min<=0))
	{
		return -1;
	}

	histNorm = (double *)malloc(histLen*sizeof(double));

	sum = 0.0;
	for (i=0;i<histLen;i++)
	{
		sum += (double)(hist[i]);
	}

	if (sum<1.0)
	{
		return -1;
	}

	*m = 0.0;

	for (i=0;i<histLen;i++)
	{
		histNorm[i] = (double)(hist[i])/sum;
		*m += histNorm[i]*(double)i;
	}

	tmp_k_min = k_min;
	tmp_k_max = k_max;
	prevBestFit = DBL_MAX;

	while(1)
	{
		minError = DBL_MAX;

		step = log(tmp_k_max/tmp_k_min)/9;

		for (i=0;i<10;i++)
		{
			tmpK = tmp_k_min*exp(step*(double)i);

			error = 0.0;

			sum = 0.0;

			for (j=0;j<histLen;j++)
			{
				error += (exp(gammln(tmpK+j)-gammln(tmpK)-gammln((double)j+1)+tmpK*log(tmpK/(tmpK+*m))+(double)j*log(*m/(tmpK+*m)))-histNorm[j])
					*(exp(gammln(tmpK+j)-gammln(tmpK)-gammln((double)j+1)+tmpK*log(tmpK/(tmpK+*m))+(double)j*log(*m/(tmpK+*m)))-histNorm[j]);
			}
			error = sqrt(error);

			if (error<minError)
			{
				bestFit = tmpK;
				minError = error;
			}
		}
		
		if (fabs(log(prevBestFit/bestFit))<0.01)
		{
			*k = bestFit;
			break;
		}
		
		tmp_k_min = bestFit/exp(step)<k_min?k_min:bestFit/exp(step);
		tmp_k_max = bestFit*exp(step)>k_max?k_max:bestFit*exp(step);

		prevBestFit = bestFit;
	}

	free(histNorm);

	return 1;
}

//Compute log-gamma function
double gammln(double xx)
{
	double x,y,tmp, ser;
	double cof[6] = {76.18009172947146, -86.50532032941677, 
		24.01409824083091, -1.231739572450155,
		0.1208650973866179e-2, -0.5395239384953e-5};
	int j;

	if (xx<0.0)
	{
		return -1.0;
	}

	if (xx==0.0)
	{
		return 0.0;
	}

	y=x=xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;

	for (j=0;j<=5;j++)
	{
		ser += cof[j]/++y;
	}

	return -tmp+log(2.5066282746310005*ser/x);
}