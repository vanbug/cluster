/**************************************************************************************************
math_api.h: header file: the APIs for the mathematical calculation in CCAT
Author: Xu Han (hanxu@gis.a-star.edu.sg)
Date: 15/07/2008
***************************************************************************************************/

#ifndef MATH_API_H
#define MATH_API_H
#endif

typedef struct
{
	double maxX;
	double minX;
	double maxY;
	double minY;
	double *density;
	int binNumX;
	int binNumY;
}HISTOGRAM_STRUCT_2D;

//Overlap: Caculate the overlap of two regions. Return 0 if no overlap
int Overlap(int start1, int end1, int start2, int end2);

//GetMaxPos: Get the maximum position in a profile.
int GetMaxPos(double *profile, int len);

//SmoothProfile: Smooth a profile by summing over a sliding window
void SmoothProfile(int *array, int len, int halfWinSize);

//Quicksort an array in real values.
void QuicksortF(double *a, int lo, int hi);

//Quicksort an array in integer values, in ascending order
void QuicksortI(int *a, int lo, int hi);

//BTreeSearchingF: Searching value in array, which was organized in descending order previously
int  bTreeSearchingF(double value, double *a, int lo, int hi);

//BTreeSearchingI: Searching value in array, which was organized in ascending order previously
int  bTreeSearchingI(int value, int *a, int lo, int hi);

//ComputeMean: compute mean of the array
double ComputeMedian(double *array, int len);

//fit a negative binomial distribution
int FitNegBinomDist(int *hist, int histLen, double *m, double *k, double k_min, double k_max);

//Compute log-gamma function
double gammln(double xx);