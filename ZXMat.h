#ifndef ZXMAT_H
#define ZXMAT_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include "matrix.h"



struct ZXsp1
{
	Matrix *knots;
	int number;
	Matrix *coefs;
};
typedef struct ZXsp1 stusp1;

int ZXUnique(double  A[], int nCnt);
void ZXSort(double ArrVal[], int nCnt);
void ZXdiff(double  PrimVal[], double  DiffVal[], int nCnt);
int ZXIsEmptyMat(Matrix _Mat);


#endif
