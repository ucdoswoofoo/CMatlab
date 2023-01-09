/*
 * File: spmak2.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 08-Dec-2022 23:55:19
 */

 /* Include Files */
#include "spmak2.h"

/* Function Definitions */

/*
 * Arguments    : const double knot
 *                const double coef
 *                struct0_T *sp
 * Return Type  : void
 */
void spmak2(Matrix *knot, Matrix *coef, stusp1* sp)
{
	int number;

	number      = coef->column;   //可能需要改为去第1行的元素个数
	sp->number  = number;
	sp->knots   = knot;
	sp->coefs   = coef;

}


