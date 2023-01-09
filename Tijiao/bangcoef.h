#ifndef BANGCOEF_H
#define BANGCOEF_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "ZXMat.h"
#include "fnder2.h"


extern void  bangcoef(Matrix* tknot, Matrix* talpha, Matrix* ttemp, const stusp1* sp, double Vmax, double Amax, double Jmax, double RVmax_A, double RVmax_C, double RAmax, double RJmax);

#endif




