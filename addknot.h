#ifndef ADDKNOT_H
#define ADDKNOT_H

/* Include Files */
#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ZXMat.h"
#include "matrix.h"


extern void addknot(stusp1* nsp, Matrix* superknots, const stusp1* sp, Matrix* pkd, double wchorderr, double wdataerr, double eps,Matrix* pkpar,Matrix* dchord, double rd, double reps);

#endif
