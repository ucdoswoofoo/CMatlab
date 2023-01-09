#ifndef FEATUREPOINT_H
#define FEATUREPOINT_H

	/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "ZXMat.h"


extern void   featurepoint(Matrix feature, Matrix k, Matrix d, Matrix* segwrefq, const double dmax);

#endif


