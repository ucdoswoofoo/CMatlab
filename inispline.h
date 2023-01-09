#ifndef INISPLINE_H
#define INISPLINE_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "ZXMat.h"


//function sp = inispline(segdata, segwrefq, eps)

extern void inispline(stusp1* sp, Matrix segdata, Matrix segwrefq, const double eps);

#endif



