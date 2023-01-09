#ifndef GSMOOTH2_H
#define GSMOOTH2_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "ZXMat.h"


//function[gg, g2] = gsmooth2(sp)

extern void gsmooth2(Matrix gg, Matrix g2, const stusp1* sp);

#endif