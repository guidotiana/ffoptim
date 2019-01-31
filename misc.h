/*
 * misc.h
 *
 *  Created on: Apr 22, 2009
 *      Author: guido
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "struct.h"

double Dist(struct vector a, struct vector b);
double Angle(struct vector B, struct vector C, struct vector D);
struct vector CrossProd(struct vector A, struct vector B);

