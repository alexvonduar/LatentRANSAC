#ifndef HTOOLS_H
#define HTOOLS_H

#include "MathFunctions.h"

namespace HTools
{
	void computeDataMatrix(double* data_matrix, unsigned int num_points, double* points);
	void crossprod(double *out, const double *a, const double *b, unsigned int st);
	inline bool isClockWise(double x1, double y1, double x2, double y2, double x3, double y3);
}
#endif