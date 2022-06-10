#pragma once
#include <stdlib.h>

typedef struct polinom {
	double third;
	double second;
	double first;
	double zero;
}Polinom;

double* create_array(int size);

Polinom* create_spline(double* x, double* y, int size);

void print_spline(Polinom* spline, int size);