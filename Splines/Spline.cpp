#define _CRT_SECURE_NO_WARNINGS
#include "Spline.h"
#include <stdlib.h>
#include <stdio.h>

double* create_array(int size) {
	double* array = (double*)calloc(size, sizeof(double));
	return array;
}

double* get_h(double* x, int size) {
	double* h = create_array(size - 1);
	for (int i = 0; i < size - 1; i++) {
		h[i] = x[i + 1] - x[i];
	}
	return h;
}

double* progon(double* h, double* y, int size) {
	double* gamma = create_array(size);

	double* a = create_array(size);
	double* b = create_array(size);
	double* c = create_array(size);
	double* f = create_array(size);
	for (int i = 1; i < size - 1; i++) {
		a[i] = h[i - 1] / 6.0;
		b[i] = (h[i - 1] + h[i]) / 3.0;
		c[i] = h[i + 1] / 6.0;
		f[i] = (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1];
	}
	a[1] = 0.0;
	c[size - 2] = 0.0;
	double* A = create_array(size);
	double* B = create_array(size);
	A[1] = -c[1] / b[1];
	B[1] = f[1] / b[1];
	for (int i = 2; i < size - 1; i++) {
		A[i] = (-c[i] / (b[i] + a[i] * A[i - 1]));
		B[i] = (f[i] - a[i] * B[i - 1]) / (b[i] + a[i] * A[i - 1]);
	}
	gamma[size - 1] = 0;
	gamma[size - 2] = B[size - 2];
	for (int i = size - 3; i >= 1; i--) {
		gamma[i] = A[i] * gamma[i + 1] + B[i];
	}
	gamma[0] = 0;
	free(a); free(b); free(c); free(f);
	free(A); free(B);
	return gamma;
}

Polinom* get_spline_coefficients(Polinom* spline, double* x, double* y, double* h, double* gamma, int size) {
	for (int i = 0; i < size - 1; i++) {
		spline[i].third = (gamma[i + 1] - gamma[i]) / (6 * h[i]);
		spline[i].second = 3 * (x[i + 1] * gamma[i] - gamma[i + 1] * x[i]) / (6 * h[i]);
		spline[i].first = (y[i + 1] - y[i]) / h[i] + \
			(gamma[i] * (h[i] * h[i] - 3 * x[i + 1] * x[i + 1]) - gamma[i + 1] * (h[i] * h[i] * x[i] - 3 * x[i] * x[i])) / (6 * h[i]);
		spline[i].zero = (y[i] * x[i + 1] - y[i + 1] * x[i]) / h[i] + \
			(gamma[i] * (-h[i] * h[i] * x[i + 1] + x[i + 1] * x[i + 1] * x[i + 1]) + gamma[i + 1] * (h[i] * h[i] * x[i] - x[i] * x[i] * x[i])) / (6 * h[i]);
	}
	return spline;
}

Polinom* create_spline(double* x, double* y, int size) {
	double* h = get_h(x, size);
	double* gamma = progon(h, y, size);

	Polinom* spline = (Polinom*)malloc(size * sizeof(Polinom));
	spline = get_spline_coefficients(spline, x, y, h, gamma, size);

	free(h); free(gamma);
	return spline;
}

void print_spline_i(Polinom spline_i) {
	printf("%lfx^3+%lfx^2+%lfx+%lf\n", spline_i.third, spline_i.second, spline_i.first, spline_i.zero);
}

void print_spline(Polinom* spline, int size) {
	printf("Got a spline\n");
	for (int i = 0; i < size - 1; i++) {
		printf("[x_%d;x_%d]:  ",i,i+1);
		print_spline_i(spline[i]);
	}
}