#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "Spline.h"
#include <stdlib.h>
#include "Intersectionpoint.h"
#include <stdio.h>

void get_array_values(double* array, int size) {
	for (int i = 0; i < size; i++)
		scanf("%lf", &array[i]);
}

int main() {
	printf("Welcome to the Spline search program!\n");
	printf("Please input the count of a first array of points\n");
	int size1;
	scanf("%d", &size1);
	double* X1 = create_array(size1);
	double* Y1 = create_array(size1);
	printf("Then input the array of X separated by a space and press Enter\n");
	get_array_values(X1, size1);
	printf("Then input the array of Y separated by a space and press Enter\n");
	get_array_values(Y1, size1);
	printf("Now input the count of a second array of points\n");
	int size2;
	scanf("%d", &size2);
	double* X2 = create_array(size2);
	double* Y2 = create_array(size2);
	printf("Then input the array of X separated by a space and press Enter\n");
	get_array_values(X2, size2);
	printf("Then input the array of Y separated by a space and press Enter\n");
	get_array_values(Y2, size2);

	Polinom* spline1 = create_spline(X1, Y1, size1);
	Polinom* spline2 = create_spline(X2, Y2, size2);

	print_spline(spline1, size1);
	print_spline(spline2, size2);

	find_intersection_point(X1, size1, X2, size2, spline1, spline2);

	free(X1); free(Y1); free(X2); free(Y2);

	free(spline1); free(spline2);
}