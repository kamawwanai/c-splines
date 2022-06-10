#define _CRT_SECURE_NO_WARNINGS
#include "Intersectionpoint.h"
#include "stdlib.h"
#include "Spline.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

double function(Polinom f, double x) {
	double res = f.third * x * x * x + f.second * x * x + f.first * x + f.zero;
	return res;
}

bool method_horda(double a, double b, double eps, Polinom f, Polinom ff) {
	double c = a - (function(f, a) / (function(f, b) - function(f, a)) * (b - a));
	int it = 1;
	while (abs(function(f, c) >= eps)) {
		if (it > 10000) {
			return false;
		}
		it++;
		if (function(f, a) * function(f, c) < 0)
			b = c;
		else
			a = c;
		c= a - (function(f, a) / (function(f, b) - function(f, a)) * (b - a));
	}
	printf("Find intersection point:(%lf,%lf)", c, function(ff, c));
	return true;
}

void qsortRecursive(double* mas, int size) {
    int i = 0;
    int j = size - 1;
    int mid = mas[size / 2];
    do {
        while (mas[i] < mid) {
            i++;
        }
        while (mas[j] > mid) {
            j--;
        }
        if (i <= j) {
            int tmp = mas[i];
            mas[i] = mas[j];
            mas[j] = tmp;
            i++;
            j--;
        }
    } while (i <= j);
    if (j > 0) {
        qsortRecursive(mas, j + 1);
    }
    if (i < size) {
        qsortRecursive(&mas[i], size - i);
    }
}

bool if_in_mas(double* mas, double size, double element) {
    for (int i = 0; i < size; i++) {
        if (mas[i] == element)
            return true;
    }
    return false;
}

void find_intersection_point(double* X1, int size1, double* X2, int size2, Polinom* spline1, Polinom* spline2) {
    double* mas = (double*)calloc(size1 + size2, sizeof(double));
    int size_mas = 0;
    for (int i = 0; i < size1; i++) {
        if (X1[i] < X2[0]||X1[i]>X2[size2-1])
            continue;
        mas[i] = X1[i];
        size_mas++;
    }
    for (int j = 0; j < size2; j++) {
        bool flag = true;
        for (int i = 0; i < size1; i++) {
            if (X2[j] == mas[i]||X2[j]>X1[size1-1]||X2[j]<X1[0]) {
                flag = false;
                break;
            }
        }
        if (flag) {
            mas[size_mas] = X2[j];
            size_mas++;
        }
    }
    qsortRecursive(mas, size_mas);
    if (size_mas == 0) {
        printf("The intersection point cannot be found exactly\n");
        return;
    }
    for (int i = 0; i<size_mas-1; i++) {
        double left = mas[i];
        double right = mas[i + 1];
        int k = 0, k1 = 1; //X1
        int l = 0, l1 = 1; //X2
        while (X1[k]<left || X1[k]>right || X1[k1]<left || X1[k1]>right) {
            k++;
            k1++;
        }
        while (X1[l]<left || X1[l]>right || X1[l1]<left || X1[l1]>right) {
            l++;
            l1++;
        }
        Polinom func;
        func.third = spline1[k].third - spline2[l].third;
        func.second = spline1[k].second - spline2[l].second;
        func.first = spline1[k].first - spline2[l].first;
        func.zero = spline1[k].zero - spline2[l].zero;
        if (method_horda(left, right, 0.001, func, spline1[k])) {
            free(mas);
            return;
        }
    }
    free(mas);
    printf("The intersection point cannot be found exactly\n");
    return;
}