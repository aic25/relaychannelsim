#include "const-function-declerations.h"
void MatrixInverse_CLxCL(double(*a)[CL][2]){
	int i, j, k;
	double p[2], q[2], power;
	double temp[CL][CL][2];

	for (k = 0; k < CL; k++){
		p[0] = a[k][k][0];
		p[1] = a[k][k][1];
		power = sqr(p[0]) + sqr(p[1]);
		a[k][k][0] = 1.0;
		a[k][k][1] = 0.0;

		for (j = 0; j < CL; j++){
			temp[k][j][0] = (a[k][j][0] * p[0] + a[k][j][1] * p[1]) / power;
			temp[k][j][1] = (a[k][j][1] * p[0] - a[k][j][0] * p[1]) / power;
			a[k][j][0] = temp[k][j][0];
			a[k][j][1] = temp[k][j][1];
		}

		for (i = 0; i < CL; i++){
			if (i != k){
				q[0] = a[i][k][0];
				q[1] = a[i][k][1];
				a[i][k][0] = 0.0;
				a[i][k][1] = 0.0;

				for (j = 0; j < CL; j++){
					temp[i][j][0] = q[0] * a[k][j][0] - q[1] * a[k][j][1];
					temp[i][j][1] = q[1] * a[k][j][0] + q[0] * a[k][j][1];
					a[i][j][0] -= temp[i][j][0];
					a[i][j][1] -= temp[i][j][1];
				}
			}
		}
	}
}

void MatrixHermite_CLxCL(double(*m1)[CL][2], double(*m2)[CL][2]){
	int tx, rx;

	for (tx = 0; tx < CL; tx++){
		for (rx = 0; rx < CL; rx++){
			m2[tx][rx][0] = m1[rx][tx][0];
			m2[tx][rx][1] = -m1[rx][tx][1];
		}
	}
}

void MatrixTranspose_CLxCL(double(*m1)[CL][2], double(*m2)[CL][2]){
	int tx, rx;

	for (tx = 0; tx < CL; tx++){
		for (rx = 0; rx < CL; rx++){
			m2[tx][rx][0] = m1[rx][tx][0];
			m2[tx][rx][1] = m1[rx][tx][1];
		}
	}
}

void MatrixAddMatrixToMatrix_CLxCL(double(*a)[CL][2], double(*b)[CL][2], double(*c)[CL][2]){
	int i, j, k;

	for (i = 0; i < CL; i++){
		for (j = 0; j < CL; j++){
			for (k = 0; k < 2; k++) {
				c[i][j][k] = a[i][j][k] + b[i][j][k];
			}
		}
	}

}

void MatrixSubMatrixToMatrix_CL(double(*a)[CL][2], double(*b)[CL][2], double(*c)[CL][2]){
	int i, j, k;

	for (i = 0; i < CL; i++)
	for (j = 0; j < CL; j++)
	for (k = 0; k < 2; k++)
		c[i][j][k] = a[i][j][k] - b[i][j][k];

}

void MatrixMulMatrixToMatrix_CLxCL(double(*a)[CL][2], double(*b)[CL][2], double(*c)[CL][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < CL; i++)
	for (j = 0; j < CL; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < CL; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
	}

}

void MatrixMulMatrixToMatrix_CLxCLXCLxCL(double(*a)[CL][2], double(*b)[CL][2], double(*c)[CL][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < CL; i++)
	for (j = 0; j < CL; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < CL; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
	}
}


void MatrixMulVectorToVector_CLxCL(double(*a)[CL][2], double(*v1)[2], double(*v2)[2]){
	int i, j;
	double sum[2];


	for (i = 0; i < CL; i++) {
		sum[0] = 0.0;
		sum[1] = 0.0;
		for (j = 0; j < CL; j++) {
			sum[0] += (a[i][j][0] * v1[j][0] - a[i][j][1] * v1[j][1]);
			sum[1] += (a[i][j][0] * v1[j][1] + a[i][j][1] * v1[j][0]);
		}
		v2[i][0] = sum[0];
		v2[i][1] = sum[1];
	}
}

void ScalarMulMatrixToMatrix_CLxCL(double *s, double(*a)[CL][2], double(*c)[CL][2]){
	int i, j;

	for (i = 0; i < CL; i++)
	for (j = 0; j < CL; j++) {
		c[i][j][0] = s[0] * a[i][j][0] - s[1] * a[i][j][1];
		c[i][j][1] = s[0] * a[i][j][1] + s[1] * a[i][j][0];
	}
}

void VectorMulVectorToScalar_CL(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < CL; i++) {
		s[0] += (v1[i][0] * v2[i][0] - v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] + v1[i][1] * v2[i][0]);
	}
}

void VectorMulVectorToScalarH_CL(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < CL; i++) {
		s[0] += (v1[i][0] * v2[i][0] + v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] - v1[i][1] * v2[i][0]);
	}
}


void VectorAddVectorToVector_CL(double(*v1)[2], double(*v2)[2], double(*v3)[2]){
	int i;

	for (i = 0; i < CL; i++) {
		v3[i][0] = v1[i][0] + v2[i][0];
		v3[i][1] = v1[i][1] + v2[i][1];
	}
}

void VectorSubVectorToVector_CL(double(*v1)[2], double(*v2)[2], double(*v3)[2]){
	int i;

	for (i = 0; i < CL; i++) {
		v3[i][0] = v1[i][0] - v2[i][0];
		v3[i][1] = v1[i][1] - v2[i][1];
	}
}

void VectorMulMatrixToVector_CLxCL(double(*v1)[2], double(*a)[CL][2], double(*v2)[2]){
	int i, j;

	for (i = 0; i < CL; i++) {
		v2[i][0] = 0.0;
		v2[i][1] = 0.0;
		for (j = 0; j < CL; j++) {
			v2[i][0] += v1[j][0] * a[j][i][0] - v1[j][1] * a[j][i][1];
			v2[i][1] += v1[j][0] * a[j][i][1] + v1[j][1] * a[j][i][0];
		}
	}

}

void VectorMulVectorToMatrix_CLxCL(double(*v1)[2], double(*v2)[2], double(*a)[CL][2]){
	int i, j;

	for (i = 0; i < CL; i++)
	for (j = 0; j < CL; j++) {
		a[i][j][0] = v1[i][0] * v2[j][0] - v1[i][1] * v2[j][1];
		a[i][j][1] = v1[i][0] * v2[j][1] + v1[i][1] * v2[j][0];
	}
}

void VectorConjugate_CL(double(*v1)[2], double(*v2)[2]){
	int j;

	for (j = 0; j < CL; j++) {
		v2[j][0] = v1[j][0];
		v2[j][1] = -v1[j][1];
	}
}

void ScalarMulVectorToVector_CL(double *s, double(*v1)[2], double(*v2)[2]){
	int i;

	for (i = 0; i < CL; i++) {
		v2[i][0] = s[0] * v1[i][0] - s[1] * v1[i][1];
		v2[i][1] = s[0] * v1[i][1] + s[1] * v1[i][0];
	}
}

void VectorDistance_CL(double(*v1)[2], double(*v2)[2], double s){	/* Doesnt work */
	int i;
	double ith_distance[CL];
	s = 0;

	for (i = 0; i < CL; i++) {
		ith_distance[i] = ScalarDistance(v1[i], v2[i]);
		s += ith_distance[i];
	}

}

double MagnitudeOfDifference_CL(double(*v1)[2], double(*v2)[2]){
	double Difference[CL][2];
	double res;
	int i;
	res = 0;
	for ( i = 0; i < CL; i++)
	{
		Difference[i][0] = v1[i][0] - v2[i][0];
		Difference[i][1] = v1[i][1] - v2[i][1];
		res += sqrt(sqr(Difference[i][0]) + sqr(Difference[i][1]));
	}
	return res;
}

void VectorEqualVector_CL(double(*v1)[2],double(*v2)[2]){
	int i;
	for (i = 0; i < CL; i++){
		v2[i][0] = v1[i][0];
		v2[i][1] = v1[i][1];
	}
}
void ScalarMulScalar(double(*m1), double(*m2), double(*r)){
	r[0] = m1[0] * m2[0] - m1[1] * m2[1];
	r[1] = m1[1] * m2[0] + m1[0] * m2[1];
}

void ScalarMulScalar_H(double(*m1), double(*m2), double(*r)){
	r[0] = m1[0] * m2[0] + m1[1] * m2[1];
	r[1] = m1[1] * m2[0] - m1[0] * m2[1];
}

void sirala(double dizi[], int boy, int icerik[]){
	int i, j;
	double temp;

	for (i = 0; i < boy; i++){
		icerik[i] = i;
	}
	for (i = 0; i < boy; i++){
		for (j = 0; j<boy - 1; j++){
			if (dizi[j]>dizi[j + 1]){
				temp = dizi[j + 1];
				dizi[j + 1] = dizi[j];
				dizi[j] = temp;
				temp = icerik[j + 1];
				icerik[j + 1] = icerik[j];
				icerik[j] = temp;
			}
		}
	}
}