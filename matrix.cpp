#include "const.h"
void MatrixInverse(double(*a)[Tx][2]){
	int i, j, k;
	double p[2], q[2], power;
	double temp[Tx][Tx][2];

	for (k = 0; k < Tx; k++){
		p[0] = a[k][k][0];
		p[1] = a[k][k][1];
		power = sqr(p[0]) + sqr(p[1]);
		a[k][k][0] = 1.0;
		a[k][k][1] = 0.0;

		for (j = 0; j < Tx; j++){
			temp[k][j][0] = (a[k][j][0] * p[0] + a[k][j][1] * p[1]) / power;
			temp[k][j][1] = (a[k][j][1] * p[0] - a[k][j][0] * p[1]) / power;
			a[k][j][0] = temp[k][j][0];
			a[k][j][1] = temp[k][j][1];
		}

		for (i = 0; i < Tx; i++){
			if (i != k){
				q[0] = a[i][k][0];
				q[1] = a[i][k][1];
				a[i][k][0] = 0.0;
				a[i][k][1] = 0.0;

				for (j = 0; j < Tx; j++){
					temp[i][j][0] = q[0] * a[k][j][0] - q[1] * a[k][j][1];
					temp[i][j][1] = q[1] * a[k][j][0] + q[0] * a[k][j][1];
					a[i][j][0] -= temp[i][j][0];
					a[i][j][1] -= temp[i][j][1];
				}
			}
		}
	}
}

void MatrixHermite(double(*m1)[Tx][2], double(*m2)[Rx][2]){
	int tx, rx;

	for (tx = 0; tx < Tx; tx++){
		for (rx = 0; rx < Rx; rx++){
			m2[tx][rx][0] = m1[rx][tx][0];
			m2[tx][rx][1] = -m1[rx][tx][1];
		}
	}
}

void MatrixTranspose(double(*m1)[Tx][2], double(*m2)[Rx][2]){
	int tx, rx;

	for (tx = 0; tx < Tx; tx++){
		for (rx = 0; rx < Rx; rx++){
			m2[tx][rx][0] = m1[rx][tx][0];
			m2[tx][rx][1] = m1[rx][tx][1];
		}
	}
}


void VectorMulVectorToScalar2(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < M; i++) {
		s[0] += (v1[i][0] * v2[i][0] - v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] + v1[i][1] * v2[i][0]);
	}

}

void MatrixAddMatrixToMatrix(double(*a)[Tx][2], double(*b)[Tx][2], double(*c)[Tx][2]){
	int i, j, k;

	for (i = 0; i < Tx; i++){
		for (j = 0; j < Tx; j++){
			for (k = 0; k < 2; k++) {
				c[i][j][k] = a[i][j][k] + b[i][j][k];
			}
		}
	}

}

void MatrixSubMatrixToMatrix(double(*a)[Rx][2], double(*b)[Rx][2], double(*c)[Rx][2]){
	int i, j, k;

	for (i = 0; i < Rx; i++)
	for (j = 0; j < Rx; j++)
	for (k = 0; k < 2; k++)
		c[i][j][k] = a[i][j][k] - b[i][j][k];

}

void MatrixMulMatrixToMatrix(double(*a)[Rx][2], double(*b)[Tx][2], double(*c)[Tx][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < Tx; i++)
	for (j = 0; j < Tx; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < Rx; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
	}

}

void MatrixMulMatrixToMatrix2(double(*a)[Tx][2], double(*b)[Rx][2], double(*c)[Rx][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < Tx; i++)
	for (j = 0; j < Rx; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < Tx; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
	}
}


void MatrixMulVectorToVector(double(*a)[Rx][2], double(*v1)[2], double(*v2)[2]){
	int i, j;
	double sum[2];


	for (i = 0; i < Rx; i++) {
		sum[0] = 0.0;
		sum[1] = 0.0;
		for (j = 0; j < Rx; j++) {
			sum[0] += (a[i][j][0] * v1[j][0] - a[i][j][1] * v1[j][1]);
			sum[1] += (a[i][j][0] * v1[j][1] + a[i][j][1] * v1[j][0]);
		}
		v2[i][0] = sum[0];
		v2[i][1] = sum[1];
	}
}

void MatrixMulVectorToVector2(double(*a)[Tx][2], double(*v1)[2], double(*v2)[2]){
	int i, j;
	double sum[2];

	for (i = 0; i < Tx; i++) {
		sum[0] = 0.0;
		sum[1] = 0.0;
		for (j = 0; j < Tx; j++) {
			sum[0] += (a[i][j][0] * v1[j][0] - a[i][j][1] * v1[j][1]);
			sum[1] += (a[i][j][0] * v1[j][1] + a[i][j][1] * v1[j][0]);
		}
		v2[i][0] = sum[0];
		v2[i][1] = sum[1];
	}
}

void ScalarMulMatrixToMatrix(double *s, double(*a)[Rx][2], double(*c)[Rx][2]){
	int i, j;

	for (i = 0; i < Rx; i++)
	for (j = 0; j < Rx; j++) {
		c[i][j][0] = s[0] * a[i][j][0] - s[1] * a[i][j][1];
		c[i][j][1] = s[0] * a[i][j][1] + s[1] * a[i][j][0];
	}
}

void VectorMulVectorToScalar(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < Rx; i++) {
		s[0] += (v1[i][0] * v2[i][0] - v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] + v1[i][1] * v2[i][0]);
	}
}

void VectorMulVectorToScalarH(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < Rx; i++) {
		s[0] += (v1[i][0] * v2[i][0] + v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] - v1[i][1] * v2[i][0]);
	}
}


void VectorAddVectorToVector(double(*v1)[2], double(*v2)[2], double(*v3)[2]){
	int i;

	for (i = 0; i < Rx; i++) {
		v3[i][0] = v1[i][0] + v2[i][0];
		v3[i][1] = v1[i][1] + v2[i][1];
	}
}

void VectorSubVectorToVector(double(*v1)[2], double(*v2)[2], double(*v3)[2]){
	int i;

	for (i = 0; i < Rx; i++) {
		v3[i][0] = v1[i][0] - v2[i][0];
		v3[i][1] = v1[i][1] - v2[i][1];
	}
}

void VectorMulMatrixToVector(double(*v1)[2], double(*a)[Rx][2], double(*v2)[2]){
	int i, j;

	for (i = 0; i < Rx; i++) {
		v2[i][0] = 0.0;
		v2[i][1] = 0.0;
		for (j = 0; j < Rx; j++) {
			v2[i][0] += v1[j][0] * a[j][i][0] - v1[j][1] * a[j][i][1];
			v2[i][1] += v1[j][0] * a[j][i][1] + v1[j][1] * a[j][i][0];
		}
	}

}

void VectorMulVectorToMatrix(double(*v1)[2], double(*v2)[2], double(*a)[Rx][2]){
	int i, j;

	for (i = 0; i < Rx; i++)
	for (j = 0; j < Rx; j++) {
		a[i][j][0] = v1[i][0] * v2[j][0] - v1[i][1] * v2[j][1];
		a[i][j][1] = v1[i][0] * v2[j][1] + v1[i][1] * v2[j][0];
	}
}

void VectorConjugate(double(*v1)[2], double(*v2)[2]){
	int j;

	for (j = 0; j < Rx; j++) {
		v2[j][0] = v1[j][0];
		v2[j][1] = -v1[j][1];
	}
}

void VectorConjugate2(double(*v1)[2], double(*v2)[2]){
	int j;

	for (j = 0; j < M; j++) {
		v2[j][0] = v1[j][0];
		v2[j][1] = -v1[j][1];
	}
}
void ScalarMulVectorToVector(double *s, double(*v1)[2], double(*v2)[2]){
	int i;

	for (i = 0; i < Rx; i++) {
		v2[i][0] = s[0] * v1[i][0] - s[1] * v1[i][1];
		v2[i][1] = s[0] * v1[i][1] + s[1] * v1[i][0];
	}
}

void ScalarDivScalarToScalar(double *s1, double *s2, double *s3){
	double pow;

	pow = s2[0] * s2[0] + s2[1] * s2[1];
	if (pow == 0.0)
		printf("CANNOT DIVIDE\n");
	s3[0] = (s1[0] * s2[0] + s1[1] * s2[1]) / pow;
	s3[1] = (s1[1] * s2[0] - s1[0] * s2[1]) / pow;

}

double ScalarDistance(double *s1, double *s2){
	double s3;
	s3 = sqrt(pow((s1[0] - s2[0]), 2) + pow((s1[1] - s2[1]), 2));
	return s3;
}

void VectorDistance(double(*v1)[2], double(*v2)[2], double s){	/* Doesnt work */
	int i;
	double ith_distance[Rx];
	s = 0;

	for (i = 0; i < Rx; i++) {
		ith_distance[i] = ScalarDistance(v1[i], v2[i]);
		s += ith_distance[i];
	}

}
