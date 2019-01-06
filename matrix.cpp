#include "const-function-declerations.h"			   
void MatrixInverse_BSxBS(double(*a)[BS][2]){
	int i, j, k;
	double p[2], q[2], power;
	double temp[BS][BS][2];

	for (k = 0; k < BS; k++){
		p[0] = a[k][k][0];
		p[1] = a[k][k][1];
		power = sqr(p[0]) + sqr(p[1]);
		a[k][k][0] = 1.0;
		a[k][k][1] = 0.0;

		for (j = 0; j < BS; j++){
			temp[k][j][0] = (a[k][j][0] * p[0] + a[k][j][1] * p[1]) / power;
			temp[k][j][1] = (a[k][j][1] * p[0] - a[k][j][0] * p[1]) / power;
			a[k][j][0] = temp[k][j][0];
			a[k][j][1] = temp[k][j][1];
		}

		for (i = 0; i < BS; i++){
			if (i != k){
				q[0] = a[i][k][0];
				q[1] = a[i][k][1];
				a[i][k][0] = 0.0;
				a[i][k][1] = 0.0;

				for (j = 0; j < BS; j++){
					temp[i][j][0] = q[0] * a[k][j][0] - q[1] * a[k][j][1];
					temp[i][j][1] = q[1] * a[k][j][0] + q[0] * a[k][j][1];
					a[i][j][0] -= temp[i][j][0];
					a[i][j][1] -= temp[i][j][1];
				}
			}
		}
	}
}

void MatrixHermite_TSxBS(double(*m1)[BS][2], double(*m2)[TS][2]){
	int tx, rx;

	for (tx = 0; tx < BS; tx++){
		for (rx = 0; rx < TS; rx++){
			m2[tx][rx][0] = m1[rx][tx][0];
			m2[tx][rx][1] = -m1[rx][tx][1];
		}
	}
}

void MatrixTranspose_TSxBS(double(*m1)[BS][2], double(*m2)[TS][2]){
	int tx, rx;

	for (tx = 0; tx < BS; tx++){
		for (rx = 0; rx < TS; rx++){
			m2[tx][rx][0] = m1[rx][tx][0];
			m2[tx][rx][1] = m1[rx][tx][1];
		}
	}
}


void VectorMulVectorToScalar_M(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < M; i++) {
		s[0] += (v1[i][0] * v2[i][0] - v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] + v1[i][1] * v2[i][0]);
	}

}

void MatrixAddMatrixToMatrix_BSxBS(double(*a)[BS][2], double(*b)[BS][2], double(*c)[BS][2]){
	int i, j, k;

	for (i = 0; i < BS; i++){
		for (j = 0; j < BS; j++){
			for (k = 0; k < 2; k++) {
				c[i][j][k] = a[i][j][k] + b[i][j][k];
			}
		}
	}

}

void MatrixSubMatrixToMatrix_TS(double(*a)[TS][2], double(*b)[TS][2], double(*c)[TS][2]){
	int i, j, k;

	for (i = 0; i < TS; i++)
	for (j = 0; j < TS; j++)
	for (k = 0; k < 2; k++)
		c[i][j][k] = a[i][j][k] - b[i][j][k];

}

void MatrixMulMatrixToMatrix_BSxTS(double(*a)[TS][2], double(*b)[BS][2], double(*c)[BS][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < BS; i++)
	for (j = 0; j < BS; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < TS; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
	}

}

void MatrixMulMatrixToMatrix_BSxBSXBSxTS(double(*a)[BS][2], double(*b)[TS][2], double(*c)[TS][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < BS; i++)
	for (j = 0; j < TS; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < BS; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
	}
}

void MatrixMulMatrixToMatrix_BSxBSXBSxCL(double(*a)[BS][2], double(*b)[CL][2], double(*c)[CL][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < BS; i++)
	for (j = 0; j < CL; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < BS; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
	}
}

void MatrixMulVectorToVector_TSxTS(double(*a)[TS][2], double(*v1)[2], double(*v2)[2]){
	int i, j;
	double sum[2];


	for (i = 0; i < TS; i++) {
		sum[0] = 0.0;
		sum[1] = 0.0;
		for (j = 0; j < TS; j++) {
			sum[0] += (a[i][j][0] * v1[j][0] - a[i][j][1] * v1[j][1]);
			sum[1] += (a[i][j][0] * v1[j][1] + a[i][j][1] * v1[j][0]);
		}
		v2[i][0] = sum[0];
		v2[i][1] = sum[1];
	}
}

void MatrixMulVectorToVector_BSxBS(double(*a)[BS][2], double(*v1)[2], double(*v2)[2]){
	int i, j;
	double sum[2];

	for (i = 0; i < BS; i++) {
		sum[0] = 0.0;
		sum[1] = 0.0;
		for (j = 0; j < BS; j++) {
			sum[0] += (a[i][j][0] * v1[j][0] - a[i][j][1] * v1[j][1]);
			sum[1] += (a[i][j][0] * v1[j][1] + a[i][j][1] * v1[j][0]);
		}
		v2[i][0] = sum[0];
		v2[i][1] = sum[1];
	}
}

void ScalarMulMatrixToMatrix_TSxTS(double *s, double(*a)[TS][2], double(*c)[TS][2]){
	int i, j;

	for (i = 0; i < TS; i++)
	for (j = 0; j < TS; j++) {
		c[i][j][0] = s[0] * a[i][j][0] - s[1] * a[i][j][1];
		c[i][j][1] = s[0] * a[i][j][1] + s[1] * a[i][j][0];
	}
}

void VectorMulVectorToScalar_TS(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < TS; i++) {
		s[0] += (v1[i][0] * v2[i][0] - v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] + v1[i][1] * v2[i][0]);
	}
}

void VectorMulVectorToScalarH_TS(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < TS; i++) {
		s[0] += (v1[i][0] * v2[i][0] + v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] - v1[i][1] * v2[i][0]);
	}
}

void VectorMulVectorToScalarH_Np(double(*v1)[2], double(*v2)[2], double *s) {
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < Np; i++) {
		s[0] += (v1[i][0] * v2[i][0] + v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] - v1[i][1] * v2[i][0]);
	}
}


void VectorAddVectorToVector_TS(double(*v1)[2], double(*v2)[2], double(*v3)[2]){
	int i;

	for (i = 0; i < TS; i++) {
		v3[i][0] = v1[i][0] + v2[i][0];
		v3[i][1] = v1[i][1] + v2[i][1];
	}
}

void VectorSubVectorToVector_TS(double(*v1)[2], double(*v2)[2], double(*v3)[2]){
	int i;

	for (i = 0; i < TS; i++) {
		v3[i][0] = v1[i][0] - v2[i][0];
		v3[i][1] = v1[i][1] - v2[i][1];
	}
}

void VectorMulMatrixToVector_TSxTS(double(*v1)[2], double(*a)[TS][2], double(*v2)[2]){
	int i, j;

	for (i = 0; i < TS; i++) {
		v2[i][0] = 0.0;
		v2[i][1] = 0.0;
		for (j = 0; j < TS; j++) {
			v2[i][0] += v1[j][0] * a[j][i][0] - v1[j][1] * a[j][i][1];
			v2[i][1] += v1[j][0] * a[j][i][1] + v1[j][1] * a[j][i][0];
		}
	}

}

void VectorMulVectorToMatrix_TSxTS(double(*v1)[2], double(*v2)[2], double(*a)[TS][2]){
	int i, j;

	for (i = 0; i < TS; i++)
	for (j = 0; j < TS; j++) {
		a[i][j][0] = v1[i][0] * v2[j][0] - v1[i][1] * v2[j][1];
		a[i][j][1] = v1[i][0] * v2[j][1] + v1[i][1] * v2[j][0];
	}
}

void VectorConjugate_TS(double(*v1)[2], double(*v2)[2]){
	int j;

	for (j = 0; j < TS; j++) {
		v2[j][0] = v1[j][0];
		v2[j][1] = -v1[j][1];
	}
}

void VectorConjugate_M(double(*v1)[2], double(*v2)[2]){
	int j;

	for (j = 0; j < M; j++) {
		v2[j][0] = v1[j][0];
		v2[j][1] = -v1[j][1];
	}
}

void ScalarMulVectorToVector_TS(double *s, double(*v1)[2], double(*v2)[2]){
	int i;

	for (i = 0; i < TS; i++) {
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

void VectorDistance_TS(double(*v1)[2], double(*v2)[2], double s){	/* Doesnt work */
	int i;
	double ith_distance[TS];
	s = 0;

	for (i = 0; i < TS; i++) {
		ith_distance[i] = ScalarDistance(v1[i], v2[i]);
		s += ith_distance[i];
	}

}
