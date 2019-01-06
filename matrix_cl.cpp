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

void MatrixInverse_4x4(double(*a)[4][2]){
	int i, j, k;
	double p[2], q[2], power;
	double temp[4][4][2];

	for (k = 0; k < 4; k++){
		p[0] = a[k][k][0];
		p[1] = a[k][k][1];
		power = sqr(p[0]) + sqr(p[1]);
		a[k][k][0] = 1.0;
		a[k][k][1] = 0.0;

		for (j = 0; j < 4; j++){
			temp[k][j][0] = (a[k][j][0] * p[0] + a[k][j][1] * p[1]) / power;
			temp[k][j][1] = (a[k][j][1] * p[0] - a[k][j][0] * p[1]) / power;
			a[k][j][0] = temp[k][j][0];
			a[k][j][1] = temp[k][j][1];
		}

		for (i = 0; i < 4; i++){
			if (i != k){
				q[0] = a[i][k][0];
				q[1] = a[i][k][1];
				a[i][k][0] = 0.0;
				a[i][k][1] = 0.0;

				for (j = 0; j < 4; j++){
					temp[i][j][0] = q[0] * a[k][j][0] - q[1] * a[k][j][1];
					temp[i][j][1] = q[1] * a[k][j][0] + q[0] * a[k][j][1];
					a[i][j][0] -= temp[i][j][0];
					a[i][j][1] -= temp[i][j][1];
				}
			}
		}
	}
}

void MatrixHermite_4x4(double(*m1)[4][2], double(*m2)[4][2]){
	int tx, rx;

	for (tx = 0; tx < 4; tx++){
		for (rx = 0; rx < 4; rx++){
			m2[tx][rx][0] = m1[rx][tx][0];
			m2[tx][rx][1] = -m1[rx][tx][1];
		}
	}
}

void MatrixTranspose_4x4(double(*m1)[4][2], double(*m2)[4][2]){
	int tx, rx;

	for (tx = 0; tx < 4; tx++){
		for (rx = 0; rx < 4; rx++){
			m2[tx][rx][0] = m1[rx][tx][0];
			m2[tx][rx][1] = m1[rx][tx][1];
		}
	}
}

void MatrixAddMatrixToMatrix_4x4(double(*a)[4][2], double(*b)[4][2], double(*c)[4][2]){
	int i, j, k;

	for (i = 0; i < 4; i++){
		for (j = 0; j < 4; j++){
			for (k = 0; k < 2; k++) {
				c[i][j][k] = a[i][j][k] + b[i][j][k];
			}
		}
	}

}

void MatrixSubMatrixToMatrix_4(double(*a)[4][2], double(*b)[4][2], double(*c)[4][2]){
	int i, j, k;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 2; k++)
				c[i][j][k] = a[i][j][k] - b[i][j][k];

}

void MatrixMulMatrixToMatrix_4x4(double(*a)[4][2], double(*b)[4][2], double(*c)[4][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < 4; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
		}

}

void MatrixMulMatrixToMatrix_4x4X4x4(double(*a)[4][2], double(*b)[4][2], double(*c)[4][2]){
	int i, j, k;
	double sum[2];

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
		sum[0] = 0.0;
		sum[1] = 0.0;

		for (k = 0; k < 4; k++) {
			sum[0] += (a[i][k][0] * b[k][j][0] - a[i][k][1] * b[k][j][1]);
			sum[1] += (a[i][k][0] * b[k][j][1] + a[i][k][1] * b[k][j][0]);
		}
		c[i][j][0] = sum[0];
		c[i][j][1] = sum[1];
		}
}


void MatrixMulVectorToVector_4x4(double(*a)[4][2], double(*v1)[2], double(*v2)[2]){
	int i, j;
	double sum[2];


	for (i = 0; i < 4; i++) {
		sum[0] = 0.0;
		sum[1] = 0.0;
		for (j = 0; j < 4; j++) {
			sum[0] += (a[i][j][0] * v1[j][0] - a[i][j][1] * v1[j][1]);
			sum[1] += (a[i][j][0] * v1[j][1] + a[i][j][1] * v1[j][0]);
		}
		v2[i][0] = sum[0];
		v2[i][1] = sum[1];
	}
}

void ScalarMulMatrixToMatrix_4x4(double *s, double(*a)[4][2], double(*c)[4][2]){
	int i, j;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
		c[i][j][0] = s[0] * a[i][j][0] - s[1] * a[i][j][1];
		c[i][j][1] = s[0] * a[i][j][1] + s[1] * a[i][j][0];
		}
}

void VectorMulVectorToScalar_4(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < 4; i++) {
		s[0] += (v1[i][0] * v2[i][0] - v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] + v1[i][1] * v2[i][0]);
	}
}

void VectorMulVectorToScalarH_4(double(*v1)[2], double(*v2)[2], double *s){
	int i;

	s[0] = 0.0;
	s[1] = 0.0;
	for (i = 0; i < 4; i++) {
		s[0] += (v1[i][0] * v2[i][0] + v1[i][1] * v2[i][1]);
		s[1] += (v1[i][0] * v2[i][1] - v1[i][1] * v2[i][0]);
	}
}


void VectorAddVectorToVector_4(double(*v1)[2], double(*v2)[2], double(*v3)[2]){
	int i;

	for (i = 0; i < 4; i++) {
		v3[i][0] = v1[i][0] + v2[i][0];
		v3[i][1] = v1[i][1] + v2[i][1];
	}
}

void VectorSubVectorToVector_4(double(*v1)[2], double(*v2)[2], double(*v3)[2]){
	int i;

	for (i = 0; i < 4; i++) {
		v3[i][0] = v1[i][0] - v2[i][0];
		v3[i][1] = v1[i][1] - v2[i][1];
	}
}

void VectorMulMatrixToVector_4x4(double(*v1)[2], double(*a)[4][2], double(*v2)[2]){
	int i, j;

	for (i = 0; i < 4; i++) {
		v2[i][0] = 0.0;
		v2[i][1] = 0.0;
		for (j = 0; j < 4; j++) {
			v2[i][0] += v1[j][0] * a[j][i][0] - v1[j][1] * a[j][i][1];
			v2[i][1] += v1[j][0] * a[j][i][1] + v1[j][1] * a[j][i][0];
		}
	}

}

void VectorMulVectorToMatrix_4x4(double(*v1)[2], double(*v2)[2], double(*a)[4][2]){
	int i, j;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
		a[i][j][0] = v1[i][0] * v2[j][0] - v1[i][1] * v2[j][1];
		a[i][j][1] = v1[i][0] * v2[j][1] + v1[i][1] * v2[j][0];
		}
}

void VectorConjugate_4(double(*v1)[2], double(*v2)[2]){
	int j;

	for (j = 0; j < 4; j++) {
		v2[j][0] = v1[j][0];
		v2[j][1] = -v1[j][1];
	}
}

void ScalarMulVectorToVector_4(double *s, double(*v1)[2], double(*v2)[2]){
	int i;

	for (i = 0; i < 4; i++) {
		v2[i][0] = s[0] * v1[i][0] - s[1] * v1[i][1];
		v2[i][1] = s[0] * v1[i][1] + s[1] * v1[i][0];
	}
}

void VectorDistance_4(double(*v1)[2], double(*v2)[2], double s){	/* Doesnt work */
	int i;
	double ith_distance[4];
	s = 0;

	for (i = 0; i < 4; i++) {
		ith_distance[i] = ScalarDistance(v1[i], v2[i]);
		s += ith_distance[i];
	}

}

double MagnitudeOfDifference_4(double(*v1)[2], double(*v2)[2]){
	double Difference[4][2];
	double res;
	int i;
	res = 0;
	for (i = 0; i < 4; i++)
	{
		Difference[i][0] = v1[i][0] - v2[i][0];
		Difference[i][1] = v1[i][1] - v2[i][1];
		res += sqrt(sqr(Difference[i][0]) + sqr(Difference[i][1]));
	}
	return res;
}

void VectorEqualVector_4(double(*v1)[2], double(*v2)[2]){
	int i;
	for (i = 0; i < 4; i++){
		v2[i][0] = v1[i][0];
		v2[i][1] = v1[i][1];
	}
}


