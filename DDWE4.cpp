#include "const-function-declerations.h"		 

extern double EbNo, CNR;
extern int loop;
extern double channel_gain[N][2], p_in[N];
extern errno_t err;	/* Error vector */
extern double P_for_DDWE[N][CL][CL][2], j[CL][2];
extern double OutPflag;
double OutPflag2 = 1;

const std::string currentDateTime() {
	// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}

void DDWE_Update(
	int f_index,
	int(*Decision),
	int R_index,
	double(*received_vector)[2],
	double(*reference_symbol),
	double(*reference_symbol_bs3),
	double(*old_weights)[2],
	double(*new_weights)[2],
	double(*H)[BS][2],
	double(*e)){

	/*if (f_index == 5){
		cout << "user 0 ref" << reference_symbol[0] << "+j" << reference_symbol[1] << endl;
		cout << "user 1 ref" << reference_symbol_bs3[0] << "+j" << reference_symbol_bs3[1] << endl;
		getchar();
	}	  */

#ifndef DDWE_VARIABLES
	Matrix3cf TempoH;
	Vector3f EigValue_Differences;
	ComplexEigenSolver<Matrix3cf> TempoHEig;   
	double
		bir[2] = { 1, 0 },
		sifir[CL][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 },
		//e[2],
		temp[2],
		temp_2[2],
		received_symbol[2],
		tempvec[CL][2],
		wh[CL][2],
		tempvec2[CL][2],
		oldw_store[CL][2],
		x[CL][2],
		xh[CL][2],
		lpx[CL][2],
		lP_for_DDWE[N][CL][CL][2],
		tempmat[CL][CL][2],
		s_qpsk[4][2] = { { OneBySqrt2, OneBySqrt2 }, { OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, OneBySqrt2 } };
	int s3_index, s1_index, i;
#endif // !DDWE_VARIABLES						   

#ifndef UPDATE_WEIGHTS
	for (i = 0; i < TS; i++)
	{
		x[i][0] = received_vector[i][0];
		x[i][1] = received_vector[i][1];
	}
	x[CL - 1][0] = reference_symbol_bs3[0];
	x[CL - 1][1] = reference_symbol_bs3[1];


	VectorConjugate_CL(x, xh);
	ScalarMulMatrixToMatrix_CLxCL(ff_inverse, P_for_DDWE[f_index], lP_for_DDWE[f_index]);
	MatrixMulVectorToVector_CLxCL(lP_for_DDWE[f_index], x, lpx);	// nom
	VectorMulVectorToScalar_CL(xh, lpx, temp);
	temp[0] += bir[0];	// denom
	ScalarDivScalarToScalar(bir, temp, temp_2);
	ScalarMulVectorToVector_CL(temp_2, lpx, j);	// eq. 1

	VectorEqualVector_CL(old_weights, oldw_store);

	VectorConjugate_CL(old_weights, wh);
	VectorMulVectorToScalar_CL(x, wh, received_symbol);

	e[0] = reference_symbol[0] - received_symbol[0];
	e[1] = -(reference_symbol[1] - received_symbol[1]);

	//if (sqr(e[f_index][index][0]) + sqr(e[f_index][index][1]) < 1)
	//{
	ScalarMulVectorToVector_CL(e, j, tempvec);
	VectorAddVectorToVector_CL(old_weights, tempvec, tempvec2);
	VectorAddVectorToVector_CL(tempvec2, sifir, new_weights);

	VectorMulMatrixToVector_CLxCL(xh, lP_for_DDWE[f_index], tempvec);
	VectorMulVectorToMatrix_CLxCL(j, tempvec, tempmat);
	MatrixSubMatrixToMatrix_CL(lP_for_DDWE[f_index], tempmat, P_for_DDWE[f_index]);
	//}

#endif // !UPDATE_WEIGHTS

#ifdef SD							
	double
		power[2],
		gH[CL][2],
		g_i[BS][2],
		power_i[2],
		g_iH[BS][2],
		R[BS][2],
		no = pow(10.0, -CNR / 10.0) / CODING_RATE;
	VectorMulMatrixToVector_TSxTS(new_weights, H, R);
	g_i[0][0] = 0.0; g_i[0][1] = 0.0;
	for (size_t i = 1; i < M; i++){g_i[i][0] = R[i][0]; g_i[i][1] = R[i][1];}
	VectorConjugate_M(g_i, g_iH);
	VectorMulVectorToScalar_M(g_iH, g_i, power_i);	/* power_i = g_i^H x g_i */
	VectorConjugate_TS(new_weights, gH);
	VectorMulVectorToScalar_TS(gH, new_weights, power);
	channel_gain[f_index][0] = R[0][0];
	channel_gain[f_index][1] = -R[0][1];
	p_in[f_index] = no * power[0] + power_i[0];	/* Incoming signal power. */
#endif // SD			
}	/* DDWE_Update :: end */