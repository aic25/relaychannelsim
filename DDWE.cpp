#include "const.h"

extern double CNR;
extern int loop;
extern double channel_gain[N][2], p_in[N];
extern errno_t err;	/* Error vector */
extern double P_for_DDWE[N][CL][CL][2], j[CL][2];

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
	double(*received_vector)[2],
	double(*reference_symbol),
	double(*old_weights)[2],
	double(*new_weights)[2],
	double(*H)[BS][2],
	double(*e),
	double(*reference_symbol_bs3),
	int f_index){

#ifndef DDWE_VARIABLES
	double
		power[2],
		gH[CL][2],
		g_i[BS][2],
		power_i[2],
		g_iH[BS][2];
	double
		no = pow(10.0, -CNR / 10.0) / CODING_RATE,
		temp[2],
		temp2[2],
		temp3[2],
		R[BS][2],
		tempvec[CL][2],
		tempvec2[CL][2],
		tempvec3[CL][2],
		tempvec4[CL][2],
		received_symbol[2],
		tempmat[CL][CL][2],
		tempmat2[CL][CL][2],
		s_qpsk[4][2] = { { OneBySqrt2, OneBySqrt2 }, { OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, OneBySqrt2 } };
	int s3_index, s1_index;
#endif // !DDWE_VARIABLES


#ifndef UPDATE_WEIGHTS
	for (size_t i = 0; i < TS; i++)
	{		 
		tempvec3[i][0] = received_vector[i][0];
		tempvec3[i][1] = received_vector[i][1];
	}
	tempvec3[CL - 1][0] = reference_symbol_bs3[0];
	tempvec3[CL - 1][1] = reference_symbol_bs3[1];

	VectorMulVectorToScalar_CL(tempvec3, old_weights, received_symbol);
	e[0] = reference_symbol[0] - received_symbol[0];
	e[1] = reference_symbol[1] - received_symbol[1];

	//e[0] = s_qpsk[Decision[0]][0] - received_symbol[0];
	//e[1] = s_qpsk[Decision[0]][1] - received_symbol[1];

	VectorConjugate_CL(tempvec3, tempvec);	/* y^C_m[n] */
	MatrixMulVectorToVector_CLxCL(P_for_DDWE[f_index], tempvec, tempvec2);	/* P[n-1] x y^C_m[n] */
	VectorMulVectorToScalar_CL(tempvec3, tempvec2, temp);	/* y^T_m[n] x P[n-1] x y^C_m[n] */
	temp[0] += FORGETTING_FACTOR;	/* lamda + y^T_m[n] x P[n-1] x y^C_m[n] */
	temp2[0] = 1.0; temp2[1] = 0.0;	/* {}^-1 */
	ScalarDivScalarToScalar(temp2, temp, temp3); /* {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */
	ScalarMulVectorToVector_CL(temp3, tempvec2, j);	/*  y^H_m[n] x {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */

	temp2[0] = FF_INVERSE; temp2[1] = 0.0;	/* lamda^-1 */
	ScalarMulMatrixToMatrix_CLxCL(temp2, P_for_DDWE[f_index], tempmat);	/* lamda^-1 x P[n-1] */
	VectorMulMatrixToVector_CLxCL(tempvec3, tempmat, tempvec);	/* y^T_m[n] x lamda^-1 x P[n-1] */
	VectorMulVectorToMatrix_CLxCL(j, tempvec, tempmat2);	/* g[n] x y^T_m[n] x lamda^-1 x P[n-1] */
	MatrixSubMatrixToMatrix_CL(tempmat, tempmat2, P_for_DDWE[f_index]);	/* lamda^-1 x P[n-1] - g[n] x y^T_m[n] x lamda^-1 x P[n-1] */

	ScalarMulVectorToVector_CL(e, j, tempvec4);
	VectorAddVectorToVector_CL(old_weights, tempvec4, new_weights);
#endif // !UPDATE_WEIGHTS

#ifndef UPDATE_POWER
#ifdef SD
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
#endif // !POWER_UPDATE
}	/* DDWE_Update :: end */