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
	int f_index){

#ifndef DDWE_VARIABLES
	double
		power[2],
		gH[TS][2],
		g_i[BS][2],
		power_i[2],
		g_iH[BS][2];
	double
		no = pow(10.0, -CNR / 10.0) / CODING_RATE,
		temp[2],
		temp2[2],
		temp3[2],
		R[BS][2],
		tempvec[TS][2],
		tempvec2[TS][2],
		received_symbol[2],
		tempmat[TS][TS][2],
		tempmat2[TS][TS][2];
#endif // !DDWE_VARIABLES


#ifndef UPDATE_WEIGHTS
	VectorMulVectorToScalar_TS(received_vector, old_weights, received_symbol);
	e[0] = reference_symbol[0] - received_symbol[0];
	e[1] = reference_symbol[1] - received_symbol[1];

	VectorConjugate_TS(received_vector, tempvec);	/* y^C_m[n] */
	MatrixMulVectorToVector_TSxTS(P_for_DDWE[f_index], tempvec, tempvec2);	/* P[n-1] x y^C_m[n] */
	VectorMulVectorToScalar_TS(received_vector, tempvec2, temp);	/* y^T_m[n] x P[n-1] x y^C_m[n] */
	temp[0] += FORGETTING_FACTOR;	/* lamda + y^T_m[n] x P[n-1] x y^C_m[n] */
	temp2[0] = 1.0; temp2[1] = 0.0;	/* {}^-1 */
	ScalarDivScalarToScalar(temp2, temp, temp3); /* {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */
	ScalarMulVectorToVector_TS(temp3, tempvec2, j);	/*  y^H_m[n] x {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */

	temp2[0] = FF_INVERSE; temp2[1] = 0.0;	/* lamda^-1 */
	ScalarMulMatrixToMatrix_TSxTS(temp2, P_for_DDWE[f_index], tempmat);	/* lamda^-1 x P[n-1] */
	VectorMulMatrixToVector_TSxTS(received_vector, tempmat, tempvec);	/* y^T_m[n] x lamda^-1 x P[n-1] */
	VectorMulVectorToMatrix_TSxTS(j, tempvec, tempmat2);	/* g[n] x y^T_m[n] x lamda^-1 x P[n-1] */
	MatrixSubMatrixToMatrix_TS(tempmat, tempmat2, P_for_DDWE[f_index]);	/* lamda^-1 x P[n-1] - g[n] x y^T_m[n] x lamda^-1 x P[n-1] */

	ScalarMulVectorToVector_TS(e, j, received_vector);
	VectorAddVectorToVector_TS(old_weights, received_vector, new_weights);
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