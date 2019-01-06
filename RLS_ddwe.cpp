#include "const.h"

extern double CNR, EbNo;
extern int loop;
extern double channel_gain[N][2], p_in[N], channel_power[N][2];
extern errno_t err;	/* Error vector */
extern double P_for_DDWE2[N][FILTER_ORDER + 1][FILTER_ORDER + 1][2], j[FILTER_ORDER + 1][2];

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
		g[TS][2],
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
		received_symbol[2],
		//j[FILTER_ORDER + 1][2],
		tempvec[FILTER_ORDER + 1][2],
		tempvec2[FILTER_ORDER + 1][2],
		tempvec3[FILTER_ORDER + 1][2],
		tempmat[FILTER_ORDER + 1][FILTER_ORDER + 1][2],
		tempmat2[FILTER_ORDER + 1][FILTER_ORDER + 1][2];
#endif // !DDWE_VARIABLES


#ifndef UPDATE_WEIGHTS

	VectorMulVectorToScalar(received_vector, old_weights, received_symbol);
	e[0] = reference_symbol[0] - received_symbol[0];
	e[1] = reference_symbol[1] - received_symbol[1];

	VectorConjugate(received_vector, tempvec2);	/* y^C_m[n] */
	MatrixMulVectorToVector(P_for_DDWE2[f_index], tempvec2, tempvec3);	/* P[n-1] x y^C_m[n] */
	VectorMulVectorToScalar(received_vector, tempvec3, temp);	/* y^T_m[n] x P[n-1] x y^C_m[n] */
	temp[0] += FORGETTING_FACTOR;	/* lamda + y^T_m[n] x P[n-1] x y^C_m[n] */
	temp2[0] = 1.0; temp2[1] = 0.0;	/* {}^-1 */
	ScalarDivScalarToScalar(temp2, temp, temp3); /* {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */
	ScalarMulVectorToVector(temp3, tempvec3, j);	/*  y^H_m[n] x {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */

	temp2[0] = FF_INVERSE; temp2[1] = 0.0;	/* lamda^-1 */
	ScalarMulMatrixToMatrix(temp2, P_for_DDWE2[f_index], tempmat);	/* lamda^-1 x P[n-1] */
	VectorMulMatrixToVector(received_vector, tempmat, tempvec2);	/* y^T_m[n] x lamda^-1 x P[n-1] */
	VectorMulVectorToMatrix(j, tempvec2, tempmat2);	/* g[n] x y^T_m[n] x lamda^-1 x P[n-1] */
	MatrixSubMatrixToMatrix(tempmat, tempmat2, P_for_DDWE2[f_index]);	/* lamda^-1 x P[n-1] - g[n] x y^T_m[n] x lamda^-1 x P[n-1] */

	ScalarMulVectorToVector(e, j, received_vector);
	VectorAddVectorToVector(old_weights, received_vector, new_weights);
#endif // !UPDATE_WEIGHTS

#ifndef POWER_UPDATE
#ifdef SD
	VectorMulMatrixToVector(new_weights, H, R);
	g_i[0][0] = 0.0; g_i[0][1] = 0.0; g_i[1][0] = R[1][0]; g_i[1][1] = R[1][1];
	VectorConjugate2(g_i, g_iH);
	VectorMulVectorToScalar2(g_iH, g_i, power_i);	/* power_i = g_i^H x g_i */
	VectorConjugate(new_weights, gH);
	VectorMulVectorToScalar(gH, new_weights, power);
	channel_gain[f_index][0] = R[0][0];
	channel_gain[f_index][1] = -R[0][1];
	p_in[f_index] = no * power[0] + power_i[0];	/* Incoming signal power. */
#endif // SD						  
#endif // !POWER_UPDATE
}	/* DDWE_Update :: end */