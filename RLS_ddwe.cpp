#include "const-function-declerations.h"		 

extern double CNR;
extern int loop;
extern double channel_gain[N][2], p_in[N];
extern errno_t err;	/* Error vector */
extern double P_for_DDWE[N][TS][TS][2], j[TS][2];

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
	int i;
	double
		power[2],
		gH[TS][2],
		g_i[BS][2],
		power_i[2],
		g_iH[BS][2];
	double
		no = pow(10.0, -CNR / 10.0) / CODING_RATE,
		no_loss = pow(10.0, -(CNR - LOSS) / 10.0) / CODING_RATE,
		bir[2] = { 1, 0 },
		sifir[TS][2] = { { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 },
		temp[2],
		temp2[2],
		temp3[2],
		R[BS][2],
		x[TS][2],
		xh[TS][2],
		lpx[TS][2],
		lP_for_DDWE[N][TS][TS][2],
		tempvec[TS][2],
		tempvech[TS][2],
		tempvec2[TS][2],
		received_symbol[2],
		tempmat[TS][TS][2],
		tempmat2[TS][TS][2];
#endif // !DDWE_VARIABLES


#ifndef UPDATE_WEIGHTS
	for (i = 0; i < TS; i++)
	{
		x[i][0] = received_vector[i][0];
		x[i][1] = received_vector[i][1];
	}

	VectorConjugate_TS(x, xh);
	ScalarMulMatrixToMatrix_TSxTS(ff_inverse, P_for_DDWE[f_index], lP_for_DDWE[f_index]);
	MatrixMulVectorToVector_TSxTS(lP_for_DDWE[f_index], x, lpx);	// nom
	VectorMulVectorToScalar_TS(xh, lpx, temp2);
	temp2[0] += bir[0];	// denom
	ScalarDivScalarToScalar(bir, temp2, temp3);
	ScalarMulVectorToVector_TS(temp3, lpx, j);	// eq. 1

	VectorConjugate_TS(old_weights, tempvech);
	VectorMulVectorToScalar_TS(x, tempvech, received_symbol);

	e[0] = reference_symbol[0] - received_symbol[0];
	e[1] = -(reference_symbol[1] - received_symbol[1]);

	ScalarMulVectorToVector_TS(e, j, tempvec);
	VectorAddVectorToVector_TS(old_weights, tempvec, tempvec2);
	VectorAddVectorToVector_TS(tempvec2, sifir, new_weights);

	VectorMulMatrixToVector_TSxTS(xh, lP_for_DDWE[f_index], tempvec);
	VectorMulVectorToMatrix_TSxTS(j, tempvec, tempmat);
	MatrixSubMatrixToMatrix_TS(lP_for_DDWE[f_index], tempmat, P_for_DDWE[f_index]);
#endif // !UPDATE_WEIGHTS

#ifndef UPDATE_POWER
#ifdef SD
	VectorMulMatrixToVector_TSxTS(new_weights, H, R);
	g_i[0][0] = 0.0; g_i[0][1] = 0.0;
	for (size_t i = 1; i < M; i++){ g_i[i][0] = R[i][0]; g_i[i][1] = R[i][1]; }
	VectorConjugate_M(g_i, g_iH);
	VectorMulVectorToScalar_M(g_iH, g_i, power_i);	/* power_i = g_i^H x g_i */

	/*-----skip*/
	VectorConjugate_TS(new_weights, gH);
	VectorMulVectorToScalar_TS(gH, new_weights, power);
	/*------------skip*/
	power[0] = new_weights[0][0] * new_weights[0][0] + new_weights[0][1] * new_weights[0][1];
	power[1] = new_weights[1][0] * new_weights[1][0] + new_weights[1][1] * new_weights[1][1];

	channel_gain[f_index][0] = R[0][0];
	channel_gain[f_index][1] = -R[0][1];
	p_in[f_index] = (no_loss)*power[0] + (no)*power[1] + power_i[0];	/* Incoming signal power. */
#endif // SD						  
#endif // !POWER_UPDATE
}	/* DDWE_Update :: end */