#include "const.h"

extern double CNR;
extern int loop;
extern double channel_gain[N][2], p_in[N];
extern errno_t err;	/* Error vector */
extern double P_for_DDWE[N][CL][CL][2], j[CL][2];
extern double OutPflag;

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
	double(*e)[NUMBER_OF_STEPS][2],
	double(*reference_symbol_bs3),
	int f_index,
	int index){

#ifndef DDWE_VARIABLES
	static double
		pf_counter;
	double
		power[2],
		gH[CL][2],
		g_i[BS][2],
		power_i[2],
		g_iH[BS][2],
		no = pow(10.0, -CNR / 10.0) / CODING_RATE,
		bir[2] = { 1, 0 },
		sifir[CL][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 },
		temp[2],
		temp2[2],
		temp3[2],
		received_symbol[2],
		R[BS][2],
		tempvec[CL][2],
		tempvech[CL][2],
		tempvec2[CL][2],
		tempvec3[CL][2],
		x[CL][2],
		xh[CL][2],
		lpx[CL][2],
		lP_for_DDWE[N][CL][CL][2],
		tempmat[CL][CL][2],
		tempmat2[CL][CL][2],
		s_qpsk[4][2] = { { OneBySqrt2, OneBySqrt2 }, { OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, OneBySqrt2 } };
	int s3_index, s1_index;
#endif // !DDWE_VARIABLES

#ifndef UPDATE_WEIGHTS
	for (size_t i = 0; i < TS; i++)
	{
		x[i][0] = received_vector[i][0];
		x[i][1] = received_vector[i][1];
	}
	x[CL - 1][0] = reference_symbol_bs3[0];
	x[CL - 1][1] = reference_symbol_bs3[1];


	VectorConjugate_CL(x, xh);
	ScalarMulMatrixToMatrix_CLxCL(ff_inverse, P_for_DDWE[f_index], lP_for_DDWE[f_index]);
	MatrixMulVectorToVector_CLxCL(lP_for_DDWE[f_index], x, lpx);	// nom
	VectorMulVectorToScalar_CL(xh, lpx, temp2);
	temp2[0] += bir[0];	// denom
	ScalarDivScalarToScalar(bir, temp2, temp3);
	ScalarMulVectorToVector_CL(temp3, lpx, j);	// eq. 1

	VectorConjugate_CL(old_weights, tempvech);
	VectorMulVectorToScalar_CL(x, tempvech, received_symbol);
	//VectorMulVectorToScalar_CL(xh, old_weights, received_symbol);
	/*temp[0] = 0.0;
	temp[1] = 0.0;
	for (ts = 0; ts < CL; ts++)
	{
	temp[0] += W[f_index][0][ts][0] * tempvec[ts][0] - W[f_index][0][ts][1] * tempvec[ts][1];
	temp[1] += W[f_index][0][ts][0] * tempvec[ts][1] + W[f_index][0][ts][1] * tempvec[ts][0];
	} */

	e[f_index][index][0] = reference_symbol[0] - received_symbol[0];
	e[f_index][index][1] = -(reference_symbol[1] - received_symbol[1]);

	if ((sqr(e[f_index][index][0]) + sqr(e[f_index][index][1]) > 50) || (pf_counter > 0)){
		if (pf_counter == 0) { pf_counter = 5; cout << "if1:" << pf_counter; }
		if (pf_counter > 0) { pf_counter--; cout << "if2" << pf_counter; }
		cout << endl;
		cout << "Counter: " << pf_counter << endl;
		cout << "f_index: " << f_index << "\tindex: " << index << endl;
		cout << received_symbol[0] << "\t" << received_symbol[1] << endl;
		cout << reference_symbol[0] << "\t" << reference_symbol[1] << endl;
		cout << reference_symbol_bs3[0] << "\t" << reference_symbol_bs3[1] << endl;
		cout << "Received vector:" << endl;
		cout << received_vector[0][0] << "\t" << received_vector[0][1] << "\t" << sqr(received_vector[0][1]) + sqr(received_vector[0][0]) << endl;
		cout << received_vector[1][0] << "\t" << received_vector[1][1] << "\t" << sqr(received_vector[1][1]) + sqr(received_vector[1][0]) << endl;
		cout << "Old weights:" << endl;
		cout << old_weights[0][0] << "\t" << old_weights[0][1] << endl;
		cout << old_weights[1][0] << "\t" << old_weights[1][1] << endl;
		cout << old_weights[2][0] << "\t" << old_weights[2][1] << endl;

		}		
	//if (sqr(e[f_index][index][0]) + sqr(e[f_index][index][1]) < 1)
	//{
		ScalarMulVectorToVector_CL(e[f_index][index], j, tempvec);
		VectorAddVectorToVector_CL(old_weights, tempvec, tempvec2);
		VectorAddVectorToVector_CL(tempvec2, sifir, new_weights);

		VectorMulMatrixToVector_CLxCL(xh, lP_for_DDWE[f_index], tempvec);
		VectorMulVectorToMatrix_CLxCL(j, tempvec, tempmat);
		MatrixSubMatrixToMatrix_CL(lP_for_DDWE[f_index], tempmat, P_for_DDWE[f_index]);	
	//}


	if ((sqr(e[f_index][index][0]) + sqr(e[f_index][index][1]) > 50) || (pf_counter > 0)){
			cout << "New weights:" << endl;
			cout << new_weights[0][0] << "\t" << new_weights[0][1] << endl;
			cout << new_weights[1][0] << "\t" << new_weights[1][1] << endl;
			cout << new_weights[2][0] << "\t" << new_weights[2][1] << endl;
			OutPflag = loop;
			getchar();
			}	
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