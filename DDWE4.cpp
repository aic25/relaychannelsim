#include "const-function-declerations.h"		 

extern double EbNo, CNR;
extern int loop;
extern double channel_gain[N][2], p_in[N];
extern errno_t err;	/* Error vector */
extern double sigma_old[N][2];
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

void DF_CombCoeff(
	int symbol_index,
	int f_index,
	int(*Decision),
	int R_index,
	int S_index,
	double(*received_vector)[2],
	double(*reference_symbol),
	double(*reference_symbol_bs3),
	double(*old_weights)[2],
	double(*new_weights)[2],
	double(*H)[BS][2],
	double(*e)){

#ifndef DDWE_VARIABLES
	Matrix3cf TempoH;
	Vector3f EigValue_Differences;
	ComplexEigenSolver<Matrix3cf> TempoHEig;
	double
		bir[2] = { 1, 0 },
		sifir[CL][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 },
		u[2],
		d[2],
		k[2],
		wu[2],
		ke[2],
		p_denom[2],
		p_nom[2],
		ffxs_old[2],
		sigma_new[2],
		/*temp[2],
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
		tempmat[CL][CL][2],	 */
		s_qpsk[4][2] = { { OneBySqrt2, OneBySqrt2 }, { OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, OneBySqrt2 } };
	int s3_index, s1_index, i;
#endif // !DDWE_VARIABLES						   

#ifndef UPDATE_WEIGHTS
	if ((f_index == 5) && (DBG == ON))cout << "step: " << symbol_index << endl;
	if ((f_index == 5) && (DBG == ON))cout << "y[0]: " << received_vector[0][0] << " " << received_vector[0][1] << endl;
	if ((f_index == 5) && (DBG == ON))cout << "y[1]: " << received_vector[1][0] << " " << received_vector[1][1] << endl;
	if (S_index==0){
		u[0] = received_vector[1][0] -
			(
			H[1][R_index][0] * reference_symbol[0] -
			H[1][R_index][1] * reference_symbol[1]
			) -
			(
			H[1][(R_index + 1) % 2][0] * reference_symbol_bs3[0] -
			H[1][(R_index + 1) % 2][1] * reference_symbol_bs3[1]
			);
		u[1] = received_vector[1][1] -
			(
			H[1][R_index][0] * reference_symbol[1] +
			H[1][R_index][1] * reference_symbol[0]
			) -
			(
			H[1][(R_index + 1) % 2][0] * reference_symbol_bs3[1] +
			H[1][(R_index + 1) % 2][1] * reference_symbol_bs3[0]
			);
		if ((f_index == 5) && (DBG == ON))cout << "u: " << u[0] << " " << u[1] << endl;
		d[0] = received_vector[0][0] -
			(
			H[0][R_index][0] * reference_symbol[0] -
			H[0][R_index][1] * reference_symbol[1]
			) -
			(
			H[0][(R_index + 1) % 2][0] * reference_symbol_bs3[0] -
			H[0][(R_index + 1) % 2][1] * reference_symbol_bs3[1]
			);
		d[1] = received_vector[0][1] -
			(
			H[0][R_index][0] * reference_symbol[1] +
			H[0][R_index][1] * reference_symbol[0]
			) -
			(
			H[0][(R_index + 1) % 2][0] * reference_symbol_bs3[1] +
			H[0][(R_index + 1) % 2][1] * reference_symbol_bs3[0]
			);
	}
	else{ 
		d[0] = received_vector[1][0] -
			(
			H[1][R_index][0] * reference_symbol[0] -
			H[1][R_index][1] * reference_symbol[1]
			) -
			(
			H[1][(R_index + 1) % 2][0] * reference_symbol_bs3[0] -
			H[1][(R_index + 1) % 2][1] * reference_symbol_bs3[1]
			);
		d[1] = received_vector[1][1] -
			(
			H[1][R_index][0] * reference_symbol[1] +
			H[1][R_index][1] * reference_symbol[0]
			) -
			(
			H[1][(R_index + 1) % 2][0] * reference_symbol_bs3[1] +
			H[1][(R_index + 1) % 2][1] * reference_symbol_bs3[0]
			);
		if ((f_index == 5) && (DBG == ON))cout << "u: " << u[0] << " " << u[1] << endl;
		u[0] = received_vector[0][0] -
			(
			H[0][R_index][0] * reference_symbol[0] -
			H[0][R_index][1] * reference_symbol[1]
			) -
			(
			H[0][(R_index + 1) % 2][0] * reference_symbol_bs3[0] -
			H[0][(R_index + 1) % 2][1] * reference_symbol_bs3[1]
			);
		u[1] = received_vector[0][1] -
			(
			H[0][R_index][0] * reference_symbol[1] +
			H[0][R_index][1] * reference_symbol[0]
			) -
			(
			H[0][(R_index + 1) % 2][0] * reference_symbol_bs3[1] +
			H[0][(R_index + 1) % 2][1] * reference_symbol_bs3[0]
			);
	}
	if ((f_index == 5) && (DBG == ON))cout << "d: " << d[0] << " " << d[1] << endl;
	if ((f_index == 5) && (DBG == ON))cout << "sigma_o: " << sigma_old[f_index][0] << " " << sigma_old[f_index][1] << endl;
	ScalarMulScalar(ff, sigma_old[f_index], ffxs_old);
	if ((f_index == 5) && (DBG == ON))cout << "ffxs_old: " << ffxs_old[0] << " " << ffxs_old[1] << endl;
	sigma_new[0] = ffxs_old[0] + (sqr(u[0]) + sqr(u[1]));
	sigma_new[1] = ffxs_old[1];

	if ((f_index == 5) && (DBG == ON))cout << "sigma_n: " << sigma_new[0] << " " << sigma_new[1] << endl;
	//if ((f_index == 5)&&(DBG==ON))cout << "div ddwe" << u[0] << "+j" << u[1] << " / " << sigma_new[0] << "+j" << sigma_new[1] << " = " << k[0] << "+j" << k[1] << endl; 
	ScalarDivScalarToScalar(u, sigma_new, k);
	if ((f_index == 5) && (DBG == ON))cout << "div ddwe" << u[0] << "+j" << u[1] << " / " << sigma_new[0] << "+j" << sigma_new[1] << " = " << k[0] << "+j" << k[1] << endl;
	if ((f_index == 5) && (DBG == ON))cout << "k: " << k[0] << " " << k[1] << endl;

	ScalarMulScalar_H(u, old_weights[0], wu);
	if ((f_index == 5) && (DBG == ON))cout << "wu: " << wu[0] << " " << wu[1] << endl;
	e[0] = d[0] - wu[0];
	e[1] = d[1] - wu[1];
	if ((f_index == 5) && (DBG == ON))cout << "e: " << e[0] << " " << e[1] << endl;
	if ((f_index == 5) && (DBG == ON))cout << "|e|^2: " << sqr(e[0]) + sqr(e[1]) << endl;
	if ((f_index == 5) && (DBG == ON))getchar();
	ScalarMulScalar_H(k, e, ke);
	new_weights[0][0] = old_weights[0][0] + ke[0];
	new_weights[0][1] = old_weights[0][1] + ke[1];

	sigma_old[f_index][0] = sigma_new[0];
	sigma_old[f_index][1] = sigma_new[1];

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
}	/* DF_CombCoeff :: end */