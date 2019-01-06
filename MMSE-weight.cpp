#include "const-function-declerations.h"		
extern double CNR, channel_power[N][2];
double channel_gain[N][2], p_in[N];

void MMSE_Weight_Matrix(double(*H)[TS][BS][2], double(*W)[BS][CL][2]) {
	double no, hbhrH3[2], hbhbH3[2];
	int f_index;
	std::complex<double> h_m1r3, h_m1b3, h_m1r3H, h_m1b3H, w;
	/* This one looks like numerical CNR calculation, while taking into account the coding rate */
	//no = pow(10.0, -CNR / 10.0) / CODING_RATE;	/* Noise power. if you assume one for signal power... */
	for (f_index = 0; f_index < N; f_index++) {
		/*ScalarMulScalar_H(H[f_index][0][2], H[f_index][1][2], hbhrH3);
		ScalarMulScalar_H(H[f_index][0][2], H[f_index][0][2], hbhbH3);
		hbhbH3[0] += no;
		ScalarDivScalarToScalar(hbhrH3, hbhbH3, W[f_index][0][0]);
		W[f_index][0][0][0] = -W[f_index][0][0][0];
		W[f_index][0][0][1] = -W[f_index][0][0][1];
		REMOVE COMMENT WHEN NOT THEORETICAL			
		ScalarMulScalar_H(H[f_index][1][2], H[f_index][0][2], hbhrH3);
		ScalarMulScalar_H(H[f_index][1][2], H[f_index][1][2], hbhbH3);
		hbhbH3[0] = hbhbH3[0] * (1 + pow(10.0, -5));
		hbhbH3[1] = hbhbH3[1] * (1 + pow(10.0, -5));
		hbhbH3[0] += no;
		ScalarDivScalarToScalar(hbhrH3, hbhbH3, W[f_index][0][0]);
		W[f_index][0][0][0] = -W[f_index][0][0][0];
		W[f_index][0][0][1] = -W[f_index][0][0][1];	 */ 		


		h_m1b3 = complex<double>(H[f_index][0][2][0], H[f_index][0][2][1]);
		h_m1b3H = complex<double>(H[f_index][0][2][0], -H[f_index][0][2][1]);
		h_m1r3 = complex<double>(H[f_index][1][2][0], H[f_index][1][2][1]);
		h_m1r3H = complex<double>(H[f_index][1][2][0], -H[f_index][1][2][1]);
		w = -(h_m1r3*h_m1b3H) / ((h_m1r3H*h_m1r3).real()*(1 + pow(10.0, -5)) + pow(10.0, -CNR / 10.0));
		W[f_index][0][0][0] = w.real();
		W[f_index][0][0][1] = w.imag();  
		/*h_m1b3 = complex<double>(H[f_index][0][2][0], H[f_index][0][2][1]);
		h_m1b3H = complex<double>(H[f_index][0][2][0], -H[f_index][0][2][1]);
		h_m1r3 = complex<double>(H[f_index][1][2][0], H[f_index][1][2][1]);
		h_m1r3H = complex<double>(H[f_index][1][2][0], -H[f_index][1][2][1]);
		w = -(h_m1b3*h_m1r3H) / ((h_m1b3*h_m1b3H) + pow(10.0, -CNR / 10.0));
		W[f_index][0][0][0] = w.real();
		W[f_index][0][0][1] = w.imag();	*/  

		/*cout << "hr - w*hb:"
		<< pow(abs(complex<double>(H[f_index][1][2][0], H[f_index][1][2][1])
		- complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1])
		*complex<double>(H[f_index][0][2][0], H[f_index][0][2][1])), 2) << endl;
		cout << "hb - w*hr:"
		<< pow(abs(complex<double>(H[f_index][0][2][0], H[f_index][0][2][1])
		- complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1])
		*complex<double>(H[f_index][1][2][0], H[f_index][1][2][1])), 2) << endl;
		cout << "hr + w*hb:"
		<< pow(abs(complex<double>(H[f_index][1][2][0], H[f_index][1][2][1])
		+ complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1])
		*complex<double>(H[f_index][0][2][0], H[f_index][0][2][1])), 2) << endl;
		cout << "hb + w*hr:"
		<< pow(abs(complex<double>(H[f_index][0][2][0], H[f_index][0][2][1])
		+ complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1])
		*complex<double>(H[f_index][1][2][0], H[f_index][1][2][1])), 2) << endl;
		getchar(); 	*/
	}
}