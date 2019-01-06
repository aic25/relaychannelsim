#include "const-function-declerations.h"		
extern int R_index[N];
extern double CNR, channel_power[N][2];
double channel_gain[N][2], p_in[N];

void MMSE_Weight_Matrix(double(*H)[TS][BS][2], double(*W)[BS][CL][2]) {
	double no, absw3_I[2], absw3_II[2];// , Ruu[3][3][2], Ruu_inv[3][3][2], Rud[3][2], Ryy[2][2][2], HH[3][2][2];
	int f_index, i, j;
	MatrixXcf h(2, 3), hc(2, 3);
	Matrix2cf Ryy;
	Matrix3cf Ruu;
	Vector3cf Rud, w;

	/* This one looks like numerical CNR calculation, while taking into account the coding rate */
	no = pow(10.0, -CNR / 10.0) / CODING_RATE;	/* Noise power. if you assume one for signal power... */

	for (f_index = 0; f_index < N; f_index++) {	 
#if THEORETICAL==OFF
		/*normal*/

		for (i = 0; i < 2; i++)
			for (j = 0; j < 3; j++) {
				h(i, j) = complex<double>(H[f_index][i][j][0], H[f_index][i][j][1]);
			}
		hc = h.conjugate();
		Ryy = h*h.adjoint();
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++) {
				if (i == j) Ruu(i, j) = Ryy(i, j) + (float)no;
				else Ruu(i, j) = Ryy(i, j);
			}

		/* R_index = 0 */
		Ruu(0, 2) = h(0, 1); Ruu(1, 2) = h(1, 1);
		Ruu(2, 0) = hc(0, 1); Ruu(2, 1) = hc(1, 1);
		Ruu(2, 2) = complex<double>(1, 0);
		Rud(0) = h(0, 0);
		Rud(1) = h(1, 0);
		Rud(2) = 0;
		w = Ruu.inverse()*Rud;
		for (i = 0; i < 3; i++) {
			W[f_index][0][i][0] = w(i).real();
			W[f_index][0][i][1] = w(i).imag();
			//cout << "w_" << i << ": " << w(i) << endl;
		}

#else
		/* R_index = 0 */
		/*trial	*/
		for (i = 0; i < 2; i++)
			for (j = 0; j < 3; j++) {
				h(i, j) = complex<double>(H[f_index][i][j][0], H[f_index][i][j][1]);
			}
		hc = h.conjugate();
		Ryy = h*h.adjoint();		
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++) {
				if (i == j) Ruu(i+1, j+1) = Ryy(i, j) + (float)no;
				else Ruu(i+1, j+1) = Ryy(i, j);
			}

		
		Ruu(2, 0) = -h(1, 1); Ruu(1, 0) = -h(0, 1);
		Ruu(0, 2) = -hc(1, 1); Ruu(0, 1) = -hc(0, 1);
		Ruu(0, 0) = complex<double>(1, 0);	 
		Rud(0) = 0;
		Rud(1) = h(0, 0);
		Rud(2) = h(1, 0);
		w = Ruu.inverse()*Rud;
		for (i = 0; i < 3; i++) {
			W[f_index][0][i][0] = w(i).real();
			W[f_index][0][i][1] = w(i).imag();
			//cout << "w_" << i << ": " << w(i) << endl;
		}								  
#endif

#ifdef MMSE_WSEL
		/* R_index = 1 */
		Ruu(0, 2) = h(0, 0); Ruu(1, 2) = h(1, 0);
		Ruu(2, 0) = hc(0, 0); Ruu(2, 1) = hc(1, 0);
		Ruu(2, 2) = complex<double>(1, 0);
		Rud(0) = h(0, 1);
		Rud(1) = h(1, 1);
		Rud(2) = 0;
		w = Ruu.inverse()*Rud;
		for (i = 0; i < 3; i++) {
			W[f_index][1][i][0] = w(i).real();
			W[f_index][1][i][1] = w(i).imag();
		}

		VectorMulVectorToScalarH_CL(W[f_index][0], W[f_index][0], absw3_I);
		VectorMulVectorToScalarH_CL(W[f_index][1], W[f_index][1], absw3_II);
		if (absw3_I[0] < absw3_II[0]) { R_index[f_index] = 0; }
		else { R_index[f_index] = 1; }
#endif // MMSE_WSEL
	}
}
