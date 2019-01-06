#include "const-function-declerations.h"			
extern errno_t err;

void WHadamard_initialize(double(*transmitted_signal)[BS][BURST][2], double(*wh_sequence)[Np][N][2], double(*transmitted_signal_actual)[BURST][2]){

	if (pow(2, Np_POWER2) != Np || Np > 256){ printf("Make Np power of 2 and less than 300.\nAlso check if 2^Np_POWER2 = Np.\n"); getchar(); exit(1); }
	double code[Np][Np], codei[Np][Np];
	int h_init, v_init, h, v, Pow, i, j, ts, bs, step, matrix_quarter, h_index, v_index, f_index, m, n;
	int index;
	double time_signal[BS][Nf][2];
	double FFT_signal[BS][Nf][2];
	double IFFT_signal[BS][Nf][2];
	/*  i, j: auxiliary indexes
		h: for horizontal
		v: for vertical
		m: preamble index
		n: frequency sampling index
		ts: time slot
		bs: base station
		*/

	for (i = 0; i < Np; i++){ for (j = 0; j < Np; j++){ code[i][j] = 0; } }

	code[0][0] = (double)OneBySqrt2;

	for (step = 0; step < Np_POWER2; step++){
		Pow = intPow(2, step);
		for (matrix_quarter = 0; matrix_quarter < 4; matrix_quarter++){
			h_init = matrix_quarter % 2;
			v_init = matrix_quarter > 1 ? 1 : 0;
			h_init = h_init*Pow;
			v_init = v_init*Pow;
			for (h_index = 0; h_index < Pow; h_index++){
				for (v_index = 0; v_index < Pow; v_index++){
					h = h_init + h_index;
					v = v_init + v_index;
					code[h][v] = (matrix_quarter == 3) ? -code[h_index][v_index] : code[h_index][v_index];
				}
			}
		}
	}	/* end :: Hadamard matrix generation */

	for (h = 0; h < Np; h++){ for (v = 0; v < Np; v++){ codei[h][v] = code[h][Np - v - 1]; } } /* complex part. real part is orthogonal all by itself(meaning even if u dont add complex part it will still be orthogonal). */

	for (f_index = 0; f_index < N; f_index++){	/* Only use first two sequences. same for all sub-carriers. */
		for (m = 0; m < Np; m++){
			wh_sequence[0][m][f_index][0] = code[WH_CS1][m];	/* wh_sequence[transmitter index][symbol index][subcarrier index] */
			wh_sequence[0][m][f_index][1] = codei[WH_CS1][m];
			wh_sequence[1][m][f_index][0] = code[WH_CS2][m];
			wh_sequence[1][m][f_index][1] = codei[WH_CS2][m];
			wh_sequence[2][m][f_index][0] = code[WH_CS3][m];
			wh_sequence[2][m][f_index][1] = codei[WH_CS3][m];
		}
	}	/* end :: chose codes to use */

	/* ~~~~~~~~~~~~~~Sending~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	for (m = 0; m < Np; m++){//1	
		for (n = 0; n < Nf; n++) { for (bs = 0; bs < BS; bs++)	{ FFT_signal[bs][n][0] = 0.0; FFT_signal[bs][n][1] = 0.0; } }

		for (n = 0; n < N / 2; n++) {	/* ?Why */
			for (bs = 0; bs < BS; bs++)	{
				FFT_signal[bs][n][0] = wh_sequence[bs][m][n][0];
				FFT_signal[bs][n][1] = wh_sequence[bs][m][n][1];
			}
		}

		for (n = N / 2; n < N; n++) {	/* ?Why */
			for (bs = 0; bs < BS; bs++)	{
				FFT_signal[bs][n + Nf - N][0] = wh_sequence[bs][m][n][0];
				FFT_signal[bs][n + Nf - N][1] = wh_sequence[bs][m][n][1];
			}
		}

		for (bs = 0; bs < BS; bs++)	{
			(Fourier_MOD == ON) ? IFFT(FFT_signal[bs], IFFT_signal[bs]) : IDFT(FFT_signal[bs], IFFT_signal[bs]);
		}

		for (n = 0; n < Nf; n++){
			for (bs = 0; bs < BS; bs++)	{
				time_signal[bs][n][0] = IFFT_signal[bs][n][0];
				time_signal[bs][n][1] = IFFT_signal[bs][n][1];
			}
		}

		index = m * SAMPLEN;

		/**/for (i = 0; i < GI; i++){
			for (bs = 0; bs < BS; bs++)	{
				for (ts = 0; ts < TS; ts++)	{
					transmitted_signal[ts][bs][i + index][0] = time_signal[bs][Nf - GI + i][0];
					transmitted_signal[ts][bs][i + index][1] = time_signal[bs][Nf - GI + i][1];
				}
				transmitted_signal_actual[bs][i + index][0] = time_signal[bs][Nf - GI + i][0];
				transmitted_signal_actual[bs][i + index][1] = time_signal[bs][Nf - GI + i][1];
			}
		}
		for (i = 0; i < Nf; i++){
			for (bs = 0; bs < BS; bs++)	{
				for (ts = 0; ts < TS; ts++)	{
					transmitted_signal[ts][bs][i + GI + index][0] = time_signal[bs][i][0];
					transmitted_signal[ts][bs][i + GI + index][1] = time_signal[bs][i][1];
				}
				transmitted_signal_actual[bs][i + GI + index][0] = time_signal[bs][i][0];
				transmitted_signal_actual[bs][i + GI + index][1] = time_signal[bs][i][1];
			}
		}
	}
}	/* end :: WH_initialize */

void WH_channel_estimation(double(*DFT_signal)[TS][N][2], double(*H)[TS][BS][2], double(*wh_sequence)[Np][N][2]){
	double temp[2], whsum;
	temp[0] = temp[1] = 0; whsum = 0;
	for (int tx = 0; tx < BS; tx++){
		for (int rx = 0; rx < TS; rx++){
			for (int f_index = 0; f_index < N; f_index++){
				for (int s_index = 0; s_index < Np; s_index++){	/* This 0 and 1 below is for Q-I */
					temp[0] += DFT_signal[s_index][rx][f_index][0] * wh_sequence[tx][s_index][f_index][0] + DFT_signal[s_index][rx][f_index][1] * wh_sequence[tx][s_index][f_index][1];
					temp[1] += DFT_signal[s_index][rx][f_index][1] * wh_sequence[tx][s_index][f_index][0] - DFT_signal[s_index][rx][f_index][0] * wh_sequence[tx][s_index][f_index][1];
					whsum += sqr(wh_sequence[tx][s_index][f_index][0]) + sqr(wh_sequence[tx][s_index][f_index][1]);
				}
				H[f_index][rx][tx][0] = temp[0] / whsum;
				H[f_index][rx][tx][1] = temp[1] / whsum;
				temp[0] = 0; temp[1] = 0; whsum = 0;
			}
		}
	}
}

/* For integer power calculation: */
int intPow(int x, int p)
{
	if (p == 0) return 1;
	if (p == 1) return x;

	int tmp = intPow(x, p / 2);
	if (p % 2 == 0) return tmp * tmp;
	else return x * tmp * tmp;
}