#include "const.h"

extern errno_t err;

void WHadamard_initialize(double(*transmitted_signal)[BS][BURST][2], double(*wh_sequence)[Np][N]){

	if (pow(2, Np_POWER2) != Np || Np > 256){
		printf("Make Np power of 2 and less than 300.\nAlso check if 2^Np_POWER2 = Np.\n");
		getchar();
		exit(1);
	}
	int code[Np][Np];
	int h_init, v_init, h, v, Pow, m, n;
	for (size_t i = 0; i < Np; i++){
		for (size_t j = 0; j < Np; j++){
			code[i][j] = 0;
		}
	}
	code[0][0] = 1;
	for (int step = 0; step < Np_POWER2; step++){
		Pow = intPow(2, step);
		for (int matrix_quarter = 0; matrix_quarter < 4; matrix_quarter++){
			h_init = matrix_quarter % 2;
			v_init = matrix_quarter > 1 ? 1 : 0;
			h_init = h_init*Pow;
			v_init = v_init*Pow;
			for (int h_index = 0; h_index < Pow; h_index++){
				for (int v_index = 0; v_index < Pow; v_index++){
					h = h_init + h_index;
					v = v_init + v_index;
					code[h][v] = (matrix_quarter == 3) ? -code[h_index][v_index] : code[h_index][v_index];
				}
			}
		}
	}	/* end :: Hadamard matrix generation */
	for (int f_index = 0; f_index < N; f_index++){	/* Only use first two sequences. same for all sub-carriers. */
		for (int m = 0; m < Np; m++){
			wh_sequence[0][m][f_index] = code[WH_CS1][m];	/* wh_sequence[transmitter index][symbol index][subcarrier index] */
			wh_sequence[1][m][f_index] = code[WH_CS2][m];
		}
	}	/* end :: chose codes to use */

	int index;
	double time_signal[Nf][2], time_signal2[Nf][2];
	double FFT_signal[Nf][2], FFT_signal2[Nf][2];
	double IFFT_signal[Nf][2], IFFT_signal2[Nf][2];

	/* ~~~~~~~~~~~~~~Sending~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	for (int m = 0; m < Np; m++){//1	
		for (int n = 0; n < Nf; n++) {
			FFT_signal2[n][0] = 0.0;
			FFT_signal[n][1] = 0.0;
			FFT_signal2[n][1] = 0.0;
			FFT_signal[n][0] = 0.0;
		}

		for (int n = 0; n < N / 2; n++) {	/* ?Why */
			FFT_signal2[n][0] = wh_sequence[0][m][n];
			FFT_signal2[n][1] = 0;
			FFT_signal[n][0] = wh_sequence[1][m][n];
			FFT_signal[n][1] = 0;
		}

		for (int n = N / 2; n < N; n++) {	/* ?Why */
			FFT_signal2[n + Nf - N][0] = wh_sequence[0][m][n];
			FFT_signal2[n + Nf - N][1] = 0;
			FFT_signal[n + Nf - N][0] = wh_sequence[1][m][n];
			FFT_signal[n + Nf - N][1] = 0;
		}

		IFFT(FFT_signal2, IFFT_signal2);	/* Swapped two parts in frequency domain and go back to time domain */
		IFFT(FFT_signal, IFFT_signal);

		for (int n = 0; n < Nf; n++){
			time_signal[n][0] = IFFT_signal[n][0];
			time_signal[n][1] = IFFT_signal[n][1];
			time_signal2[n][0] = IFFT_signal2[n][0];
			time_signal2[n][1] = IFFT_signal2[n][1];
		}

		index = m * SAMPLEN;

		for (int i = 0; i < GI; i++){
			transmitted_signal[0][0][i + index][0] = time_signal2[Nf - GI + i][0];
			transmitted_signal[0][0][i + index][1] = time_signal2[Nf - GI + i][1];

			transmitted_signal[0][1][i + index][0] = time_signal[Nf - GI + i][0];
			transmitted_signal[0][1][i + index][1] = time_signal[Nf - GI + i][1];


			transmitted_signal[1][0][i + index][0] = time_signal2[Nf - GI + i][0];
			transmitted_signal[1][0][i + index][1] = time_signal2[Nf - GI + i][1];

			transmitted_signal[1][1][i + index][0] = time_signal[Nf - GI + i][0];
			transmitted_signal[1][1][i + index][1] = time_signal[Nf - GI + i][1];

		}
		for (int i = 0; i < Nf; i++){
			transmitted_signal[0][0][i + GI + index][0] = time_signal2[i][0];
			transmitted_signal[0][0][i + GI + index][1] = time_signal2[i][1];

			transmitted_signal[0][1][i + GI + index][0] = time_signal[i][0];
			transmitted_signal[0][1][i + GI + index][1] = time_signal[i][1];


			transmitted_signal[1][0][i + GI + index][0] = time_signal2[i][0];
			transmitted_signal[1][0][i + GI + index][1] = time_signal2[i][1];

			transmitted_signal[1][1][i + GI + index][0] = time_signal[i][0];
			transmitted_signal[1][1][i + GI + index][1] = time_signal[i][1];
		}
	}
}	/* end :: WH_initialize */

void WH_channel_estimation(double(*DFT_signal)[TS][N][2], double(*H)[TS][BS][2], double(*wh_sequence)[Np][N]){
	double temp[2], whsum;
	temp[0] = temp[1] = 0; whsum = 0;
	for (int tx = 0; tx < BS; tx++){
		for (int rx = 0; rx < TS; rx++){
			for (int f_index = 0; f_index < N; f_index++){
				for (int s_index = 0; s_index < Np; s_index++){	/* This 0 and 1 below is for Q-I */
					temp[0] += DFT_signal[s_index][rx][f_index][0] * wh_sequence[tx][s_index][f_index];
					temp[1] += DFT_signal[s_index][rx][f_index][1] * wh_sequence[tx][s_index][f_index];
					whsum += sqr(wh_sequence[tx][s_index][f_index]);
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