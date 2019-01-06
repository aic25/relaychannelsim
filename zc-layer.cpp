#include "const-function-declerations.h"
double zc_signal_a[BS][Np][2], zc_signal_b[UTx][N][2], freq_zc_signal[N][2];
extern double wl[TS][BS][PATH][2];

void Zadoff_Chu_layer_initialize(){	/* Think of it like a matrix on frequency and time axises. */
	int i, j, u, utx;
	double temp;

	/* Using Eulers formula */
	for (u = 0; u < BS; u++){	/* BS = 2 */
		for (i = 0; i < Np; i++){	/* Np = 2 */
			temp = (double)-2.0 * PI * (u + 1) * i / Np;	/* Recall two different sequences are just phase shifted versions of each other */
			zc_signal_a[u][i][0] = cos(temp);
			zc_signal_a[u][i][1] = sin(temp);
		}
	}
	for (utx = 0; utx < UTx; utx++){	/* UTx = 1 */	/* I wanted to say but those a and b does not seem like sequences for different users. The loop numbers are decided differently */
		for (i = 0; i < N; i++){
			j = utx * -N / 4 + i;
			temp = (double)-PI * j * j * 11.0 / N;
			zc_signal_b[utx][i][0] = cos(temp);
			zc_signal_b[utx][i][1] = sin(temp);
		}
	}
}	/* Those a and b could be for vertical and horizontal. */

void Zadoff_Chu_generator_layer1(double(*transmitted_signal)[BS][BURST][2]){

	int i, Nl, utx, n, index, u, m, tx;	
	double zc_signal_ab[N][2];
	double time_signal[Nf][2];
	double FFT_signal[N][2], FFT_signal2[Nf][2];
	double IFFT_signal[Nf][2];

	Nl = Np * Nf;	/* Total length of training data, what bothers me is that he used N up there */

	for (i = 0; i < N; i++){

		zc_signal_ab[i][0] = zc_signal_b[0][i][0];
		zc_signal_ab[i][1] = zc_signal_b[0][i][1];

	}

	FFT_N(zc_signal_ab, FFT_signal);	/* ?Why :: maybe can produce the sequence in time domain and using the property of being able to keep CAZAC property in both domains*/

	for (n = 0; n < N; n++) {

		freq_zc_signal[n][0] = FFT_signal[(n + N / 2) % N][0];	/* ? Why this swap */
		freq_zc_signal[n][1] = FFT_signal[(n + N / 2) % N][1];	/* Global variable */

	}

	for (u = 0; u < BS; u++){	/* Loop over base stations */
		for (utx = 0; utx < UTx; utx++){	/* Loop over antennas */
			tx = UTx * u + utx;	/* ?Seems more like a stream index counter */
			for (m = 0; m < Np; m++){//1	/* Loop over training OFDM symbols. I guess the loop is over time, which means for one subcarrier :: ? */
				for (i = 0; i < N; i++){	/* Remember training signal blocks? this could be matrixwise convolution */

					zc_signal_ab[i][0] = zc_signal_a[u][m][0] * zc_signal_b[utx][i][0] - zc_signal_a[u][m][1] * zc_signal_b[utx][i][1];
					zc_signal_ab[i][1] = zc_signal_a[u][m][0] * zc_signal_b[utx][i][1] + zc_signal_a[u][m][1] * zc_signal_b[utx][i][0];

				}

				FFT_N(zc_signal_ab, FFT_signal);	/* Switch to frequency domain */

				/* WH :: copy below */

				for (n = 0; n < Nf; n++) {
					FFT_signal2[n][0] = 0.0;
					FFT_signal2[n][1] = 0.0;
				}

				for (n = 0; n < N / 2; n++) {	/* ?Why */
					FFT_signal2[n][0] = FFT_signal[n][0];
					FFT_signal2[n][1] = FFT_signal[n][1];
				}
				/* Wtih Nf = N, FFT_signal2 = FFT_signal */

				for (n = N / 2; n < N; n++) {	/* ?Why */
					FFT_signal2[n + Nf - N][0] = FFT_signal[n][0];
					FFT_signal2[n + Nf - N][1] = FFT_signal[n][1];
				}

				//IFFT(FFT_signal2, IFFT_signal);	/* Swapped two parts in frequency domain and go back to time domain */
				IDFT(FFT_signal2, IFFT_signal);

				for (n = 0; n < Nf; n++){
					time_signal[n][0] = IFFT_signal[n][0];
					time_signal[n][1] = IFFT_signal[n][1];
				}

				/******************************* GI insert *******************************/
				index = m * SAMPLEN;

				for (i = 0; i < GI; i++){
					transmitted_signal[0][tx][i + index][0] = time_signal[Nf - GI + i][0];
					transmitted_signal[0][tx][i + index][1] = time_signal[Nf - GI + i][1];
					
				}
				for (i = 0; i < Nf; i++){
					transmitted_signal[0][tx][i + GI + index][0] = time_signal[i][0];
					transmitted_signal[0][tx][i + GI + index][1] = time_signal[i][1];
				}
			}
		}
	}
}

void Zadoff_Chu_generator_layer2(double(*transmitted_signal)[BS][BURST][2]){//4
	int i, Nl, utx, n, index, u, m, tx;
	double time_signal[Nf][2];
	double zc_signal_ab[N][2];
	double FFT_signal[N][2], FFT_signal2[Nf][2];
	double IFFT_signal[Nf][2];

	Nl = Np * Nf;	/* Not being used. */

	for (i = 0; i < N; i++){
		zc_signal_ab[i][0] = zc_signal_b[0][i][0];
		zc_signal_ab[i][1] = zc_signal_b[0][i][1];
	}
	FFT_N(zc_signal_ab, FFT_signal);
	for (n = 0; n < N; n++) {	/* What is this for? */
		freq_zc_signal[n][0] = FFT_signal[(n + N / 2) % N][0];
		freq_zc_signal[n][1] = FFT_signal[(n + N / 2) % N][1];
	}
	for (u = 0; u < BS; u++){
		for (utx = 0; utx < UTx; utx++){
			tx = UTx * u + utx;
			for (m = 0; m < Np; m++){
				for (i = 0; i < N; i++){
					zc_signal_ab[i][0] = zc_signal_a[u][m][0] * zc_signal_b[utx][i][0] - zc_signal_a[u][m][1] * zc_signal_b[utx][i][1];
					zc_signal_ab[i][1] = zc_signal_a[u][m][0] * zc_signal_b[utx][i][1] + zc_signal_a[u][m][1] * zc_signal_b[utx][i][0];
				}

				FFT_N(zc_signal_ab, FFT_signal);

				for (n = 0; n < Nf; n++) {
					FFT_signal2[n][0] = 0.0;
					FFT_signal2[n][1] = 0.0;
				}

				for (n = 0; n < N / 2; n++) {
					FFT_signal2[n][0] = FFT_signal[n][0];
					FFT_signal2[n][1] = FFT_signal[n][1];
				}
				for (n = N / 2; n < N; n++) {
					FFT_signal2[n + Nf - N][0] = FFT_signal[n][0];
					FFT_signal2[n + Nf - N][1] = FFT_signal[n][1];
				}

				IFFT(FFT_signal2, IFFT_signal); /* IFFT */

				for (n = 0; n < Nf; n++){
					time_signal[n][0] = IFFT_signal[n][0];
					time_signal[n][1] = IFFT_signal[n][1];
				}
				/******************************* GI insert *******************************/
				/* This second version is only for the second time slot. Being identical to the first above means relay also knows the sequence. */
				index = m * SAMPLEN;
				for (i = 0; i < GI; i++){
					transmitted_signal[1][tx][i + index][0] = time_signal[Nf - GI + i][0];
					transmitted_signal[1][tx][i + index][1] = time_signal[Nf - GI + i][1];

				}
				for (i = 0; i < Nf; i++){
					transmitted_signal[1][tx][i + GI + index][0] = time_signal[i][0];
					transmitted_signal[1][tx][i + GI + index][1] = time_signal[i][1];
				}
			}
		}
	}
}




