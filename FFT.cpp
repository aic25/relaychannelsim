#include "const-function-declerations.h"	  
double Wf[Nf / 2][2];
double Wf_N[N / 2][2];

void FFT_initialize(){
	int i;
	double phase;

	Wf[0][0] = 1.0;
	Wf[0][1] = 0.0;
	for (i = 1; i < Nf / 2; i++){
		phase = (double)i*2.0*PI / Nf;
		Wf[i][0] = cos(phase);
		Wf[i][1] = sin(phase);
	}

	Wf_N[0][0] = 1.0;
	Wf_N[0][1] = 0.0;
	for (i = 1; i < N / 2; i++){
		phase = (double)i*2.0*PI / N;
		Wf_N[i][0] = cos(phase);
		Wf_N[i][1] = sin(phase);
	}
}

void IFFT(double(*frequency_signal)[2], double(*time_signal)[2]){
	int i, j, j1, j2, it, iter, k, kk, xp, xp2, xp3, carrier, carrier1, carrier2, m;
	double temp[2], temp1[2], temp2[2];

	for (iter = 0, i = Nf;; iter++){
		if ((i >>= 1) == 0) break;
	}

	xp2 = Nf;
	xp3 = 1;
	for (it = 0; it < iter; it++){
		xp = xp2;
		xp2 = (xp2 >> 1);
		if (it != 0){
			xp3 = (xp3 << 1);
		}
		for (k = 0; k < xp2; k++){
			if (k != 0){
				carrier = k*xp3;
			}
			i = k - xp;
			for (j = xp; j <= Nf; j += xp){
				carrier1 = j + i;
				carrier2 = carrier1 + xp2;
				for (kk = 0; kk < 2; kk++){
					temp1[kk] = frequency_signal[carrier1][kk];
					temp2[kk] = frequency_signal[carrier2][kk];

					temp[kk] = temp1[kk] - temp2[kk];
				}
				frequency_signal[carrier1][0] = temp1[0] + temp2[0];
				frequency_signal[carrier1][1] = temp1[1] + temp2[1];

				if (k == 0){
					frequency_signal[carrier2][0] = temp[0];
					frequency_signal[carrier2][1] = temp[1];
				}
				else {
					frequency_signal[carrier2][0] = temp[0] * Wf[carrier][0] - temp[1] * Wf[carrier][1];
					frequency_signal[carrier2][1] = temp[0] * Wf[carrier][1] + temp[1] * Wf[carrier][0];
				}
			}
		}
	}
	/********Bit Reversal to {0,1, ...... ,Nf} **********/
	j1 = Nf / 2;
	j2 = Nf - 1;
	j = 1;

	for (i = 1; i < j2; i++){
		if (i < j){
			carrier1 = i - 1;
			carrier2 = j - 1;
			for (kk = 0; kk < 2; kk++){
				temp[kk] = frequency_signal[carrier1][kk];
				frequency_signal[carrier1][kk] = frequency_signal[carrier2][kk];
				frequency_signal[carrier2][kk] = temp[kk];
			}
		}
		k = j1;

		while (k < j){
			j -= k;
			k >>= 1;
		}

		j += k;
	}

	for (m = 0; m < Nf; m++){
		time_signal[m][0] = sqrt((double)1.0 / Nf) * frequency_signal[m][0];
		time_signal[m][1] = sqrt((double)1.0 / Nf) * frequency_signal[m][1];
	}
}


void FFT(double(*time_signal)[2], double(*frequency_signal)[2]){
	int i, j, j1, j2, it, iter, k, kk, xp, xp2, xp3, carrier, carrier1, carrier2, carrier_index;
	double temp[2], temp1[2], temp2[2];

	for (iter = 0, i = Nf;; iter++){
		if ((i >>= 1) == 0) break;
	}

	xp2 = Nf;
	xp3 = 1;
	for (it = 0; it < iter; it++){
		xp = xp2;
		xp2 = (xp2 >> 1);
		if (it != 0){
			xp3 = (xp3 << 1);
		}
		for (k = 0; k < xp2; k++){
			if (k != 0){
				carrier = k*xp3;
			}
			i = k - xp;
			for (j = xp; j <= Nf; j += xp){
				carrier1 = j + i;
				carrier2 = carrier1 + xp2;
				for (kk = 0; kk < 2; kk++){
					temp1[kk] = time_signal[carrier1][kk];
					temp2[kk] = time_signal[carrier2][kk];

					temp[kk] = temp1[kk] - temp2[kk];
				}
				time_signal[carrier1][0] = temp1[0] + temp2[0];
				time_signal[carrier1][1] = temp1[1] + temp2[1];

				if (k == 0){
					time_signal[carrier2][0] = temp[0];
					time_signal[carrier2][1] = temp[1];
				}
				else {
					time_signal[carrier2][0] = temp[0] * Wf[carrier][0] + temp[1] * Wf[carrier][1];
					time_signal[carrier2][1] = -temp[0] * Wf[carrier][1] + temp[1] * Wf[carrier][0];
				}
			}
		}
	}
	/********Bit Reversal to {0,1, ...... ,Nf} **********/
	j1 = Nf / 2;
	j2 = Nf - 1;
	j = 1;

	for (i = 1; i < j2; i++){
		if (i < j){
			carrier1 = i - 1;
			carrier2 = j - 1;
			for (kk = 0; kk < 2; kk++){
				temp[kk] = time_signal[carrier1][kk];
				time_signal[carrier1][kk] = time_signal[carrier2][kk];
				time_signal[carrier2][kk] = temp[kk];
			}
		}
		k = j1;

		while (k < j){
			j -= k;
			k >>= 1;
		}

		j += k;
	}

	for (carrier_index = 0; carrier_index < Nf; carrier_index++){
		for (kk = 0; kk < 2; kk++){
			frequency_signal[carrier_index][kk] = sqrt((double)1.0 / Nf) * time_signal[carrier_index][kk];
		}
	}
}

void IFFT_N(double(*frequency_signal)[2], double(*time_signal)[2]){
	int i, j, j1, j2, it, iter, k, kk, xp, xp2, xp3, carrier, carrier1, carrier2, m;
	double temp[2], temp1[2], temp2[2];

	for (iter = 0, i = N;; iter++){
		if ((i >>= 1) == 0) break;
	}

	xp2 = N;
	xp3 = 1;
	for (it = 0; it < iter; it++){
		xp = xp2;
		xp2 = (xp2 >> 1);
		if (it != 0){
			xp3 = (xp3 << 1);
		}
		for (k = 0; k < xp2; k++){
			if (k != 0){
				carrier = k*xp3;
			}
			i = k - xp;
			for (j = xp; j <= N; j += xp){
				carrier1 = j + i;
				carrier2 = carrier1 + xp2;
				for (kk = 0; kk < 2; kk++){
					temp1[kk] = frequency_signal[carrier1][kk];
					temp2[kk] = frequency_signal[carrier2][kk];

					temp[kk] = temp1[kk] - temp2[kk];
				}
				frequency_signal[carrier1][0] = temp1[0] + temp2[0];
				frequency_signal[carrier1][1] = temp1[1] + temp2[1];

				if (k == 0){
					frequency_signal[carrier2][0] = temp[0];
					frequency_signal[carrier2][1] = temp[1];
				}
				else {
					frequency_signal[carrier2][0] = temp[0] * Wf_N[carrier][0] - temp[1] * Wf_N[carrier][1];
					frequency_signal[carrier2][1] = temp[0] * Wf_N[carrier][1] + temp[1] * Wf_N[carrier][0];
				}
			}
		}
	}
	/********Bit Reversal to {0,1, ...... ,N} **********/
	j1 = N / 2;
	j2 = N - 1;
	j = 1;

	for (i = 1; i < j2; i++){
		if (i < j){
			carrier1 = i - 1;
			carrier2 = j - 1;
			for (kk = 0; kk < 2; kk++){
				temp[kk] = frequency_signal[carrier1][kk];
				frequency_signal[carrier1][kk] = frequency_signal[carrier2][kk];
				frequency_signal[carrier2][kk] = temp[kk];
			}
		}
		k = j1;

		while (k < j){
			j -= k;
			k >>= 1;
		}

		j += k;
	}

	for (m = 0; m < N; m++){
		time_signal[m][0] = sqrt((double)1.0 / N) * frequency_signal[m][0];
		time_signal[m][1] = sqrt((double)1.0 / N) * frequency_signal[m][1];
	}
}

void FFT_N(double(*time_signal)[2], double(*frequency_signal)[2]){ //for signal at zc channel estimation
	int i, j, j1, j2, it, iter, k, kk, xp, xp2, xp3, carrier, carrier1, carrier2, carrier_index;
	double temp[2], temp1[2], temp2[2];

	for (iter = 0, i = N;; iter++){
		if ((i >>= 1) == 0) break;
	}

	xp2 = N;
	xp3 = 1;
	for (it = 0; it < iter; it++){
		xp = xp2;
		xp2 = (xp2 >> 1);
		if (it != 0){
			xp3 = (xp3 << 1);
		}
		for (k = 0; k < xp2; k++){
			if (k != 0){
				carrier = k*xp3;
			}
			i = k - xp;
			for (j = xp; j <= N; j += xp){
				carrier1 = j + i;
				carrier2 = carrier1 + xp2;
				for (kk = 0; kk < 2; kk++){
					temp1[kk] = time_signal[carrier1][kk];
					temp2[kk] = time_signal[carrier2][kk];

					temp[kk] = temp1[kk] - temp2[kk];
				}
				time_signal[carrier1][0] = temp1[0] + temp2[0];
				time_signal[carrier1][1] = temp1[1] + temp2[1];

				if (k == 0){
					time_signal[carrier2][0] = temp[0];
					time_signal[carrier2][1] = temp[1];
				}
				else {
					time_signal[carrier2][0] = temp[0] * Wf_N[carrier][0] + temp[1] * Wf_N[carrier][1];
					time_signal[carrier2][1] = -temp[0] * Wf_N[carrier][1] + temp[1] * Wf_N[carrier][0];
				}
			}
		}
	}
	/********Bit Reversal to {0,1, ...... ,N} **********/
	j1 = N / 2;
	j2 = N - 1;
	j = 1;

	for (i = 1; i < j2; i++){
		if (i < j){
			carrier1 = i - 1;
			carrier2 = j - 1;
			for (kk = 0; kk < 2; kk++){
				temp[kk] = time_signal[carrier1][kk];
				time_signal[carrier1][kk] = time_signal[carrier2][kk];
				time_signal[carrier2][kk] = temp[kk];
			}
		}
		k = j1;
		while (k < j){
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	for (carrier_index = 0; carrier_index < N; carrier_index++){
		for (kk = 0; kk < 2; kk++){
			frequency_signal[carrier_index][kk] = sqrt((double)1.0 / N) * time_signal[carrier_index][kk];
		}
	}
}
