#include "const.h"

extern double
wl[TS][BS][PATH*DELAY + RS_PATH*RS_DELAY - 1][2];
extern double
zc_signal_a[BS][Np][2],
zc_signal_b[UTx][N][2],
freq_zc_signal[N][2];

int timing_delta[TS];

void ideal_channel_estimation4all_ch(double(*H)[TS][BS][2]){
	int ts, n, bs;
	double Ht[Nf][2], ht[Nf][2], coef;

	coef = sqrt((double)Nf);
	for (bs = 0; bs < BS; bs++){
		for (ts = 0; ts < TS; ts++){
			if ((ts == 0) || RS_DF == ON)
			{
				for (n = 0; n < Nf; n++){
#ifdef MULTIPATH
					if ((n % (DELAY) == 0) && (n / DELAY) < PATH){
						ht[n][0] = wl[ts][bs][n / DELAY][0];
						ht[n][1] = wl[ts][bs][n / DELAY][1];
					}
					else{
						ht[n][0] = ht[n][1] = 0.0;
					}
#else	   
					if ((n % (DELAY) == 0) && (n / DELAY) < 1){
						ht[n][0] = wl[ts][bs][n / DELAY][0];
						ht[n][1] = wl[ts][bs][n / DELAY][1];
					}
#endif
				}
			}
			else {
				for (n = 0; n < Nf; n++){
#ifdef MULTIPATH
					if ((n % (1) == 0) && (n / 1) < (PATH*DELAY + RS_PATH*RS_DELAY - 1)){
						ht[n][0] = wl[ts][bs][n / 1][0];
						ht[n][1] = wl[ts][bs][n / 1][1];
					}
					else{
						ht[n][0] = ht[n][1] = 0.0;
					}
#else	   
					if ((n % (DELAY) == 0) && (n / DELAY) < 1){
						ht[n][0] = wl[ts][bs][n / DELAY][0];
						ht[n][1] = wl[ts][bs][n / DELAY][1];
					}
#endif
				}
			}

			FFT(ht, Ht);

			for (n = 0; n < N; n++){
				H[n][ts][bs][0] = coef * Ht[(n + Nf - N / 2) % Nf][0];	/* ? */
				H[n][ts][bs][1] = coef * Ht[(n + Nf - N / 2) % Nf][1];
			}
		}
	}
}

void channel_estimation4ZC_4all_ch(double(*received_signal)[BURST][2], double(*H)[TS][BS][2]){
	int data_number, ts, n, tx, utx, bs;
	double DFT_signal[TS][N][2], correlation[N][2], Ht[N][2], coef;
	double Ht2[N][2], ht[N][2], ht2[N][2];

	coef = 1.0 / Np;
	for (ts = 0; ts < TS; ts++)
	{
		for (bs = 0; bs < BS; bs++)
		{
			for (n = 0; n < N; n++)
				correlation[n][0] = correlation[n][1] = 0.0;
			for (data_number = 0; data_number < Np; data_number++)
			{
				DFT(ts, data_number, received_signal, DFT_signal);
				for (n = 0; n < N; n++)
				{
					correlation[n][0] += coef * (DFT_signal[ts][n][0] * zc_signal_a[bs][data_number][0] - DFT_signal[ts][n][1] * -zc_signal_a[bs][data_number][1]);
					correlation[n][1] += coef * (DFT_signal[ts][n][0] * -zc_signal_a[bs][data_number][1] + DFT_signal[ts][n][1] * zc_signal_a[bs][data_number][0]);
				}
			}
			for (n = 0; n < N; n++)
			{
				Ht[n][0] = correlation[n][0] * freq_zc_signal[n][0] - correlation[n][1] * -freq_zc_signal[n][1];
				Ht[n][1] = correlation[n][0] * -freq_zc_signal[n][1] + correlation[n][1] * freq_zc_signal[n][0];
			}

			for (n = 0; n < N; n++)
			{
				for (n = 0; n < N; n++)
				{
					Ht2[n][0] = Ht[(n + N / 2) % N][0];
					Ht2[n][1] = Ht[(n + N / 2) % N][1];
				}
				IFFT_N(Ht2, ht);
			}

			for (utx = 0; utx < UTx; utx++)
			{
				tx = UTx * bs + utx;

				for (n = 0; n < N; n++)
				{
#if(DC_CUT==OFF)	/* DC_CUT is not defined in const.h!! */
					if (n < TAPLN / 2)
					{
#endif					
						ht2[n][0] = ht[(n + utx*N / 4 + N) % N][0];
						ht2[n][1] = ht[(n + utx*N / 4 + N) % N][1];

#if(DC_CUT==OFF)
					}
					else
					{
						ht2[n][0] = ht2[n][1] = 0.0;
					}
#endif
				}

				FFT_N(ht2, Ht2);
				for (n = 0; n < N; n++)
				{
					H[n][ts][tx][0] = Ht2[(n + N / 2) % N][0];
					H[n][ts][tx][1] = Ht2[(n + N / 2) % N][1];
				}
			}
		}
	}
}