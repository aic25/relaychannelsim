#include "const.h"			
extern double CNR;
extern double signal_power, noise_power;

double wl[TS][BS][PATH*DELAY + RS_PATH*RS_DELAY - 1][2];	/* Used for pulling out channel coefficients. */

/* received signal = transmitted_sigvec */
void multipath_propagation(double(*transmitted_sigmat)[BS][BURST][2], double(*received_signal)[BURST][2]){
	int time, path_index, index, iq, ts, bs;
	double Div, K;
	double temp[2];
	double w[PATH][2];	/* w? */
	double w_extended[PATH*DELAY][2];	/* w? */
	double w_r[RS_PATH][2];	/* w? */
	double w_r_extended[RS_PATH*RS_DELAY][2];	/* w? */
	double w_rb[PATH + RS_PATH - 1][2];	/* w? */
	double w_rb_extended[PATH*DELAY + RS_PATH*RS_DELAY - 1][2];	/* w? */
	double rs_rcvd_sig[BS][BURST][2];
	double noise[2], noise_amplitude;
	noise_amplitude = sqrt(pow(10.0, -RLY_CNR / 10.0) / CODING_RATE);

	for (ts = 0; ts < TS; ts++){
		for (time = 0; time < BURST; time++) {
			for (iq = 0; iq < 2; iq++) received_signal[ts][time][iq] = 0.0;
		}
	}

	for (ts = 0; ts < TS; ts++){
		if (ts == 0){
			for (bs = 0; bs < BS; bs++){
				for (time = 0; time < BURST; time++) {
#ifdef PROPAGATION
					//fading_process(time, w, (ts == 1) ? 0 : (0 - LOSS));	/* Most likely generates path coefficients. Time is the location on the long data stream. */
					fading_process(time, w, (0 - LOSS));

					for (path_index = 0; path_index < PATH + RS_PATH - 1; path_index++) {	/* To be used in other parts of the code to drag channel coefficents */
						if (path_index < PATH)
						{
							wl[ts][bs][path_index][0] = w[path_index][0];	/* kk was preferred here */
							wl[ts][bs][path_index][1] = w[path_index][1];
						}
						else{
							wl[ts][bs][path_index][0] = 0;	/* kk was preferred here */
							wl[ts][bs][path_index][1] = 0;
						}
					}
#ifdef MULTIPATH
					for (path_index = 0; path_index < PATH; path_index++) {	/* Signal goes thru channel */
						index = time - path_index * DELAY;

						if (index >= 0) {
							for (iq = 0; iq < 2; iq++) temp[iq] = transmitted_sigmat[ts][bs][index][iq];	/* Why use temp here? */

							received_signal[ts][time][0] += w[path_index][0] * temp[0] - w[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							received_signal[ts][time][1] += w[path_index][1] * temp[0] + w[path_index][0] * temp[1];
						}
					}
#else		  
					received_signal[ts][time][0] += PATH * w[0][0] * transmitted_sigmat[ts][bs][time][0] - w[0][1] * transmitted_sigmat[ts][bs][time][1];
					received_signal[ts][time][1] += PATH * w[0][1] * transmitted_sigmat[ts][bs][time][0] + w[0][0] * transmitted_sigmat[ts][bs][time][1];
#endif			
#else				   
					received_signal[ts][time][0] = transmitted_sigmat[ts][0][time][0];
					received_signal[ts][time][1] = transmitted_sigmat[ts][0][time][1];
#endif
					signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
				}
			}
		}
		else if (ts == 1){
#ifdef TWO_STG_RLY
			for (bs = 0; bs < BS; bs++){
				for (time = 0; time < BURST; time++) {
#ifdef PROPAGATION
					fading_process(time, w, 0);	/* Most likely generates path coefficients. Time is the location on the long data stream. */

					w[0][0] = 1;
					w[0][1] = 0;
					Div = pow(2, 8);

					for (path_index = 0; path_index < PATH; path_index++) {	/* To be used in other parts of the code to drag channel coefficents */
						if (path_index > 0){ w[path_index][0] = w[path_index - 1][0] / Div; w[path_index][1] = w[path_index - 1][1] / Div; }
						wl[ts][bs][path_index][0] = w[path_index][0];	/* kk was preferred here */
						wl[ts][bs][path_index][1] = w[path_index][1];
					}
#ifdef MULTIPATH
					for (path_index = 0; path_index < PATH; path_index++) {	/* Signal goes thru channel */
						index = time - path_index * DELAY;

						if (index >= 0) {
							for (iq = 0; iq < 2; iq++) temp[iq] = transmitted_sigmat[ts][bs][index][iq];	/* Why use temp here? */

							rs_rcvd_sig[bs][time][0] += w[path_index][0] * temp[0] - w[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							rs_rcvd_sig[bs][time][1] += w[path_index][1] * temp[0] + w[path_index][0] * temp[1];
						}
					}
#else		  
					received_signal[ts][time][0] += PATH * w[0][0] * transmitted_sigmat[ts][bs][time][0] - w[0][1] * transmitted_sigmat[ts][bs][time][1];
					received_signal[ts][time][1] += PATH * w[0][1] * transmitted_sigmat[ts][bs][time][0] + w[0][0] * transmitted_sigmat[ts][bs][time][1];
#endif			
#else				   
					received_signal[ts][time][0] = transmitted_sigmat[ts][0][time][0];
					received_signal[ts][time][1] = transmitted_sigmat[ts][0][time][1];
#endif
					signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
				}
			}
			for (bs = 0; bs < BS; bs++){
				for (time = 0; time < BURST; time++) {
#ifdef PROPAGATION
					fading_process(time, w, 0);	/* Most likely generates path coefficients. Time is the location on the long data stream. */

					for (path_index = 0; path_index < PATH; path_index++) {	/* To be used in other parts of the code to drag channel coefficents */
						wl[ts][bs][path_index][0] = w[path_index][0];	/* kk was preferred here */
						wl[ts][bs][path_index][1] = w[path_index][1];
					}
#ifdef MULTIPATH
					for (path_index = 0; path_index < PATH; path_index++) {	/* Signal goes thru channel */
						index = time - path_index * DELAY;

						if (index >= 0) {
							for (iq = 0; iq < 2; iq++) temp[iq] = rs_rcvd_sig[bs][index][iq];	/* Why use temp here? */

							received_signal[ts][time][0] += w[path_index][0] * temp[0] - w[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							received_signal[ts][time][1] += w[path_index][1] * temp[0] + w[path_index][0] * temp[1];
						}
					}
#else		  
					received_signal[ts][time][0] += PATH * w[0][0] * transmitted_sigmat[ts][bs][time][0] - w[0][1] * transmitted_sigmat[ts][bs][time][1];
					received_signal[ts][time][1] += PATH * w[0][1] * transmitted_sigmat[ts][bs][time][0] + w[0][0] * transmitted_sigmat[ts][bs][time][1];
#endif			
#else				   
					received_signal[ts][time][0] = transmitted_sigmat[ts][0][time][0];
					received_signal[ts][time][1] = transmitted_sigmat[ts][0][time][1];
#endif
					signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
				}
			}
#else
			for (bs = 0; bs < BS; bs++){
				for (time = 0; time < BURST; time++) {
#ifdef PROPAGATION
					fading_process(time, w, 0);	/* Most likely generates path coefficients. Time is the location on the long data stream. */
					Relay_fading_process(time, w_r);
					//K = Relay_LOS_exponentialPATH(w_r);
					zerotap(w, w_extended, PATH, DELAY);
					zerotap(w_r, w_r_extended, RS_PATH, RS_DELAY);
					convolve(w_extended, PATH*DELAY, w_r_extended, RS_PATH*RS_DELAY, w_rb_extended);
					//convolve(w, PATH, w_r, RS_PATH, w_rb);

					/*printSignal("w", w, PATH);
					printSignal("w_r", w_r, RS_PATH);
					//printSignal("w_rb", w_rb, PATH + RS_PATH - 1);
					printSignal("w_ex", w_extended, PATH*DELAY);
					printSignal("w_r_ex", w_r_extended, RS_PATH*RS_DELAY);
					printSignal("w_rb_ex", w_rb_extended, PATH*DELAY + RS_PATH*RS_DELAY - 1);
					getchar();*/

					//w[0][0] = 1;
					//w[0][1] = 0;
					//Div = pow(2, 8);

					for (path_index = 0; path_index < PATH*DELAY + RS_PATH*RS_DELAY - 1; path_index++) {	/* To be used in other parts of the code to drag channel coefficents */
						//if (path_index > 0){ w[path_index][0] = w[path_index - 1][0] / Div; w[path_index][1] = w[path_index - 1][1] / Div; }
						wl[ts][bs][path_index][0] = w_rb_extended[path_index][0];	/* kk was preferred here */
						wl[ts][bs][path_index][1] = w_rb_extended[path_index][1];
					}
#ifdef MULTIPATH
					for (path_index = 0; path_index < PATH*DELAY + RS_PATH*RS_DELAY - 1; path_index++) {	/* Signal goes thru channel */
						//for (path_index = 0; path_index < PATH*DELAY; path_index++) {	/* Signal goes thru channel */
						//index = time - path_index * RS_DELAY;
						index = time - path_index;

						if (index >= 0) {
							for (iq = 0; iq < 2; iq++) temp[iq] = transmitted_sigmat[ts][bs][index][iq];	/* Why use temp here? */

							received_signal[ts][time][0] += w_rb_extended[path_index][0] * temp[0] - w_rb_extended[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							received_signal[ts][time][1] += w_rb_extended[path_index][1] * temp[0] + w_rb_extended[path_index][0] * temp[1];
							//received_signal[ts][time][0] += w_extended[path_index][0] * temp[0] - w_extended[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							//received_signal[ts][time][1] += w_extended[path_index][1] * temp[0] + w_extended[path_index][0] * temp[1];
						}
					}
					signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
					noise_generator(noise_amplitude, noise, 0);
					for (path_index = 0; path_index < PATH; path_index++) {	/* Signal goes thru channel */
						index = time - path_index * DELAY;

						if (index >= 0) {
							for (iq = 0; iq < 2; iq++) temp[iq] = noise[iq];	/* Why use temp here? */

							received_signal[ts][time][0] += w[path_index][0] * temp[0] - w[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							received_signal[ts][time][1] += w[path_index][1] * temp[0] + w[path_index][0] * temp[1];
							noise_power += sqr(w[path_index][0] * temp[0] - w[path_index][1] * temp[1]) + sqr(w[path_index][1] * temp[0] + w[path_index][0] * temp[1]);	/* why use w? this w is not combining vector, it is the channel */
						}
					}
#else		  
					received_signal[ts][time][0] += PATH * w[0][0] * transmitted_sigmat[ts][bs][time][0] - w[0][1] * transmitted_sigmat[ts][bs][time][1];
					received_signal[ts][time][1] += PATH * w[0][1] * transmitted_sigmat[ts][bs][time][0] + w[0][0] * transmitted_sigmat[ts][bs][time][1];
#endif			
#else				   
					received_signal[ts][time][0] = transmitted_sigmat[ts][0][time][0];
					received_signal[ts][time][1] = transmitted_sigmat[ts][0][time][1];
#endif
					//					signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
				}
			}
#endif // 2STG_RLY 
		}
	}
}

void fading_process(int time, double(*w)[2], double CNR){

	int	i, kk;
	double phase, x;
	static double
		ftemp,
		omega[PATH][RAY_COMPONENT],
		p0[PATH][RAY_COMPONENT],
		a[PATH][RAY_COMPONENT],
		power[PATH][RAY_COMPONENT],
		phase0;

	if (time == 0) {	/* Generate the channel at time instance 0, means constant over transmission duration. */
		ftemp = 2.0 * PI / SAMPLING_RATE;	/* 2xpixf = 2xpix1/T */

		for (i = 0; i < PATH; i++){
			for (kk = 0; kk < RAY_COMPONENT; kk++) {

				omega[i][kk] = ftemp * FD * cos(2.0 * PI * (double)rand() / RAND_MAX);
				p0[i][kk] = 2.0 * PI * (double)rand() / RAND_MAX;

				/* Channel profile (recall Proakis) */
				if (EXPMODEL == ON)
					//exponential_PATH(i, kk, power);
					exponential_PATH_rlystyle(i, kk, power);
				//Relay_exponential_PATH(i, kk, power);
				else
					samelevel_PATH(i, kk, power);

				x = (double)rand() / RAND_MAX;	if (x < 1.0e-6) x = 1.0e-6;
				a[i][kk] = sqrt(power[i][kk] * -log(x) * (pow(10.0, CNR / 10.0) / CODING_RATE));
			}
		}
		phase0 = 2.0 * PI * (double)rand() / RAND_MAX; /* Bu da ne? */
	}

	for (i = 0; i < PATH; i++) {
		w[i][0] = 0.0;
		w[i][1] = 0.0;

		for (kk = 0; kk < RAY_COMPONENT; kk++) {
			phase = (double)omega[i][kk] * time + p0[i][kk];
			w[i][0] += a[i][kk] * cos(phase);
			w[i][1] += a[i][kk] * sin(phase);
		}
	}
}


void samelevel_PATH(int i, int kk, double(*power)[RAY_COMPONENT]){
	power[i][kk] = (double) 1.0 / (RAY_COMPONENT * PATH);
}

void exponential_PATH(int i, int kk, double(*power)[RAY_COMPONENT]){
	int d;
	double path[PATH], sum;

	sum = 0.0;
	for (d = 0; d < PATH; d++){
		path[d] = exp((-DUR * log(10.0) * ((double)d)) / (10.0 * (double)(PATH - 1)));
		sum += path[d];
	}
	power[i][kk] = (double)path[i] / (((double)RAY_COMPONENT) * sum);
}

