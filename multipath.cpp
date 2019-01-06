#include "const.h"			
extern double CNR;
extern double signal_power;

double wl[TS][BS][PATH][2];	/* Used for pulling out channel coefficients. */

/* received signal = transmitted_sigvec */
void multipath_propagation(double(*transmitted_sigmat)[BS][BURST][2], double(*received_signal)[BURST][2]){
	int time, path_index, index, iq, ts, bs;
	double temp[2];
	double w[PATH][2];	/* w? */

	for (ts = 0; ts < TS; ts++){
		for (time = 0; time < BURST; time++) {
			for (iq = 0; iq < 2; iq++) received_signal[ts][time][iq] = 0.0;
		}
	}

	for (ts = 0; ts < TS; ts++){
		for (bs = 0; bs < BS; bs++){
			for (time = 0; time < BURST; time++) {
				fading_process(time, w, (ts == 1) ? CNR : (CNR - LOSS));	/* Most likely generates path coefficients. Time is the location on the long data stream. */

				for (path_index = 0; path_index < PATH; path_index++) {	/* To be used in other parts of the code to drag channel coefficents */
					wl[ts][bs][path_index][0] = w[path_index][0];	/* kk was preferred here */
					wl[ts][bs][path_index][1] = w[path_index][1];
				}				 

				for (path_index = 0; path_index < PATH; path_index++) {	/* Signal goes thru channel */
					index = time - path_index * DELAY;

					if (index >= 0) {
						for (iq = 0; iq < 2; iq++) temp[iq] = transmitted_sigmat[ts][bs][index][iq];	/* Why use temp here? */

						received_signal[ts][time][0] += w[path_index][0] * temp[0] - w[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
						received_signal[ts][time][1] += w[path_index][1] * temp[0] + w[path_index][0] * temp[1];
					}
				}
				signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
			}
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
				if (EXPMODEL == 1)
					exponential_PATH(i, kk, power);
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

