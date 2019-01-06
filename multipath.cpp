#include "const-function-declerations.h"		
extern double CNR;
extern double signal_power, noise_power;
extern double path_ms[PATH][RAY_COMPONENT], sum_ms;

double wl[TS][BS][PATH*DELAY + RS_PATH*RS_DELAY - 1][2];	/* Used for pulling out channel coefficients. */

/* received signal = transmitted_sigvec */
void multipath_propagation(double(*transmitted_sigmat)[BS][BURST][2], double(*received_signal)[BURST][2]){
	int time, path_index, index, iq, ts, bs, i;
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
		if ((ts == 0) || (RS_DF == ON)){
			for (bs = 0; bs < BS; bs++){
				for (time = 0; time < BURST; time++) {
#ifdef PROPAGATION
					fading_process(time, w, (ts == 0) ? (0 - LOSS) : 0);

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
					for (path_index = 0; path_index < PATH; path_index++) {	/* Signal goes thru channel */
						index = time - path_index * DELAY;

						if (index >= 0) {
							for (iq = 0; iq < 2; iq++) temp[iq] = transmitted_sigmat[ts][bs][index][iq];	/* Why use temp here? */
							received_signal[ts][time][0] += w[path_index][0] * temp[0] - w[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							received_signal[ts][time][1] += w[path_index][1] * temp[0] + w[path_index][0] * temp[1];
						}
					}
#else				   
					received_signal[ts][time][0] = transmitted_sigmat[ts][0][time][0];
					received_signal[ts][time][1] = transmitted_sigmat[ts][0][time][1];
#endif
				}

				for (time = 0; time < BURST; time++) {
					signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
				}
			}
		}
		else if ((ts == 1) && (RS_DF == OFF)) {
			for (bs = 0; bs < BS; bs++){
				for (time = 0; time < BURST; time++) {
#ifdef PROPAGATION
					fading_process(time, w, 0);	/* Most likely generates path coefficients. Time is the location on the long data stream. */
					Relay_fading_process(time, w_r);	 /* cnr is 0 inside. */
					//K = Relay_LOS_exponentialPATH(w_r);
					zerotap(w, w_extended, PATH, DELAY);
					zerotap(w_r, w_r_extended, RS_PATH, RS_DELAY);
					//convolve(w, PATH, w_r, PATH, w_rb);
					convolve(w_extended, PATH*DELAY, w_r_extended, RS_PATH*RS_DELAY, w_rb_extended);
					for (path_index = 0; path_index < PATH*DELAY + RS_PATH*RS_DELAY - 1; path_index++) {	/* To be used in other parts of the code to drag channel coefficents */
						wl[ts][bs][path_index][0] = w_rb_extended[path_index][0];	/* kk was preferred here */
						wl[ts][bs][path_index][1] = w_rb_extended[path_index][1];
					}
					for (path_index = 0; path_index < PATH*DELAY + RS_PATH*RS_DELAY - 1; path_index++) {	/* Signal goes thru channel */
						index = time - path_index; /* delay is taken as 1 */

						if (index >= 0) {
							for (iq = 0; iq < 2; iq++) temp[iq] = transmitted_sigmat[ts][bs][index][iq];	/* Why use temp here? */
							received_signal[ts][time][0] += w_rb_extended[path_index][0] * temp[0] - w_rb_extended[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							received_signal[ts][time][1] += w_rb_extended[path_index][1] * temp[0] + w_rb_extended[path_index][0] * temp[1];
							//received_signal[ts][time][0] += w_extended[path_index][0] * temp[0] - w_extended[path_index][1] * temp[1];	/* why use w? this w is not combining vector, it is the channel */
							//received_signal[ts][time][1] += w_extended[path_index][1] * temp[0] + w_extended[path_index][0] * temp[1];
						}
					}
					noise_generator(noise_amplitude, noise, 0);
					for (path_index = 0; path_index < PATH; path_index++) {
						index = time - path_index * DELAY;

						if (index >= 0) {
							for (iq = 0; iq < 2; iq++) temp[iq] = noise[iq];

							received_signal[ts][time][0] += w[path_index][0] * temp[0] - w[path_index][1] * temp[1];
							received_signal[ts][time][1] += w[path_index][1] * temp[0] + w[path_index][0] * temp[1];
							//noise_power += sqr(w[path_index][0] * temp[0] - w[path_index][1] * temp[1]) + sqr(w[path_index][1] * temp[0] + w[path_index][0] * temp[1]);	
							noise_power += sqr(w[path_index][0] * temp[0] - w[path_index][1] * temp[1]) + sqr(w[path_index][1] * temp[0] + w[path_index][0] * temp[1]);
							signal_power -= sqr(w[path_index][0] * temp[0] - w[path_index][1] * temp[1]) + sqr(w[path_index][1] * temp[0] + w[path_index][0] * temp[1]);
						}
					}
#else				   
					received_signal[ts][time][0] = transmitted_sigmat[ts][0][time][0];
					received_signal[ts][time][1] = transmitted_sigmat[ts][0][time][1];
#endif
					//					signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
				}

				for (time = 0; time < BURST; time++) {
					signal_power += sqr(received_signal[ts][time][0]) + sqr(received_signal[ts][time][1]);
				}
			}
		}
		else { cout << "Something wrong at multipath.cpp." << endl;	getchar();	exit(3); }
	}
}

void fading_process(int time, double(*w)[2], double CNR){
	int	path_index, component_index;
	double ftemp, phase, x, power[PATH][RAY_COMPONENT];
	static double
		omega[PATH][RAY_COMPONENT],
		p0[PATH][RAY_COMPONENT],
		a[PATH][RAY_COMPONENT];

	if (time == 0) {	/* Generate the channel at time instance 0, means constant over transmission duration. */
		ftemp = 2.0 * PI / SAMPLING_RATE;	/* 2xpixf = 2xpix1/T */
		for (path_index = 0; path_index < PATH; path_index++){
			for (component_index = 0; component_index < RAY_COMPONENT; component_index++) {
				omega[path_index][component_index] = ftemp * FD * cos(2.0 * PI * (double)rand() / RAND_MAX);
				p0[path_index][component_index] = 2.0 * PI * (double)rand() / RAND_MAX;
				power[path_index][component_index] = path_ms[path_index][component_index];
				x = (double)rand() / RAND_MAX;
				if (x < 1.0e-6) x = 1.0e-6;
				a[path_index][component_index] = sqrt(power[path_index][component_index] * (((path_index == 0) && (LOS == ON)) ? 1 : (-log(x))) * (pow(10.0, CNR / 10.0) / CODING_RATE));
			}
		}
	}
	for (path_index = 0; path_index < PATH; path_index++) {
		w[path_index][0] = 0.0;
		w[path_index][1] = 0.0;

		for (component_index = 0; component_index < ((path_index == 0 && LOS == ON) ? 1 : RAY_COMPONENT); component_index++){
			phase = (double)omega[path_index][component_index] * time + p0[path_index][component_index];
			w[path_index][0] += a[path_index][component_index] * cos(phase);
			w[path_index][1] += a[path_index][component_index] * sin(phase);
		}
	}
}

