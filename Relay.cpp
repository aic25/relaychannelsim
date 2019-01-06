#include "const.h"
extern double CNR, BS_signal_power[BS];
double h_r[BS][PATH][2];
void  relay(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps],int bs){

#if SWITCH == 0
	int k;
	for (k = 0; k < BURST; k++){
		relay_out[k][0] = relay_in[k][0];
		relay_out[k][1] = relay_in[k][1];
	}

#elif SWITCH == 1
#ifdef DF
	Relay_DF(relay_in, relay_out, bit, bs);
#else
	Relay_AF(relay_in, relay_out, bit, bs);
#endif // DF	 
#endif
}

void Relay_DF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs){
	double relay_received_signal[BURST][2];

#ifdef RLY_PROPOGATION
	Relay_propagation(relay_in, relay_received_signal, bs);	/* Signal goes thru channel and reaches the relay element */
#else
	for (size_t k = 0; k < BURST; k++){
		relay_received_signal[k][0] = relay_in[k][0];
		relay_received_signal[k][1] = relay_in[k][1];
	}
#endif // RLY_PROPOGATION

#ifdef RLY_NOISE
#ifdef RLY_PROPOGATION
	noise_addition(0, relay_received_signal, relay_received_signal);	/* given the CNR add noise */
#else
	noise_addition(RLY_CNR, relay_received_signal, relay_received_signal);
#endif // RLY_PROPOGATION
#endif // RLY_NOISE

	Relay_receiver(relay_received_signal, bit, bs);	/* HD decoder. */
	Relay_transmitter(bit, relay_out);	/* Encode, modulate etc. */
}

void Relay_AF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs){
	double RS_signal_power[BS], beta;

	RS_signal_power[bs] = 0;

#ifdef RLY_PROPOGATION
	Relay_propagation(relay_in, relay_out, bs);	/* Signal goes thru channel and reaches the relay element */
#else
	for (size_t k = 0; k < BURST; k++){
		relay_out[k][0] = relay_in[k][0];
		relay_out[k][1] = relay_in[k][1];
	}
#endif // RLY_PROPOGATION

#ifdef RLY_NOISE
#ifdef RLY_PROPOGATION
	noise_addition(0, relay_out, relay_out);	/* given the CNR add noise */
#else
	noise_addition(RLY_CNR, relay_out, relay_out);
#endif // RLY_PROPOGATION
#endif // RLY_NOISE				

	/* Beta is calculated */
	for (size_t k = 0; k < BURST; k++){
		RS_signal_power[bs] += sqr(relay_out[k][0]) + sqr(relay_out[k][1]);
	}
	RS_signal_power[bs] /= BURST;

	beta = sqrt(BS_signal_power[bs] / RS_signal_power[bs]);

	for (size_t k = 0; k < BURST; k++){
		relay_out[k][0] = beta*relay_out[k][0];
		relay_out[k][1] = beta*relay_out[k][1];
	}
}
/************************* Relay Propagation and Noise addition *******************************/

void Relay_propagation(double(*transmitted_signal)[2], double(*received_signal)[2], int bs){
	int time, path_index, index, i, j;
	double K;

	// cout << "tx sig at relay_prop: " << transmitted_signal[0][0] << " " << transmitted_signal[0][1] << endl;

	for (time = 0; time < BURST; time++) {
		Relay_fading_process(time, h_r[bs]);	/* Channel generator */

		received_signal[time][0] = 0.0;
		received_signal[time][1] = 0.0;

		for (path_index = 0; path_index < RS_PATH; path_index++) {
			index = time - path_index * RS_DELAY;
			if (0 <= index && index < BURST) {
				received_signal[time][0] += h_r[bs][path_index][0] * transmitted_signal[index][0] - h_r[bs][path_index][1] * transmitted_signal[index][1];
				received_signal[time][1] += h_r[bs][path_index][1] * transmitted_signal[index][0] + h_r[bs][path_index][0] * transmitted_signal[index][1];
			}
		}
	}
}

void Relay_fading_process(int time, double(*w)[2]){
	int	path_index, kk;
	double ftemp, phase, x, power[RS_PATH][RS_RAY_COMPONENT];
	static double omega[RS_PATH][RS_RAY_COMPONENT], p0[RS_PATH][RS_RAY_COMPONENT], a[RS_PATH][RS_RAY_COMPONENT];

	if (time == 0) {
		ftemp = 2.0 * PI / SAMPLING_RATE;
		for (path_index = 0; path_index < RS_PATH; path_index++){
			for (kk = 0; kk < RS_RAY_COMPONENT; kk++) {
				omega[path_index][kk] = ftemp * FDr * cos(2.0 * PI * (double)rand() / RAND_MAX);
				p0[path_index][kk] = 2.0 * PI * (double)rand() / RAND_MAX;
				if (RLY_EXPMODEL == ON) {
					Relay_exponential_PATH(path_index, kk, power);
				}
				else{
					Relay_samelevel_PATH(path_index, kk, power);
				}
				x = (double)rand() / RAND_MAX;
				if (x < 1.0e-6) x = 1.0e-6;
				a[path_index][kk] = sqrt(power[path_index][kk] * (((path_index == 0) && (RLY_LOS == ON)) ? 1 : (-log(x))) * (pow(10.0, 0 / 10.0) / CODING_RATE));
			}
		}
	}
	for (path_index = 0; path_index < RS_PATH; path_index++) {
		w[path_index][0] = 0.0;
		w[path_index][1] = 0.0;

		for (kk = 0; kk < RS_RAY_COMPONENT; kk++) {
			phase = omega[path_index][kk] * (double)time + p0[path_index][kk];

			w[path_index][0] += a[path_index][kk] * cos(phase);
			w[path_index][1] += a[path_index][kk] * sin(phase);
		}
	}
}

void Relay_samelevel_PATH(int i, int kk, double(*power)[RS_RAY_COMPONENT]){
	power[i][kk] = (double) 1.0 / (RS_RAY_COMPONENT * RS_PATH);
}

void Relay_exponential_PATH(int i, int kk, double(*power)[RS_RAY_COMPONENT]){
	int path_index;
	double path[RS_PATH], sum;

	sum = 0.0;
	for (path_index = 0; path_index < RS_PATH; path_index++){
		path[path_index] = exp((-((double)path_index)) / ((double)(RS_PATH)*RS_TAU));
		sum += path[path_index];
	}
	power[i][kk] = (double)path[i] / (((double)RS_RAY_COMPONENT) * sum);
}
void exponential_PATH_rlystyle(int i, int kk, double(*power)[RAY_COMPONENT]){
	int path_index;
	double path[PATH], sum;

	sum = 0.0;
	for (path_index = 0; path_index < PATH; path_index++){
		path[path_index] = exp((-((double)path_index)) / ((double)(PATH)*TAU));
		sum += path[path_index];
	}
	power[i][kk] = (double)path[i] / (((double)RAY_COMPONENT) * sum);
}




