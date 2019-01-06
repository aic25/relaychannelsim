#include "const-function-declerations.h"		
extern double CNR, BS_signal_power[BS];
double h_r[BS][RS_PATH][2];
void  relay(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs){

#if ((SWITCH == 0)||(RS_DF==OFF))
	int k;
	for (k = 0; k < BURST; k++){
		relay_out[k][0] = relay_in[k][0];
		relay_out[k][1] = relay_in[k][1];
	}	  
#elif SWITCH == 1
	Relay_DF(relay_in, relay_out, bit, bs);
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
	noise_addition(RLY_CNR, relay_received_signal, relay_received_signal, 0);
#endif // RLY_NOISE

	Relay_receiver(relay_received_signal, bit, bs);	/* HD decoder. */
	Relay_transmitter(bit, relay_out);	/* Encode, modulate etc. */
}

/************************* Relay Propagation and Noise addition *******************************/

void Relay_propagation(double(*transmitted_signal)[2], double(*received_signal)[2], int bs){
	int time, path_index, index, i, j;
	double K;
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
	int	path_index, component_index;
	double ftemp, phase, x, power[RS_PATH][RS_RAY_COMPONENT];
	static double
		omega[RS_PATH][RS_RAY_COMPONENT],
		p0[RS_PATH][RS_RAY_COMPONENT],
		a[RS_PATH][RS_RAY_COMPONENT];

	if (time == 0) {
		ftemp = 2.0 * PI / SAMPLING_RATE;
		for (path_index = 0; path_index < RS_PATH; path_index++){
			for (component_index = 0; component_index < RS_RAY_COMPONENT; component_index++) {
				omega[path_index][component_index] = ftemp * FDr * cos(2.0 * PI * (double)rand() / RAND_MAX);
				p0[path_index][component_index] = 2.0 * PI * (double)rand() / RAND_MAX;

				if (RLY_EXPMODEL == ON) { Relay_exponential_PATH(path_index, component_index, power); }
				else{ Relay_samelevel_PATH(path_index, component_index, power); }

				/*Relay_exponential_PATH(path_index, component_index, power);
				cout << "Exp path: " << endl << power[path_index][component_index] << endl;
				Relay_samelevel_PATH(path_index, component_index, power);
				cout << "Same lvl path: " << endl << power[path_index][component_index] << endl;
				getchar();			   */

				x = (double)rand() / RAND_MAX;
				if (x < 1.0e-6) x = 1.0e-6;
				a[path_index][component_index] = sqrt(power[path_index][component_index] * (((path_index == 0) && (RLY_LOS == ON)) ? 1 : (-log(x))) * (pow(10.0, 0 / 10.0) / CODING_RATE));
				//a[path_index][component_index] = sqrt(power[path_index][component_index] * (-log(x)) * (pow(10.0, 0 / 10.0) / CODING_RATE));
			}
		}
	}
	for (path_index = 0; path_index < RS_PATH; path_index++) {
		w[path_index][0] = 0.0;
		w[path_index][1] = 0.0;

		for (component_index = 0; component_index < RS_RAY_COMPONENT; component_index++) {
			phase = omega[path_index][component_index] * (double)time + p0[path_index][component_index];
			w[path_index][0] += a[path_index][component_index] * cos(phase);
			w[path_index][1] += a[path_index][component_index] * sin(phase);
		}
	}
}

void Relay_samelevel_PATH(int i, int kk, double(*power)[RS_RAY_COMPONENT]){
	power[i][kk] = (double) 1.0 / (RS_RAY_COMPONENT * RS_PATH);
}

void Relay_exponential_PATH(int i, int kk, double(*power)[RS_RAY_COMPONENT]){
	int path_index;
	double temp;
	double path[RS_PATH], sum;

	sum = 0.0;
	for (path_index = 0; path_index < RS_PATH; path_index++){
		temp = (double)path_index;
		path[path_index] = (exp(-temp / (double)RS_TAU) > pow(10, -6)) ? exp(-temp / (double)RS_TAU) : 0;
		sum += path[path_index];
	}
	power[i][kk] = (double)path[i] / (((double)RS_RAY_COMPONENT) * sum);
}