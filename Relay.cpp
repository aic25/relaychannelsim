#include "const-function-declerations.h"		
extern double CNR, BS_signal_power[BS];
double h_r[BS][RS_PATH][2];
extern double path_rs[RS_PATH][RS_RAY_COMPONENT], sum_rs;
extern double CNR, EbNo;
extern int loop;

void  relay(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs){

#if((SWITCH == 0)||(RS_DF==OFF))
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
	double K, pow_xmit, pow_rcvd;
	pow_xmit = 0;
	pow_rcvd = 0;

	// cout << "tx sig at relay_prop: " << transmitted_signal[0][0] << " " << transmitted_signal[0][1] << endl;

	for (time = 0; time < BURST; time++) {
		Relay_fading_process(time, h_r[bs]);	/* Channel generator */
		//VectorNormalize(h_r[bs], RS_PATH);

		//if (time == 0){cout << "h: " << bs << endl; for (i = 0; i < RS_PATH; i++){ cout << h_r[bs][i][0] << " " << h_r[bs][i][1] << endl; }	}
		//if (time == 0){ cout << "h: " << bs << endl; for (i = 0; i < RS_PATH; i++){ cout << sqrt(sqr(h_r[bs][i][0]) + sqr(h_r[bs][i][1])) << endl; } }

		received_signal[time][0] = 0.0;
		received_signal[time][1] = 0.0;

		for (path_index = 0; path_index < RS_PATH; path_index++) {
			index = time - path_index * RS_DELAY;
			if (0 <= index && index < BURST) {
				received_signal[time][0] += h_r[bs][path_index][0] * transmitted_signal[index][0] - h_r[bs][path_index][1] * transmitted_signal[index][1];
				received_signal[time][1] += h_r[bs][path_index][1] * transmitted_signal[index][0] + h_r[bs][path_index][0] * transmitted_signal[index][1];
			}
		}

		/*if (time>Np*SAMPLEN){
			pow_xmit += sqrt(sqr(transmitted_signal[time][0]) + sqr(transmitted_signal[time][1])) / (double)(Nd*SAMPLEN);
			pow_rcvd += sqrt(sqr(received_signal[time][0]) + sqr(received_signal[time][1])) / (double)(Nd*SAMPLEN);
		} */

	}

	/*-----------debug---------------*/
	/*
	if (bs == 0){
		cout << "Bs1:  Trans. Sig Pow.: " << pow_xmit << " Rcvd. Sig. Pow.: " << pow_rcvd << endl;
	}
	else{
		cout << "Bs2:  Trans. Sig Pow.: " << pow_xmit << " Rcvd. Sig. Pow.: " << pow_rcvd << endl;
	}
	
	FILE *xmit_sign, *rcvd_sign;
	errno_t err;
	if (bs == 0){
		cout << "Bs1:  Trans. Sig Pow.: " << pow_xmit << " Rcvd. Sig. Pow.: " << pow_rcvd << endl;
		char pout_fname[50] = "xmit_signal_bs1_";
		char pout_fname_2[50] = "rcvd_signal_bs1_";
		char extention[5] = ".dat";
		char addition[10];
		_itoa(loop, addition, 10);
		std::strcat(pout_fname, addition);
		std::strcat(pout_fname_2, addition);
		std::strcat(pout_fname, extention);
		std::strcat(pout_fname_2, extention);
		if ((err = fopen_s(&xmit_sign, pout_fname, "w+")) != 0){ cout << "Error opening xmit_signal.dat!" << endl; getchar(); exit(1); }
		if ((err = fopen_s(&rcvd_sign, pout_fname_2, "w+")) != 0){ cout << "Error opening rcvd_signal.dat!" << endl; getchar(); exit(1); }
		for (time = Np*SAMPLEN; time < BURST; time++){
			fprintf(xmit_sign, "%f\t%f\n", transmitted_signal[time][0], transmitted_signal[time][1]);
			fprintf(rcvd_sign, "%f\t%f\n", received_signal[time][0], received_signal[time][1]);
		}
		fclose(xmit_sign);
		fclose(rcvd_sign);
	}
	else{
		cout << "Bs2:  Trans. Sig Pow.: " << pow_xmit  << " Rcvd. Sig. Pow.: " << pow_rcvd  << endl;
		char pout_fname[50] = "xmit_signal_bs2_";
		char pout_fname_2[50] = "rcvd_signal_bs2_";
		char extention[5] = ".dat";
		char addition[10];
		_itoa(loop, addition, 10);
		std::strcat(pout_fname, addition);
		std::strcat(pout_fname_2, addition);
		std::strcat(pout_fname, extention);
		std::strcat(pout_fname_2, extention);
		if ((err = fopen_s(&xmit_sign, pout_fname, "w+")) != 0){ cout << "Error opening xmit_signal.dat!" << endl; getchar(); exit(1); }
		if ((err = fopen_s(&rcvd_sign, pout_fname_2, "w+")) != 0){ cout << "Error opening rcvd_signal.dat!" << endl; getchar(); exit(1); }
		for (time = Np*SAMPLEN; time < BURST; time++){
			fprintf(xmit_sign, "%f\t%f\n", transmitted_signal[time][0], transmitted_signal[time][1]);
			fprintf(rcvd_sign, "%f\t%f\n", received_signal[time][0], received_signal[time][1]);
		}
		fclose(xmit_sign);
		fclose(rcvd_sign);
	}
	*/
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
				power[path_index][component_index] = path_rs[path_index][component_index];
				//cout << path_rs[path_index] << endl; getchar();
				x = (double)rand() / RAND_MAX; if (x < 1.0e-6) x = 1.0e-6;
				//a[path_index][component_index] = sqrt(power[path_index][component_index] * (((path_index == 0) && (RLY_LOS == ON)) ? 1 : (-log(x))) * (1 / (double)CODING_RATE));
				a[path_index][component_index] = sqrt(power[path_index][component_index] * (((path_index == 0) && (component_index==0) && (RLY_LOS == ON)) ? 1 : (-log(x))) * (pow(10.0, 0 / 10.0) / CODING_RATE));
				//a[path_index][component_index] = sqrt(power[path_index][component_index] * (-log(x)) * (pow(10.0, 0 / 10.0) / CODING_RATE));		 
				//if (component_index == 0){ if (path_index == 0)cout << "a: " << endl; cout << a[path_index][0] << endl; }
			}
		}
	}
	for (path_index = 0; path_index < RS_PATH; path_index++) {
		w[path_index][0] = 0.0;
		w[path_index][1] = 0.0;

		for (component_index = 0; component_index < (path_index == 0 ? (RLY_LOS == ON ? 1 : RS_RAY_COMPONENT) : RS_RAY_COMPONENT); component_index++){
			phase = omega[path_index][component_index] * (double)time + p0[path_index][component_index];
			w[path_index][0] += a[path_index][component_index] * cos(phase);
			w[path_index][1] += a[path_index][component_index] * sin(phase);
		}			 
		//VectorNormalize(w, RS_PATH);
		//if (time == 0){	if (path_index == 0)cout << "w: " << endl;	cout << w[path_index][0] << " " << w[path_index][1] << endl;}
	}
}
