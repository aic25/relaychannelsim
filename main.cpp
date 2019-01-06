#include "const.h"
double CNR, EbNo;
int loop;
int maxloop;
/* TS : time slot index. BS is for BS when TS = 0, and it is for relay when TS = 1. */
double wh_sequence[BS][Np][N][2];
double transmitted_sigmat[TS][BS][BURST][2];	/* This combines all transmitted singal into matrix form. */
double mse_sum[LOOPN][NUMBER_OF_STEPS];	/* To print out average MSE */
extern errno_t err;	/* Error vector */

void main(){
	double transmitted_signal[BS][BURST][2];	/* Signal of BS1 */
	double relay_tx_signal[BS][BURST][2];	/* Signal relayed from RL1 */
	double transmitted_sigvec[TS][BURST][2];	/* Received signal at MS1 after multipath propogation. */
	int transmitted_bit[BS][Nd][N_dbps];
	int received_bit[Nd][N_dbps];
	int bs, ts, rs, rs_bit[BS][Nd][N_dbps];


	/*~~~~~~~~~Initialize~begin~~~~~~~~~~~~~~~~~~~~~~~~*/
	setting_print();	/* Lots of printf */
	srand((unsigned)time(NULL));
	FFT_initialize();	/* Construct an FFT matrix */
	state_model_convolutional();	/* Function and parameters of convolutional encoder */

	transmitter_initialize(transmitted_sigmat);	/* Set transmitted signal vector to zero */
	transmitters_initialize(transmitted_signal);
#ifdef WH_PREAMBLE	/* for 2 user */
	WHadamard_initialize(transmitted_sigmat, wh_sequence, transmitted_signal);	/* For walsh hadamart sequence */
#else
	Zadoff_Chu_layer_initialize();	/* Most likely: generate two ZC sequences. */
	Zadoff_Chu_generator_layer1(transmitted_sigmat); /* These two is for TS. */
	Zadoff_Chu_generator_layer2(transmitted_sigmat);
#endif // WH
	/*~~~~~~~~initialize end~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	for (EbNo = EbNo_FROM; EbNo >= EbNo_TO; EbNo -= 5.0){	/* begin :: EbNo */
		CNR = (double)EbNo + 3.0 * MODULATION / 2.0;
#ifdef ADPTV_LOOP
		maxloop = intPow((int)EbNo, 2) + 10;
		//maxloop = 5;
#else
		maxloop = LOOPN;
#endif // ADPTV_LOOP
		for (loop = 0; loop < maxloop; loop++){	/* begin :: loop */
			for (rs = 0; rs < BS; rs++)
			{
				bs = rs;
				// cout << "tx sig at main: " << transmitted_signal[rs][0][0] << " " << transmitted_signal[rs][0][1] << endl;
				// cout << "tx sig at main: " << transmitted_signal[bs][0][0] << " " << transmitted_signal[bs][0][1] << endl;
				transmitter(transmitted_bit[bs], transmitted_signal[bs], bs);
				// cout << "tx sig at main: " << transmitted_signal[rs][0][0] << " " << transmitted_signal[rs][0][1] << endl;
				relay(transmitted_signal[rs], relay_tx_signal[rs], rs_bit[rs], rs);	/* I see no difference */
			}
			transmitted_signal_Matrix(transmitted_signal, relay_tx_signal, transmitted_sigmat);	/* Combine those for into the fifth */
			multipath_propagation(transmitted_sigmat, transmitted_sigvec);	/* Signal propogates thru channel. We also need the signal which gone thru an ideal channel*/
#ifdef NOISE
			for (ts = 0; ts < TS; ts++)
			{
#ifdef PROPAGATION
				noise_addition(CNR, transmitted_sigvec[ts], transmitted_sigvec[ts]);	/* Add noise to first? time slot. CNR=0 bcs it is included in propogation */
#endif
			}
#endif // NOISE
			receiver(transmitted_sigvec, received_bit);
			BER_RsPlusBs(loop, rs_bit[0], transmitted_bit[0], received_bit);
		}	/* end :: loop */
	}	/* end :: EbNo */
}	/* end :: main() */