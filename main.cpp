#include "const.h"
double CNR, EbNo;
int loop;
int maxloop;
/* TS : time slot index, BS is for BS when TS = 0, and it is for relay when TS = 1. */
double wh_sequence[BS][Np][N];
double transmitted_sigmat[TS][BS][BURST][2];	/* This combines all transmitted singal into matrix form. */
double mse_sum[LOOPN][NUMBER_OF_STEPS];	/* To print out average MSE */
extern errno_t err;	/* Error vector */

void main(){
	double transmitted_signal1[BURST][2];	/* Signal of BS1 */
	double transmitted_signal2[BURST][2];	/* Signal of BS2 */

	double relay1_tx_signal[BURST][2];	/* Signal relayed from RL1 */
	double relay2_tx_signal[BURST][2];	/* Signal relayed from RL2 */

	double transmitted_sigvec[TS][BURST][2];	/* Received signal at MS1 after multipath propogation. */

	int transmitted_bit1[Nd][N_dbps];
	int transmitted_bit2[Nd][N_dbps];
	int received_bit[Nd][N_dbps];

	int bit[Nd][N_dbps];
	int bit2[Nd][N_dbps];


	/*~~~~~~~~~Initialize~begin~~~~~~~~~~~~~~~~~~~~~~~~*/
	setting_print();	/* Lots of printf */
	srand((unsigned)time(NULL));
	FFT_initialize();	/* Construct an FFT matrix */
	state_model_convolutional();	/* Function and parameters of convolutional encoder */

	transmitter_initialize(transmitted_sigmat);	/* Set transmitted signal vector to zero */
#ifdef WH_PREAMBLE	/* for 2 user */
	WHadamard_initialize(transmitted_sigmat, wh_sequence);	/* For walsh hadamart sequence */
#else
	Zadoff_Chu_layer_initialize();	/* Most likely: generate two ZC sequences. */
	Zadoff_Chu_generator_layer1(transmitted_sigmat);
	Zadoff_Chu_generator_layer2(transmitted_sigmat);	/* Only difference is the transmitter(? isnt that for time slot) index. One is for tx = 0 one is for tx = 1. */
#endif // WH
	/*~~~~~~~~initialize end~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	for (EbNo = EbNo_FROM; EbNo >= EbNo_TO; EbNo -= 5.0){	/* begin :: EbNo */
		CNR = (double)EbNo + 3.0 * MODULATION / 2.0;
#ifdef ADPTV_LOOP
		maxloop = intPow((int)EbNo, 2) + 10;
#else
		maxloop = LOOPN;
#endif // ADPTV_LOOP
		for (loop = 0; loop < maxloop; loop++){	/* begin :: loop */
			transmitter1(transmitted_bit1, transmitted_signal1);	/* This transmitter1 and 2 are the same functions. */
			transmitter2(transmitted_bit2, transmitted_signal2);	/* Generate, and modulate the data into 10 OFDM symbols  */
			relay1(transmitted_signal1, relay1_tx_signal, bit);	/* I see no difference */
			relay2(transmitted_signal2, relay2_tx_signal, bit2);
			transmitted_signal_Matrix(transmitted_signal1, transmitted_signal2,
				relay1_tx_signal, relay2_tx_signal, transmitted_sigmat);	/* Combine those for into the fifth */
			multipath_propagation(transmitted_sigmat, transmitted_sigvec);	/* Signal propogates thru channel. We also need the signal which gone thru an ideal channel*/
#ifdef NOISE
			noise_addition((0), transmitted_sigvec[0], transmitted_sigvec[0]);	/* Add noise to first? time slot */
			noise_addition(0, transmitted_sigvec[1], transmitted_sigvec[1]);	/* Add noise to second? time slot */
#endif // NOISE
			receiver(transmitted_sigvec, received_bit);
			BER(loop, transmitted_bit1, received_bit);
		}	/* end :: loop */	  
	}	/* end :: EbNo */
}	/* end :: main() */