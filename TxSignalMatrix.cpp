#include "const-function-declerations.h"			  
void transmitted_signal_Matrix(
	double(*transmitted_signal)[BURST][2],
	double(*tx_relay_signal)[BURST][2],
	double(*transmitted_sigmat)[BS][BURST][2]
	){

	int i, ts, tx;
	for (tx = 0; tx < BS; tx++){
		for (i = (Np * SAMPLEN); i < BURST; i++){	/* Apply for data, not for parity symbols. */
			if (tx == 0){
				transmitted_sigmat[0][tx][i][0] = transmitted_signal[tx][i][0];
				transmitted_sigmat[0][tx][i][1] = transmitted_signal[tx][i][1];
			}
			if (tx > 0){
#if (INTERFERENCE==OFF)
				transmitted_sigmat[0][tx][i][0]=  0.0;
				transmitted_sigmat[0][tx][i][1]=  0.0;
#else
				transmitted_sigmat[0][tx][i][0] = transmitted_signal[tx][i][0];
				transmitted_sigmat[0][tx][i][1] = transmitted_signal[tx][i][1];
#endif
			}
			if (tx == 0){
				transmitted_sigmat[1][tx][i][0] = tx_relay_signal[tx][i][0];
				transmitted_sigmat[1][tx][i][1] = tx_relay_signal[tx][i][1];
			}
			if (tx > 0){
#if (INTERFERENCE==OFF)
				transmitted_sigmat[1][tx][i][0]= 0.0; 
				transmitted_sigmat[1][tx][i][1]= 0.0;
#else
				transmitted_sigmat[1][tx][i][0] = tx_relay_signal[tx][i][0];
				transmitted_sigmat[1][tx][i][1] = tx_relay_signal[tx][i][1];
#endif
			}
		}
	}
}