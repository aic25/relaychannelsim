#include "const.h"

void transmitted_signal_Matrix(
	double(*transmitted_signal)[2],
	double(*transmitted_signal2)[2],
	double(*tx_relay_signal)[2],
	double(*tx_relay_signal2)[2],
	double(*transmitted_sigmat)[Tx][BURST][2]
	){

	int i, ts, tx;
	for (ts = 0; ts < Rx; ts++){
		for (tx = 0; tx < Tx; tx++){
			for (i = (Np * SAMPLEN); i < BURST; i++){	/* Apply for data, not for parity symbols. */

				if (ts == 0 && tx == 0){
					transmitted_sigmat[ts][tx][i][0] = transmitted_signal[i][0];
					transmitted_sigmat[ts][tx][i][1] = transmitted_signal[i][1];
				}
				else if (ts == 0 && tx == 1){
#if (INTERFERENCE==OFF)
					transmitted_sigmat[ts][tx][i][0]=  0.0;
					transmitted_sigmat[ts][tx][i][1]=  0.0;
#endif
#if (INTERFERENCE==ON)
					transmitted_sigmat[ts][tx][i][0] = transmitted_signal2[i][0];
					transmitted_sigmat[ts][tx][i][1] = transmitted_signal2[i][1];
#endif
				}
				else if (ts == 1 && tx == 0){
					transmitted_sigmat[ts][tx][i][0] = tx_relay_signal[i][0];
					transmitted_sigmat[ts][tx][i][1] = tx_relay_signal[i][1];
				}
				else if(ts == 1 && tx == 1){
#if (INTERFERENCE==OFF)
					transmitted_sigmat[ts][tx][i][0]= 0.0; 
					transmitted_sigmat[ts][tx][i][1]= 0.0;
#endif
#if (INTERFERENCE==ON)
					transmitted_sigmat[ts][tx][i][0] = tx_relay_signal2[i][0];
					transmitted_sigmat[ts][tx][i][1] = tx_relay_signal2[i][1];
#endif
				}
				else if (ts == 2 && tx == 0){
					transmitted_sigmat[ts][tx][i][0] = transmitted_signal[i][0];
					transmitted_sigmat[ts][tx][i][1] = transmitted_signal[i][1];
				}
				else if (ts == 2 && tx == 1){
#if (INTERFERENCE==OFF)
					transmitted_sigmat[ts][tx][i][0] = 0.0;
					transmitted_sigmat[ts][tx][i][1] = 0.0;
#endif
#if (INTERFERENCE==ON)
					transmitted_sigmat[ts][tx][i][0] = transmitted_signal2[i][0];
					transmitted_sigmat[ts][tx][i][1] = transmitted_signal2[i][1];
#endif
				}
			}
		}
	}
}