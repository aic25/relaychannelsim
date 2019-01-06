#include "const.h"
double her[N][TS][2], channel_power[N][2];
extern double CNR;
extern int XI;
extern double channel_gain[N][2], p_in[N];

void mrc_weight(double(*H)[TS][BS][2]){	/* Here, w is selected as h; a kind of matched filter. */
	double h[N][TS][2], Hher[BS][TS][2], R[BS][BS][2], power[2], poweri[2], no, no_loss;

	h_desired(H, h);	/* Choose the channel of desired user */
	VectorHermite(h, her);	/* her = h^H */
	VectorMulVectorToscalar(her, h, channel_power);	/* channel_power = h^H x h */
#ifdef SD
	no = pow(10.0, -CNR / 10.0) / CODING_RATE;
	no_loss = pow(10.0, -(CNR - (double)XI) / 10.0) / CODING_RATE;
	for (size_t f_index = 0; f_index < N; f_index++)
	{
		power[0] = her[f_index][0][0] * her[f_index][0][0] + her[f_index][0][1] * her[f_index][0][1];
		power[1] = her[f_index][1][0] * her[f_index][1][0] + her[f_index][1][1] * her[f_index][1][1];
		MatrixHermite_TSxBS(H[f_index], Hher);
		MatrixMulMatrixToMatrix_BSxTS(Hher, H[f_index], R);
		poweri[0] = 0;
		for (size_t i = 1; i < M; i++)
		{
			poweri[0] += sqr(R[0][i][0]) + sqr(R[0][i][1]);	/* 0.5=no */
		}
		channel_gain[f_index][0] = R[0][0][0];
		channel_gain[f_index][1] = -R[0][0][1];
		p_in[f_index] = (no_loss)*power[0] + (no)*power[1] + poweri[0];
	}
#endif // SD

}


void h_desired(double(*H)[TS][BS][2], double(*h_m)[TS][2]){
	int ts, i;
	for (ts = 0; ts < TS; ts++)
	{
		for (i = 0; i < N; i++)
		{
			h_m[i][ts][0] = H[i][ts][0][0];
			h_m[i][ts][1] = H[i][ts][0][1];
		}
	}
}

void VectorHermite(double(*h)[TS][2], double(*her)[TS][2]){
	int i, ts;

	for (i = 0; i < N; i++)
	{
		for (ts = 0; ts < TS; ts++)
		{
			her[i][ts][0] = h[i][ts][0];
			her[i][ts][1] = -h[i][ts][1];
		}
	}
}


void VectorMulVectorToscalar(double(*her)[TS][2], double(*h)[TS][2], double(*channel_power)[2]){
	int ts, n;
	double sum[2];
	for (n = 0; n < N; n++)
	{
		sum[0] = 0;
		sum[1] = 0;
		for (ts = 0; ts < TS; ts++)
		{
			sum[0] += her[n][ts][0] * h[n][ts][0] - her[n][ts][1] * h[n][ts][1];
			sum[1] += her[n][ts][0] * h[n][ts][1] + her[n][ts][1] * h[n][ts][0];
		}
		channel_power[n][0] = sum[0];
		channel_power[n][1] = sum[1];
	}
}








