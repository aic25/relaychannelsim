#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>
#include <string.h>
#include <string>
#include <iostream>
using namespace std;

/*---------------General Adjustable Parameters--------------------------------*/
#define RLS						/* MRC - MMSE :: SD FOR MRC unavailable :: There is a normalization difference, as far as I can notice. */
#ifdef RLS
//#define DDWE
#endif
#define HD						/* HD - SD	; only HD 							 */
#define DF						/* AF - DF				 */
#define	CUNNING					/* CUNNING - ZC - Wh :: if CUNNING "ideal channel estimation" else "estimation using XX sequence"	*/
#define LOSS			10		/* Power difference between signal from bs and signal from relay (stronger) 0 - 5 - 10 - 15	*/
#define NOISE					/* Define noise for MS. */
#define INTERFERENCE	ON		/* ON - OFF								*/
#define MODULATION		QPSK	/* QPSK - QAM16 - QAM64					*/
#define EbNo_FROM		30		/* Define limits of SNR main loop */
#define EbNo_TO			0		/* for SNR loop  */
#define SWITCH			ON		/* ON - OFF :: OFF => s_r = s_b	*/
#define RLY_NOISE				/* Define noise for RS */
#define RLY_PROPOGATION			/* Dont undefine noise and propogation at the same time. use switch instead */
#define RLY_CNR			30
#define DBG						/* If defined, activates print commands in fopen functions. */
/*--  RLS  --------------------------------------------------------*/
#define PREAMBLE_LENGTH		16
#define MAX_DATE			21							/* 12 + 9 Get date to be used in a file name */
#define SIGMA_INVERSE		(1/(1))
#define FILTER_ORDER		(TS-1)						/* RLS filter order. */
#define NUMBER_OF_STEPS		(Np+Nd)						/* Max number of steps for algorithm to converge. It cost two symbol duration for MMSE to estimate the channel. */
#define FORGETTING_FACTOR	1							/* lambda */
#define FF_INVERSE			(1/FORGETTING_FACTOR)		  
#ifdef RLS	 
#define Np			16
#define Np_POWER2	4
#else								/* For MRC and MMSE */		
#define Np			2	
#define Np_POWER2	1				/* Np < 1/5*Nd ( p ?= parity :: there are two OFDM symbols for training)*/
#endif // RLS						/* for Np = 2->1, 4->2, 8->3, 16->4, 32->5, 64->6, 128->7, 256->8 */ 
/*----------------------CHANNEL------------------------------------*/
//#define DC_CUT				// on-off. idk what this is for
#define EXPMODEL		OFF		/* ON - OFF											*/
#define DUR				10.0	/* ? Used in exponential distribution				*/
#define FD				0.5		/* doppler shift */
#define FDr				0.0
#define	PATH			8		/* Number of paths									*/
#define DELAY			8		/* Uniform delay between signal reflections */
#define	RAY_COMPONENT	8		/* Number of simultaneous signal components			*/
#define SAMPLING_RATE	10.0e6	/* 10MHz :: sampling freq = 1/T_s (1/s) Ts = 0.1 us */
#define ADPTV_LOOP				/* Could cause problems if EbNo_to is negative. Actually it is dynamic loop. Number of loops is a function of SNR. */
#define	LOOPN			(2 * EbNo_FROM + 10)	/* Monte Carlo loop */
/*------------------------------------------------------------------*/
#define	CUNNING_Relay			/* ? if this is defined, relay gets an ideal channel estimation. */
/*--  Walsh-Hadamard-----------------------------------------------*/
#define WH_CS1 0
#define WH_CS2 1
#define WH_CS3 2
#ifndef ZC
#define WH_PREAMBLE				/* Use WH sequence instead of ZC */  
#endif // !ZC
#ifdef WH_PREAMBLE
#ifndef CUNNING
#define WH_CHEST				/* Use WH sequence to estimate channel */  
#endif // !CUNNING
#endif // WH_PREAMBLE 
/*------------------------------------------------------------------*/
#define BURST			(PACKETN*SAMPLEN)	/* Length of data to be processed */
#define BERFILE			"STATS.dat"			/* Do not forget to modfify */
#define DATA_COUNT_MODE						/* Related to some printing option in initialize. nothing else */
/*---------------OFDM Parameters------------------------------------*/
#define N		512					/* Generic for 512, 128, and 32. Effective number of subcarriers (most likely means excluding GI). Changing this could cause stack overflow because of the interleaver.	*/
#define Nf		(4*N)				/* ? FFT size :: duble the size of data, is this how it supposed to be? */
#define	GI		72					/* ? GI/Nf ~= 0.07 (what is this ratio? is not that supposed to be N)	*/
#define SAMPLEN	(Nf+GI)				/* Number of samples for one packet										*/
#define Nd		100					/* Number of data packets			*/
#define PACKETN	(Np+Nd)				/* Number of packets per transmission									*/
#define SYMBOL	(PACKETN*SAMPLEN)	/* Number of symbols per transmission									*/
#define N_cbps	(N*MODULATION)		/* Number of bits per packet											*/
#define TAIL	6					/* Known bits at the packet end for viterbi decoder						*/
#define STATEN	64					/* Number of states for viterbi decoder	*/
/*--------------Cellular Environment-------------------------------*/
#define BS		3			/* Number of base stations (BS) and relays										*/
#define CL		BS			/* Include H3/H1. Combiner Length */
#define KI		2			/* Known Receivers */
#define BTx		1			/* Number of BS and relay antennas												*/
#define Tx		(BS*BTx)	/* Total number of transmitting BS antennas (what of relay?). I doubt the code is designed for MIMO.	*/
#define TS		2			/* time slot (TS)																*/
#define UTx		1			/* Number of desired user antennas												*/
#define Rx		(TS*UTx)	/* Total signaling dimension for desired user									*/
#define M		(BS*BTx)	/* Number of streams on BS side													*/
#define TAPLN	GI			/* ? ((PATH-1)*DELAY+1) :: GI=72‚È‚Ì‚É(8-1)*8+1=57‚Å‚·‚¯‚Ç						*/
/*---------------Standard Coefficients------------------------------*/
#define ON			1
#define OFF			0
#define	PI			3.141592654
#define	OneBySqrt2	0.707106781
#define	OneBySqrt10	0.3162277660
#define OneBySqrt42	0.1543033499
#define QPSK		2
#define QAM16		4
#define QAM64		6
#define	sqr(x)		((x)*(x))

/*
  1 OFDM symbol duration =  (1024 + 72) * Ts =   109.6e-6 seconds
  Total 1 Packet Length Duration with 10 symbols = 109.6e-6  * 10 symbols = 1096e-6 seconds
  */
/*--  CODING  -----------------------------------------------------*/
#define CODING_RATE		0.50
#define N_dbps			(N_cbps/2)	/* Coded bps and Data bps */



/*----------------------------------------------------*/
/*--                 Prototype       ?              --*/
/*----------------------------------------------------*/
short double2short(double x, short q);
double short2double(short x, short q);
/*-----------------------------------------------------------------*/
/*--initialize-----------------------------------------------------*/
void initialization(double(*noise_power), double(*signal_power), double(*));					/* Set noise_power and signal_power to zero							*/
void setting_print();						/* Open a file and write necessary data into it.					*/
void transmitters_initialize(double(*transmitted_signal)[BURST][2]);	/* Reset the transmitted_signal vector								*/
void transmitter_initialize(double(*)[BS][BURST][2]);	/* Reset the transmitted_signal[TS][BS][BURST][COMPLEX] */
/* Additional print codes.											*/
/*--  encODER.c  --------------------------------------------------*/
void convolutional_encoder(int(*bit), int(*code));	/* Prameters defined here */
void interleaver(int(*input), int(*output));
/*--  tRaNSmITteR1.c  ---------------------------------------------*/
void transmitter(int(*bit)[N_dbps], double(*)[2], int bs, double(*BS_signal_power));	/* Prepare OFDM symbol using functions just below														*/
void info_generator(int(*data_bit));				/* Generate random bits. Tail bits set to zero.															*/
void data_insert(int nd, int(*bit)[N_cbps], double(*frequency_signal)[Nf][2]);	/* Assign codewords to constellation points :: Defined in frequency domain. */
void data_insert_regenerated(int nd, int(*bit), double(*frequency_signal)[2]);
void data_GI_insert(double(*)[2], double(*)[2]);	/* Add cyclic prefix.
/*--  mUlTIPAth.cpp  ----------------------------------------------*/
void multipath_propagation(double(*transmitted_sigmat)[BS][BURST][2], double(*Received_Signal)[BURST][2], double CNR, double(*), double(*));	/* Propogates OFDM symbol over multipath channel.											*/
void fading_process(int time, double(*w)[2], double CNR);							/* Calculate fading coefficient using multiple simultaneously received signal components.										*/
void samelevel_PATH(int i, int kk, double(*power)[RAY_COMPONENT]);		/* Normalize the fading coefficient																								*/
void exponential_PATH(int i, int kk, double(*power)[RAY_COMPONENT]);	/* I guess used to calculate fading coefficient using exponential distribution													*/
/*--  TxSignalMatrix.cpp.cpp  -------------------------------------*/	/* There is a serious THING in this code, for loops without curly brackets */
void transmitted_signal_Matrix(double(*)[BURST][2], double(*)[BURST][2], double(*)[BS][BURST][2]);	/* transmitted_signal[TS][BS][BURST][COMPLEX] Transmission between BSs and MSs. */
/*--  fFT.Cpp  ----------------------------------------------------*/
void FFT_initialize();						/* Reset FFT matrix */
void IFFT(double(*)[2], double(*)[2]);		/* Size = Nf		*/
void FFT(double(*)[2], double(*)[2]);		/* Size = Nf		*/
void IFFT_N(double(*)[2], double(*)[2]);	/* Size = N			*/
void FFT_N(double(*)[2], double(*)[2]);		/* Size = N			*/
/*--  zC_LAYEr.CPp  -----------------------------------------------*/
void Zadoff_Chu_layer_initialize();	/* Define zc_signal_a and zc_signal_b with unknown functions */
void Zadoff_Chu_generator_layer1(double(*)[BS][BURST][2]);	/* Layer one and layer two seems identical */
void Zadoff_Chu_generator_layer2(double(*)[BS][BURST][2]);	/**/
/*--  EStIMatIon.cPP ----------------------------------------------*/
void ideal_channel_estimation4all_ch(double(*H)[TS][BS][2]);	/* This one estimates ideally */
void ideal_channel_estimation4all_ch_MS2(double(*H)[TS][BS][2]);	/* This one estimates ideally */
void channel_estimation4ZC_4all_ch(double(*Received_Signal)[BURST][2], double(*H)[TS][BS][2]);	/* ? */
void time2freq4channelZC2(double(*h)[2], double(*H)[2]);	/* This function does not exist */
/*--  nOIsE.cpp ---------------------------------------------------*/
void noise_addition(double SNR, double(*input_signal)[2], double(*output_signal)[2], double(*noise_power));	/*  */
void noise_generator(double n_amplitude, double(*noise), double(*noise_power));	/*  */
/*--  MMSE - dEtecTOR.cpp  ----------------------------------------*/
void CoherentDetection4SD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][CL][2], double(*diLLR)[N_cbps]);
void CoherentDetection4HD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][TS][2], int(*diLLR)[Nd][N_cbps], int(*Decision));
/*--  MMSE - weight.cpp  ------------------------------------------*/
void MMSE_Weight_Matrix(double(*)[TS][BS][2], double(*)[BS][TS][2]);
void MMSE_Weight_Compare(double(*)[TS][BS][2], double(*)[BS][TS][2]);
/*--  RLS.cpp  ----------------------------------------------------*/
void RLS_initial_weigth(int, double(*)[TS][N][2], double(*)[TS][N][2], double(*)[BS][TS][2], double(*)[TS][BS][2], double(*wh_sequence)[Np][N], double(*e)[NUMBER_OF_STEPS][2]);
void DDWE_Update(double(*received_vector)[2], double(*reference_symbol), double(*old_weights)[2], double(*new_weights)[2], double(*H)[BS][2], double(*e)[NUMBER_OF_STEPS][2], double(*reference_symbol_bs3), int f_index, int);
const std::string currentDateTime();
/*--  decoder.cpp  ------------------------------------------------*/
void state_model_convolutional();
void SD_Viterbi_decoder(double(*decode_signal), int(*received_bit));
void deinterleaver4SD(double(*input), double(*output));
void HD_Viterbi_decoder(int(*decode), int(*received_bit));
void deinterleaver4HD(int(*input), int(*output));
/*--  receiver.cpp  -----------------------------------------------*/
void receiver(double(*Received_Signal)[BURST][2], int(*rbit)[Nd][N_dbps], double(*wh_sequence)[Np][N]);
void interference_detection(double(*Received_Signal)[BURST][2]);
void DFT(int, int, double(*)[BURST][2], double(*)[N][2]);
/*--  BER.cpp  ----------------------------------------------------*/
void BER(int loop, int maxloop, double CNR, int(*tbit)[Nd][N_dbps], int(*rbit)[Nd][N_dbps], double noise_power, double signal_power, double interference_power);
/*-------------------- RELAY 1 -------------------------------------*/
void relay(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs, double(*BS_signal_power));
void Relay_DF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs);
void Relay_AF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs, double(*BS_signal_power));
void Relay_propagation(double(*transmitted_signal)[2], double(*Received_Signal)[2], int bs);
void Relay_fading_process(int time, double(*w)[2]);
void Relay_receiver(double(*Received_Signal)[2], int(*rbit)[N_dbps], int bs);
void Relay_ideal_channel_estimation(double(*H_r)[2], int bs);
void Relay_DFT(int data_number, double(*Received_Signal)[2], double(*DFT_signal)[2]);
void Relay_CoherentDetection(int data_number, double(*DFT_signal)[2], double(*H_r)[2], int(*decode)[N_cbps]);
void Relay_transmitter(int(*bit)[N_dbps], double(*transmitted_signal)[2]);
/*--------------------MRC.cpp---------------------------------------*/
void mrc_weight(double(*)[TS][BS][2]);
void h_desired(double(*)[TS][BS][2], double(*)[TS][2]);
void VectorHermite(double(*)[TS][2], double(*)[TS][2]);
void VectorMulVectorToscalar(double(*a)[TS][2], double(*b)[TS][2], double(*c)[2]);
/*---------------Walsh-Hadamard--------------------------------------*/
void WHadamard_initialize(double(*)[BS][BURST][2], double(*)[Np][N]);
void C_WHadamard_initialize(double(*)[BS][BURST][2], double(*)[Np][N][2]);
void WH_channel_estimation(double(*DFT_signal)[TS][N][2], double(*H)[TS][BS][2], double(*)[Np][N]);
int intPow(int x, int p);
/*--  matrix_cl.cpp  -------------------------------------------------*/
void MatrixInverse_CLxCL(double(*a)[CL][2]);
void MatrixHermite_CLxCL(double(*m1)[CL][2], double(*m2)[CL][2]);
void MatrixTranspose_CLxCL(double(*m1)[CL][2], double(*m2)[CL][2]);
void MatrixAddMatrixToMatrix_CLxCL(double(*a)[CL][2], double(*b)[CL][2], double(*c)[CL][2]);
void MatrixSubMatrixToMatrix_CL(double(*)[CL][2], double(*)[CL][2], double(*)[CL][2]);
void MatrixMulMatrixToMatrix_CLxCL(double(*a)[CL][2], double(*b)[CL][2], double(*c)[CL][2]);
void MatrixMulMatrixToMatrix_CLxCLXCLxCL(double(*a)[CL][2], double(*b)[CL][2], double(*c)[CL][2]);
void MatrixMulVectorToVector_CLxCL(double(*)[CL][2], double(*)[2], double(*)[2]);
void MatrixMulVectorToVector_CLxCL(double(*)[CL][2], double(*)[2], double(*)[2]);
void ScalarMulMatrixToMatrix_CLxCL(double *, double(*)[CL][2], double(*)[CL][2]);
void VectorMulVectorToScalar_CL(double(*)[2], double(*)[2], double *);
void VectorMulVectorToScalarH_CL(double(*)[2], double(*)[2], double *);
void VectorAddVectorToVector_CL(double(*)[2], double(*)[2], double(*)[2]);
void VectorSubVectorToVector_CL(double(*)[2], double(*)[2], double(*)[2]);
void VectorMulMatrixToVector_CLxCL(double(*)[2], double(*)[CL][2], double(*)[2]);
void VectorMulVectorToMatrix_CLxCL(double(*)[2], double(*)[2], double(*)[CL][2]);
void VectorConjugate_CL(double(*)[2], double(*)[2]);
void ScalarMulVectorToVector_CL(double *, double(*)[2], double(*)[2]);
void VectorDistance_CL(double(*)[2], double(*)[2], double);	/* Complex */
/*--  matrix.cpp  -------------------------------------------------*/
void MatrixInverse_BSxBS(double(*a)[BS][2]);
void MatrixHermite_TSxBS(double(*m1)[BS][2], double(*m2)[TS][2]);
void MatrixTranspose_TSxBS(double(*m1)[BS][2], double(*m2)[TS][2]);
void VectorMulVectorToScalar_M(double(*v1)[2], double(*v2)[2], double *s);
void VectorConjugate_M(double(*v1)[2], double(*v2)[2]);
void MatrixAddMatrixToMatrix_BSxBS(double(*a)[BS][2], double(*b)[BS][2], double(*c)[BS][2]);
void MatrixSubMatrixToMatrix_TS(double(*)[TS][2], double(*)[TS][2], double(*)[TS][2]);
void MatrixMulMatrixToMatrix_BSxTS(double(*a)[TS][2], double(*b)[BS][2], double(*c)[BS][2]);
void MatrixMulMatrixToMatrix_BSxBSXBSxTS(double(*a)[BS][2], double(*b)[TS][2], double(*c)[TS][2]);
void MatrixMulVectorToVector_TSxTS(double(*)[TS][2], double(*)[2], double(*)[2]);
void MatrixMulVectorToVector_BSxBS(double(*)[BS][2], double(*)[2], double(*)[2]);
void ScalarMulMatrixToMatrix_TSxTS(double *, double(*)[TS][2], double(*)[TS][2]);
void VectorMulVectorToScalar_TS(double(*)[2], double(*)[2], double *);
void VectorMulVectorToScalarH_TS(double(*)[2], double(*)[2], double *);
void VectorAddVectorToVector_TS(double(*)[2], double(*)[2], double(*)[2]);
void VectorSubVectorToVector_TS(double(*)[2], double(*)[2], double(*)[2]);
void VectorMulMatrixToVector_TSxTS(double(*)[2], double(*)[TS][2], double(*)[2]);
void VectorMulVectorToMatrix_TSxTS(double(*)[2], double(*)[2], double(*)[TS][2]);
void VectorConjugate_TS(double(*)[2], double(*)[2]);
void ScalarMulVectorToVector_TS(double *, double(*)[2], double(*)[2]);
void ScalarDivScalarToScalar(double *, double *, double *);
double ScalarDistance(double *, double *);
void VectorDistance_TS(double(*)[2], double(*)[2], double);	/* Complex */
/*------state machine for Viterbi algorithm -----------------------*/
struct StateInfo{
	int output0;	/* output of encoder1 */
	int output1;	/* output of encoder2 */
	int next_state;
	int pre_state;
};