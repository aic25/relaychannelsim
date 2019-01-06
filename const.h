#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>
#include <string.h>
#include <string>

/*---------------General Adjustable Parameters--------------------------------*/
#define MMSE						/* MRC - MMSE :: SD FOR MRC unavailable :: There is a normalization difference, as far as I can notice. */
#ifdef RLS
#define DDWE
#endif
#define SD						/* HD - SD								*/
#define DF						/* AF - DF				 */
#define	ZC						/* CUNNING - ZC - Wh :: if CUNNING "ideal channel estimation" else "estimation using XX sequence"	*/
#define LOSS			0		/* Power difference between signal from bs and signal from relay (stronger) 0 - 5 - 10 - 15	*/
#define NOISE					/* Define noise for MS. */
#define RLY_NOISE				/* Define noise for RS */
#define RLY_PROPOGATION			/* Dont undefine noise and propogation at the same time. use switch instead */
#define INTERFERENCE	ON		/* ON - OFF								*/
#define MODULATION		QPSK	/* QPSK - QAM16 - QAM64					*/
#define SWITCH			ON		/* ON - OFF :: OFF => relay re-transmits the exact received signal	*/
#define RLY_CNR			20
#define EbNo_FROM		30		/* Define limits of SNR main loop */
#define EbNo_TO			0		/* for SNR loop  */
#define DBG						/* If defined, activates print commands in fopen functions. */
/*--  RLS  --------------------------------------------------------*/
#define PREAMBLE_LENGTH		4
#define CH_STEP_ADJUST		(Np - PREAMBLE_LENGTH)
#define MAX_DATE			21							/* 12 + 9 Get date to be used in a file name */
#define SIGMA_INVERSE		(1/(pow(10.0, 0 / 10.0)))
#define FILTER_ORDER		(TS-1)						/* RLS filter order. */
#define NUMBER_OF_STEPS		(Np+Nd)						/* Max number of steps for algorithm to converge. It cost two symbol duration for MMSE to estimate the channel. */
#define FORGETTING_FACTOR	1.0							/* lambda */
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
#define FD				0.0		/* doppler shift */
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
#define WH_CS1 (0)
#define WH_CS2 (1)
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
#define Nf		(2*N)				/* ? FFT size :: duble the size of data, is this how it supposed to be? */
#define	GI		72					/* ? GI/Nf ~= 0.07 (what is this ratio? is not that supposed to be N)	*/
#define SAMPLEN	(Nf+GI)				/* Number of samples for one packet										*/
#define Nd		100					/* Number of data packets			*/
#define PACKETN	(Np+Nd)				/* Number of packets per transmission									*/
#define SYMBOL	(PACKETN*SAMPLEN)	/* Number of symbols per transmission									*/
#define N_cbps	(N*MODULATION)		/* Number of bits per packet											*/
#define TAIL	6					/* Known bits at the packet end for viterbi decoder						*/
#define STATEN	64					/* Number of states for viterbi decoder	*/
/*--------------Cellular Environment-------------------------------*/
#define BS		2			/* Number of base stations (BS) and relays										*/
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
void initialization(void);					/* Set noise_power and signal_power to zero							*/
void setting_print();						/* Open a file and write necessary data into it.					*/
void transmitter1_initialize(double(*)[2]);	/* Reset the transmitted_signal vector								*/
void transmitter2_initialize(double(*)[2]);	/* Same as transmitter1_initialize									*/
void transmitter_initialize(double(*)[BS][BURST][2]);	/* Reset the transmitted_signal[TS][BS][BURST][COMPLEX] */
/* Additional print codes.											*/
/*--  encODER.c  --------------------------------------------------*/
void convolutional_encoder(int(*bit), int(*code));	/* Prameters defined here */
void interleaver(int(*input), int(*output));
/*--  tRaNSmITteR1.c  ---------------------------------------------*/
void transmitter1(int(*bit)[N_dbps], double(*)[2]);	/* Prepare OFDM symbol using functions just below														*/
void info_generator(int(*data_bit));				/* Generate random bits. Tail bits set to zero.															*/
void data_insert(int nd, int(*bit)[N_cbps], double(*frequency_signal)[Nf][2]);	/* Assign codewords to constellation points :: Defined in frequency domain. */
void data_insert_trsym(int nd, int(*bit), double(*frequency_signal)[2]);
void data_GI_insert(double(*)[2], double(*)[2]);	/* Add cyclic prefix.																					*/
/*--  tRaNSmITteR2.c  ---------------------------------------------*/
void transmitter2(int(*bit)[N_dbps], double(*)[2]);	/* Same function */
/*--  mUlTIPAth.cpp  ----------------------------------------------*/
void multipath_propagation(double(*transmitted_sigmat)[BS][BURST][2], double(*received_signal)[BURST][2]);	/* Propogates OFDM symbol over multipath channel.											*/
void fading_process(int time, double(*w)[2], double CNR);							/* Calculate fading coefficient using multiple simultaneously received signal components.										*/
void samelevel_PATH(int i, int kk, double(*power)[RAY_COMPONENT]);		/* Normalize the fading coefficient																								*/
void exponential_PATH(int i, int kk, double(*power)[RAY_COMPONENT]);	/* I guess used to calculate fading coefficient using exponential distribution													*/
/*--  TxSignalMatrix.cpp.cpp  -------------------------------------*/	/* There is a serious THING in this code, for loops without curly brackets */
void transmitted_signal_Matrix(double(*)[2], double(*)[2], double(*)[2], double(*)[2], double(*)[BS][BURST][2]);	/* transmitted_signal[TS][BS][BURST][COMPLEX] Transmission between BSs and MSs. */
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
void channel_estimation4ZC_4all_ch(double(*received_signal)[BURST][2], double(*H)[TS][BS][2]);	/* ? */
void time2freq4channelZC2(double(*h)[2], double(*H)[2]);	/* This function does not exist */
/*--  nOIsE.cpp ---------------------------------------------------*/
void noise_addition(double SNR, double(*input_signal)[2], double(*output_signal)[2]);	/*  */
void noise_generator(double n_amplitude, double(*noise));	/*  */
/*--  MMSE - dEtecTOR.cpp  ----------------------------------------*/
void CoherentDetection4SD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][TS][2], double(*diLLR)[N_cbps]);
void CoherentDetection4HD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][TS][2], int(*diLLR)[N_cbps]);
/*--  MMSE - weight.cpp  ------------------------------------------*/
void MMSE_Weight_Matrix(double(*)[TS][BS][2], double(*)[BS][TS][2]);
void MMSE_Weight_Compare(double(*)[TS][BS][2], double(*)[BS][TS][2]);
/*--  RLS.cpp  ----------------------------------------------------*/
void RLS_initial_weigth(int, double(*)[TS][N][2], double(*)[TS][N][2], double(*)[BS][TS][2], double(*)[TS][BS][2]);
void DDWE_Update(double(*received_vector)[2], double(*reference_symbol), double(*old_weights)[2], double(*new_weights)[2], double(*H)[BS][2], double(*e), int f_index);
std::string get_date(void);
const std::string currentDateTime();
/*--  decoder.cpp  ------------------------------------------------*/
void state_model_convolutional();
void SD_Viterbi_decoder(double(*decode_signal), int(*received_bit));
void deinterleaver4SD(double(*input), double(*output));
void HD_Viterbi_decoder(int(*decode), int(*received_bit));
void deinterleaver4HD(int(*input), int(*output));
/*--  receiver.cpp  -----------------------------------------------*/
void receiver(double(*received_signal)[BURST][2], int(*rbit)[N_dbps]);
void DFT(int, int, double(*)[BURST][2], double(*)[N][2]);
/*--  BER.cpp  ----------------------------------------------------*/
void BER(int loop, int(*tbit)[N_dbps], int(*rbit)[N_dbps]);				
/*-------------------- RELAY 1 -------------------------------------*/
void relay1(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]);
void Relay1_DF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]);
void Relay1_AF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]);
void Relay1_propagation(double(*transmitted_signal)[2], double(*received_signal)[2]);
void Relay_fading_process(int time, double(*w)[2]);
void Relay1_receiver(double(*received_signal)[2], int(*rbit)[N_dbps]);
void Relay_ideal_channel_estimation(double(*H_r)[2]);
void Relay_DFT(int data_number, double(*received_signal)[2], double(*DFT_signal)[2]);
void Relay_CoherentDetection(int data_number, double(*DFT_signal)[2], double(*H_r)[2], int(*decode)[N_cbps]);
void Relay1_transmitter(int(*bit)[N_dbps], double(*transmitted_signal)[2]);
void Relay1_receiver_af(double(*received_signal)[2], double(*transmitted_signal)[2]);
void Relay2_receiver_af(double(*received_signal)[2], double(*transmitted_signal)[2]);
/*--------------------- RELAY2 -------------------------------------*/
void relay2(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]);
void Relay2_DF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]);
void Relay2_AF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]);
void Relay2_propagation(double(*transmitted_signal)[2], double(*received_signal)[2]);
void Relay2_receiver(double(*received_signal)[2], int(*rbit)[N_dbps]);
void Relay2_ideal_channel_estimation(double(*H_r)[2]);
void Relay2_transmitter(int(*bit)[N_dbps], double(*transmitted_signal)[2]);
/*--------------------MRC.cpp---------------------------------------*/
void mrc_weight(double(*)[TS][BS][2]);
void h_desired(double(*)[TS][BS][2], double(*)[TS][2]);
void VectorHermite(double(*)[TS][2], double(*)[TS][2]);
void VectorMulVectorToscalar(double(*a)[TS][2], double(*b)[TS][2], double(*c)[2]);
/*---------------Walsh-Hadamard--------------------------------------*/
void WHadamard_initialize(double(*)[BS][BURST][2], double(*)[Np][N]);
void WH_channel_estimation(double(*DFT_signal)[TS][N][2], double(*H)[TS][BS][2], double(*)[Np][N]);
int intPow(int x, int p);
/*--  matrix.cpp  -------------------------------------------------*/
void MatrixInverse(double(*a)[BS][2]);
void MatrixHermite(double(*m1)[BS][2], double(*m2)[TS][2]);
void MatrixTranspose(double(*m1)[BS][2], double(*m2)[TS][2]);
void VectorMulVectorToScalar2(double(*v1)[2], double(*v2)[2], double *s);
void VectorConjugate2(double(*v1)[2], double(*v2)[2]);
void MatrixAddMatrixToMatrix(double(*a)[BS][2], double(*b)[BS][2], double(*c)[BS][2]);
void MatrixSubMatrixToMatrix(double(*)[TS][2], double(*)[TS][2], double(*)[TS][2]);
void MatrixMulMatrixToMatrix(double(*a)[TS][2], double(*b)[BS][2], double(*c)[BS][2]);
void MatrixMulMatrixToMatrix2(double(*a)[BS][2], double(*b)[TS][2], double(*c)[TS][2]);
void MatrixMulVectorToVector(double(*)[TS][2], double(*)[2], double(*)[2]);
void MatrixMulVectorToVector2(double(*)[BS][2], double(*)[2], double(*)[2]);
void ScalarMulMatrixToMatrix(double *, double(*)[TS][2], double(*)[TS][2]);
void VectorMulVectorToScalar(double(*)[2], double(*)[2], double *);
void VectorMulVectorToScalarH(double(*)[2], double(*)[2], double *);
void VectorAddVectorToVector(double(*)[2], double(*)[2], double(*)[2]);
void VectorSubVectorToVector(double(*)[2], double(*)[2], double(*)[2]);
void VectorMulMatrixToVector(double(*)[2], double(*)[TS][2], double(*)[2]);
void VectorMulVectorToMatrix(double(*)[2], double(*)[2], double(*)[TS][2]);
void VectorConjugate(double(*)[2], double(*)[2]);
void ScalarMulVectorToVector(double *, double(*)[2], double(*)[2]);
void ScalarDivScalarToScalar(double *, double *, double *);
double ScalarDistance(double *, double *);
void VectorDistance(double(*)[2], double(*)[2], double);	/* Complex */  
/*------state machine for Viterbi algorithm -----------------------*/
struct StateInfo{
	int output0;	/* output of encoder1 */
	int output1;	/* output of encoder2 */
	int next_state;
	int pre_state;
};