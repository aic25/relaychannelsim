
#include "const.h"

/*----------------------------------------------------*/
/*--                 Prototype       ?              --*/
/*----------------------------------------------------*/
short double2short(double x, short q);
double short2double(short x, short q);
/*-----------------------------------------------------------------*/
/*--initialize-----------------------------------------------------*/		
void init_channel_models();
void initialization(void);					/* Set noise_power and signal_power to zero							*/
void setting_print();						/* Open a file and write necessary data into it.					*/
void transmitters_initialize(double(*transmitted_signal)[BURST][2]);	/* Reset the transmitted_signal vector								*/
void transmitter_initialize(double(*)[BS][BURST][2]);	/* Reset the transmitted_signal[TS][BS][BURST][COMPLEX] */
void printSignal(const char* Name, double Signal[/* SignalLen */][2], size_t SignalLen);
/* Additional print codes.											*/
/*--  encODER.c  --------------------------------------------------*/
void convolutional_encoder(int(*bit), int(*code));	/* Prameters defined here */
void interleaver(int(*input), int(*output));
/*--  tRaNSmITteR1.c  ---------------------------------------------*/
void transmitter(int(*bit)[N_dbps], double(*)[2], int bs);	/* Prepare OFDM symbol using functions just below														*/
void info_generator(int(*data_bit));				/* Generate random bits. Tail bits set to zero.															*/
void data_insert(int nd, int(*bit)[N_cbps], double(*frequency_signal)[Nf][2]);	/* Assign codewords to constellation points :: Defined in frequency domain. */
void data_insert_regenerated(int nd, int(*bit), double(*frequency_signal)[2]);
void data_GI_insert(double(*)[2], double(*)[2]);	/* Add cyclic prefix.
													/*--  mUlTIPAth.cpp  ----------------------------------------------*/
void multipath_propagation(double(*transmitted_sigmat)[BS][BURST][2], double(*received_signal)[BURST][2]);	/* Propogates OFDM symbol over multipath channel.											*/
void fading_process(int time, double(*w)[2], double CNR);							/* Calculate fading coefficient using multiple simultaneously received signal components.										*/
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
void channel_estimation4ZC_4all_ch(double(*received_signal)[BURST][2], double(*H)[TS][BS][2]);	/* ? */
void time2freq4channelZC2(double(*h)[2], double(*H)[2]);	/* This function does not exist */
/*--  nOIsE.cpp ---------------------------------------------------*/
void noise_addition(double SNR, double(*input_signal)[2], double(*output_signal)[2], int count_for_power);	/*  */
void noise_generator(double n_amplitude, double(*noise), int count_for_power);	/*  */
/*--  MMSE - dEtecTOR.cpp  ----------------------------------------*/
void CoherentDetection4SD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][TS][2], double(*diLLR)[N_cbps]);
void CoherentDetection4HD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][CL][2], double(*H)[TS][BS][2], int(*diLLR)[Nd][N_cbps], int(*Decision), int(*) ,int(*));
/*--  MMSE - weight.cpp  ------------------------------------------*/
void MMSE_Weight_Matrix(double(*)[TS][BS][2], double(*)[BS][TS][2]);
void MMSE_Weight_Compare(double(*)[TS][BS][2], double(*)[BS][TS][2]);
/*--  RLS.cpp  ----------------------------------------------------*/
void RLS_CombCoeff(int, int, int, double(*)[TS][N][2], double(*)[TS][N][2], double(*)[BS][CL][2], double(*)[TS][BS][2], double(*e)[NUMBER_OF_STEPS][2]);
void DF_CombCoeff(int, int f_index, int(*Decision), int R_index, int, double(*received_vector)[2], double(*reference_symbol), double(*reference_symbol_bs3), double(*old_weights)[2], double(*new_weights)[2], double(*H)[BS][2], double(*e));
void RLS_Error(double e_rls[N][NUMBER_OF_STEPS][2]);
const std::string currentDateTime();
/*--  decoder.cpp  ------------------------------------------------*/
void state_model_convolutional();
void SD_Viterbi_decoder(double(*decode_signal), int(*received_bit));
void deinterleaver4SD(double(*input), double(*output));
void HD_Viterbi_decoder(int(*decode), int(*received_bit));
void deinterleaver4HD(int(*input), int(*output));
/*--  receiver.cpp  -----------------------------------------------*/
void receiver(double(*received_signal)[BURST][2], int(*rbit)[Nd][N_dbps]);
void DFT(int, int, double(*)[BURST][2], double(*)[N][2]);
/*--  BER.cpp  ----------------------------------------------------*/
void BER_RsPlusBs(int loop, int(*rsbit)[Nd][N_dbps], int(*tbit)[Nd][N_dbps], int(*rbit)[Nd][N_dbps]);
/*-------------------- RELAY 1 -------------------------------------*/
void relay(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs);
void Relay_DF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs);
void Relay_propagation(double(*transmitted_signal)[2], double(*received_signal)[2], int bs);
void Relay_fading_process(int time, double(*w)[2]);
void Relay_receiver(double(*received_signal)[2], int(*rbit)[N_dbps], int bs);
void Relay_ideal_channel_estimation(double(*H_r)[2], int bs);
void Relay_DFT(int data_number, double(*received_signal)[2], double(*DFT_signal)[2]);
void Relay_CoherentDetection(int data_number, double(*DFT_signal)[2], double(*H_r)[2], int(*decode)[N_cbps]);
void Relay_transmitter(int(*bit)[N_dbps], double(*transmitted_signal)[2]);
/*--------------------MRC.cpp---------------------------------------*/
void mrc_weight(double(*)[TS][BS][2]);
void h_desired(double(*)[TS][BS][2], double(*)[TS][2]);
void VectorHermite(double(*)[TS][2], double(*)[TS][2]);
void VectorMulVectorToscalar(double(*a)[TS][2], double(*b)[TS][2], double(*c)[2]);
/*---------------Walsh-Hadamard--------------------------------------*/
void WHadamard_initialize(double(*)[BS][BURST][2], double(*)[Np][N][2], double(*transmitted_signal_actual)[BURST][2]);
void WH_channel_estimation(double(*DFT_signal)[TS][N][2], double(*H)[TS][BS][2], double(*)[Np][N][2]);
void WH_chest_update(double(*DFT_signal)[2], double(*H)[BS][2], double detected_signal[2], int f_index);
int intPow(int x, int p);
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
void MatrixMulMatrixToMatrix_BSxBSXBSxCL(double(*a)[BS][2], double(*b)[CL][2], double(*c)[CL][2]);
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
void convolve(const double Signal[/* SignalLen */][2], size_t SignalLen,
	const double Kernel[/* KernelLen */][2], size_t KernelLen,
	double Result[/* SignalLen + KernelLen - 1 */][2]);
void zerotap(const double Signal[/* SignalLen */][2],
	double SignalExtended[/* SignalLen * SignalDelay */][2],
	size_t SignalLen, size_t SignalDelay);
void ScalarMulScalar(double(*m1), double(*m2), double(*r));
void ScalarMulScalar_H(double(*m1), double(*m2), double(*r));
void MatrixEquate(double(*Matrix_Bir)[3][2], const double Matrix_Iki[2][3][2]);
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
void MatrixMulVectorToVector_3x2(double(*a)[3][2], double(*v1)[2], double(*v2)[2]);
void ScalarMulMatrixToMatrix_CLxCL(double *, double(*)[CL][2], double(*)[CL][2]);
void VectorMulVectorToScalar_CL(double(*)[2], double(*)[2], double *);
void VectorMulVectorToScalarH_CL(double(*)[2], double(*)[2], double *);
void VectorAddVectorToVector_CL(double(*)[2], double(*)[2], double(*)[2]);
void VectorSubVectorToVector_CL(double(*)[2], double(*)[2], double(*)[2]);
void VectorMulMatrixToVector_CLxCL(double(*)[2], double(*)[CL][2], double(*)[2]);
void VectorMulMatrixToVector_TSxCL(double(*)[2], double(*)[CL][2], double(*)[2]);
void VectorMulVectorToMatrix_CLxCL(double(*)[2], double(*)[2], double(*)[CL][2]);
void VectorConjugate_CL(double(*)[2], double(*)[2]);
void ScalarMulVectorToVector_CL(double *, double(*)[2], double(*)[2]);
void VectorDistance_CL(double(*)[2], double(*)[2], double);	/* Complex */
double MagnitudeOfDifference_CL(double(*v1)[2], double(*v2)[2]);
void VectorEqualVector_CL(double(*v1)[2], double(*v2)[2]);	   
void MatrixInverse_4x4(double(*a)[4][2]);
void MatrixHermite_4x4(double(*m1)[4][2], double(*m2)[4][2]);
void MatrixTranspose_4x4(double(*m1)[4][2], double(*m2)[4][2]);
void MatrixAddMatrixToMatrix_4x4(double(*a)[4][2], double(*b)[4][2], double(*c)[4][2]);
void MatrixSubMatrixToMatrix_4(double(*)[4][2], double(*)[4][2], double(*)[4][2]);
void MatrixMulMatrixToMatrix_4x4(double(*a)[4][2], double(*b)[4][2], double(*c)[4][2]);
void MatrixMulMatrixToMatrix_4x4X4x4(double(*a)[4][2], double(*b)[4][2], double(*c)[4][2]);
void MatrixMulVectorToVector_4x4(double(*)[4][2], double(*)[2], double(*)[2]);
void MatrixMulVectorToVector_4x4(double(*)[4][2], double(*)[2], double(*)[2]);
void ScalarMulMatrixToMatrix_4x4(double *, double(*)[4][2], double(*)[4][2]);
void VectorMulVectorToScalar_4(double(*)[2], double(*)[2], double *);
void VectorMulVectorToScalarH_4(double(*)[2], double(*)[2], double *);
void VectorAddVectorToVector_4(double(*)[2], double(*)[2], double(*)[2]);
void VectorSubVectorToVector_4(double(*)[2], double(*)[2], double(*)[2]);
void VectorMulMatrixToVector_4x4(double(*)[2], double(*)[4][2], double(*)[2]);
void VectorMulVectorToMatrix_4x4(double(*)[2], double(*)[2], double(*)[4][2]);
void VectorConjugate_4(double(*)[2], double(*)[2]);
void ScalarMulVectorToVector_4(double *, double(*)[2], double(*)[2]);
void VectorDistance_4(double(*)[2], double(*)[2], double);	/* Complex */
double MagnitudeOfDifference_4(double(*v1)[2], double(*v2)[2]);
void VectorEqualVector_4(double(*v1)[2], double(*v2)[2]);
/*------state machine for Viterbi algorithm -----------------------*/
struct StateInfo{
	int output0;	/* output of encoder1 */
	int output1;	/* output of encoder2 */
	int next_state;
	int pre_state;
};
/*--------dFt.cPp-----------*/
void DFT_initialize();
void IDFT(double(*input)[2], double(*output)[2]);
void DFT(double(*input)[2], double(*output)[2]);