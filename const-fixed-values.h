#define Fourier_MOD			OFF	

#define FILTER_ORDER		(TS-1)									/* RLS filter order. */
#define NUMBER_OF_STEPS		(Np+Nd)									/* Max number of steps for algorithm to converge. It cost two symbol duration for MMSE to estimate the channel. */
#define FORGETTING_FACTOR	0.990										/* lambda */
#define FF_INVERSE			(1/FORGETTING_FACTOR)		  	

#define	RAY_COMPONENT	8											/* Number of simultaneous signal components			*/
#define SAMPLING_RATE	10.0e6										/* 10MHz :: sampling freq = 1/T_s (1/s) Ts = 0.1 us */

#define DBG															/* If defined, activates print commands in fopen functions. */
#define CORRECTION_TRESHOLD 0.01
#define MODULATION		QPSK										/* QPSK - QAM16 - QAM64					*/

#define EXPMODEL		ON											/* ON - OFF											*/

#define ADPTV_LOOP													/* Could cause problems if EbNo_to is negative. Actually it is dynamic loop. Number of loops is a function of SNR. */
#define	LOOPN			(2 * EbNo_FROM + 10)						/* Monte Carlo loop */

#define LOS				OFF
#define RLY_LOS			ON
#define RLY_NOISE													/* Define noise for RS */
#define RLY_PROPOGATION												/* Dont undefine noise and propogation at the same time. use switch instead */

#define DF															/* AF - DF				 */
#define HD															/* HD - SD								*/


#define RLY_EXPMODEL		ON

#define RS_DUR				(1)										/* ? Used in exponential distribution				*/

//#define RS_TAU			(2.8994/RS_PATH)						// 0.6229 0.9139 1.4793 2.8994 for 6 3 0 -3		
 
/*--  Walsh-Hadamard-----------------------------------------------*/
#define WH_CS1 (0)													
#define WH_CS2 (1)
#define WH_CS3 (2)
#ifndef ZC
#define WH_PREAMBLE				/* Use WH sequence instead of ZC */  
#endif // !ZC
#ifdef WH_PREAMBLE
#ifndef CUNNING
#define WH_CHEST				/* Use WH sequence to estimate channel */  
#endif // !CUNNING
#endif // WH_PREAMBLE 
/*------------------------------------------------------------------*/
#define BURST			(PACKETN*SAMPLEN)							/* Length of data to be processed */
#define BERFILE			"STATS.dat"									/* Do not forget to modfify */
/*---------------OFDM Parameters------------------------------------*/
#define SAMPLEN	(Nf+GI)												/* Number of samples for one packet										*/
#define Nf		(N)													/* ? FFT size :: duble the size of data, is this how it supposed to be? */
#define PACKETN	(Np+Nd)												/* Number of packets per transmission									*/
#define SYMBOL	(PACKETN*SAMPLEN)									/* Number of symbols per transmission									*/
#define N_cbps	(N*MODULATION)										/* Number of bits per packet											*/
#define TAIL	6													/* Known bits at the packet end for viterbi decoder						*/
#define TAPLN	GI			/* ? ((PATH-1)*DELAY+1) :: GI=72‚È‚Ì‚É(8-1)*8+1=57‚Å‚·‚¯‚Ç						*/

/*---------------Standard Coefficients------------------------------*/
#define ON			1
#define OFF			0
#define	PI			3.141592654
#define	OneBySqrt2	0.707106781
#define	OneBySqrt10	0.3162277660
#define OneBySqrt42	0.1543033499	 
#define CONST_E		2.71828
#define QPSK		2
#define QAM16		4
#define QAM64		6
#define	sqr(x)		((x)*(x))
#define ELEMENT_COUNT(X) (sizeof(X) / sizeof((X)[0]))

/*
1 OFDM symbol duration =  (1024 + 72) * Ts =   109.6e-6 seconds
Total 1 Packet Length Duration with 10 symbols = 109.6e-6  * 10 symbols = 1096e-6 seconds
*/

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

