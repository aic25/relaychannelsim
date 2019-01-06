#include <iostream>		  
#include <Eigen/Eigenvalues>
using namespace Eigen;
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>
#include <string.h>
#include <string>	   
#include "const-fixed-values.h"		
using namespace std;

/*---------------General Adjustable Parameters--------------------------------*/
//#define MMSE														/* MRC - MMSE :: SD FOR MRC unavailable :: There is a normalization difference, as far as I can notice. */
#define RLS
//#define DDWE			/*- no define--*/							/*define only if rls is defined.*/

#define CUNNING
//#define Wh 														/* CUNNING - ZC - Wh :: if CUNNING "ideal channel estimation" else "estimation using XX sequence"	*/
//#define ZC

#define LOSS				10										/* Power difference between signal from bs and signal from relay (stronger) 0 - 5 - 10 - 15	*/
#define NOISE														/* Define noise for MS. */
#define INTERFERENCE		ON										/* ON - OFF								*/
#define PROPAGATION
#define MULTIPATH
#define EbNo_FROM			35										/* Define limits of SNR main loop */
#define EbNo_TO				0										/* for SNR loop  */
/*--  RLS  --------------------------------------------------------*/
#define Np					4											/* for Np = 2->1, 4->2, 8->3, 16->4, 32->5, 64->6, 128->7, 256->8 */
#define Np_POWER2			2											/* Np < 1/5*Nd ( p ?= parity :: there are two OFDM symbols for training)*/					 
#define PREAMBLE_LENGTH		4
#define Metric_Order		1
/*----------------------CHANNEL------------------------------------*/																//#define DC_CUT				// on-off. idk what this is for
#define DUR					10.0									/* ? Used in exponential distribution				*/
#define FD					0.0										/* doppler shift */
#define	PATH				15										/* Number of paths									*/
#define DELAY				1										/* Uniform delay between signal reflections */
/*----------------------------Relay Parameters----------------------------------------*/
#define RLY_CNR				30
#define RS_PATH				10
#define RS_DELAY			1										
#define RS_RAY_COMPONENT	8																	 
#define FDr					0.0											 
#define RS_TAU				(0.42)
#define TAU					(1.44)
/*---------OFDM and Viterbi---------*/
#define N					128										/* Generic for 512, 128, and 32. Effective number of subcarriers (most likely means excluding GI). Changing this could cause stack overflow because of the interleaver.	*/
#define	GI					16										/* ? GI/Nf ~= 0.07 (what is this ratio? is not that supposed to be N)	*/
#define Nd					100										/* Number of data packets			*/
#define STATEN				64										/* Number of states for viterbi decoder	*/
/* Both is off for AF. Both is on for DF */
#define RS_DF				ON										/* Both is off for AF. Both is on for DF */
#define SWITCH				ON										/* ON - OFF :: OFF => relay re-transmits the exact received signal	*/		 

