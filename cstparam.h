#ifndef _CSTPARAM_
#define _CSTPARAM_

//
// This part is for Causal Set Theory analysis
//----------------------------------------------------------------------------------------
#define N_AGRID 50000
char BASEFOLDER[1000];
int SEEDNUMBER;
double AlphaOmOmmax[21][3], a_omega_InterP[21][N_AGRID], a_InterP[N_AGRID], a_omega_interpolated[N_AGRID];
//AlphaOmOmmax - alpha, OmegaM, MaxOmegaM

//
// Supernova Data
//----------------------------------------------------------------------------------------
#define N_SUPERNOVA 1048
double supernovaStr[N_SUPERNOVA][6]; 

#endif
