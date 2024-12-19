#include<stdio.h>
#include<stdlib.h>


double run_plc(char filename[],int rank1)
{

  char command[1000];
  double A[20],B[40];
  

  A[0] = 46.054;  //A_cib_217
  A[1] = -1.3;   //cib_index
  A[2] = 0.656701;    //xi_sz_cib
  A[3] = 7.0838;    //A_sz
  A[4] = 248.212;  //ps_A_100_100
  A[5] = 50.6799;  //ps_A_143_143
  A[6] = 53.3096;  //ps_A_143_217
  A[7] = 121.861;  //ps_A_217_217
  A[8] = 0.00384693;    //ksz_norm

  A[9] = 8.80249;    //gal545_A_100
  A[10] = 11.0066;   //gal545_A_143
  A[11] = 20.1589;  //gal545_A_143_217
  A[12] = 95.4504;  //gal545_A_217
  
  A[13] = 1.0;   //A_sbpx_100_100_TT
  A[14] = 1.0;   //A_sbpx_143_143_TT
  A[15] = 1.0;   //A_sbpx_143_217_TT
  A[16] = 1.0;   //A_sbpx_217_217_TT

  A[17] = 0.999742; //calib_100T
  A[18] = 0.998187; //calib_217T
  A[19] = 1.00044;   //A_planck



  B[0] = 46.054;  //A_cib_217
  B[1] = -1.3;   //cib_index
  B[2] = 0.656701;    //xi_sz_cib
  B[3] = 7.0838;    //A_sz
  B[4] = 248.212;  //ps_A_100_100
  B[5] = 50.6799;  //ps_A_143_143
  B[6] = 53.3096;  //ps_A_143_217
  B[7] = 121.861;  //ps_A_217_217
  B[8] = 0.00384693;    //ksz_norm

  B[9] = 8.80249;    //gal545_A_100
  B[10] = 11.0066;   //gal545_A_143
  B[11] = 20.1589;   //gal545_A_143_217
  B[12] = 95.4504;   //gal545_A_217
  

  B[13] = 0.055;    // galf_EE_A_100
  B[14] = 0.04;     // galf_EE_A_100_143
  B[15] = 0.094;    // galf_EE_A_100_217
  B[16] = 0.086;    // galf_EE_A_143
  B[17] = 0.21;     // galf_EE_A_143_217
  B[18] = 0.7;      // galf_EE_A_217
  B[19] = -2.4;     // galf_EE_index
  B[20] = 0.113825; // galf_TE_A_100
  B[21] = 0.134564; // galf_TE_A_100_143
  B[22] = 0.478734; // galf_TE_A_100_217
  B[23] = 0.224849; // galf_TE_A_143
  B[24] = 0.664945; // galf_TE_A_143_217
  B[25] = 2.08249;  // galf_TE_A_217
  B[26] = -2.4;     // gal  f_TE_index

  B[27] = 1; // A_cnoise_e2e_100_100_EE = 1
  B[28] = 1; // A_cnoise_e2e_143_143_EE = 1
  B[29] = 1; // A_cnoise_e2e_217_217_EE = 1
  B[30] = 1; // A_sbpx_100_100_TT = 1
  B[31] = 1; // A_sbpx_143_143_TT = 1
  B[32] = 1; // A_sbpx_143_217_TT = 1
  B[33] = 1; // A_sbpx_217_217_TT = 1
  B[34] = 1; // A_sbpx_100_100_EE = 1
  B[35] = 1; // A_sbpx_100_143_EE = 1
  B[36] = 1; // A_sbpx_100_217_EE = 1
  B[37] = 1; // A_sbpx_143_143_EE = 1
  B[38] = 1; // A_sbpx_143_217_EE = 1
  B[39] = 1; // A_sbpx_217_217_EE = 1


  B[40] = 0.999742; //calib_100T
  B[41] = 0.998187; //calib_217T

  B[42] = 1.021; //calib_100P = 1.021
  B[43] = 0.966; //calib_143P = 0.966
  B[44] = 1.04;  //calib_217P = 1.04
  B[45] = 1.0;  //A_pol = 1

  B[46] = 1.00044;   //A_planck


  FILE *readfile, *TTwrite, *TEBwrite, *TTTEEEwrite;
  char writeTT[200], writeTEB[200], writeTTTEEE[200]; //, basefolder[150];

  // int seed;
  // FILE *fseed;
  // fseed = fopen("seed.in","r");
  // fscanf(fseed,"%d",&seed);
  // fclose(fseed);

  // sprintf(basefolder,"/rds/general/user/sdas5/home/EverpresentLambda/SCoPE/seed_%d", seed);
  printf("\n%s",filename);
  readfile = fopen(filename,"r");
	
  sprintf(writeTTTEEE,"%s/TTTEEE_%d.d",BASEFOLDER,rank1);
  // printf("\n%s",writeTTTEEE);

  TTTEEEwrite = fopen(writeTTTEEE,"w");

  int i,temp;
  double TTd,TEd,EEd,BBd,Td;
  double TTtemp[2508]={0.0},TEtemp[2508]={0.0},EEtemp[2508]={0.0},BBtemp[2508]={0.0};


  for(i=2;i<=2508;i++)
  {
     fscanf(readfile,"%d %lf %lf %lf %lf",&temp,&TTd,&EEd,&BBd,&TEd);
     TTtemp[i] = TTd/(double)(temp*(temp+1)/2)*3.14159256;
     TEtemp[i] = TEd/(double)(temp*(temp+1)/2)*3.14159256;
     EEtemp[i] = EEd/(double)(temp*(temp+1)/2)*3.14159256;
     BBtemp[i] = BBd/(double)(temp*(temp+1)/2)*3.14159256;
     }


  for(i=0;i<=2508;i++)
     fprintf(TTTEEEwrite,"%e\n",TTtemp[i]);
  for(i=0;i<=2508;i++)
     fprintf(TTTEEEwrite,"%e\n",EEtemp[i]);
  for(i=0;i<=2508;i++)
     fprintf(TTTEEEwrite,"%e\n",TEtemp[i]);

  for(i=0;i<47;i++)
     fprintf(TTTEEEwrite,"%e\n",B[i]);

  fclose(readfile);
  fclose(TTTEEEwrite);

  sprintf(command,"rm %s",filename);
  
  sprintf(command,"/rds/general/user/sdas5/home/EverpresentLambda/plancklikelihood/code/plc_3.0/plc-3.01/bin/clik_example_C /rds/general/user/sdas5/home/EverpresentLambda/plancklikelihood/baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik %s > clickdump.out",writeTTTEEE);
  system(command);


  sprintf(command,"%s.like",writeTTTEEE);

  // printf("I am my command :  %s",command);

  FILE *fpchk;
l1k:
  fpchk = fopen(command,"r");
  if(fpchk == NULL)
    goto l1k;


  double result1;
  fscanf(fpchk,"%lf",&result1);
  fclose(fpchk);

  printf("Chisquare is (From runplc): %lf\n",result1);
  
  return -2*(result1);
}