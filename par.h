#ifndef _PAR_H_
#define _PAR_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "cstparam.h"
// #include "CST_eq_state.h"

double rombint(double (*func)(double, double, double), double a, double b, double tol, double h0, double omegam) // trapezoid approx to integral of func in interval [a,b]
{
    int MAXJ = 5;
    int MAXITER = 30, nint, i, j, k, jmax;
    double g[MAXJ + 1], h, gmax, g0, fourj, g1, error;

    h = 0.5 * (b - a);
    gmax = h * (func(a, h0, omegam) + func(b, h0, omegam)); // first trapezoid approx to integral of func in interval [a,b]
    g[0] = gmax;
    nint = 1;
    error = 1.e20;


    for (i = 0; !((i > MAXITER) || ((i > 5) && (fabs(error) < tol))); i++)
    {
        /*     Calculate next trapezoidal rule approximation to integral. */
        g0 = 0.0;
        for (k = 1; k <= nint; k++)
            g0 += func(a + (k + k - 1) * h, h0, omegam);

        g0 = 0.5 * g[0] + h * g0;
        h = 0.5 * h;
        nint *= 2;
        jmax = (i < MAXJ) ? i : MAXJ;
        fourj = 1.;

        for (j = 1; j <= jmax; j++)
        /* Use Richardson extrapolation. */
        {
            fourj *= 4.0;
            g1 = g0 + (g0 - g[j - 1]) / (fourj - 1.0);
            g[j - 1] = g0;
            g0 = g1;
        }

        /*  if (fabs(g0) > tol) */
        if (fabs(g0) > tol) // 1.e-30)
            error = 1.0 - gmax / g0;

        else
            error = gmax;
        gmax = g0;
        g[jmax] = g0;
    }

    if ((i > MAXITER) && (fabs(error) > tol))
        printf("rombint failed to converge; integral=%g, error=%g\n", g0, error);

    return g0;
}

double omegaDyn(double a)
{
  int imin = 0, imax = N_AGRID-1;
  do{
    if(a < a_InterP[(imin+imax)/2]) imax = (imin+imax)/2;
    else imin = (imin+imax)/2;
  } while( (imax-imax) > 1 );
  int i = imin;
  double val = (a_omega_interpolated[i+1] - a_omega_interpolated[i])/(a_InterP[i+1] - a_InterP[i])*(a - a_InterP[i]) + a_omega_interpolated[i];
  return val * pow(a, 4);
  }

double CSTdtauda(double a, double h0, double omegam)
{
    double tcmb = 2.7254;
    double grhog = (1.4952e-13) * pow(tcmb, 4); // 8*pi*G/(3c^2)*4*3*sigma_B/c^3 T^4           ,this is 8pi.G/(c^2) *(4sigma_B*T_present^4/c)/c^2 in Mpc^-2
    double grhom = (3.3379e-11) * h0 * h0;      // 3/c^2 = 3.3379e-11                          ,this is 3H0^2/c^2 in Mpc^-2

    double warned = 0;
    double dtauda1;
    double grho2; // Its a temporary variable for storing the total density

    grho2 = grhom * omegam * a + grhog + grhom * (1 - omegam) * omegaDyn(a);

    //if (grho2 < 0)
    //{
    //    warned = 1; // hit a warning
    //    wcount = wcount + 1;
    //    throw grho2;
    //}

    dtauda1 = sqrt(3.0 / grho2); 
    return (dtauda1);
}

double getAlpha(double OmegaM)
{
  if(OmegaM > AlphaOmOmmax[0][1])
    return AlphaOmOmmax[0][0];
  for( int i=0; i < 22; i++ )
  {
    //printf("\n-- %e %e",AlphaOmOmmax[i][0],AlphaOmOmmax[i][1]);
    if((AlphaOmOmmax[i][1] > OmegaM) && (AlphaOmOmmax[i+1][1] < OmegaM))
    {
      if(AlphaOmOmmax[i+1][1] > 9999.9)
        return AlphaOmOmmax[i][0];
      return AlphaOmOmmax[i][0] + (AlphaOmOmmax[i+1][0] - AlphaOmOmmax[i][0])/(AlphaOmOmmax[i+1][1] - AlphaOmOmmax[i][1])*( OmegaM - AlphaOmOmmax[i][1] );
      }
    }
  }

double get_lowboundMatter()
{
  double OmegaM_lowbound = 1.0;
  for( int i=0; i < 22; i++ )
    if(OmegaM_lowbound > AlphaOmOmmax[i][1])
      OmegaM_lowbound = AlphaOmOmmax[i][1];
  return OmegaM_lowbound;
}

double supernovaChi2(double h0, double omegam)
{
  double Chi2run1, Chi2run = 0;
  double a, b, c;

  for (int j = 0; j < N_SUPERNOVA; j++)
  {
    a = 1.0 /( supernovaStr[j][1] + 1.0 );
    b = 1;
    if( j > 0 ) b = 1.0 /( supernovaStr[j-1][1] + 1.0 );
    if( j > 0 ) c = supernovaStr[j-1][5];

    supernovaStr[j][5] = c + rombint(CSTdtauda, a, b, 1.0e-7, h0, omegam);
    supernovaStr[j][5] = supernovaStr[j][5] / supernovaStr[j][1];        
    supernovaStr[j][5] = log(supernovaStr[j][5])*5/log(10) - 19.35 + 25.0;
    Chi2run += (supernovaStr[j][5] - supernovaStr[j][2])*(supernovaStr[j][5] - supernovaStr[j][2])/supernovaStr[j][3]/supernovaStr[j][3];
    }
  return Chi2run;
  }

int param(int mid,double task[])
{
  double As;
  As = exp(task[5])*1.0e-10;
  char param[150],omegafile[1000]; //, basefolder[150];
  double Hub = task[2];
  double Omegam = task[0]/Hub/Hub;
  double Omegal = 1.0 - Omegam;
  double f,aomega;
  char command[400];
  int commind, size;
  
  FILE *fpaomega;
  FILE *fpwrite;

  sprintf(param,"%s/param_%d.ini",BASEFOLDER,mid);
  printf("\n%s",param);
  
  int alphaInt;
  double alphaFloat;
  double alpha = getAlpha(Omegam);

  alphaInt = (int)(alpha*1000);
  alphaFloat = alpha*1000 - alphaInt;

  fprintf(stdout, "\nFromslave(%d) par.h :Omega, alpha: %e %e max:%e %e %e %d %e (%e %e)\n", mid, task[0], Hub, get_lowboundMatter(), Omegam, alpha, alphaInt, alphaFloat, AlphaOmOmmax[alphaInt-1][1],AlphaOmOmmax[alphaInt][1]);
  fflush(stdout);
  printf("\nThis is test: %d, %s",mid,BASEFOLDER);
A133:  
  printf("\nThis is test: %d",mid);
  sprintf(omegafile,"%s/a_omega_%d.d",BASEFOLDER,mid);
  printf("%s",omegafile); //exit(1);
  fpaomega = fopen(omegafile,"w");
  for(int i = 0; i < N_AGRID; i++)
  { 
    aomega = a_omega_InterP[alphaInt-1][i] * (1 - alphaFloat) + a_omega_InterP[alphaInt][i] * alphaFloat;

    if( Omegam + (1-Omegam)*aomega*a_InterP[i]*a_InterP[i]*a_InterP[i] < 0)
      aomega = -Omegam/a_InterP[i]/a_InterP[i]/a_InterP[i]/(1-Omegam) +.001/a_InterP[i]/a_InterP[i]/(1-Omegam);

    a_omega_interpolated[i] = aomega;
    fprintf(fpaomega,"%e %e\n",a_InterP[i], aomega);
    }
  fclose(fpaomega);

  fpwrite =fopen(param,"w+"); 
  fprintf(fpwrite,"output_root = test_%d\n",mid);
  fprintf(fpwrite,"get_scalar_cls = T\n");
  fprintf(fpwrite,"get_vector_cls = F\n");
  fprintf(fpwrite,"get_tensor_cls = F\n");
  fprintf(fpwrite,"get_transfer   = F\n");
  fprintf(fpwrite,"do_lensing     = T\n");
  fprintf(fpwrite,"do_nonlinear = 0\n");
  fprintf(fpwrite,"l_max_scalar      = 2700\n");
  fprintf(fpwrite,"l_max_tensor      = 1500\n");
  fprintf(fpwrite,"k_eta_max_tensor  = 5000\n");
  fprintf(fpwrite,"use_physical   = T\n");
  fprintf(fpwrite,"ombh2          = %e\n",task[1]);
  fprintf(fpwrite,"omch2          = %e\n",Omegam*Hub*Hub-task[1]);
  fprintf(fpwrite,"omnuh2         = 0\n");
  fprintf(fpwrite,"omk            = 0\n");
  fprintf(fpwrite,"hubble         = %e\n",100*task[2]);
  fprintf(fpwrite,"#w              = -1.0\n");
  fprintf(fpwrite,"usew0wa = F\n");
  fprintf(fpwrite,"#if usew0wa = T, read in w_0, w_a \n");
  fprintf(fpwrite,"#w0             = -1.0\n");
  fprintf(fpwrite,"#wa             =  0.0\n");
  fprintf(fpwrite,"use_tabulated_w = T\n");
  fprintf(fpwrite,"#if usew0wa = F, read (a,w) from the following user-supplied file\n");
  fprintf(fpwrite,"wafile = a_omega_%d.d\n",mid);
  fprintf(fpwrite,"#constant comoving sound speed of the dark energy (1=quintessence)\n");
  fprintf(fpwrite,"cs2_lam        =  1.0\n");
  fprintf(fpwrite,"temp_cmb           = 2.726\n");
  fprintf(fpwrite,"helium_fraction    = 0.24\n");
  fprintf(fpwrite,"massless_neutrinos = 3.046\n");
  fprintf(fpwrite,"massive_neutrinos  = 0\n");
  fprintf(fpwrite,"nu_mass_eigenstates = 1\n");
  fprintf(fpwrite,"nu_mass_degeneracies = 0  \n");
  fprintf(fpwrite,"nu_mass_fractions = 1\n");
  fprintf(fpwrite,"initial_power_num         = 1\n");
  fprintf(fpwrite,"pivot_scalar              = 0.05\n");
  fprintf(fpwrite,"pivot_tensor              = 0.002\n");
  fprintf(fpwrite,"scalar_amp(1)             = %e\n",As);
  fprintf(fpwrite,"scalar_spectral_index(1)  = %e\n",task[4]);
  fprintf(fpwrite,"scalar_nrun(1)            = 0\n");
  fprintf(fpwrite,"tensor_spectral_index(1)  = 0\n");
  fprintf(fpwrite,"initial_ratio(1)          = 1\n");
  fprintf(fpwrite,"reionization         = T\n");
  fprintf(fpwrite,"re_use_optical_depth = T\n");
  fprintf(fpwrite,"re_optical_depth     = %e\n",task[3]);
  fprintf(fpwrite,"re_redshift          = 11\n");
  fprintf(fpwrite,"re_delta_redshift    = 1.5\n");
  fprintf(fpwrite,"re_ionization_frac   = -1\n");
  fprintf(fpwrite,"RECFAST_fudge = 1.14\n");
  fprintf(fpwrite,"RECFAST_fudge_He = 0.86\n");
  fprintf(fpwrite,"RECFAST_Heswitch = 6\n");
  fprintf(fpwrite,"RECFAST_Hswitch  = T\n");
  fprintf(fpwrite,"initial_condition   = 1\n");
  fprintf(fpwrite,"initial_vector = -1 0 0 0 0\n");
  fprintf(fpwrite,"vector_mode = 0\n");
  fprintf(fpwrite,"COBE_normalize = F\n");
  fprintf(fpwrite,"CMB_outputscale = 7.4311e12\n");
  fprintf(fpwrite,"transfer_high_precision = F\n");
  fprintf(fpwrite,"transfer_kmax           = 2\n");
  fprintf(fpwrite,"transfer_k_per_logint   = 0\n");
  fprintf(fpwrite,"transfer_num_redshifts  = 1\n");
  fprintf(fpwrite,"transfer_interp_matterpower = T\n");
  fprintf(fpwrite,"transfer_redshift(1)    = 0\n");
  fprintf(fpwrite,"transfer_filename(1)    = transfer_out.dat\n");
  fprintf(fpwrite,"transfer_matterpower(1) = matterpower.dat\n");
  fprintf(fpwrite,"scalar_output_file = scalCls.dat\n");
  fprintf(fpwrite,"vector_output_file = vecCls.dat\n");
  fprintf(fpwrite,"tensor_output_file = tensCls.dat\n");
  fprintf(fpwrite,"total_output_file  = totCls.dat\n");
  fprintf(fpwrite,"lensed_output_file = lensedCls.dat\n");
  fprintf(fpwrite,"lensed_total_output_file  =lensedtotCls.dat\n");
  fprintf(fpwrite,"lens_potential_output_file = lenspotentialCls.dat\n");
  fprintf(fpwrite,"FITS_filename      = scalCls.fits\n");
  fprintf(fpwrite,"do_lensing_bispectrum = F\n");
  fprintf(fpwrite,"do_primordial_bispectrum = F\n");
  fprintf(fpwrite,"bispectrum_nfields = 1\n");
  fprintf(fpwrite,"bispectrum_slice_base_L = 0\n");
  fprintf(fpwrite,"bispectrum_ndelta=3\n");
  fprintf(fpwrite,"bispectrum_delta(1)=0\n");
  fprintf(fpwrite,"bispectrum_delta(2)=2\n");
  fprintf(fpwrite,"bispectrum_delta(3)=4\n");
  fprintf(fpwrite,"bispectrum_do_fisher= F\n");
  fprintf(fpwrite,"bispectrum_fisher_noise=0\n");
  fprintf(fpwrite,"bispectrum_fisher_noise_pol=0\n");
  fprintf(fpwrite,"bispectrum_fisher_fwhm_arcmin=7\n");
  fprintf(fpwrite,"bispectrum_full_output_file=\n");
  fprintf(fpwrite,"bispectrum_full_output_sparse=F\n");
  fprintf(fpwrite,"bispectrum_export_alpha_beta=F\n");
  fprintf(fpwrite,"feedback_level = 1\n");
  fprintf(fpwrite,"lensing_method = 1\n");
  fprintf(fpwrite,"accurate_BB = F\n");
  fprintf(fpwrite,"massive_nu_approx = 1\n");
  fprintf(fpwrite,"accurate_polarization   = T\n");
  fprintf(fpwrite,"accurate_reionization   = T\n");
  fprintf(fpwrite,"do_tensor_neutrinos     = T\n");
  fprintf(fpwrite,"do_late_rad_truncation   = T\n");
  fprintf(fpwrite,"number_of_threads       = 0\n");
  fprintf(fpwrite,"high_accuracy_default=F\n");
  fprintf(fpwrite,"accuracy_boost          = 1\n");
  fprintf(fpwrite,"l_accuracy_boost        = 1\n");
  fprintf(fpwrite,"l_sample_boost          = 1\n");
  fclose(fpwrite);



  sprintf(command,"rm %s/test_%d_*",BASEFOLDER,mid);
  system(command);

  sprintf(command,"rm %s/cambdump.out",BASEFOLDER);
  system(command);

  // Sometimes a_omega files are not getting written. No idea why. Trying some fixing
  // ---------------------------------------------------------------------------------

  sprintf(omegafile,"%s/a_omega_%d.d",BASEFOLDER,mid);
  fpaomega = fopen(omegafile,"r");
  if (NULL != fpaomega) {
  fseek (fpaomega, 0, SEEK_END);
  size = ftell(fpaomega);
  fclose(fpaomega);
  if(size != 0)
  {
    commind = 3993;  
    sprintf(command,"/rds/general/user/sdas5/home/EverpresentLambda/cosmos/camb_de_lzz0/camb %s/param_%d.ini >cambdump.out",BASEFOLDER,mid);
    }

    }
  else 
  {
    printf("\nMyid (%d): Size 0, Why no idea",mid);
    goto A133;
    }
    

  commind = system(command);


  FILE *fpchk;
  int tempi = 0,tempiii;
  double tempfloat,tempfloat1;



  sprintf(command,"%s/test_%d_lensedCls.dat",BASEFOLDER,mid);
  tempiii = 0;
l2i:
  fpchk = fopen(command,"r+");
  tempiii++;
  if(fpchk == NULL)
  {
    //sleep(1);
    if(tempiii<51)
       goto l2i;
  }
  fscanf(fpchk,"%d %lf %lf %lf",&tempiii,&tempfloat,&tempfloat1,&tempfloat1);
  

  if(isnan(tempfloat))
  {
     tempi = 53;
     }
  fclose(fpchk);



  if(tempi>50)
    return 0;
  else
    return 1; 

}

#endif
