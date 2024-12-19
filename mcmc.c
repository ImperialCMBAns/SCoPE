#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<mpi.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_math.h>





#include "cstparam.h"
#include "par.h"
#include "runplc.h"

//
// CMBAns header files
// ---------------------------------------------------------------------

#define NOPARAM 6			// The number of parameters in our model
#define MAX_TASKARRAY_SIZE 50		// this is the maximum task array size
#define NR_END 1                        // Numerical Recipr function  
#define FREE_ARG char*

//
// The number of chains; 1-9 chains are supported, if you want more,
// you'll have to slightly modify the code. No big problem, though. 
//----------------------------------------------------------------
#define CHAINS 3
static short int Astier = 0;
short int PARAMETERS = NOPARAM;	// The number of parameters in our model

//
// The size of the array, i.e. all possible values 
// (parameters, likelihoods, derived parameters) that we wish 
// to write to a file. Automatically set by the program   
//-------------------------------------------------------------------------------------
unsigned short int TASKARRAY_SIZE;  
unsigned short int MULTIPURPOSEPOS;	// will be set automatically 
unsigned short int TAKEPOS,PROBPOS;
unsigned short int ADAPTIVEPOS;
unsigned int LOGLIKEPOS;
int RANDPOS;
//
// Other information (derived parameters such as sigma 8 etc) can be stored from this
// position on.  MULTIPURPOSE_REQUIRED is the number of variables that you can
// store in the chain on top of the parameters and likelihoods. 
//-------------------------------------------------------------------------------------
static unsigned int MULTIPURPOSE_REQUIRED=1; 

//
//  SDSS Flag
//-------------------------------------------------------------------------------------
static short int SDSS = 0;
static short int TwodF = 0;
static short int SDSS_BAO = 0;  	// Baryon acoustic peak 
static short int SDSSLRG = 0;   	// Luminous reg galaxies
static short int LYA_MCDONALD = 0;  	// Lyman alpha code of Pat McDonald

//
//  Other Flags
//-------------------------------------------------------------------------------------
static short int WMAP7 = 1; 		// WMAP 3-year yes or no ?
static short int BOOMERANG03 = 0; 	// Boomerang from 2003 flight
static short int VSA = 0; 
static short int ACBAR = 0;  
static short int CBI = 0;
static short int Riess06 = 0;

//
//  Needed to keep track of   parameters  needed for experiments
//-------------------------------------------------------------------------------------
short int AstierParameterPos=0;
short int WMAP5ParameterPos=0;
short int SDSSLRGParameterPos=0;

//
// Start using the estimated covariance matrix for steps after BEGINCOVUPDATE points have
// been computed in a chain 
//-------------------------------------------------------------------------------------
static unsigned int BEGINCOVUPDATE=50;
static unsigned int MAXCHAINLENGTH=50000;

//
// Some constants needed for multiprocessor application
//-------------------------------------------------------------------------------------
static int GIVEMETASK=1;   // < MPI flag indicating slave()'s wish to receive a task from master()
static int TAKERESULT=2;   // < MPI flag indicating slave()'s request to master() to take result
static int TAKETASK=3;     // < MPI flag indicating master()'s request to slave to take task.

//
// For the adaptive stepsize algoritm, variable step length; these numbers are
// rather arbitrary, but seem to work well
// Please note, that after FREEZE_IN, the covariance matrix without the 
// step factor will be used. Hence, the average stay at a certain parameter point
// is approx 3 
//-------------------------------------------------------------------------------------
static double INCREASE_STEP=1.15; //1.15
static double DECREASE_STEP=0.9; //0.9
static int  HIGH_STEP_BOUND=5;
static int  LOW_STEP_BOUND=3;


// For the covariance updater; recompute covariance matrix after UPDATE_TIME steps 
const int UPDATE_TIME=5;

// Start using the estimated covariance matrix for steps after BEGINCOVUPDATE points 
// have been computed in a chain 
// const unsigned int BEGINCOVUPDATE=100;

// For the adaptive stepsize: Do not compute covariance Matrix with more than MAXPREVIOUSPOINTS*/

const unsigned int MAXPREVIOUSPOINTS=5000;


// Minimum number of points for R-statistic calculation
const  unsigned int RMIN_POINTS = 30;  // minimum number of points for R-statistic calculation

//
// only if FREEZE_IN=true: If R-statistic for all parameters less than RBREAK and
// number of points more than MIN_SIZE_FOR_FREEZE_IN then freeze-in. 
// If you want to be really conservative, you can set RBREAK=1.1
//---------------------------------------------------------------------------------------
const double RBREAK=1.2;
const unsigned int MIN_SIZE_FOR_FREEZE_IN=100;

int SLAVEPARCHAIN = 1;
extern void test_likelihood_(int *inputflag,double *cl_in_tt,double *cl_in_te,double *cl_in_ee, double *cl_in_bb);

double drfactor;






//
//  A task is just a double array.  The benefit of having
//  this structure anyhow is that we can very easily store tasks
//  in lists and all that.
//  Each double variable corresponds to either a parameter,
//  a likelihood, some multipurpose stuff or just no information
//  at all. 
//  I find it convenient to leave [0..7] for parameters
//  [8..14] is for likelihoods and after that come multipurpose 
//  stuff. If you only need 4 parameters, there will be garbage
//  in [5..7], but as you don't have to look at it, who cares, if 
//  you don't ?
//---------------------------------------------------------------------------------------

typedef struct Task {            // The weigh of this point in the chain,   
  double  f[MAX_TASKARRAY_SIZE]; // i.e. how long does it rest at the point
  int Multiplicity; 	  	 //
                       	         // The number of times a different model has been 
                                 // simulated until the step has been taken. 
                                 // This is Multiplicity minus the number of times the
                                 // step proposal crossed the boundaries of parameter Space. 
                                 // We take this number as a meassure to increase or 
  int ReallyInvestigated;        // decrease step sizes
} Task;                          //-----------------------------------------------------------


Task *free_Task()
{
  Task *tsk;
  tsk = (Task *)malloc(sizeof(Task));
  for(unsigned int i=0;i<MAX_TASKARRAY_SIZE;i++)
    tsk->f[i]=0.0;
  tsk->Multiplicity = 0;
  tsk->ReallyInvestigated=0;
  return tsk;
}

void Task_copy(Task *tska,Task *tskb)
{
  for(unsigned int i=0;i<MAX_TASKARRAY_SIZE;i++)
    tska->f[i]=tskb->f[i];
  tska->Multiplicity = tskb->Multiplicity;
  tska->ReallyInvestigated = tskb->ReallyInvestigated; 
}

//
//  Small class to compute a "rolling" average over values.
//  The constructer takes as an argument the maximum number
//  of values. By calling push(x), the value x will be stored and
//  if the number of values exceeds the limit size, the one stored
//  earliest will be erased. Mathematically speaking:
//
//  average = (\sum_0^size x_i ) / min(size,#values stored)
//  and at each push(), x_i = x_i+1 
//  
//  The average is updated at each call to average() by updating
//  the sum of all values. Each 100 calls to average(), it is computed
//  from scratch to limit elimination of significant figures
// -------------------------------------------------------------------------------------------

typedef struct RollingAverage {
  unsigned int Size; 		// the maximum number of values
  unsigned int Number;  	// the number of values
  double Sum; 			// the sum of all values
  double *y; 			// the values
  unsigned int idx;  		// current position for push()
  int PerformFullCounter;       // Every 100 calls to average(), 
                                //we re-compute the average from scratch
}RollingAverage;

void Rolling_Average_push(RollingAverage *RA, double x) // push x into the average-array
{
  RA->Sum -= RA->y[RA->idx];                            // subtract the value we will overwrite
  RA->Sum += x;  
  RA->y[RA->idx++] = x;                                 // store the new value
  if (RA->idx == RA->Size) RA->idx = 0;                 // roll-over
  RA->Number++;  
  if (RA->Number >= RA->Size) RA->Number = RA->Size;    // at most we have Size numbers in the array
}  

void RollingAverage_clear();  	                        // clear the array

double RollingAverage_average(RollingAverage *RA)        // compute the average
{
  if (RA->Number == 0)
    {
      return 0;
    }
  if (RA->PerformFullCounter++ < 10) 
    return RA->Sum/RA->Number;

  RA->PerformFullCounter = 0;                           // reset 
  double Sum = 0.0;
  for (unsigned int i = 0; i < RA->Number; i++) {       // at most to number, which might later be size
    Sum += RA->y[i];
  }
  return RA->Sum/RA->Number;  
}

RollingAverage *new_RollingAverage(int Size)
{
  RollingAverage *RA;
  RA = (RollingAverage *)malloc(sizeof(RollingAverage));
  RA->y=(double *)malloc(Size*sizeof(double));
  RA->Size = Size;
  RA->idx = 0;
  if(!RA)
    {
      printf("Error in allocating Rolling average.");
      exit(1);
    }
  return RA;	
}

typedef struct MultiGaussian
{
  double Scale;   // bool Lock;
                  // void generateAndTransform();
	          // The row number of the covariance matrix and number of random variables sought.
  unsigned int SIZE;
  double **MasterMatrix;       // Encodes the eigenvectors/transformation matrix
  double *eigenvalues;	       // The eigenvalues
  double *generatedValues;     // The random values, generated anew if throwDice() is called
  
  double *lbounds;	       // Upper and lower bounds are encoded here
  double *hbounds;

  double *randomq;
  
  double *center;	       // The mean of the gaussian to sample from		
}MultiGaussian;

MultiGaussian *new_MultiGaussian(int Size)
{
  MultiGaussian *MG;
  MG = (MultiGaussian *)malloc(sizeof(MultiGaussian));
  MG->SIZE = Size;
  MG->MasterMatrix = (double **)malloc(Size*sizeof(double *));
  for (unsigned int i=0;i<Size;i++)
    MG->MasterMatrix[i] = (double *)malloc(Size*sizeof(double)); 
  
  MG->eigenvalues = (double *)malloc(Size*sizeof(double));
  MG->generatedValues = (double *)malloc(Size*sizeof(double));
  MG->lbounds  = (double *)malloc(Size*sizeof(double));
  MG->hbounds  = (double *)malloc(Size*sizeof(double));
  MG->center   = (double *)malloc(Size*sizeof(double));
  MG->randomq   = (double *)malloc(Size*sizeof(double));
  return MG;
}

//	
// Set bounds for the random variables. This needs to be called before any other
//  function, or unpredictive behaviour will result. In most applications, you will
//  probably set this only once.
//  \param  lowBound a vector with lower Bounds for the random variables
//  \param highBound a vector with upper Bounds for the random variables
//-----------------------------------------------------------------------------------------------------  
void MultiGaussian_setBounds(MultiGaussian *MG, double lowBound[], double highBound[]) 
{
  for(unsigned int i=0; i <MG->SIZE; i++ ) 
    {
      MG->lbounds[i]=lowBound[i];
      MG->hbounds[i]=highBound[i];
    }
}

void printInfo(MultiGaussian *MG) 
{ 
  printf("************************************\n");
  printf("MultiGaussian::printInfo:\n");
  printf("Eigenvalues (we used scale = %d):\n",MG->Scale);
  for(unsigned int i=0; i< MG->SIZE; i++)
    {
      printf("%e  --> sigma:  %e  ",MG->eigenvalues[i],sqrt(MG->eigenvalues[i]));
      printf("eigenvector w/o scale: %e  ",MG->eigenvalues[i]/MG->Scale);
      printf("  --> sigma: %e\n",sqrt(MG->eigenvalues[i]/MG->Scale));
    }
  
  printf("\nEigenvectors:\n");
  for(unsigned int i=0; i < MG->SIZE; i++)
    {
    for(unsigned int j=0; j < MG->SIZE ; j++)
      {
	printf("%e  ",MG->MasterMatrix[i][j]) ;
      }
    printf("\n");
    }
  printf("***********************************\n");
}

void nrerror(char s[])
{
  printf("%s",s);
  exit(1);
}

double **convert_matrix(double *a,long nrl,long nrh,long ncl,long nch)
// allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
// declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
// and ncol=nch-ncl+1. The routine should be called with the address
// &a[0][0] as the first argument.
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;
  
  // allocate pointers to rows 
  m=(double **) malloc((unsigned int) ((nrow+NR_END)*sizeof(double *)));
  if (!m)	nrerror("allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;
  
  // set pointers to rows 
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  // return pointer to array of pointers to rows 
  return m;
}

void free_convert_matrix(double **b,long nrl,long nrh,long ncl,long nch)
{                     	
  free((FREE_ARG) (b+nrl-NR_END));  // free a matrix allocated by convert_matrix() 
}

double *vector(long nl,long nh)     // allocate a float vector with subscript range v[nl..nh] 
{
  double *v;
  
  v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

void free_vector(double *v,long nl,long nh)
{                                  // free a float vector allocated with vector()
  free((FREE_ARG) (v+nl-NR_END));
}

double *new_double(unsigned short int n)
{
  double *v;
  v = (double *)malloc(n*sizeof(double));
  return v;
}

double fabs(double take)
{
  if(take<0)
    return -take;
  return take;
}

double posRnd(double max) {
  double rand_max = 1.0/RAND_MAX;

  return  rand() * rand_max *max;
}

double ran1(float x)
{
  double rand_max = 1.0/RAND_MAX;
  return 2*(rand()*rand_max - 0.5)*x;
}

double gasdev(double mean, double std)
{
  double rsq,v1,v2;
  do{
    v1 = ran1(1.0);
    v2 = ran1(1.0);
    rsq = v1*v1 + v2*v2;
  }while(rsq>=1.0 || rsq == 0.0);

  double fac = v1*sqrt(-2.0*log(rsq)/rsq);
  return fac*std+mean;
}

void generateRandom(MultiGaussian *MG)
{
  double *y;
  y = (double *)malloc(MG->SIZE * sizeof(double));
  
  for(unsigned int i=0;i<MG->SIZE;i++)
    y[i]=gasdev(0.0,sqrt(MG->eigenvalues[i]));
  
  for(unsigned int i = 0; i<MG->SIZE;i++)
    {
      MG->randomq[i]=y[i];
      MG->generatedValues[i]=0;
      for(unsigned int j=0;j<MG->SIZE;j++)
	{
	  MG->generatedValues[i] += MG->MasterMatrix[i][j]*y[j];
	}
    }
}

int throwDice(Task chain,Task *next,MultiGaussian *MG,int scale)
{ 
  for(unsigned int i=0;i<MG->SIZE;i++)
    {
      MG->center[i]=chain.f[i];
    }

  generateRandom(MG);


  for(unsigned int i=0;i<PARAMETERS;i++)
    {
      next->f[i]= MG->generatedValues[i]*pow(drfactor,scale)+MG->center[i];
      next->f[RANDPOS+i] = MG->randomq[i]*pow(drfactor,scale);
    }

  for(unsigned int i=0;i<MG->SIZE;i++)
    {
      if(MG->generatedValues[i]*pow(drfactor,scale) + MG->center[i] < MG->lbounds[i])      return 0;
      if(MG->generatedValues[i]*pow(drfactor,scale) + MG->center[i] > MG->hbounds[i])      return 0;
    }

  if( ( next->f[0] / (next->f[2] * next->f[2])) < get_lowboundMatter() )    return 0;   // Omega_m should not be less than the lowerbound of Omega_m
  if( ( next->f[0] / (next->f[2] * next->f[2])) > 1.0 )    return 0;

  return 1;
}

//
//  Generate a set of eigenvectors and eigenvalues from the covariance matrix.
//  Thus, if you whish to sample from a new distribution with covariance matrix S,
//  you need to call this function again. This function uses the eigenvalue and
//  eigenvector routines of the GNU Scientific Library, so if you need this for
//  large matrices it is probably best to change this code and use LAPACK or
//  something similar.
//  param covarianceMatrix the covariance matrix of the distribution you want to draw from
// -------------------------------------------------------------------------------------------------------------

void generateEigenvectors(MultiGaussian *MG, double **covarianceMatrix,double scale) 
{ 
  MG->Scale = scale;
  double *data;
  data = (double *)malloc((MG->SIZE)*(MG->SIZE)*sizeof(double));        // reformat [][] into [] for creating matrix
  int k=0;
  for(unsigned int j=0; j< MG->SIZE; j++)               		// copy into tempArray
    for(unsigned int i=0; i< MG->SIZE; i++) 
      {
	data[k]=covarianceMatrix[i][j];
	k++;
      }

  gsl_matrix_view m = gsl_matrix_view_array(data,MG->SIZE,MG->SIZE);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(MG->SIZE);
  
  gsl_vector *eval = gsl_vector_alloc(MG->SIZE);
  gsl_matrix *evec = gsl_matrix_alloc(MG->SIZE,MG->SIZE);

  gsl_eigen_symmv(&m.matrix,eval,evec,w);
  
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

  for(unsigned short int i=0;i < MG->SIZE;i++)
    {
      MG->eigenvalues[i] = gsl_vector_get(eval,i);
      for(unsigned short int j = 0; j < MG->SIZE; j++)
	MG->MasterMatrix[i][j] = gsl_matrix_get(evec,i,j);
    }

  for(unsigned short int i=0;i < MG->SIZE; i++)
    MG->eigenvalues[i] *= scale;

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  free(data);
}

unsigned short int ADAPTIVE = 1;
unsigned short int FREEZE_IN = 1;

void printExtInfo(MultiGaussian *MG)
{
  printf("\n******************************************************************************\n");
  printf("EXTENDED information :\n");
  printf("Size : %d\n",MG->SIZE);
  printf("Bounds:\n");
  for(unsigned int i=0;i<MG->SIZE;i++)
    printf("%d) %e  %e  %e\n",i,MG->lbounds[i],MG->hbounds[i],MG->generatedValues[i]);
  
  printf("\n******************************************************************************\n");
  printf("\n\n");
  for(unsigned int i = 0;i<MG->SIZE;i++)
    {
      printf("%e  :  ",MG->eigenvalues[i]);
      for(unsigned int j=0;j<MG->SIZE;j++)
	printf("  %.2e",MG->MasterMatrix[i][j]);
      printf("\n");
    }
  printf("\n******************************************************************************\n");   
}

int myslave_rank(int ii)
{
  return ii*SLAVEPARCHAIN;
}

int middleman_slave_rank(int ii,int rank)
{
  return rank + ii;
}

int middleman(int rank)
{
   if(rank%SLAVEPARCHAIN == 1) return 1;
   return 0;
}

int mymiddlemann(int rank)
{
   return ((rank-1)/SLAVEPARCHAIN)*SLAVEPARCHAIN+1;
}

double normd(double a[],int size)
{
   double norm = 0.0;
   for(int ii=0;ii<size;ii++)
     norm += a[ii]*a[ii];
   return norm;
}


int isAllawedOmegaM(double OmegaM)
{
  double MaxOmega;
  for( int i=0; i < 22; i++ )
  {
    if((AlphaOmOmmax[i][1] > OmegaM) && (AlphaOmOmmax[i+1][1] < OmegaM))
    {
      if(AlphaOmOmmax[i+1][1] > 9999.0) return 0;
      MaxOmega =  AlphaOmOmmax[i][2] + (AlphaOmOmmax[i+1][2] - AlphaOmOmmax[i][2])/(AlphaOmOmmax[i+1][1] - AlphaOmOmmax[i][1])*( AlphaOmOmmax[i+1][1] - OmegaM );
      if( MaxOmega > OmegaM ) return 1;
      else return 0;
      }
    }
}

//
// The master is responsible for sending tasks to the slaves and collecting 
// the information from the slaves after the calculation. In other words, 
// this is the coordinator. The slaves don't know about each others existence. 
// The master also monitors the output. There is only one master, but an 
// arbitrary number of slaves.
//-------------------------------------------------------------------------------------

int master(double **startnum,unsigned short int Restart) 
{
  unsigned int chainSize[CHAINS];              // chain size
  unsigned int chainTotalPerformed[CHAINS];    // total number of performed models
  unsigned int chainBack[CHAINS]={0};          // int count[CHAINS];
  unsigned int checkSize=0;                    // Last time when GR statistics was calculated
  Task chain[CHAINS][MAXCHAINLENGTH];
  unsigned int covarianceUndateFlag[CHAINS]={0};

  // For dynamical stepsize
  int steps_since_update[CHAINS]={0};
  double EntireFactor[CHAINS];
  RollingAverage *roll[CHAINS];               // rolling average of EntireFactor for after burn-in
  
  double *sigma[CHAINS];
  double **covMatrix[CHAINS]; 
  for(unsigned int ii = 0;ii<CHAINS;ii++)
  {
    sigma[ii] = new_double(TASKARRAY_SIZE);
    covMatrix[ii]=(double **)malloc(PARAMETERS*sizeof(double*));
    for(unsigned int jj = 0;jj<PARAMETERS;jj++)
	    covMatrix[ii][jj]=(double *)malloc(PARAMETERS*sizeof(double));
    }
  
  
  Task *next[CHAINS];  
  double AverageMultiplicity[CHAINS];         // just for our information

  for (unsigned int i=0; i < CHAINS; i++)     //initialize to 0 
    AverageMultiplicity[i]=0;
  
  MultiGaussian *mGauss[CHAINS];              // multivariate gaussian random number generator
  
  double *lowbound;
  double *highbound;
  

  unsigned int break_at;
  unsigned int N;

  // now let us calculate the mean of each parameter
  double *y[CHAINS];
  double *dist_y;     // distribution mean
  double *B;          // variance between chains
  double *W;          // variance within chains
  double *R;          // monitoring parameter	

  double *average;    //[PARAMETERS];	
  double **cov;       //[PARAMETERS][PARAMETERS];	 

  lowbound = new_double(PARAMETERS);
  highbound = new_double(PARAMETERS);
  dist_y = new_double(PARAMETERS);
  B = new_double(PARAMETERS);
  W = new_double(PARAMETERS);
  R = new_double(PARAMETERS);
  average = new_double(PARAMETERS);


  for(unsigned int ii = 0;ii<CHAINS; ii++)
    y[ii] = new_double(PARAMETERS);

  cov = (double **)malloc(PARAMETERS*sizeof(double*));
  for(unsigned int ii = 0;ii<PARAMETERS; ii++)
    cov[ii] = new_double(PARAMETERS);
  

  //int myseed;
  //FILE *omMax;
  //omMax = fopen("/rds/general/user/sdas5/home/EverpresentLambda/SCoPE/seed_10892/omMax.d","r");
  //fscanf(omMax,"%d %lf",&myseed,&lowbound[0]);
  //fclose(omMax);

  //
  //set the (flat)  priors, i.e. parameter boundaries
  //------------------------------------------------------------------
  //lowbound[0]=0.08;  highbound[0] = 0.24;    // omega_mh2 0.05 - 0.24
  
  // lowbound[0] = get_lowboundMatter(); highbound[0] = 0.50;  // Omega_m
  
  lowbound[0] = get_lowboundMatter()*.55*.55; highbound[0] = 0.24;
  if(lowbound[0]<0.06) lowbound[0] = 0.06;

  lowbound[1]=0.010; highbound[1] = 0.03;    // omega_bh2 -0.03
  lowbound[2]=0.40;  highbound[2] = 0.75;    // alpha 0.0 - 0.03  // h  //50 - 85
  lowbound[3]=0.01;  highbound[3] = 0.2;     // optdlss 0.0 -0.3
  lowbound[4]=0.8;   highbound[4] = 1.2;     // n  -0.14
  lowbound[5]=2.5;   highbound[5] = 3.5;     // ln (10^10 A_s) - 2\tau 
//  lowbound[6]=-2.0;   highbound[6] = 0.0;      // w0
//  lowbound[7]=0.5;   highbound[7] = 1.0;     // cs   
//  lowbound[6]=-2.0;  highbound[6] = 1.0;     // w0
//  lowbound[7]=-6.0; highbound[7] = 2.0;     // wa 

  // Must check what is the RandMax defined in your stdlib.h. It must be > 5*10^4
  srand ( (unsigned)time ( NULL ) );   // Generate a random seed
                                       //  Just keep on, if you have more parameters
                                       //  this sets the starting points and stepsizes for each chain for the first step


  int intelligentstartnum = 1;
  for (unsigned int i = 0; i < CHAINS; i++) 
  {
      roll[i] = (RollingAverage *)malloc(sizeof(RollingAverage));
      roll[i] = new_RollingAverage(500);
      
      Task t;
      for (unsigned int k=0; k < TASKARRAY_SIZE; k++) 
	      t.f[k]=0;                      // for safety reasons, initialize (otherwise garbage in first line)
      
      mGauss[i]=new_MultiGaussian(PARAMETERS);
      MultiGaussian_setBounds(mGauss[i],lowbound,highbound);
      
      //
      //Generate Starting points (inside the prior)
      //---------------------------------------------------------------------------------------------

      if(!intelligentstartnum)
      {
      for(unsigned int j=0; j<PARAMETERS; j++) 
	      t.f[j] = lowbound[j] + posRnd(highbound[j]-lowbound[j]);
      }
      else 
      {
      for(unsigned int j=0; j<PARAMETERS; j++)
        t.f[j] = startnum[i][j]; 
      }
      //
      // the sigmas; these need to be chosen with care, which is why we use dynamic stepsize
      // the sigmas only matter in the beginning (or if you de-select ADAPTIVE) 
      //----------------------------------------------------------------------------------------------
      //sigma[i][0] = 0.0007;  // omega_mh2
      sigma[i][0] = 0.0005;
      sigma[i][1] = 0.00007; // omega_bh2
      sigma[i][2] = 0.003; // h
      sigma[i][3] = 0.002;  // optdlss
      sigma[i][4] = 0.0022;  // n
      sigma[i][5] = 0.008;  // ln (10^10 A_s) - 2\tau  
      // sigma[i][6] = 0.1;  // w0
      // sigma[i][7] = 0.1;  // wa      

      
      //
      // More sigmas, if you have more parameters
      //---------------------------------------------------------------------------------------------------
      
      for(unsigned int k=0; k< PARAMETERS; k++)        //fill with sigmas
        for(unsigned int j=0; j< PARAMETERS; j++)
          {
            if (j==k) 
              covMatrix[i][j][j]=sigma[i][j]*sigma[i][j];
            else 
              covMatrix[i][j][k]=0.0;
          }
      
      generateEigenvectors(mGauss[i],covMatrix[i],1.0);       // .....Just for Debugging purpose
            
      steps_since_update[i]=0; 
      EntireFactor[i]  = 1.0;
      
      //
      // we set the loglike for this starting point to be really small, this will wash away in burn in
      //---------------------------------------------------------------------------------------------------------------
      t.f[LOGLIKEPOS] = -1e100;
      t.Multiplicity = 0;
      t.ReallyInvestigated = 0;
      if(Restart == 0) 
      	chainBack[i]++;
      
      Task_copy(&chain[i][chainBack[i]-1],&t);

      chainSize[i] = 0;    
      chainTotalPerformed[i]=0;
      printf("Myrank is :%d %d",i,myslave_rank(i));

      for(unsigned int k=0; k < PARAMETERS ;k++)  //TASKARRAY_SIZE
        t.f[PROBPOS+k] = posRnd(1.0);
      t.f[PROBPOS+PARAMETERS] = posRnd(1.0);

      printf("\nMyrank :\n");
      for(int ii=0;ii<TASKARRAY_SIZE;ii++)
         printf("%e  ",t.f[ii]);

      MPI_Send(t.f, TASKARRAY_SIZE, MPI_DOUBLE,myslave_rank(i)+1,TAKETASK,MPI_COMM_WORLD);
    }  

  // 
  // The chain output will be written to these files
  //------------------------------------------------------------------------------------------------------
  FILE *data[CHAINS];
  FILE *head[CHAINS]; 	       	// for information: this is the model currently at the top of the chain
  FILE *investigated[CHAINS];	// all models, whether taken or not taken and their chi^2's


  fprintf(stdout,"\nTHis is a function \n\n");
  fflush(stdout);

  //
  // Change the string "base" below, if you like to get your results into a different 
  // (existing) directory
  //-------------------------------------------------------------------------------------------------------
  
  char montecarlo_name[]="montecarlo_chain";
  char head_name[]= "head";
  char investigated_name[]="investigated";
  char covfile_name[]="covfile";
  char inttochar[10][10];
  // printf("\n\n\n22222222222222222222dfgkjldfjgsklhjgkdjfhgksdjhfgkdjgdskfjghdskjfghlsdfkgh");
   
  //
  // If we have restarted, we need to read in and set some stuff
  //--------------------------------------------------------------------------------------------------
   if (Restart == 1) 
    {
      for (unsigned int k =0; k < CHAINS; k++) 
	sprintf(inttochar[k],"c_%d.dat",k+1);
  
       for (unsigned int k =0; k < CHAINS; k++) 
	{
	  FILE *readcov;
	  readcov = fopen(inttochar[k],"r+");
	  if(readcov == NULL)
	    {
	      printf("Covariance file corrupt or non-existent. Cannot restart. : %s\n",inttochar[k]);
	      exit(1);
	    }
	  Task t;
	  for (unsigned int j=0; j < TASKARRAY_SIZE; j++) 
	    {
	      fscanf(readcov, "%lf",&t.f[j]);
	    }
	  
	  fscanf(readcov,"%d",&t.Multiplicity);
	  fscanf(readcov,"%d",&t.ReallyInvestigated);
	  
	  for(unsigned int j=0; j < PARAMETERS; j++)
	    for (unsigned int m=0; m < PARAMETERS; m++)
	      fscanf(readcov,"%lf",&covMatrix[k][j][m]);
	  

          Task_copy(&chain[k][chainBack[k]-1],&t);

          chainSize[k] = 0;
          chainTotalPerformed[k]=0;


	  fscanf(readcov,"%d",&readcov);
	  double OptimalFactor = 2.4 / sqrt((double)PARAMETERS);
	  generateEigenvectors(mGauss[k],covMatrix[k],OptimalFactor*OptimalFactor);
          EntireFactor[k] = OptimalFactor;
	  ADAPTIVE=0;
	  FREEZE_IN=0;
	}
    }
  

  //
  // End of reading in to restart
  //---------------------------------------------------------------------------  

  char montecarlo_name1[500];
  char head_name1[500];
  char investigated_name1[500];
  char progress_name1[500];
  char gelmanRubin_name1[500];
  char covarianceMatrices_name1[500];
  for (unsigned int k=1; k<= CHAINS; k++) 
    {
      sprintf(montecarlo_name1,"%s/%s_%d.d",BASEFOLDER,montecarlo_name,k);
      sprintf(head_name1,"%s/%s_%d.d",BASEFOLDER,head_name,k);
      sprintf(investigated_name1,"%s/%s_%d.d",BASEFOLDER,investigated_name,k);
      data[k-1] = fopen(montecarlo_name1,"w"); 
      head[k-1] = fopen(head_name1,"w");
      investigated[k-1] = fopen(investigated_name1,"w");
    }  


  printf("\n******************************************************\n");
  FILE *progress,*gelmanRubin,*covarianceMatrices;


  sprintf(progress_name1,"%s/progress.txt",BASEFOLDER);
  sprintf(gelmanRubin_name1,"%s/gelmanRubin.txt",BASEFOLDER);
  sprintf(covarianceMatrices_name1,"%s/covMatrix.txt",BASEFOLDER);  
  progress = fopen(progress_name1,"w");   	      // a summary of information will be written to this file 
  gelmanRubin = fopen(gelmanRubin_name1,"w");        // Gelman-Rubin statistical output into this file
  covarianceMatrices = fopen(covarianceMatrices_name1,"w");   // covariance output for this file
  
  if(Restart == 1) 
    fprintf(progress,"\n\n::::::::::::::::::: RESTARTED FROM HERE ::::::::::::::::::::::::\n\n");

  short int work = 1;  			//as long as its value is true, the chains will run.  
  MPI_Status  status;
  double result[MAX_TASKARRAY_SIZE]; 	// should be larger than TASKARRAY_SIZE

  //
  // Now that we have initialized everything, we go into an 
  // eternal loop: We wait for a message from the slaves (either
  // GIMMETASK or TAKERESULT
  // and send out work and receive the   results.
  // The master will also calculate the Gelman-Rubin statistic,
  // which takes up most of it's size
  //-------------------------------------------------------------------------------------------------
  int mult[CHAINS][SLAVEPARCHAIN];
  int tempi = 0;
  int take = 0;
  unsigned int min_size=100*100, min_size_performed=100*100;
  int loopstop=0;
  double randq[NOPARAM],randq1[NOPARAM],randold;
  int localpeakadjust;
  do{
    printf("\nIN Master reached : #STOP3100");
    MPI_Recv(result, TASKARRAY_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE,  MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    printf("\nIN Master reached : #STOP3100a");

    if (status.MPI_TAG == GIVEMETASK) 
    {
      int i = (status.MPI_SOURCE - 1)/SLAVEPARCHAIN;               // which slave sent this?
      randold = 0.0;                                               // Just initiallize in case        
      localpeakadjust = 1;
    
      for(int ichains=0;ichains<SLAVEPARCHAIN;ichains++)
      {
        //
        // generate a new set of parameters from the transition kernel and send to slave
        //--------------------------------------------------------------------------------------------
        mult[i][ichains] = 0; 
        next[i] = free_Task();                            
        for (;;) 
        {          
	        if(throwDice(chain[i][chainBack[i]-localpeakadjust],next[i],mGauss[i],ichains))
	      {
		break;                                 // forever (if not break)
	      }                                // if succesful step, quit
        

        //fflush(stdout);

//	    chain[i][chainBack[i]-1].Multiplicity++;   // Hence we have to increase the weight of current point
            mult[i][ichains]++; 
//	    if(chain[i][chainBack[i]-1].Multiplicity > 10)      //Hence we have to increase the weight of current point
            if(mult[i][ichains] > 10)
	        printf("Stack here with Multiplicity : %d\n",chain[i][chainBack[i]-1].Multiplicity);
	    
//	    if(chain[i][chainBack[i]-1].Multiplicity > 3000)   
          if(mult[i][ichains] > 9000)
          {
            printf("\n\nSomething has gone wrong...\n Check the covariance matrix .. ");
            exit(1);
            }
          }
        
        printf("\nMYY: (%d) ",i);
        for(int iii=0;iii<6;iii++)
            printf("%e ",next[i]->f[iii]);
        printf("\n");

        // printf("After generating random values.. step a :%d\n",ichains);	
	      // just for information: this file is always one ahead 
        for(unsigned int k=0; k < PARAMETERS ;k++)  //TASKARRAY_SIZE
        {
          fprintf(head[i],"(%e %e)\t",next[i]->f[k],chain[i][chainBack[i]-localpeakadjust].f[k]);
            if(ichains == 0)
            {
              randq[k] = next[i]->f[RANDPOS+k]; 
              }  
          }

        if(ichains == 0)
        {
           randold = normd(randq,PARAMETERS); 
           next[i]->f[RANDPOS] = randold;
           next[i]->f[RANDPOS + 1] = 0.0;
         }        
        else
        {
           for(int ij = 0;ij < PARAMETERS;ij++)
             randq1[ij] = next[i]->f[PROBPOS+ij] - randq[ij];
           next[i]->f[RANDPOS] = normd(randq1,PARAMETERS);
           next[i]->f[RANDPOS + 1] = randold;
        } 

        next[i]->f[PROBPOS] = posRnd(1.0);
        next[i]->f[ADAPTIVEPOS] = (double)ADAPTIVE; 
        next[i]->f[TAKEPOS] = (double)ADAPTIVE;
	
	      fprintf(head[i],"\n");
        
        fprintf(stdout,"Sending to slave : %d %d",status.MPI_SOURCE+ichains,status.MPI_SOURCE+ichains);
        fflush(stdout);
 
        printf("I am master, trying to send to : %d",status.MPI_SOURCE+ichains);
	      MPI_Send(next[i]->f, TASKARRAY_SIZE, MPI_DOUBLE,status.MPI_SOURCE+ichains, TAKETASK, MPI_COMM_WORLD);
        }
      }
  
    if (status.MPI_TAG == TAKERESULT)                                           // the slave sent us a result
    {
	    int i = (status.MPI_SOURCE-1)/SLAVEPARCHAIN;                            // which slave sent this?
      fprintf(stdout,"\n This is a test (%d)  %d  %d test",i,(int)result[TAKEPOS],i+1);
      fflush(stdout);

      //
      // First, keep track of the model, no matter whether we take or do not take the step
      // store this in a file called "investigatedxyz.d"
      //----------------------------------------------------------------------------------------------------
      for(unsigned int k=0; k<= PARAMETERS;k++)
        fprintf(investigated[i],"%e  ",result[k]);
	

	    // we have the result, now we have to decide whether or not to take the step
      take = (int)result[TAKEPOS];

      if(take ==1)
        chainTotalPerformed[i] += (int)(result[PROBPOS]);
      else
        chainTotalPerformed[i] += SLAVEPARCHAIN; 

      // Count chain multiplicity
      for(int fooj = 0;fooj<(int)(result[PROBPOS]);fooj++)
        chain[i][chainBack[i]-1].Multiplicity += mult[i][fooj];

	    double loglike1 = chain[i][chainBack[i]-1].f[LOGLIKEPOS];        // get the Likelihood of the previous parameter set
	
      fprintf(investigated[i],"%e  %e ",loglike1,result[LOGLIKEPOS+1]);
      fprintf(investigated[i]," %d\n",i);
	
       
	    if( ADAPTIVE == 1 ) 
      {
	      Rolling_Average_push(roll[i],pow(drfactor,(int)result[PROBPOS]-1));             // no matter if take or not, note EntireFactor for average
        EntireFactor[i] = pow(drfactor,(int)result[PROBPOS]-1);
        }
	
      //
      // we know now whether or not to take the step. 
      // ----------------------------------------------------------------------------------------
      if (take == 1)
	    {	  
        // first, write the last one in the chain to data file
        for(unsigned int k=0; k< PARAMETERS;k++) 
        fprintf(data[i],"%e  ",chain[i][chainBack[i]-1].f[k]);

        fprintf(data[i],"%e  %e  ",chain[i][chainBack[i]-1].f[LOGLIKEPOS],chain[i][chainBack[i]-1].f[LOGLIKEPOS+1]);
        fprintf(data[i],"%d \n",chain[i][chainBack[i]-1].Multiplicity); // Multiplicity, how often is this point in the chain?
        
        
        // ADAPTIVE: increase stepsize, if we take steps  too often
        if(ADAPTIVE == 1 )
          if (chain[i][chainBack[i]-1].ReallyInvestigated < LOW_STEP_BOUND) 
        if (EntireFactor[i] < 3) 
        { 
          EntireFactor[i] *= INCREASE_STEP;
          printf("\nEntire Factor : %e",EntireFactor[i]);
          generateEigenvectors(mGauss[i],covMatrix[i],EntireFactor[i]*EntireFactor[i]);      
          }	 
        // end of ADAPTIVE increase
	    
	    
	    Task t;	    
	    for (unsigned int k=0; k < TASKARRAY_SIZE; k++) 
	      chain[i][chainBack[i]].f[k] = result[k];
	    chainBack[i]++;
	    steps_since_update[i]++;                     // a taken step contributes towards steps_since_update counting
	    chainSize[i]++;

            chain[i][chainBack[i]-1].Multiplicity += (int)(result[PROBPOS]);              // Enhance multiplicty of the old one again
            chain[i][chainBack[i]-1].ReallyInvestigated += (int)(result[PROBPOS]);       // we did a simulation and we didn't take the step

	  }  
	
	else 
	  {
	    chain[i][chainBack[i]-1].Multiplicity += SLAVEPARCHAIN;  	         // Enhance multiplicty of the old one again
	    chain[i][chainBack[i]-1].ReallyInvestigated += SLAVEPARCHAIN;       // we did a simulation and we didn't take the step
	  }

  //
  // We should not start covariance Update if the MCMC steps are not started yet
  // -----------------------------------------------------------------------------------
  if(chainSize[i] >= BEGINCOVUPDATE)
  {
    if(fabs(chain[i][chainBack[i]-1].f[LOGLIKEPOS] - chain[i][chainBack[i]/2].f[LOGLIKEPOS])<500.0)
      covarianceUndateFlag[i] = 1;
    else
      covarianceUndateFlag[i] = 0;
    }

  //
  //  Check if we can kill this job
  //-------------------------------------------    
  {  
    min_size = 10000; 
    min_size_performed = 100000;

    for (unsigned int ii = 0; ii < CHAINS; ii++) 
      if ( chainSize[ii] < min_size) 
        min_size = chainSize[ii];
        
    for (unsigned int ii = 0; ii < CHAINS; ii++) 
      if ( chainTotalPerformed[ii]< min_size_performed) 
        min_size_performed = chainTotalPerformed[ii];

    double min_chisquare = 1e8;;
    if((min_size>170) || (min_size_performed > 2000))
    {
      for (unsigned int ii = 0; ii < CHAINS; ii++)
        if(min_chisquare > fabs(chain[ii][chainBack[ii]-1].f[LOGLIKEPOS]))
          min_chisquare = fabs(chain[ii][chainBack[ii]-1].f[LOGLIKEPOS]);
      
      if(min_chisquare > 3000.0)
      {
        char mcmc_command[1000];
        printf("\n\nVery High Chi2:%e\n\n",min_chisquare);
        //sprintf(mcmc_command,"mv $PWD ${PWD%/*}/more_3000/");
        //system(mcmc_command);
        exit(1);
        }
      }
    }

  //
	// Adaptivly setting stepsize according to covariance matrix
  // Here, we estimate the covariance
	// -----------------------------------------------------------------------------------

	printf("\nIN Master reached : #STOP3101: %d",loopstop);

	if(ADAPTIVE == 1 && steps_since_update[i] >= UPDATE_TIME && chainSize[i] >= BEGINCOVUPDATE && covarianceUndateFlag[i] == 1) 
	  {
      //  drfactor = 0.7;
	    int TotalSize =0;
	    
	    //	
	    // Here, we determine the number of distinct points in the chain
	    // that will be taken into account to estimate the covariance matrix
	    //------------------------------------------------------------------------------------------
	    unsigned int covSize[CHAINS];  	 
	    if (chainSize[i] < 2*BEGINCOVUPDATE) 	// early on
	      covSize[i] = chainSize[i]/2;	  					 
	    else 														
	      covSize[i] = chainSize[i] - BEGINCOVUPDATE;  // all but the first few begincovupdate points
	    
	    covSize[i] = (covSize[i] < MAXPREVIOUSPOINTS)?covSize[i]:MAXPREVIOUSPOINTS;  // at most MAXPREVIOUSPOINTS

	    //
	    // Initialize
	    // --------------------------------------------------------------------------
	    for(unsigned int k=0; k < PARAMETERS; k++)
	    { 
        average[k]=0.0;
        for(unsigned int j=0; j<PARAMETERS; j++) 
          cov[k][j]=0.0;
	      }

	   for(unsigned int ii=0;ii<CHAINS;ii++)
	   {
	     for(unsigned int n=(chainSize[i]-covSize[ii]); n < chainSize[i]; n++) 
	      {
          for (unsigned int k=0; k< PARAMETERS; k++) 
            average[k]+= chain[ii][n].f[k] *chain[ii][n].Multiplicity;
          TotalSize += chain[ii][n].Multiplicity;                 // total weight of all points
          }	    
	    	}
	    
	    for(unsigned int k=0; k < PARAMETERS; k++)			//normalize
	    {
		    average[k]=average[k]/((double)TotalSize);
	      }
      
      for(unsigned int ii=0;ii<CHAINS;ii++)
	    {
	    for (unsigned int n = chainSize[ii]-covSize[ii]; n  < chainSize[ii]; n++) 
	      { 							// run over all points
        for (unsigned int k =0; k < PARAMETERS; k++) 
          for(unsigned int j=0; j<PARAMETERS;j++)
            cov[k][j] += (chain[ii][n].f[k] - average[k])*(chain[ii][n].f[j] - average[j]) * chain[ii][n].Multiplicity;
	      }  
	    }

	    for(unsigned int k=0; k < PARAMETERS; k++)                  // normalize
	    {
	      for(unsigned int j=0; j < PARAMETERS; j++) 
        {
          cov[k][j]/=(TotalSize-1.0);
          }
	      }
	    
      for(unsigned int ii=0;ii<CHAINS;ii++)
	      for(unsigned int k=0; k< PARAMETERS; k++)
	        for (unsigned int j=0; j< PARAMETERS; j++)
		        covMatrix[ii][k][j]=cov[k][j];
	    
      for(unsigned int ii=0;ii<CHAINS;ii++)
      {
        generateEigenvectors(mGauss[ii],cov,EntireFactor[ii]*EntireFactor[ii]);
        }

	    // output information, in text file as well as binary format for re-starting
	    fprintf(covarianceMatrices,"Chain: %d Step: %d Points used: %d\n",i+1,chainSize[i],covSize[i]);
	    
	    
	    
	    for(unsigned int j=0; j < PARAMETERS; j++)
	    {
        for (unsigned int k=0; k < PARAMETERS; k++)
            fprintf(covarianceMatrices,"%e  ",cov[j][k]);
        fprintf(covarianceMatrices,"\n");
	      }

	    
	    fprintf(covarianceMatrices,"\nEntire factor: %e",EntireFactor[i]);
	    fprintf(covarianceMatrices,"\n***********************************************************\n");
	    
            for(unsigned int ii=0;ii<CHAINS;ii++)
	       steps_since_update[ii]=0;
	  }
	// END OF ADAPTIVE  COVARIANCE
	
	//
	// In addition and to help the adaptive covariance,
	// we use an adaptive step size multiplicator, called EntireFactor
	// Whenever we take a step, we consider increasing the step size,
	// as frequent approval corresponds to small stepsizes. 
	// This is performed several lines above when storing the new step
	// in the chain. Here, we consider the opposite. If we hang around
	// for some time at the same spot, we decrease the step size..
	//-----------------------------------------------------------------------------------------------
	
	// ADAPTIVE: decrease if we didn't take for a long time
	if (ADAPTIVE == 1) 
	  if (chain[i][chainBack[i]-1].ReallyInvestigated > HIGH_STEP_BOUND) 
	    if (EntireFactor[i] > 1e-1) 
	      { 
		EntireFactor[i] *=  DECREASE_STEP;
		generateEigenvectors(mGauss[i],covMatrix[i],EntireFactor[i]*EntireFactor[i]);
	      }	
	// END ADATPIVE decrease
	
  printf("\nIN Master reached : #STOP3102: %d",loopstop);
	//
	// now, we re-evaluate 
	// first, let us find out how long the shortest chain is 
	//-----------------------------------------------------------------------------------------------
	min_size=1000*1000;
	for (unsigned int ii = 0; ii < CHAINS; ii++) 
	  if ( chainSize[ii] < min_size) 
	    min_size = chainSize[ii];
	
	// output of progress information
	fprintf(progress,"performed/Chainsize: ");
	for (unsigned int ii=0; ii<CHAINS; ii++)
	  fprintf(progress,"%d/%d  ",chainTotalPerformed[ii],chainSize[ii]);  
	
	fprintf(progress,"min_size: %d    Multiplicity [Really]: ",min_size);
	for (unsigned int ii=0; ii< CHAINS; ii++)
	  fprintf(progress,"%d[%d]  ",chain[ii][chainBack[ii]-1].Multiplicity,chain[ii][chainBack[ii]-1].ReallyInvestigated);
  
  fprintf(progress,"Adaptive: ");
  	for (unsigned int ii=0; ii< CHAINS; ii++)
	  fprintf(progress,"%d  ", covarianceUndateFlag[ii]);
	fprintf(progress,"\n");
	
	//
	// Start with the Gelman and Rubin(1992) statistic;
	// compare variances within the chain with the variances between the chains
	// R[k] should be smaller than 1.1 for convergence
	// If you would like to run a synthetic distribution for checks, this statistics
	// will slow you down. You can speed things up if you replace the if
	// statement below by something like:
	//        if ( (min_size) > RMIN_POINTS && Miscmath::posRnd(1.0) > 0.95) { 
	// hence, only about every 20 times, the statistics is re-calculated.
	// ---------------------------------------------------------------------------------------------

	if ((min_size > RMIN_POINTS) && (posRnd(1.0)>0.95))
	  {
	    break_at = min_size/2; 
	    N = min_size - break_at;

	    //
	    // Initialize 
	    // ------------------------------------------------------------
	    for (unsigned int k =  0; k < PARAMETERS; k++) 
	      {     
		dist_y[k]= 0; 
		B[k] = 0;
		W[k]=0;
		for (unsigned int ii=0; ii < CHAINS; ii++) 
		  y[ii][k]=0;
	      }
    
	    double M = (double)CHAINS;
	    for (unsigned int ii = 0; ii < CHAINS; ii++) 
	      {                                                           // for all chains
		int TotalMultiplicity = 0; 
		for (unsigned int n = break_at; n < min_size; n++) 
		  {                                                       // traverse through the chains
		    for (unsigned int k = 0; k < PARAMETERS; k++)  	  // for all parameter	      
		      y[ii][k] += chain[ii][n].f[k] *chain[ii][n].Multiplicity; 	  // add up all parameters 

		    TotalMultiplicity += chain[ii][n].Multiplicity;
		  }
		
		N = TotalMultiplicity;
		for (unsigned int k =0; k < PARAMETERS; k++)
		  {
		    y[ii][k] /= (double) N;
		    dist_y[k] += y[ii][k] / M;  // M chains hence mean is sum / CHAINS
		  }	
	      }	
		    
	    // variance between chains
	    for (unsigned int k=0; k < PARAMETERS; k++) 
	      for (unsigned int ii = 0; ii < CHAINS; ii++) 
		{
		  double t1 = y[ii][k] - dist_y[k];
		  B[k] +=t1*t1 / (M-1.0);           	// as defined in (22) (why 1/(M-1) and not M ???) 	  
		}
	    // variance within a chain

	    // once again, go back to break_at in the chains	    
	    // get W, the variance within chains
	    // run over all points
	    for (unsigned int ii = 0; ii < CHAINS; ii++) 	      // run over all chains
	      for (unsigned int n =  break_at; n  < min_size; n++)    // run over all points
		{ 
		  for (unsigned int k =0; k < PARAMETERS; k++) 
		    {   
		      double t2 = chain[ii][n].f[k] - y[ii][k];
		      W[k] += t2*t2*chain[ii][n].Multiplicity;   

		    }
		}

	    
	    // normalize W and get R
	    for (unsigned int k =0; k < PARAMETERS; k++) 
	      {
		W[k] /= (M*(N-1.0));
		R[k] = (N-1.0)/N * W[k]  + B[k]*(1.0 + 1.0/N);
		R[k] /= W[k];
	      }
	    
	    
	    // if we want FREEZE_IN, test if convergence has been reached:
	    if(FREEZE_IN == 1 && min_size > MIN_SIZE_FOR_FREEZE_IN )
	      {
		unsigned short int ConvergenceReached = 1;
		for(unsigned int m=0; m<PARAMETERS; m++) 
		  if(R[m] > RBREAK) 
		    ConvergenceReached=0;
		
		if(ConvergenceReached == 1) 
		  {
		    ADAPTIVE = 0;    //stop with adaptive stepsize
		    FREEZE_IN = 0;   //never check again (otherwise no real freeze-in)
                    drfactor = 1.0; 
		    FILE *final;
		    final = fopen("freezeInEigenvectors.txt","w");
		    
		    

		    for (unsigned int k=1; k<= CHAINS; k++)
		      {
		        fclose(data[k-1]); 
                        sprintf(montecarlo_name1,"%s_freeze_%d.d",montecarlo_name,k);
		        data[k-1] = fopen(montecarlo_name1,"w"); // full name is montecarlo_chain_1.d
		      }

		    for (unsigned int j=0; j < CHAINS; j++) 
		      {
			Task *keep;
			keep = chain[j] + chainBack[j]-1;
			// the next 2 lines write the flag to each montecarlo file 
			// indicating that all models before this have to be discarded
			// the flag is a zero weight 
			

			// Push back the last one...
			Task_copy(chain[j],keep);
			//	chainSize[j]=0;
			// USE Dunkley et. al. result for optimal sigma_t: 2.4 / sqrt(Parameters)
			double OptimalFactor = 2.4 / sqrt((double)PARAMETERS);                     // ........... Check it 
			OptimalFactor = RollingAverage_average(roll[j]);
			EntireFactor[j]=OptimalFactor;   
			if(EntireFactor[j]>1.0)
			  EntireFactor[j] = 0.2;
			fprintf(final,"CHAIN[%d] Eingenvector and value info: ",j);


			fprintf(final,"************************************\n");
			fprintf(final,"MultiGaussian::printInfo:\n");
			fprintf(final,"Eigenvalues (we used scale = %d):\n",mGauss[j]->Scale);
			for(unsigned int ii=0; ii< mGauss[j]->SIZE; ii++)
			  {
			    fprintf(final,"%e  --> sigma:  %e  ",mGauss[j]->eigenvalues[ii],sqrt(mGauss[j]->eigenvalues[ii]));
			    fprintf(final,"eigenvector w/o scale: %e  ",mGauss[j]->eigenvalues[ii]/mGauss[j]->Scale);
			    fprintf(final,"  --> sigma: %e\n",sqrt(mGauss[j]->eigenvalues[ii]/mGauss[j]->Scale));
			  }
			
			fprintf(final,"\nEigenvectors:\n");
			for(unsigned int ii=0; ii < mGauss[j]->SIZE; ii++)
			  {
			    for(unsigned int jj=0; jj < mGauss[j]->SIZE ; jj++)
			      {
				fprintf(final,"%e  ",mGauss[j]->MasterMatrix[ii][jj]) ;
			      }
			    fprintf(final,"\n");
			  }
			fprintf(final,"***********************************\n");

			printInfo(mGauss[j]);
			fprintf(final,"\n\n");		
		      }
		    
		    fclose(final);	
		    fprintf(covarianceMatrices,"Stopped adaptive stepsize at %d for all chains (FREEZE_IN=true)",min_size); 
		  }
	      }

	    
	    fprintf(progress,"Statistics: \n");
	    
	    for (unsigned int k = 0; k < PARAMETERS; k++) 
	      fprintf(progress,"%d  dist_y: %e   B: %e   W: %e     R[k]:  %e \n",k,dist_y[k],B[k],W[k],R[k]); 
	    
	    if (FREEZE_IN) 
	      fprintf(progress,"Multiplicity (and EntireFactor): "); 
	    else 
	      fprintf(progress,"Multiplicity (and frozen EntireFactor): ");
	    
	    for (unsigned int ii = 0; ii < CHAINS; ii++) 
	      fprintf(progress,"%d (%e)  ",chain[ii][chainBack[ii]-1].Multiplicity,EntireFactor[ii]);
	    fprintf(progress,"\n");
	    

	    if((min_size)!=checkSize)
	      {
		fprintf(gelmanRubin,"%d   ",break_at);
		for (unsigned int n=0 ; n < PARAMETERS; n++)
		  {
		    fprintf(gelmanRubin,"%e   ",R[n]);
		  }
		fprintf(gelmanRubin,"%d \n",min_size);
		checkSize=min_size;
	      }
	    
	  } // end if min_size > old_size
	
	if(FREEZE_IN == 0)
	{
          loopstop++;
	  printf("\n**** %d",loopstop);
	}
	unsigned int max_size=0;
	for(int ii=0;ii<CHAINS;ii++)
	  max_size = (max_size>chainBack[ii])?max_size:chainBack[ii];
	
  printf("\nIN Master reached : #STOP3103 : %d %d %d %d\n",FREEZE_IN,loopstop,max_size,MAXCHAINLENGTH);

	if((FREEZE_IN == 0 && loopstop>18000)||max_size > MAXCHAINLENGTH)
	  work = 0;
	

	fflush(data[i]);
	fflush(head[i]);
	fflush(investigated[i]);
	fflush(progress);
	fflush(gelmanRubin);
	fflush(covarianceMatrices);

      } // if TAKE_RESULT

    ++tempi;
  }  while(work);//++tempi<100); //exit(1);// end of DO it should be while(work);//

  fclose(progress);
  fclose(gelmanRubin);
  fclose(covarianceMatrices);

  for (unsigned int k=1; k<= CHAINS; k++) 
    {
      fclose(data[k-1]);
      fclose(head[k-1]);
      fclose(investigated[k-1]);
    }

  printf("Hi, I am out.So finally exiting. No mistake in the middle. Everything works perfectly. No need to worry.");
  exit(1);
}

int fact(int a)
{
  if(a<=1)return 1;
  return a*fact(a-1);
}

double mind(double x, double y)
{
  return (x<y)?x:y;
}

//
// The slaves are responsible for computing  models and its associated likelihoods. They recive the
// parameters from the Master, report back to him the results of the computation and request new jobs
//--------------------------------------------------------------------------------------------------------
int slave(int rank) 
{

  printf("\nI am from slave");
  MPI_Status  status;
  char inputfilename[30],command[200];
  double task[MAX_TASKARRAY_SIZE];
  double Chi2X;
  double ratio;
  FILE *fplens;
  int lfile;                    //We compute it if cmb required to get sigma_8 in multipurpose
  double loglike = -1.0e30;
  double Cl_TT[1200];
  double Cl_TE[1200];
  double Cl_EE[1200];
  double Cl_BB[1200];
  unsigned short int takenindex = 0,take =0;
  int magicnum1,magicnum2;
  
  double alpha12,alpha32,alpha13,l2,q1,newss;
  //
  // start an infinite loop;
  // basically, a lot of this is just setting the parameters we have received from the master
  // and running cmbans, obtaining the likelihoods and derived parameters (e.g. sigma8 ) 
  // and sending the result to the master
  //-----------------------------------------------------------------------------------------------------
  char cfplens[100];
  int cnt = 0;


  for (;;) 
    {  
  
     if(middleman(rank) == 1  || ( SLAVEPARCHAIN == 1 && middleman(rank) == 0)) 
	    MPI_Send(task,TASKARRAY_SIZE,MPI_DOUBLE,0, GIVEMETASK, MPI_COMM_WORLD);
      MPI_Recv(task,TASKARRAY_SIZE,MPI_DOUBLE,0,TAKETASK,MPI_COMM_WORLD,&status); //master's answer   
      
      takenindex = 0;
      //
      // Specify the parameters in the CMBAns
      //----------------------------------------------------------------------------------
      magicnum1 = param(rank,task);      
      
      if(magicnum1 !=0)
	    {
	      sprintf(cfplens,"%s/test_%d_lensedCls.dat",BASEFOLDER,rank);
        Chi2X = run_plc(cfplens,rank);

//	  sprintf(inputfilename,"input_%d.inp",rank);
//	  fplens = fopen(inputfilename,"w+");	
//	  fprintf(fplens,"test_%d_lensedCls.dat\noutput_%d.out",rank,rank);
//	  fclose(fplens);
//	  sprintf(command,"rm output_%d.out",rank);
//	  system(command);              
//	  sprintf(command,"/data1/student/csantanud/wmap_likelihood_v5/test<%s >wmapdump_%d.d",inputfilename,rank);
//	  system(command);
	  
//	  sprintf(inputfilename,"output_%d.out",rank);
//	  magicnum2 = 0;
//	Levfps:
//	  fplens = fopen(inputfilename,"r+");
//	  if(fplens == NULL)
//	    {
//	      sleep(1);
//	      magicnum2++;
//	      if(magicnum2 <201)
//		goto Levfps;
//	    }
//	  if(magicnum2<201)
//	    {
//	      fscanf(fplens,"%lf",&Chi2X);
//	      fclose(fplens);
//	    }
//	  else
//	    Chi2X = 10.0e10;  
	}
      else
	{
	  Chi2X = 10.0e10; 
	}
      Chi2X = -Chi2X;//-Cl_BB[0];
      
      //
      // Having computed the model, we can compare it with the data
      // using the rich number of functions available in the analyzeit-class
      //-----------------------------------------------------------------------------------------
      
      if (WMAP7 == 1) 
	    task[LOGLIKEPOS+1] = Chi2X;
         
      //
      // Now, we add up all loglike's we would like to include.
      // master() will only compare the value of task[LOGLIKEPOS]
      // for determining the likelihood. 
      // Hence, you may compute as many likelihoods as you like, and
      // store them in task[]. For as long as you don't add a log likelihood
      // to task[LOGLIKEPOS], this will not influence the run of the MCMC
      //-------------------------------------------------------------------------------

      task[LOGLIKEPOS] = 0.0;
      if (WMAP7 == 1) 			task[LOGLIKEPOS] += task[LOGLIKEPOS+1]; 
      if (CBI == 1)  			task[LOGLIKEPOS] += task[LOGLIKEPOS+10]; 
      if (VSA == 1)  			task[LOGLIKEPOS] += task[LOGLIKEPOS+11];  
      if (ACBAR == 1)  			task[LOGLIKEPOS] += task[LOGLIKEPOS+12]; 
      if (BOOMERANG03 == 1)	        task[LOGLIKEPOS] += task[LOGLIKEPOS+13]; 
      if (TwodF == 1)  			task[LOGLIKEPOS] += task[LOGLIKEPOS+14]; 
      if (SDSS == 1)  			task[LOGLIKEPOS] += task[LOGLIKEPOS+15];
      if (Riess06 == 1)  		task[LOGLIKEPOS] += task[LOGLIKEPOS+16];  
      if (Astier == 1)  		task[LOGLIKEPOS] += task[LOGLIKEPOS+17];  
      if (SDSS_BAO == 1)  		task[LOGLIKEPOS] += task[LOGLIKEPOS+18];  
      if (SDSSLRG == 1)  		task[LOGLIKEPOS] += task[LOGLIKEPOS+19];  
      if (LYA_MCDONALD == 1)	        task[LOGLIKEPOS] += task[LOGLIKEPOS+20]; 
      
      //
      // Send the whole task - line back to master()+ last minute sanity check     
      // ------------------------------------------------------------------------------
      for (unsigned int k = 0;k < TASKARRAY_SIZE; k++)  
      {
        if (isnan(task[k])) 
          {
            printf("isnan element");
            exit(1);
          }
      }

      printf("\nMy middleman Rank is: %d", middleman(rank));

      if(middleman(rank) != 1 && SLAVEPARCHAIN > 1)
	      MPI_Send(task, TASKARRAY_SIZE, MPI_DOUBLE,mymiddlemann(rank), TAKERESULT,MPI_COMM_WORLD);
      

      take = 1;    
      takenindex = 0;
      task[TAKEPOS] = 0;
     
      printf("\nMy middleman Rank is: %d %d", middleman(rank), rank);
      
      if((middleman(rank) ==1) || ( SLAVEPARCHAIN == 1 && middleman(rank) == 0))
	    {
        for(int ii = rank;ii<rank+SLAVEPARCHAIN;ii++)
          {
            take = 1; 
            if(ii!=rank  && SLAVEPARCHAIN > 1)
            {
              printf("I am waiting to receive something: %d %d %d",ii,rank, SLAVEPARCHAIN);
              MPI_Recv(task,TASKARRAY_SIZE,MPI_DOUBLE,ii,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
              printf("I received all the values: %d %d %d",ii,rank, SLAVEPARCHAIN);
              }
            if(takenindex == 0)
            {
              if(task[LOGLIKEPOS] <= loglike)
                {
                  if(ii == rank || task[ADAPTIVEPOS]== 0) 
                  {
                    ratio = 1.0/exp((loglike - task[LOGLIKEPOS])/2.0);
                    newss = task[LOGLIKEPOS]; 
                    }
                  else
                  {
                    alpha12 = mind(1,exp(0.5*(newss-loglike)));
                    alpha32 = mind(1,exp(0.5*(newss-task[LOGLIKEPOS])));
                    l2 = exp(0.5*(task[LOGLIKEPOS]-loglike));
                    // Test // In our case we have already one negative sign
                    q1 = exp(-0.5*(task[RANDPOS]-task[RANDPOS+1]));
                    alpha13 = l2*q1*(1-alpha32)/(1-alpha12);
                    ratio = alpha13;
                    if(alpha13 != alpha13)
                      ratio = 0.0;
                    }
                  if(task[PROBPOS]>ratio)
                    take=0;
                  } 
              printf("\nI am trying to send taken index");
              if(take == 1)
                {
                  takenindex = 1; 
                  loglike = task[LOGLIKEPOS];
                  task[PROBPOS] = 1.0*(ii-rank+1);
                  task[TAKEPOS] = 1;
                  MPI_Send(task, TASKARRAY_SIZE, MPI_DOUBLE,0, TAKERESULT,MPI_COMM_WORLD); 
                } 
              } 
	          }
	  
        if(takenindex == 0)    
        {  
          task[TAKEPOS] = 0.0; 
          MPI_Send(task, TASKARRAY_SIZE, MPI_DOUBLE,0, TAKERESULT,MPI_COMM_WORLD);      
          }
        
      }
    } 
}	

int slave0(int rank,int runperchain) 
{
  MPI_Status  status;
  char inputfilename[30],command[200];
  double *task;
  double Chi2X;
  double ratio;
  FILE *fplens;
  int lfile;                    //We compute it if cmb required to get sigma_8 in multipurpose
  double loglike = -1.0e30;
  double Cl_TT[1200];
  double Cl_TE[1200];
  double Cl_EE[1200];
  double Cl_BB[1200];
  unsigned short int takenindex = 0,take =0;
  int magicnum1,magicnum2;
  char cfplens[100];
  double alpha12,alpha32,alpha13,l2,q1,newss;

  task = new_double(PARAMETERS+1);
  //
  // start an infinite loop;
  // basically, a lot of this is just setting the parameters we have received from the master
  // and running cmbans, obtaining the likelihoods and derived parameters (e.g. sigma8 ) 
  // and sending the result to the master
  //-----------------------------------------------------------------------------------------------------
  printf("My rank is here %d",rank);
  int cnt = 0;
  for (int i=0;i<runperchain;i++) 
    {  
      MPI_Recv(task,PARAMETERS+1,MPI_DOUBLE,0,200,MPI_COMM_WORLD,&status); //master's answer   
         
      //
      // Specify the parameters in the CMBAns
      //----------------------------------------------------------------------------------
      magicnum1 = param(rank,task);      
      printf("\nReceived (abcdefgh) : %d %d",rank,magicnum1);

      if(magicnum1 !=0)
      {
        sprintf(cfplens,"%s/test_%d_lensedCls.dat",BASEFOLDER, rank);
	      Chi2X = run_plc(cfplens,rank);
//	sprintf(inputfilename,"input_%d.inp",rank);
//	fplens = fopen(inputfilename,"w+");	
//	fprintf(fplens,"test_%d_lensedCls.dat\noutput_%d.out",rank,rank);
//	fclose(fplens);
//	sprintf(command,"rm output_%d.out",rank);
//	system(command);              
//        printf("\nReceived (a1b1c1d1e1f1g1h) : %d %d",rank,magicnum1);
//        sprintf(command,"rm wmapdump.d",rank);
//        system(command);

//        printf("\nReceived (a2b2c2d2e2f2g2h) : %d %d",rank,magicnum1);

//	sprintf(command,"/data1/student/csantanud/wmap_likelihood_v5/test<%s >wmapdump_%d.d",inputfilename,rank);
//	system(command);

//        sprintf(command,"rm wmapdump.d",rank);
//        system(command);

//        printf("\nChisquare : %e",Chi2X);
	
//	sprintf(inputfilename,"output_%d.out",rank);
//	magicnum2 = 0;
  //    Levfps:
//	fplens = fopen(inputfilename,"r+");
//	if(fplens == NULL)
//	  {
//	    sleep(1);
//	    magicnum2++;
//	    if(magicnum2 <20)
//	      goto Levfps;
//	  }
//	if(magicnum2<20)
////	  {
//            fscanf(fplens,"%lf",&Chi2X);
  //          fclose(fplens);
//	  }
//	else
//	  Chi2X = 10.0e10;  
      }
      else
	    {
	      Chi2X = 10.0e10; 
	      }
  
      //printf("\nChisquare : %e",Chi2X);    

      //
      // Having computed the model, we can compare it with the data
      // using the rich number of functions available in the analyzeit-class
      //-----------------------------------------------------------------------------------------
      
      if (WMAP7 == 1) 
	    task[PARAMETERS] = Chi2X;
         
      printf("\nChisquare is (%d): %e",rank,Chi2X);
      //
      // Send the whole task - line back to master()+ last minute sanity check     
      // ------------------------------------------------------------------------------
      for (unsigned int k = 0;k < PARAMETERS+1; k++)  
      {
        if (isnan(task[k])) 
          {
            printf("isnan element");
            exit(1);
          }
      }

      //printf("\n\nThis is :\n\n");
      //fprintf(stdout,"\nI will now send : %e",task[PARAMETERS]); 
      fflush(stdout);
      MPI_Send(task,PARAMETERS+1, MPI_DOUBLE,0,201,MPI_COMM_WORLD); 
    }

}

double **master0(int worldsize,int runperchain)
{

  MPI_Status  status;
   double *startval[100];
   for(int i=0;i<100;i++)
     startval[i] = new_double(PARAMETERS+1);
   
  double *lowbound,*highbound;

  lowbound = new_double(PARAMETERS);
  highbound = new_double(PARAMETERS);


  //printf("\nI am master printing\n");
  //int myseed;
  //FILE *omMax;
  //omMax = fopen("/rds/general/user/sdas5/home/EverpresentLambda/SCoPE/seed_10892/omMax.d","r");
  //fscanf(omMax,"%d %lf",&myseed,&lowbound[0]);
  //fclose(omMax);
  //printf("\nI am master printing: %e\n",lowbound[0]);


  //
  //set the (flat)  priors, i.e. parameter boundaries
  //------------------------------------------------------------------
  lowbound[0] = 0.10; highbound[0] = 0.18;
  if(get_lowboundMatter()*.69*.69>lowbound[0]) lowbound[0] = get_lowboundMatter()*.69*.69;
  if( highbound[0] < lowbound[0])  highbound[0] = lowbound[0] +.02;
  // THIS is FOR Omega_m
  //lowbound[0] = (get_lowboundMatter()>0.30) ? get_lowboundMatter():0.30;  
  //highbound[0] = 0.35;                       // Omega_m 
  //if(lowbound[0]>highbound[0]) highbound[0] = lowbound[0] + 0.02;
  
  lowbound[1]=0.021; highbound[1] = 0.026;    // omega_bh2 -0.03
  lowbound[2]=0.50;  highbound[2] = 0.69;     // h  //50 - 85
  lowbound[3]=0.09;  highbound[3] = 0.12;     // optdlss 0.0 -0.3
  lowbound[4]=0.85;   highbound[4] = 1.0;     // n  -0.14
  lowbound[5]=2.9;   highbound[5] = 3.2;     // ln (10^10 A_s) - 2\tau 


  for(int i=0;i<100;i++)
  {
     for(int j=0;j<=PARAMETERS;j++)
     {
       startval[i][j] = lowbound[j] + posRnd(highbound[j]-lowbound[j]);
       printf("%e   ",startval[i][j]);  
     }
     printf("\n");  
   }

  double *startvaltemp;
  startvaltemp = new_double(PARAMETERS+1);
  int slavenumber = worldsize/runperchain;


  for(int i=0;i<runperchain;i++)
  {
    for(int j=0;j<worldsize/runperchain;j++)
        MPI_Send(startval[i+j*runperchain],PARAMETERS+1, MPI_DOUBLE,j+1,200,MPI_COMM_WORLD);
    for(int j=0;j<worldsize/runperchain;j++)
    {  
      MPI_Recv(startvaltemp,PARAMETERS+1,MPI_DOUBLE,MPI_ANY_SOURCE,201,MPI_COMM_WORLD,&status); 
      for(int k=0;k<=PARAMETERS;k++)
      {
        startval[i+(status.MPI_SOURCE-1)*runperchain][k] = startvaltemp[k];
        }
      }
    }     

  int temparrange[100];
  double temploglike[100];

  for(int i=0;i<100;i++)
  {
     temparrange[i] = i;
     temploglike[i]=startval[i][PARAMETERS];
     }
   
  // Arrange

  int temparrangeflip;
  double temploglikeflip;   

  for(int i=0;i<CHAINS;i++)
    for(int j=i;j<worldsize;j++)
    {
       if(temploglike[i]>temploglike[j])
       {
          temploglikeflip = temploglike[i];
          temparrangeflip = temparrange[i];

          temploglike[i] = temploglike[j];
          temparrange[i] = temparrange[j];
  
          temploglike[j] = temploglikeflip;
          temparrange[j] = temparrangeflip;          
       }
    }

   double **startvalfinal;
    startvalfinal=(double **)malloc(CHAINS*sizeof(double*));

   for(int i=0;i<CHAINS;i++)
     startvalfinal[i] = new_double(PARAMETERS+1);

  for(int i=0;i<CHAINS;i++)
  {
     printf("\n");
     for(int j=0;j<=PARAMETERS;j++)
     { 
         startvalfinal[i][j]=startval[temparrange[i]][j]; 
         printf("%e  ",startvalfinal[i][j]);
      }   
  }
 
  
  return startvalfinal;
  }

void setVariables() 
{
  if (Astier == 1) 
    {
      AstierParameterPos = PARAMETERS;
      PARAMETERS += 2; 			       // two more parameters for Astier
    }
  
  if (SDSSLRG == 1) 
    {
      SDSSLRGParameterPos = PARAMETERS;
      PARAMETERS +=1; 			       // bias
    }
  
  LOGLIKEPOS = PARAMETERS + 1;         	       // leave at least one blank collumn
  // at least at position 15. That means that most of the time, the likelihoods will be at exactly the same
  // position in the file. That's more convenient
  LOGLIKEPOS = (LOGLIKEPOS > 15) ? LOGLIKEPOS : 15;
  MULTIPURPOSEPOS = LOGLIKEPOS + 22;			// leave one blank collumn between

  ADAPTIVEPOS = MULTIPURPOSEPOS + 1;
  TAKEPOS = ADAPTIVEPOS + 1;
  PROBPOS = TAKEPOS + 1;
  RANDPOS = PROBPOS +1;
  TASKARRAY_SIZE = PROBPOS + PARAMETERS + MULTIPURPOSE_REQUIRED + 1;

  if (TASKARRAY_SIZE > MAX_TASKARRAY_SIZE) 
    {
      printf("Error 64 : Maximum array size crossed. Exiting with exit code 1");
      exit(1);
    }
}

int main(int argc, char* argv[]) 
{
  // printf("\nIts working here before");
  MPI_Init(&argc, &argv);
  int myrank;
  FILE *fpsupernova;

  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  unsigned short int Restart = 0;
  if (argc == 2) 
    {
      char *arg;
      arg = argv[1];
      if (arg == "-restart") 
	Restart = 1; 
      else 
	{
	  printf("There is some problem. The code will stop.");
	  exit(0);
	}
    }
  setVariables();


  int worldsize,totalrun,runperchain;
  drfactor = 0.7;
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);  

  //
  // This part is added for CST analysis
  // All the processors must do this part
  //----------------------------------------------------

  char filetoopen[1000];
  FILE *fp_filetoopen;
  int tempint_fileread;
  double tempdouble_fileread;
  printf("\nMyrank : %d",myrank);

  if(myrank == 0)
  {
    scanf("%s",&BASEFOLDER);
    scanf("%d",&SEEDNUMBER);
  }
  MPI_Bcast(BASEFOLDER, 1000, MPI_CHAR, 0,  MPI_COMM_WORLD);
  MPI_Bcast(&SEEDNUMBER, 1, MPI_INT, 0,  MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);


  printf("\nMy Basefolder is: %s\n",BASEFOLDER);

  
  for(int loopi=1;loopi<22;loopi++)
  {
    sprintf(filetoopen,"%s/OmegaM_seed_%d_%03d.d",BASEFOLDER,SEEDNUMBER,10*loopi);
    //printf("\nKKKKK: %d %s",myrank,filetoopen);
    if((fp_filetoopen = fopen(filetoopen,"r")) != NULL)
    {
      AlphaOmOmmax[loopi-1][0] = loopi*0.001;
      fscanf(fp_filetoopen,"%d %lf %lf %lf",&tempint_fileread,&tempdouble_fileread,&AlphaOmOmmax[loopi-1][1],&tempdouble_fileread);
      fclose(fp_filetoopen);
      //printf("\nYYYY:%d %d %e %e %e\n",myrank, tempint_fileread,tempdouble_fileread,AlphaOmOmmax[loopi-1][1],tempdouble_fileread);

      sprintf(filetoopen,"%s/omMax_seed_%d_%03d.d",BASEFOLDER,SEEDNUMBER,10*loopi);      
      fp_filetoopen = fopen(filetoopen,"r");
      fscanf(fp_filetoopen,"%d %lf",&tempint_fileread,&AlphaOmOmmax[loopi-1][2]);
      fclose(fp_filetoopen);      
      //printf("\nXXXX:%d %d %e %e %e\n",myrank, tempint_fileread,tempdouble_fileread,AlphaOmOmmax[loopi-1][1],tempdouble_fileread);

      sprintf(filetoopen,"%s/a_omega_seed_%d_%03d.d",BASEFOLDER,SEEDNUMBER,10*loopi);      
      fp_filetoopen = fopen(filetoopen,"r");
      for(int loopj = 0; loopj < N_AGRID; loopj++)
        fscanf(fp_filetoopen,"%lf %lf",&a_InterP[loopj],&a_omega_InterP[loopi-1][loopj]);
      fclose(fp_filetoopen); 
      
      //printf("\nZZZZ:%d %d %e %e %e\n",myrank, tempint_fileread,tempdouble_fileread,AlphaOmOmmax[loopi-1][1],tempdouble_fileread);
      }
    else{
      AlphaOmOmmax[loopi-1][0] = loopi*0.001;
      AlphaOmOmmax[loopi-1][1] = 9999999.99;  // SET Infility
      AlphaOmOmmax[loopi-1][2] = 9999999.99;  // SET Infinity
      
      for(int loopj = 0; loopj < N_AGRID; loopj++)
        a_omega_InterP[loopi-1][loopj] = 9999999.99;   // SET Infinity
      }
    }
  

    fpsupernova = fopen("/rds/general/user/sdas5/home/EverpresentLambda/SCoPE/expansion_history/supernova.csv","r");
    for(int i=0;i<N_SUPERNOVA;i++){
        fscanf(fpsupernova,"%lf %lf %lf %lf", &supernovaStr[i][0], &supernovaStr[i][1], &supernovaStr[i][2], &supernovaStr[i][3]);
        supernovaStr[i][1] = 1/(1+supernovaStr[i][0]);
    }
    fclose(fpsupernova);

    //for(int loopi=1;loopi<22;loopi++)
    //  printf("\n%d %e %e %e",myrank,AlphaOmOmmax[loopi-1][0],AlphaOmOmmax[loopi-1][1],AlphaOmOmmax[loopi-1][2]);
    
    //printf("\nMyrank : %d (%d)",myrank,SEEDNUMBER);

  if(worldsize>101)
  {
    totalrun = 100;
    runperchain = 1; 
   }
  else
  {
    runperchain = 100/(worldsize-1);
    totalrun = runperchain*(worldsize-1);
   }

  double **startnum;

  if(myrank == 0) 
    startnum = master0(totalrun,runperchain);
  else 
    slave0(myrank,runperchain);

  printf("Finally out of prechain...");

  if(myrank == 0) 
    master(startnum,Restart);
  else 
    slave(myrank);

  MPI_Finalize();

  return 0;
}
