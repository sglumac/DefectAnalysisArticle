/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Cosmin Petra @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Adjoint sensitivity example problem
 *
 * This IVP is a stiff system of 6 non-linear DAEs of index 1. The 
 * problem originates from Akzo Nobel Central research in Arnhern, 
 * The Netherlands, and describes a chemical process in which 2 
 * species are mixed, while carbon dioxide is continuously added.
 * See http://pitagora.dm.uniba.it/~testset/report/chemakzo.pdf  
 * 
 * IDAS also computes the sensitivities with respect to initial
 * conditions of the following quantity:
 *   G = int_t0^t1 y1 dt
 * The sensitivity of G is the solution of the adjoint system at t0. 
 * -----------------------------------------------------------------
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <sundials/sundials_math.h>
#include <nvector/nvector_serial.h>

/* Accessor macros */
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component */

/* Problem Constants */
#define NEQ 6
#define T0 RCONST(0.0)
#define TF RCONST(180.0)

#define RTOL  RCONST(1.0e-08)
#define ATOL  RCONST(1.0e-10)
#define RTOLB RCONST(1.0e-06)
#define ATOLB RCONST(1.0e-08)
#define RTOLQ RCONST(1.0e-10)
#define ATOLQ RCONST(1.0e-12)


#define ZERO  RCONST(0.0)
#define HALF  RCONST(0.5)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define STEPS 150

typedef struct {
  realtype k1, k2, k3, k4;
  realtype K, klA, Ks, pCO2, H;
} *UserData;

static int res(realtype t, N_Vector yy, N_Vector yd, N_Vector resval, void *userdata);

static int resB(realtype tt, 
                N_Vector yy, N_Vector yp,
                N_Vector yyB, N_Vector ypB, N_Vector rrB,
                void *user_dataB);

static int rhsQ(realtype t, N_Vector yy, N_Vector yp, 
              N_Vector qdot, void *user_data);

static void PrintOutput(realtype tfinal, N_Vector yB, N_Vector ypB);
static int check_flag(void *flagvalue, const char *funcname, int opt);


/* Main program */
int main()
{
  UserData data;
  void *mem;
  N_Vector yy, yp, rr, q;
  N_Vector yB, ypB;
  int ncheck, flag;
  realtype time;
  long int nst, nstB;
  int indexB;
    

  mem = NULL;
  yy = yp = NULL;

  printf("\nAdjoint Sensitivity Example for Akzo-Nobel Chemical Kinetics\n");
  printf("-------------------------------------------------------------\n");
  printf("Sensitivity of G = int_t0^tf (y1) dt with respect to IC.\n");
  printf("-------------------------------------------------------------\n\n");
  /* Allocate user data. */
  data = (UserData) malloc(sizeof(*data));

  /* Fill user's data with the appropriate values for coefficients. */
  data->k1 = RCONST(18.7);
  data->k2 = RCONST(0.58);
  data->k3 = RCONST(0.09);
  data->k4 = RCONST(0.42);
  data->K = RCONST(34.4);
  data->klA = RCONST(3.3);
  data->Ks = RCONST(115.83);
  data->pCO2 = RCONST(0.9);
  data->H = RCONST(737.0);

  /* Allocate N-vectors. */
  yy = N_VNew_Serial(NEQ);
  if (check_flag((void *)yy, "N_VNew_Serial", 0)) return(1);
  yp = N_VNew_Serial(NEQ);
  if (check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);

  /* Consistent IC for  y, y'. */
#define y01 0.444
#define y02 0.00123
#define y03 0.00
#define y04 0.007
#define y05 0.0
  Ith(yy,1) = RCONST(y01);
  Ith(yy,2) = RCONST(y02);
  Ith(yy,3) = RCONST(y03);
  Ith(yy,4) = RCONST(y04);
  Ith(yy,5) = RCONST(y05);
  Ith(yy,6) = data->Ks * RCONST(y01) * RCONST(y04);

  /* Get y' = - res(t0, y, 0) */
  N_VConst(ZERO, yp);

  rr = N_VNew_Serial(NEQ);
  res(T0, yy, yp, rr, data);
  N_VScale(-ONE, rr, yp);
  N_VDestroy_Serial(rr);
  
 /* Create and initialize q0 for quadratures. */
  q = N_VNew_Serial(1);
  if (check_flag((void *)q, "N_VNew_Serial", 0)) return(1);
  Ith(q,1) = ZERO;

  /* Call IDACreate and IDAInit to initialize IDA memory */
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);

  flag = IDAInit(mem, res, T0, yy, yp);
  if(check_flag(&flag, "IDAInit", 1)) return(1);


  /* Set tolerances. */
  flag = IDASStolerances(mem, RTOL, ATOL);
  if(check_flag(&flag, "IDASStolerances", 1)) return(1);

  /* Attach user data. */
  flag = IDASetUserData(mem, data);
  if(check_flag(&flag, "IDASetUser", 1)) return(1);
  
  /* Attach linear solver. */
  flag = IDADense(mem, NEQ);

  /* Initialize QUADRATURE(S). */
  flag = IDAQuadInit(mem, rhsQ, q);
  if (check_flag(&flag, "IDAQuadInit", 1)) return(1);

  /* Set tolerances and error control for quadratures. */
  flag = IDAQuadSStolerances(mem, RTOLQ, ATOLQ);
  if (check_flag(&flag, "IDAQuadSStolerances", 1)) return(1);

  flag = IDASetQuadErrCon(mem, TRUE);
  if (check_flag(&flag, "IDASetQuadErrCon", 1)) return(1);
 
  /* Prepare ADJOINT. */
  flag = IDAAdjInit(mem, STEPS, IDA_HERMITE);
  if (check_flag(&flag, "IDAAdjInit", 1)) return(1);

  /* FORWARD run. */
  printf("Forward integration ... ");
  flag = IDASolveF(mem, TF, &time, yy, yp, IDA_NORMAL, &ncheck);
  if (check_flag(&flag, "IDASolveF", 1)) return(1);

  flag = IDAGetNumSteps(mem, &nst);
  if (check_flag(&flag, "IDAGetNumSteps", 1)) return(1);

  printf("done ( nst = %ld )\n",nst);

  flag = IDAGetQuad(mem, &time, q);
  if (check_flag(&flag, "IDAGetQuad", 1)) return(1);

  printf("G:          %24.16f \n",Ith(q,1));
  printf("--------------------------------------------------------\n\n");


  /* BACKWARD run */

  /* Initialize yB */
  yB = N_VNew_Serial(NEQ);
  if (check_flag((void *)yB, "N_VNew_Serial", 0)) return(1);
  N_VConst(ZERO, yB);
  
  ypB = N_VNew_Serial(NEQ);
  if (check_flag((void *)ypB, "N_VNew_Serial", 0)) return(1);
  N_VConst(ZERO, ypB);
  Ith(ypB,1) = - ONE;

  flag = IDACreateB(mem, &indexB);
  if (check_flag(&flag, "IDACreateB", 1)) return(1);

  flag = IDAInitB(mem, indexB, resB, TF, yB, ypB);
  if (check_flag(&flag, "IDAInitB", 1)) return(1);

  flag = IDASStolerancesB(mem, indexB, RTOLB, ATOLB);
  if (check_flag(&flag, "IDASStolerancesB", 1)) return(1);

  flag = IDASetUserDataB(mem, indexB, data);
  if (check_flag(&flag, "IDASetUserDataB", 1)) return(1);

  flag = IDASetMaxNumStepsB(mem, indexB, 1000);

  flag = IDADenseB(mem, indexB, NEQ);
  if (check_flag(&flag, "IDADenseB", 1)) return(1);

  printf("Backward integration ... ");

  flag = IDASolveB(mem, T0, IDA_NORMAL);
  if (check_flag(&flag, "IDASolveB", 1)) return(1);

  IDAGetNumSteps(IDAGetAdjIDABmem(mem, indexB), &nstB);
  printf("done ( nst = %ld )\n", nstB);

  flag = IDAGetB(mem, indexB, &time, yB, ypB);
  if (check_flag(&flag, "IDAGetB", 1)) return(1);

  PrintOutput(time, yB, ypB);

  IDAFree(&mem);

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(yB);
  N_VDestroy_Serial(ypB);
  N_VDestroy_Serial(q);

  return(0);
}


static int res(realtype t, N_Vector yy, N_Vector yd, N_Vector resval, void *userdata)
{
  UserData data;
  realtype k1, k2, k3, k4;
  realtype K, klA, Ks, pCO2, H;

  realtype y1, y2, y3, y4, y5, y6;
  realtype yd1, yd2, yd3, yd4, yd5;

  realtype r1, r2, r3, r4, r5, Fin;

  data = (UserData) userdata;
  k1 = data->k1;
  k2 = data->k2;
  k3 = data->k3;
  k4 = data->k4;
  K = data->K;
  klA = data->klA;
  Ks = data->Ks;
  pCO2 = data->pCO2;
  H = data->H;

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);
  y4 = Ith(yy,4);
  y5 = Ith(yy,5);
  y6 = Ith(yy,6);

  yd1 = Ith(yd,1);
  yd2 = Ith(yd,2);
  yd3 = Ith(yd,3);
  yd4 = Ith(yd,4);
  yd5 = Ith(yd,5);

  r1 = k1 * SUNRpowerI(y1,4) * SUNRsqrt(y2);
  r2 = k2 * y3 * y4;
  r3 = k2/K * y1 * y5;
  r4 = k3 * y1 * y4 * y4;
  r5 = k4 * y6 * y6 * SUNRsqrt(y2);
  Fin = klA * ( pCO2/H - y2 );

  Ith(resval,1) = yd1 + TWO*r1 - r2 + r3 + r4;
  Ith(resval,2) = yd2 + HALF*r1 + r4 + HALF*r5 - Fin;
  Ith(resval,3) = yd3 - r1 + r2 - r3;
  Ith(resval,4) = yd4 + r2 - r3 + TWO*r4;
  Ith(resval,5) = yd5 - r2 + r3 - r5;
  Ith(resval,6) = Ks*y1*y4 - y6;

  return(0);
}

/* 
 * rhsQ routine. Computes quadrature(t,y). 
 */
 
static int rhsQ(realtype t, N_Vector yy, N_Vector yp, N_Vector qdot, void *user_data)
{
  Ith(qdot,1) = Ith(yy,1);  

  return(0);
}

#define QUARTER   RCONST(0.25)
#define FOUR      RCONST(4.0)
#define EIGHT     RCONST(8.0)

/* 
 * resB routine. Residual for adjoint system. 
 */
static int resB(realtype tt, 
                N_Vector yy, N_Vector yp,
                N_Vector yyB, N_Vector ypB, N_Vector rrB,
                void *user_dataB)
{
  UserData data;

  realtype y1, y2, y3, y4, y5, y6;

  realtype yB1, yB2, yB3, yB4, yB5, yB6;
  realtype ypB1, ypB2, ypB3, ypB4, ypB5;

  realtype k1, k2, k3, k4;
  realtype K, klA, Ks, pCO2, H;

  realtype y2tohalf, y1to3, k2overK, tmp1, tmp2;

  data = (UserData) user_dataB;
  k1 = data->k1;
  k2 = data->k2;
  k3 = data->k3;
  k4 = data->k4;
  K = data->K;
  klA = data->klA;
  Ks = data->Ks;
  pCO2 = data->pCO2;
  H = data->H;

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);
  y4 = Ith(yy,4);
  y5 = Ith(yy,5);
  y6 = Ith(yy,6);

  yB1 = Ith(yyB,1);
  yB2 = Ith(yyB,2);
  yB3 = Ith(yyB,3);
  yB4 = Ith(yyB,4);
  yB5 = Ith(yyB,5);
  yB6 = Ith(yyB,6);
  
  ypB1 = Ith(ypB,1);
  ypB2 = Ith(ypB,2);
  ypB3 = Ith(ypB,3);
  ypB4 = Ith(ypB,4);
  ypB5 = Ith(ypB,5);

  y2tohalf = sqrt(y2);
  y1to3 = y1*y1*y1;
  k2overK = k2/K;

  tmp1 = k1* y1to3 * y2tohalf; tmp2 = k3*y4*y4;
  Ith(rrB,1) = 1 +  ypB1 - (EIGHT*tmp1 + k2overK*y5 + tmp2)*yB1 
    - (TWO*tmp1+tmp2)*yB2 + (FOUR*tmp1+k2overK*y5)*yB3 
    + k2overK*y5*(yB4-yB5) - TWO*tmp2*yB4 + Ks*y4*yB6;

  tmp1 = k1 * y1*y1to3 * (y2tohalf/y2); tmp2 = k4 * y6*y6 * (y2tohalf/y2);
  Ith(rrB,2) = ypB2 - tmp1*yB1 - (QUARTER*tmp1 + QUARTER*tmp2 + klA)*yB2 
    + HALF*tmp1*yB3 + HALF*tmp2*yB5;

  Ith(rrB,3) = ypB3 + k2*y4*(yB1-yB3-yB4+yB5);

  tmp1 = k3*y1*y4; tmp2 = k2*y3;
  Ith(rrB,4) = ypB4 + (tmp2-TWO*tmp1)*yB1 - TWO*tmp1*yB2 - tmp2*yB3 
    - (tmp2+FOUR*tmp1)*yB4 + tmp2*yB5 + Ks*y1*yB6;

  Ith(rrB,5) = ypB5 - k2overK*y1*(yB1-yB3-yB4+yB5);

  Ith(rrB,6) = k4*y6*y2tohalf*(2*yB5-yB2) - yB6;


  return 0;
}


/*
 * Print results after backward integration
 */
static void PrintOutput(realtype tfinal, N_Vector yB, N_Vector ypB)
{
  printf("dG/dy0: \t%12.4e\n\t\t%12.4e\n\t\t%12.4e\n\t\t%12.4e\n\t\t%12.4e\n\t\t%12.4e\n",
         Ith(yB,1), Ith(yB,2), Ith(yB,3), Ith(yB,4), Ith(yB,5), Ith(yB,6));
  printf("--------------------------------------------------------\n\n");
}

/* 
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
