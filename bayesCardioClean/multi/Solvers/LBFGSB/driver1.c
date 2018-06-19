/*  L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License" */
/*  or "3-clause license") */
/*  Please read attached file License.txt */
#include <stdio.h>
#include <stdlib.h>
#include "lbfgsb.h"
#include "../../Amazigh/mathMatrix.h"
#include "../../Amazigh/mathVector.h"
#include "../../Amazigh/supportC.h"
#include "../../Gaussian/gaussian.h"
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798

/*
                            DRIVER 1 in Fortran 77 
    -------------------------------------------------------------- 
               SIMPLE DRIVER FOR L-BFGS-B (version 3.0) 
    -------------------------------------------------------------- 

       L-BFGS-B is a code for solving large nonlinear optimization 
            problems with simple bounds on the variables. 

       The code can also be used for unconstrained problems and is 
       as efficient for these problems as the earlier limited memory 
                         code L-BFGS. 

       This is the simplest driver in the package. It uses all the 
                   default settings of the code. 


    References: 

       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited 
       memory algorithm for bound constrained optimization'', 
       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. 

       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN 
       Subroutines for Large Scale Bound Constrained Optimization'' 
       Tech. Report, NAM-11, EECS Department, Northwestern University, 
       1994. 


         (Postscript files of these papers are available via anonymous 
          ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) 

                             *  *  * 

        March 2011   (latest revision) 
        Optimization Center at Northwestern University 
        Instituto Tecnologico Autonomo de Mexico 

        Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778: 
        L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained 
        Optimization"  (2011). To appear in  ACM Transactions on 
        Mathematical Software, 
    -------------------------------------------------------------- 
            DESCRIPTION OF THE VARIABLES IN L-BFGS-B 
    -------------------------------------------------------------- 

    n is an INTEGER variable that must be set by the user to the 
      number of variables.  It is not altered by the routine. 

    m is an INTEGER variable that must be set by the user to the 
      number of corrections used in the limited memory matrix. 
      It is not altered by the routine.  Values of m < 3  are 
      not recommended, and large values of m can result in excessive 
      computing time. The range  3 <= m <= 20 is recommended. 

    x is a DOUBLE PRECISION array of length n.  On initial entry 
      it must be set by the user to the values of the initial 
      estimate of the solution vector.  Upon successful exit, it 
      contains the values of the variables at the best point 
      found (usually an approximate solution). 

    l is a DOUBLE PRECISION array of length n that must be set by 
      the user to the values of the lower bounds on the variables. If 
      the i-th variable has no lower bound, l(i) need not be defined. 

    u is a DOUBLE PRECISION array of length n that must be set by 
      the user to the values of the upper bounds on the variables. If 
      the i-th variable has no upper bound, u(i) need not be defined. 

    nbd is an INTEGER array of dimension n that must be set by the 
      user to the type of bounds imposed on the variables: 
      nbd(i)=0 if x(i) is unbounded, 
             1 if x(i) has only a lower bound, 
             2 if x(i) has both lower and upper bounds, 
             3 if x(i) has only an upper bound. 

    f is a DOUBLE PRECISION variable.  If the routine setulb returns 
      with task(1:2)= 'FG', then f must be set by the user to 
      contain the value of the function at the point x. 

    g is a DOUBLE PRECISION array of length n.  If the routine setulb 
      returns with taskb(1:2)= 'FG', then g must be set by the user to 
      contain the components of the gradient at the point x. 

    factr is a DOUBLE PRECISION variable that must be set by the user. 
      It is a tolerance in the termination test for the algorithm. 
      The iteration will stop when 

       (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch 

      where epsmch is the machine precision which is automatically 
      generated by the code. Typical values for factr on a computer 
      with 15 digits of accuracy in double precision are: 
      factr=1.d+12 for low accuracy; 
            1.d+7  for moderate accuracy; 
            1.d+1  for extremely high accuracy. 
      The user can suppress this termination test by setting factr=0. 

    pgtol is a double precision variable. 
      On entry pgtol >= 0 is specified by the user.  The iteration 
        will stop when 

                max{|proj g_i | i = 1, ..., n} <= pgtol 

        where pg_i is the ith component of the projected gradient. 
      The user can suppress this termination test by setting pgtol=0. 

    wa is a DOUBLE PRECISION  array of length 
      (2mmax + 5)nmax + 11mmax^2 + 8mmax used as workspace. 
      This array must not be altered by the user. 

    iwa is an INTEGER  array of length 3nmax used as 
      workspace. This array must not be altered by the user. 

    task is a CHARACTER string of length 60. 
      On first entry, it must be set to 'START'. 
      On a return with task(1:2)='FG', the user must evaluate the 
        function f and gradient g at the returned value of x. 
      On a return with task(1:5)='NEW_X', an iteration of the 
        algorithm has concluded, and f and g contain f(x) and g(x) 
        respectively.  The user can decide whether to continue or stop 
        the iteration. 
      When 
        task(1:4)='CONV', the termination test in L-BFGS-B has been 
          satisfied; 
        task(1:4)='ABNO', the routine has terminated abnormally 
          without being able to satisfy the termination conditions, 
          x contains the best approximation found, 
          f and g contain f(x) and g(x) respectively; 
        task(1:5)='ERROR', the routine has detected an error in the 
          input parameters; 
      On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task 
        contains additional information that the user can print. 
      This array should not be altered unless the user wants to 
         stop the run for some reason.  See driver2 or driver3 
         for a detailed explanation on how to stop the run 
         by assigning task(1:4)='STOP' in the driver. 

    iprint is an INTEGER variable that must be set by the user. 
      It controls the frequency and type of output generated: 
       iprint<0    no output is generated; 
       iprint=0    print only one line at the last iteration; 
       0<iprint<99 print also f and |proj g| every iprint iterations; 
       iprint=99   print details of every iteration except n-vectors; 
       iprint=100  print also the changes of active set and final x; 
       iprint>100  print details of every iteration including x and g; 
      When iprint > 0, the file iterate.dat will be created to 
                       summarize the iteration. 

    csave  is a CHARACTER working array of length 60. 

    lsave is a LOGICAL working array of dimension 4. 
      On exit with task = 'NEW_X', the following information is 
        available: 
      lsave(1) = .true.  the initial x did not satisfy the bounds; 
      lsave(2) = .true.  the problem contains bounds; 
      lsave(3) = .true.  each variable has upper and lower bounds. 

    isave is an INTEGER working array of dimension 44. 
      On exit with task = 'NEW_X', it contains information that 
      the user may want to access: 
        isave(30) = the current iteration number; 
        isave(34) = the total number of function and gradient 
                        evaluations; 
        isave(36) = the number of function value or gradient 
                                 evaluations in the current iteration; 
        isave(38) = the number of free variables in the current 
                        iteration; 
        isave(39) = the number of active constraints at the current 
                        iteration; 

        see the subroutine setulb.f for a description of other 
        information contained in isave 

    dsave is a DOUBLE PRECISION working array of dimension 29. 
      On exit with task = 'NEW_X', it contains information that 
        the user may want to access: 
        dsave(2) = the value of f at the previous iteration; 
        dsave(5) = the machine precision epsmch generated by the code; 
        dsave(13) = the infinity norm of the projected gradient; 

        see the subroutine setulb.f for a description of other 
        information contained in dsave 

    -------------------------------------------------------------- 
          END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B 
    -------------------------------------------------------------- 
    */


/*
	This simple driver demonstrates how to call the L-BFGS-B code to
	solve a sample problem (the extended Rosenbrock function
	subject to bounds on the variables). The dimension n of this
	problem is variable.
	nmax is the dimension of the largest problem to be solved.
	mmax is the maximum number of limited memory corrections.
	Declare the variables needed by the code.
	A description of all these variables is given at the end of
	the driver.
	Declare a few additional variables for this sample problem.
	*/

double* code()
{


    /* Local variables */
    static double f, g[1024];
    static integer i__;
    static double l[1024];
    static integer m, n;
    static double u[1024], x[1024], t1, t2, wa[43251];
    static integer nbd[1024], iwa[3072];
/*     static char task[60]; */
    static integer taskValue;
    static integer *task=&taskValue; /* must initialize !! */
/*      http://stackoverflow.com/a/11278093/269192 */
    static double factr;
/*     static char csave[60]; */
    static integer csaveValue;
    static integer *csave=&csaveValue;
    static double dsave[29];
    static integer isave[44];
    static logical lsave[4];
    static double pgtol;
    static integer iprint;


/*     We wish to have output at every iteration. */
    iprint = 1; 
/*     iprint = 101; */
/*     We specify the tolerances in the stopping criteria. */
    factr = 1e7;
    pgtol = 1e-5;
/*     We specify the dimension n of the sample problem and the number */
/*        m of limited memory corrections stored.  (n and m should not */
/*        exceed the limits nmax and mmax respectively.) */
    n = 100;
    m = 5;
/*     We now provide nbd which defines the bounds on the variables: */
/*                    l   specifies the lower bounds, */
/*                    u   specifies the upper bounds. */
/*     First set bounds on the odd-numbered variables. */
	for (int i = 0; i < n; i++)
	{
		l[i] = -100;
		u[i] = 100;
		x[i] = 50;
	}
   

    /*     We start the iteration by initializing task. */

    *task = (integer)START;
/*     s_copy(task, "START", (ftnlen)60, (ftnlen)5); */
    /*        ------- the beginning of the loop ---------- */
L111:
    /*     This is the call to the L-BFGS-B code. */
    setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
            iprint, csave, lsave, isave, dsave);
/*     if (s_cmp(task, "FG", (ftnlen)2, (ftnlen)2) == 0) { */
    if ( IS_FG(*task) ) {
        /*        the minimization routine has returned to request the */
        /*        function f and gradient g values at the current x. */
        /*        Compute function value f for the sample problem. */
        /* Computing 2nd power */
		f = 0;
		for (int i = 0; i < n; i++)
		{
			f += x[i]* x[i];
		}
		for (int i = 0; i < n; i++)
		{
			g[i] = 2*x[i];
		}
     
        /*          go back to the minimization routine. */
        goto L111;
    }

/*     if (s_cmp(task, "NEW_X", (ftnlen)5, (ftnlen)5) == 0) { */
    if ( *task==NEW_X ) {
        goto L111;
    }
    
	for (int i = 0; i < n; i++)
		printf("%0.4f \t", x[i]);

    return x;
} /* MAIN__ */


double* solveKexpB(double** xtrainLBFGSA, double* ytrainLBFGSA, int n_pointsA, int DimensionA)
{
	int Dimension_support = DimensionA;
	int n_points_support = n_pointsA;
	double** xtrainLBFGS_support = newDouble(n_pointsA, DimensionA);
	copyMatrix(xtrainLBFGS_support, xtrainLBFGSA, n_pointsA, DimensionA);

	double* ytrainLBFGS_support = newDoubleV(n_pointsA);
	copyVector(ytrainLBFGS_support, ytrainLBFGSA, n_pointsA);
	int N = DimensionA + 1;
	double* result = newDoubleV(N);
		
	/* Local variables */
	static double f, g[1024];
	static integer i__;
	static double l[1024];
	static integer m, n;
	static double u[1024], x[1024], t1, t2, wa[43251];
	static integer nbd[1024], iwa[3072];
	/*     static char task[60]; */
	static integer taskValue;
	static integer *task = &taskValue; /* must initialize !! */
									   /*      http://stackoverflow.com/a/11278093/269192 */
	static double factr;
	/*     static char csave[60]; */
	static integer csaveValue;
	static integer *csave = &csaveValue;
	static double dsave[29];
	static integer isave[44];
	static logical lsave[4];
	static double pgtol;
	static integer iprint;


	/*     We wish to have output at every iteration. */
	iprint = 1;
	/*     iprint = 101; */
	/*     We specify the tolerances in the stopping criteria. */
	factr = 1e7;
	pgtol = 1e-5;
	/*     We specify the dimension n of the sample problem and the number */
	/*        m of limited memory corrections stored.  (n and m should not */
	/*        exceed the limits nmax and mmax respectively.) */
	n = N;
	m = 5;
	/*     We now provide nbd which defines the bounds on the variables: */
	/*                    l   specifies the lower bounds, */
	/*                    u   specifies the upper bounds. */
	/*     First set bounds on the odd-numbered variables. */
	for (int i = 0; i < N; i++)
	{
		nbd[i] = 2;
		l[i] = 0;
		u[i] = 1.0;
		x[i] = 0.0;
	}



	/*     We start the iteration by initializing task. */

	*task = (integer)START;
	/*     s_copy(task, "START", (ftnlen)60, (ftnlen)5); */
	/*        ------- the beginning of the loop ---------- */
L111:
	/*     This is the call to the L-BFGS-B code. */
	setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
		iprint, csave, lsave, isave, dsave);
	/*     if (s_cmp(task, "FG", (ftnlen)2, (ftnlen)2) == 0) { */
	if (IS_FG(*task)) {
		/*        the minimization routine has returned to request the */
		/*        function f and gradient g values at the current x. */
		/*        Compute function value f for the sample problem. */
		/* Computing 2nd power */
		double** K = KernelMV(xtrainLBFGS_support, xtrainLBFGS_support,
			n_points_support, n_points_support, Dimension_support, x, 1);

		double** id = idMatrix(n_points_support);
		double** noise = multConstant(id, 1e-12, n_points_support, n_points_support);
		double** A = addMatrix(K, noise, n_points_support, n_points_support);
		double** Ainvx = inverseMatrix(A, n_points_support);
		double* Ainvy = multVector(Ainvx, ytrainLBFGS_support, n_points_support, n_points_support);


		double yt_Ainvy = dotVector(ytrainLBFGS_support, Ainvy, n_points_support);

		double logDeterminant = logDet(A, n_points_support);
		double marginal_likelihood = -0.5*yt_Ainvy - 0.5*(logDeterminant)-0.5*log(2 * PI);
		f = -(marginal_likelihood);


		double** pK_Dimension_support = newDouble(n_points_support, n_points_support);

		double* xI = newDoubleV(Dimension_support);
		double* xJ = newDoubleV(Dimension_support);

		int indexI, indexJ;
		int k;
		for (indexI = 0; indexI < n_points_support; indexI++)
		{
			for (indexJ = 0; indexJ < n_points_support; indexJ++)
			{
				for (int i = 0; i < Dimension_support; i++)
				{
					xI[i] = xtrainLBFGS_support[indexI][i];
					xJ[i] = xtrainLBFGS_support[indexJ][i];
				}
				double* temp = subVector(xI, xJ, Dimension_support);
				double theta = x[Dimension_support];
				double coef = 0;
				for (k = 0; k < Dimension_support; k++)
				{
					coef += (temp[k] * temp[k]) / (2 * exp(x[k] * x[k]));
				}
				pK_Dimension_support[indexI][indexJ] = -(2 * theta) / (exp(theta *theta)*exp(coef));
				free(temp);
			}
		}

		double** tempTrace = multMatrix(Ainvx, pK_Dimension_support, n_points_support, n_points_support, n_points_support, n_points_support);
		double trace = traceMatrix(tempTrace, n_points_support);

		double* pdim = multVector(pK_Dimension_support, Ainvy, n_points_support, n_points_support);
		double* AinvPdim = multVector(Ainvx, pdim, n_points_support, n_points_support);
		double yAinvPdim = dotVector(ytrainLBFGS_support, AinvPdim, n_points_support);

		g[Dimension_support] = -(0.5*yAinvPdim - 0.5*trace);

		free(pdim);
		free(AinvPdim);
		freeMatrix(pK_Dimension_support, n_points_support);
		freeMatrix(tempTrace, n_points_support);
		free(xI);
		free(xJ);

		int indexTheta;
		for (indexTheta = 0; indexTheta < Dimension_support; indexTheta++)
		{
			double** pK_Theta = newDouble(n_points_support, n_points_support);

			double* xITheta = newDoubleV(Dimension_support);
			double* xJTheta = newDoubleV(Dimension_support);

			for (indexI = 0; indexI < n_points_support; indexI++)
			{
				for (indexJ = 0; indexJ < n_points_support; indexJ++)
				{
					for (int i = 0; i < Dimension_support; i++)
					{
						xITheta[i] = xtrainLBFGS_support[indexI][i];
						xJTheta[i] = xtrainLBFGS_support[indexJ][i];
					}
					double* temp = subVector(xITheta, xJTheta, Dimension_support);
					double theta = x[Dimension_support];
					double coef = 0;
					for (k = 0; k < Dimension_support; k++)
					{
						coef += (temp[k] * temp[k]) / (2 * exp(x[k] * x[k]));
					}
					pK_Theta[indexI][indexJ] = (temp[indexTheta] * temp[indexTheta] * x[indexTheta]) /
						(exp(theta *theta)*exp(coef)* exp(x[indexTheta] * x[indexTheta]));
					free(temp);
				}
			}

			double** tempTraceTheta = multMatrix(Ainvx, pK_Theta, n_points_support, n_points_support, n_points_support, n_points_support);
			double traceTheta = traceMatrix(tempTraceTheta, n_points_support);
			double* pTheta = multVector(pK_Theta, Ainvy, n_points_support, n_points_support);
			double* AinvPTheta = multVector(Ainvx, pTheta, n_points_support, n_points_support);
			double yAinvPTheta = dotVector(ytrainLBFGS_support, AinvPTheta, n_points_support);

			g[indexTheta] = -(0.5*yAinvPTheta - 0.5*traceTheta);

			freeMatrix(tempTraceTheta, n_points_support);
			freeMatrix(pK_Theta, n_points_support);
			free(xITheta);
			free(xJTheta);
			free(pTheta);
			free(AinvPTheta);
		}




		for (int j = 0; j < n_points_support; j++)
		{
			free(K[j]);
			free(id[j]);
			free(noise[j]);
			free(A[j]);
			free(Ainvx[j]);

		}

		free(K);
		free(id);
		free(noise);
		free(A);
		free(Ainvx);
		free(Ainvy);

		/*          go back to the minimization routine. */
		goto L111;
	}

	/*     if (s_cmp(task, "NEW_X", (ftnlen)5, (ftnlen)5) == 0) { */
	if (*task == NEW_X) {
		goto L111;
	}

	for (int i = 0; i < n; i++)
		printf("%0.4f \t", x[i]);

	freeMatrix(xtrainLBFGS_support, n_pointsA);
	free(ytrainLBFGS_support);

	for (int i = 0; i < N; i++)
	{
		result[i] = x[i];
	}
	return result;
}

double* solveUCBB(double** Ainv, double* ytrainLBFGSA, int n_pointsA,
	int DimensionA, double* params, double** xtrainLBFGSA,
	double* randomPoint, double kappa, int type)
{
	int Dimension_support = DimensionA;
	int n_points_support = n_pointsA;
	double **xtrainLBFGS_support = newDouble(n_pointsA, DimensionA);
	copyMatrix(xtrainLBFGS_support, xtrainLBFGSA, n_pointsA, DimensionA);

	double* ytrainLBFGS_support = newDoubleV(n_pointsA);
	copyVector(ytrainLBFGS_support, ytrainLBFGSA, n_pointsA);

	double** Ainv_support = newDouble(n_pointsA, n_pointsA);
	copyMatrix(Ainv_support, Ainv, n_pointsA, n_pointsA);

	double* param_support = newDoubleV(DimensionA + 1);
	copyVector(param_support, params, DimensionA + 1);

	double** tespoint_support = newDouble(1, DimensionA);
	int type_support = type;
	double kappa_support = kappa;

	int N = DimensionA;

	/* Local variables */
	static double f, g[1024];
	static integer i__;
	static double l[1024];
	static integer m, n;
	static double u[1024], x[1024], t1, t2, wa[43251];
	static integer nbd[1024], iwa[3072];
	/*     static char task[60]; */
	static integer taskValue;
	static integer *task = &taskValue; /* must initialize !! */
									   /*      http://stackoverflow.com/a/11278093/269192 */
	static double factr;
	/*     static char csave[60]; */
	static integer csaveValue;
	static integer *csave = &csaveValue;
	static double dsave[29];
	static integer isave[44];
	static logical lsave[4];
	static double pgtol;
	static integer iprint;


	/*     We wish to have output at every iteration. */
	iprint = 1;
	/*     iprint = 101; */
	/*     We specify the tolerances in the stopping criteria. */
	factr = 1e7;
	pgtol = 1e-5;
	/*     We specify the dimension n of the sample problem and the number */
	/*        m of limited memory corrections stored.  (n and m should not */
	/*        exceed the limits nmax and mmax respectively.) */
	n = N;
	m = 5;
	/*     We now provide nbd which defines the bounds on the variables: */
	/*                    l   specifies the lower bounds, */
	/*                    u   specifies the upper bounds. */
	for (int i = 0; i < N; i++)
	{
		nbd[i] = 2;
		l[i] = 0.0;
		u[i] = 1.0;
		x[i] = randomPoint[i];
	}



	/*     We start the iteration by initializing task. */

	*task = (integer)START;
	/*     s_copy(task, "START", (ftnlen)60, (ftnlen)5); */
	/*        ------- the beginning of the loop ---------- */
L111:
	/*     This is the call to the L-BFGS-B code. */
	setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
		iprint, csave, lsave, isave, dsave);
	/*     if (s_cmp(task, "FG", (ftnlen)2, (ftnlen)2) == 0) { */
	if (IS_FG(*task)) {
		/*        the minimization routine has returned to request the */
		/*        function f and gradient g values at the current x. */
		/*        Compute function value f for the sample problem. */
		/* Computing 2nd power */

		int i;
		
		
		for (i = 0; i < Dimension_support; i++)
		{
			tespoint_support[0][i] = x[i];
		}

		double* wi = multVector(Ainv_support, ytrainLBFGS_support, n_points_support, n_points_support);

		f = acqFunctionUCB(Ainv_support, wi, n_points_support, Dimension_support, param_support,
			xtrainLBFGS_support, tespoint_support, kappa_support, type_support);
		free(wi);
		//printf("fx  \t %0.8f \n", f);
		for (i = 0; i < Dimension_support; i++)
		{
			g[i] = getDerUCB(Ainv_support, ytrainLBFGS_support, n_points_support, Dimension_support,
				param_support, xtrainLBFGS_support, tespoint_support, type_support, i,
				kappa_support);
			//printf("g \t %d \t %0.8f \n", i, g[i]);
		}

		

		/*          go back to the minimization routine. */
		goto L111;
	}

	/*     if (s_cmp(task, "NEW_X", (ftnlen)5, (ftnlen)5) == 0) { */
	if (*task == NEW_X) {
		goto L111;
	}

	for (int i = 0; i < n; i++)
		printf("%0.4f \t", x[i]);

	double* result = newDoubleV(N);

	for (int i = 0; i < N; i++)
	{
		result[i] = x[i];
	}


	freeMatrix(xtrainLBFGS_support, n_pointsA);
	free(ytrainLBFGS_support);
	freeMatrix(Ainv_support, n_pointsA);
	free(param_support);
	freeMatrix(tespoint_support, 1);
	return result;
}

double* solveUCBBNumeric(double** Ainv, double* ytrainLBFGSA, int n_pointsA,
	int DimensionA, double* params, double** xtrainLBFGSA,
	double* randomPoint, double kappa, int type)
{
	int Dimension_support = DimensionA;
	int n_points_support = n_pointsA;
	double **xtrainLBFGS_support = newDouble(n_pointsA, DimensionA);
	copyMatrix(xtrainLBFGS_support, xtrainLBFGSA, n_pointsA, DimensionA);

	double* ytrainLBFGS_support = newDoubleV(n_pointsA);
	copyVector(ytrainLBFGS_support, ytrainLBFGSA, n_pointsA);

	double** Ainv_support = newDouble(n_pointsA, n_pointsA);
	copyMatrix(Ainv_support, Ainv, n_pointsA, n_pointsA);

	double* param_support = newDoubleV(DimensionA + 1);
	copyVector(param_support, params, DimensionA + 1);

	double** tespoint_support = newDouble(1, DimensionA);
	int type_support = type;
	double kappa_support = kappa;

	int N = DimensionA;

	/* Local variables */
	static double f, g[1024];
	static integer i__;
	static double l[1024];
	static integer m, n;
	static double u[1024], x[1024], t1, t2, wa[43251];
	static integer nbd[1024], iwa[3072];
	/*     static char task[60]; */
	static integer taskValue;
	static integer *task = &taskValue; /* must initialize !! */
									   /*      http://stackoverflow.com/a/11278093/269192 */
	static double factr;
	/*     static char csave[60]; */
	static integer csaveValue;
	static integer *csave = &csaveValue;
	static double dsave[29];
	static integer isave[44];
	static logical lsave[4];
	static double pgtol;
	static integer iprint;


	/*     We wish to have output at every iteration. */
	iprint = 1;
	/*     iprint = 101; */
	/*     We specify the tolerances in the stopping criteria. */
	factr = 1e7;
	pgtol = 1e-5;
	/*     We specify the dimension n of the sample problem and the number */
	/*        m of limited memory corrections stored.  (n and m should not */
	/*        exceed the limits nmax and mmax respectively.) */
	n = N;
	m = 5;
	/*     We now provide nbd which defines the bounds on the variables: */
	/*                    l   specifies the lower bounds, */
	/*                    u   specifies the upper bounds. */
	/*     First set bounds on the odd-numbered variables. */
	for (int i = 0; i < N; i++)
	{
		nbd[i] = 2;
		l[i] = 0;
		u[i] = 1.0;
		x[i] = randomPoint[i];
	}



	/*     We start the iteration by initializing task. */

	*task = (integer)START;
	/*     s_copy(task, "START", (ftnlen)60, (ftnlen)5); */
	/*        ------- the beginning of the loop ---------- */
L111:
	/*     This is the call to the L-BFGS-B code. */
	setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
		iprint, csave, lsave, isave, dsave);
	/*     if (s_cmp(task, "FG", (ftnlen)2, (ftnlen)2) == 0) { */
	if (IS_FG(*task)) {
		/*        the minimization routine has returned to request the */
		/*        function f and gradient g values at the current x. */
		/*        Compute function value f for the sample problem. */
		/* Computing 2nd power */
		double h = 1e-4;
		int i, j;


		for (i = 0; i < Dimension_support; i++)
		{
			tespoint_support[0][i] = x[i];
		}

		double* wi = multVector(Ainv_support, ytrainLBFGS_support, n_points_support, n_points_support);

		f = acqFunctionUCB(Ainv_support, wi, n_points_support, Dimension_support, param_support,
			xtrainLBFGS_support, tespoint_support, kappa_support, type_support);
	

		for (i = 0; i < Dimension_support; i++)
		{
			double** tempointm1 = newDouble(1, Dimension_support);
			double** tempointp1 = newDouble(1, Dimension_support);
			for (j = 0; j < Dimension_support; j++)
			{
				if (i == j)
				{
					tempointm1[0][j] = tespoint_support[0][j] - h;
					tempointp1[0][j] = tespoint_support[0][j] + h;


				}
				else
				{
					tempointm1[0][j] = tespoint_support[0][j];
					tempointp1[0][j] = tespoint_support[0][j];

				}


			}

			double fxm1 = acqFunctionUCB(Ainv_support, wi, n_points_support, Dimension_support, param_support,
				xtrainLBFGS_support, tempointm1, kappa_support, type_support);
			double fxp1 = acqFunctionUCB(Ainv_support, wi, n_points_support, Dimension_support, param_support,
				xtrainLBFGS_support, tempointp1, kappa_support, type_support);

			g[i] = (fxp1 - fxm1) / (2 * h);

			freeMatrix(tempointm1, 1);
			freeMatrix(tempointp1, 1);
		}

		free(wi);

		/*          go back to the minimization routine. */
		goto L111;
	}

	/*     if (s_cmp(task, "NEW_X", (ftnlen)5, (ftnlen)5) == 0) { */
	if (*task == NEW_X) {
		goto L111;
	}

	for (int i = 0; i < n; i++)
		printf("%0.4f \t", x[i]);

	double* result = newDoubleV(N);

	for (int i = 0; i < N; i++)
	{
		result[i] = x[i];
	}


	freeMatrix(xtrainLBFGS_support, n_pointsA);
	free(ytrainLBFGS_support);
	freeMatrix(Ainv_support, n_pointsA);
	free(param_support);
	freeMatrix(tespoint_support, 1);
	return result;
}
