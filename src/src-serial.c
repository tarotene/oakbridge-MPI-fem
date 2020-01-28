/*
//   1D Steady-State Heat Transfer 
//   FEM with Piece-wise Linear Elements
//   CG (Conjugate Gradient) Method 
//
//   d/dx(CdT/dx) + Q = 0
//   T=0@x=0  
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

int main()
{
	int NE, N, NPLU, IterMax;
	int R, Z, Q, P, DD;

	double dX, Resid, Eps, Area, QV, COND;
	double X1, X2, U1, U2, DL, Strain, Sigma, Ck;
	double QN, XL, C2, Xi, PHIa;
	double *PHI, *Rhs, *X;
	double *Diag, *AMat;
	double **W;

	int *Index, *Item, *Icelnod;
	double Kmat[2][2], Emat[2][2];

	int i, j, in1, in2, k, icel, k1, k2, jS;
	int iter;
	FILE *fp;
	double BNorm2, Rho, Rho1 = 0.0, C1, Alpha, DNorm2;
	int ierr = 1;
	int errno = 0;

	/*
// +-------+
// | INIT. |
// +-------+
*/
	fp = fopen("./input.dat", "r");
	assert(fp != NULL);
	fscanf(fp, "%d", &NE);
	fscanf(fp, "%lf %lf %lf %lf", &dX, &QV, &Area, &COND);
	fscanf(fp, "%d", &IterMax);
	fscanf(fp, "%lf", &Eps);
	fclose(fp);

	N = NE + 1;

	PHI = calloc(N, sizeof(double));
	X = calloc(N, sizeof(double));
	Diag = calloc(N, sizeof(double));
	AMat = calloc(2 * N - 2, sizeof(double));
	Rhs = calloc(N, sizeof(double));
	Index = calloc(N + 1, sizeof(int));
	Item = calloc(2 * N - 2, sizeof(int));
	Icelnod = calloc(2 * NE, sizeof(int));

	W = (double **)malloc(sizeof(double *) * 4);
	if (W == NULL)
	{
		fprintf(stderr, "Error: %s\n", strerror(errno));
		return -1;
	}
	for (i = 0; i < 4; i++)
	{
		W[i] = (double *)malloc(sizeof(double) * N);
		if (W[i] == NULL)
		{
			fprintf(stderr, "Error: %s\n", strerror(errno));
			return -1;
		}
	}

	for (i = 0; i < N; i++)
		PHI[i] = 0.0;
	for (i = 0; i < N; i++)
		Diag[i] = 0.0;
	for (i = 0; i < N; i++)
		Rhs[i] = 0.0;
	for (k = 0; k < 2 * N - 2; k++)
		AMat[k] = 0.0;
	for (i = 0; i < N; i++)
		X[i] = i * dX;
	for (icel = 0; icel < NE; icel++)
	{
		Icelnod[2 * icel] = icel;
		Icelnod[2 * icel + 1] = icel + 1;
	}

	Kmat[0][0] = +1.0;
	Kmat[0][1] = -1.0;
	Kmat[1][0] = -1.0;
	Kmat[1][1] = +1.0;
	/*
// +--------------+
// | CONNECTIVITY |
// +--------------+
*/
	for (i = 0; i < N + 1; i++)
		Index[i] = 2;
	Index[0] = 0;
	Index[1] = 1;
	Index[N] = 1;

	for (i = 0; i < N; i++)
	{
		Index[i + 1] = Index[i + 1] + Index[i];
	}

	NPLU = Index[N];

	for (i = 0; i < N; i++)
	{
		int jS = Index[i];
		if (i == 0)
		{
			Item[jS] = i + 1;
		}
		else if (i == N - 1)
		{
			Item[jS] = i - 1;
		}
		else
		{
			Item[jS] = i - 1;
			Item[jS + 1] = i + 1;
		}
	}

	/*
// +-----------------+
// | MATRIX assemble |
// +-----------------+
*/
	for (icel = 0; icel < NE; icel++)
	{
		in1 = Icelnod[2 * icel];
		in2 = Icelnod[2 * icel + 1];
		X1 = X[in1];
		X2 = X[in2];
		DL = fabs(X2 - X1);

		Ck = Area * COND / DL;
		Emat[0][0] = Ck * Kmat[0][0];
		Emat[0][1] = Ck * Kmat[0][1];
		Emat[1][0] = Ck * Kmat[1][0];
		Emat[1][1] = Ck * Kmat[1][1];

		Diag[in1] = Diag[in1] + Emat[0][0];
		Diag[in2] = Diag[in2] + Emat[1][1];

		if (icel == 0)
		{
			k1 = Index[in1];
		}
		else
		{
			k1 = Index[in1] + 1;
		}
		k2 = Index[in2];

		AMat[k1] = AMat[k1] + Emat[0][1];
		AMat[k2] = AMat[k2] + Emat[1][0];

		QN = 0.5 * QV * Area * dX;
		Rhs[in1] = Rhs[in1] + QN;
		Rhs[in2] = Rhs[in2] + QN;
	}

	/*
// +---------------------+
// | BOUNDARY conditions |
// +---------------------+
*/

	/* X=Xmin */
	i = 0;
	jS = Index[i];
	AMat[jS] = 0.0;
	Diag[i] = 1.0;
	Rhs[i] = 0.0;

	for (k = 0; k < NPLU; k++)
	{
		if (Item[k] == 0)
		{
			AMat[k] = 0.0;
		}
	}
	/*
// +---------------+
// | CG iterations |
// +---------------+
*/
	R = 0;
	Z = 1;
	Q = 1;
	P = 2;
	DD = 3;

	for (i = 0; i < N; i++)
	{
		W[DD][i] = 1.0 / Diag[i];
	}

	/*
//-- {r0}= {b} - [A]{xini} |
*/
	for (i = 0; i < N; i++)
	{
		W[R][i] = Diag[i] * PHI[i];
		for (j = Index[i]; j < Index[i + 1]; j++)
		{
			W[R][i] += AMat[j] * PHI[Item[j]];
		}
	}

	BNorm2 = 0.0;
	for (i = 0; i < N; i++)
	{
		BNorm2 += Rhs[i] * Rhs[i];
		W[R][i] = Rhs[i] - W[R][i];
	}

	clock_t start = clock();

	for (iter = 1; iter <= IterMax; iter++)
	{
		/*
//-- {z}= [Minv]{r}
*/
		for (i = 0; i < N; i++)
		{
			W[Z][i] = W[DD][i] * W[R][i];
		}

		/*
//-- RHO= {r}{z}
*/
		Rho = 0.0;
		for (i = 0; i < N; i++)
		{
			Rho += W[R][i] * W[Z][i];
		}

		/*
//-- {p} = {z} if      ITER=1  
//   BETA= RHO / RHO1  otherwise 
*/
		if (iter == 1)
		{
			for (i = 0; i < N; i++)
			{
				W[P][i] = W[Z][i];
			}
		}
		else
		{
			double Beta = Rho / Rho1;
			for (i = 0; i < N; i++)
			{
				W[P][i] = W[Z][i] + Beta * W[P][i];
			}
		}

		/*
//-- {q}= [A]{p}
*/
		for (i = 0; i < N; i++)
		{
			W[Q][i] = Diag[i] * W[P][i];
			for (j = Index[i]; j < Index[i + 1]; j++)
			{
				W[Q][i] += AMat[j] * W[P][Item[j]];
			}
		}

		/*
//-- ALPHA= RHO / {p}{q}
*/
		C1 = 0.0;
		for (i = 0; i < N; i++)
		{
			C1 += W[P][i] * W[Q][i];
		}
		Alpha = Rho / C1;

		/*
//-- {x}= {x} + ALPHA*{p}
//   {r}= {r} - ALPHA*{q}
*/
		for (i = 0; i < N; i++)
		{
			PHI[i] += Alpha * W[P][i];
			W[R][i] -= Alpha * W[Q][i];
		}

		DNorm2 = 0.0;
		for (i = 0; i < N; i++)
		{
			DNorm2 += W[R][i] * W[R][i];
		}

		Resid = sqrt(DNorm2 / BNorm2);

		// if ((iter) % 1000 == 0)
		// {
		// 	printf("%8d%s%16.6e\n",
		// 		   iter, " iters, RESID=", Resid);
		// }

		if (Resid <= Eps)
		{
			ierr = 0;
			break;
		}

		Rho1 = Rho;
	}

	clock_t end = clock();
	/* ********************************************************************
*/
	double elapsed = ((double)(end - start)) / ((double)CLOCKS_PER_SEC);

	printf("Elapsed Time: %13.3e\n", elapsed);

	fp = fopen("./out/out-serial.dat", "w");
	fprintf(fp, "%s\n", "# node id, temperature, exact temp.");
	// XL = NE * dX;
	// C2 = QV * NE * dX;

	for (i = 0; i < N; i++)
	{
		Xi = X[i];
		PHIa = (-0.5 * QV * Xi * Xi + QV * NE * dX * Xi) / COND;
		fprintf(fp, "%9d, %11.3e, %11.3e\n", i + 1, PHI[i], PHIa);
	}

	fclose(fp);

	return ierr;
}
