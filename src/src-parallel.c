#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

int main(int argc, char **argv)
{
	/* variables */
	int NE, N, NP, NPLU, IterMax, NEg, Ng, errno;
	double dX, Resid, Eps, Area, QV, COND, QN;
	double X1, X2, DL, Ck;
	double *PHI, *Rhs, *X, *Diag, *AMat;
	double *R, *Z, *Q, *P, *DD;
	int *Index, *Item, *Icelnod;
	double Kmat[2][2], Emat[2][2];

	int i, j, in1, in2, k, icel, k1, k2, jS;
	int iter, nr, neib;
	FILE *fp;

	double BNorm2, Rho, Rho1 = 0.0, C1, Alpha, Beta, DNorm2;
	int PETot, MyRank, kk, is, ir, len_s, len_r, tag;
	int NeibPETot, BufLength, NeibPE[2];

	int import_index[3], import_item[2];
	int export_index[3], export_item[2];
	double SendBuf[2], RecvBuf[2];

	double BNorm20, Rho0, C10, DNorm20;
	double StartTime, EndTime;
	int ierr = 1;

	MPI_Status *StatSend, *StatRecv;
	MPI_Request *RequestSend, *RequestRecv;

	/* initialization */

	// control data
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &PETot);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

	if (MyRank == 0)
	{
		fp = fopen("./input.dat", "r");
		assert(fp != NULL);
		fscanf(fp, "%d", &NEg);
		fscanf(fp, "%lf %lf %lf %lf", &dX, &QV, &Area, &COND);
		fscanf(fp, "%d", &IterMax);
		fscanf(fp, "%lf", &Eps);
		fclose(fp);
	}

	ierr = MPI_Bcast(&NEg, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ierr = MPI_Bcast(&IterMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ierr = MPI_Bcast(&dX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ierr = MPI_Bcast(&QV, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ierr = MPI_Bcast(&Area, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ierr = MPI_Bcast(&COND, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ierr = MPI_Bcast(&Eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// local mesh size
	Ng = NEg + 1;
	N = Ng / PETot;

	nr = Ng - N * PETot;
	if (MyRank < nr)
	{
		N++;
	}

	NE = N - 1 + 2;
	NP = N + 2;

	if (MyRank == 0)
	{
		NE = N - 1 + 1;
		NP = N + 1;
	}

	if (MyRank == PETot - 1)
	{
		NE = N - 1 + 1;
		NP = N + 1;
	}

	if (PETot == 1)
	{
		NE = N - 1;
		NP = N;
	}

	// arrays
	PHI = calloc(NP, sizeof(double));
	Diag = calloc(NP, sizeof(double));
	AMat = calloc(2 * NP - 2, sizeof(double));
	Rhs = calloc(NP, sizeof(double));
	Index = calloc(NP + 1, sizeof(int));
	Item = calloc(2 * NP - 2, sizeof(int));
	Icelnod = calloc(2 * NE, sizeof(int));

	for (i = 0; i < NP; i++)
	{
		PHI[i] = 0.0;
		Diag[i] = 0.0;
		Rhs[i] = 0.0;
	}

	for (k = 0; k < 2 * NP - 2; k++)
	{
		AMat[k] = 0.0;
	}

	for (i = 0; i < 3; i++)
	{
		import_index[i] = 0;
		export_index[i] = 0;
	}

	for (i = 0; i < 2; i++)
	{
		import_item[i] = 0;
		export_item[i] = 0;
	}

	for (icel = 0; icel < NE; icel++)
	{
		Icelnod[2 * icel] = icel;
		Icelnod[2 * icel + 1] = icel + 1;
	}

	if (PETot > 1)
	{
		if (MyRank == 0)
		{
			Icelnod[2 * (NE - 1) + 0] = N - 1;
			Icelnod[2 * (NE - 1) + 1] = N;
		}
		else if (MyRank == PETot - 1)
		{
			Icelnod[2 * (NE - 1) + 0] = N;
			Icelnod[2 * (NE - 1) + 1] = 0;
		}
		else
		{
			Icelnod[2 * (NE - 2) + 0] = N;
			Icelnod[2 * (NE - 2) + 1] = 0;
			Icelnod[2 * (NE - 1) + 0] = N - 1;
			Icelnod[2 * (NE - 1) + 1] = N + 1;
		}
	}

	Kmat[0][0] = +1.0;
	Kmat[0][1] = -1.0;
	Kmat[1][0] = -1.0;
	Kmat[1][1] = +1.0;

	// connectivity
	for (i = 0; i <= N; i++)
	{
		Index[i] = 2;
	}

	Index[0] = 0;

	Index[N + 1] = 1;
	Index[NP] = 1;

	if (MyRank == 0)
	{
		Index[1] = 1;
	}
	if (MyRank == PETot - 1)
	{
		Index[N] = 1;
	}

	for (i = 0; i < NP; i++)
	{
		Index[i + 1] = Index[i + 1] + Index[i];
	}

	NPLU = Index[NP];

	for (jS = 0; jS < 2 * NP - 2; jS++)
	{
		Item[jS] = 0;
	}

	for (i = 0; i < N; i++)
	{
		jS = Index[i];
		if (MyRank == 0 && i == 0)
		{
			Item[jS] = i + 1;
		}
		else if (MyRank == PETot - 1 && i == N - 1)
		{
			Item[jS] = i - 1;
		}
		else
		{
			Item[jS] = i - 1;
			Item[jS + 1] = i + 1;
			if (i == 0)
			{
				Item[jS] = N;
			}
			if (i == N - 1)
			{
				Item[jS + 1] = N + 1;
			}
			if (MyRank == 0 && i == N - 1)
			{
				Item[jS + 1] = N;
			}
		}
	}

	if (MyRank == 0)
	{
		Item[Index[N]] = N - 1;
	}
	else
	{
		Item[Index[N]] = 0;
	}

	if (MyRank != 0 && MyRank != PETot - 1)
	{
		Item[Index[N + 1]] = N - 1;
	}

	// communication
	NeibPETot = 2;
	if (MyRank == 0)
	{
		NeibPETot = 1;
	}
	if (MyRank == PETot - 1)
	{
		NeibPETot = 1;
	}
	if (PETot == 1)
	{
		NeibPETot = 0;
	}

	NeibPE[0] = MyRank - 1;
	NeibPE[1] = MyRank + 1;

	if (MyRank == 0)
	{
		NeibPE[0] = MyRank + 1;
	}

	if (MyRank == PETot - 1)
	{
		NeibPE[0] = MyRank - 1;
	}

	import_index[1] = 1;
	import_index[2] = 2;
	import_item[0] = N;
	import_item[1] = N + 1;

	export_index[1] = 1;
	export_index[2] = 2;
	export_item[0] = 0;
	export_item[1] = N - 1;

	if (MyRank == 0)
	{
		import_item[0] = N;
		export_item[0] = N - 1;
		// NeibPE[0] = MyRank + 1;
	}

	BufLength = 1;

	StatSend = malloc(sizeof(MPI_Status) * NeibPETot);
	StatRecv = malloc(sizeof(MPI_Status) * NeibPETot);
	RequestSend = malloc(sizeof(MPI_Request) * NeibPETot);
	RequestRecv = malloc(sizeof(MPI_Request) * NeibPETot);

	// StatSend = malloc(sizeof(MPI_Status) * 2 * NeibPETot);
	// RequestSend = malloc(sizeof(MPI_Request) * 2 * NeibPETot);

	/* matrix assemble */
	for (icel = 0; icel < NE; icel++)
	{
		in1 = Icelnod[2 * icel];
		in2 = Icelnod[2 * icel + 1];
		DL = dX;

		Ck = Area * COND / DL;
		Emat[0][0] = Ck * Kmat[0][0];
		Emat[0][1] = Ck * Kmat[0][1];
		Emat[1][0] = Ck * Kmat[1][0];
		Emat[1][1] = Ck * Kmat[1][1];

		Diag[in1] += Emat[0][0];
		Diag[in2] += Emat[1][1];

		if (MyRank == 0 && icel == 0)
		{
			k1 = Index[in1];
		}
		else
		{
			k1 = Index[in1] + 1;
		}

		k2 = Index[in2];

		AMat[k1] += Emat[0][1];
		AMat[k2] += Emat[1][0];

		QN = 0.5 * QV * Area * dX;
		Rhs[in1] += QN;
		Rhs[in2] += QN;
	}

	/* boudary conditions */

	// X = Xmin
	if (MyRank == 0)
	{
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
	}

	/* CG iterations */
	R = calloc(NP, sizeof(double));
	Z = calloc(NP, sizeof(double));
	P = calloc(NP, sizeof(double));
	Q = calloc(NP, sizeof(double));
	DD = calloc(NP, sizeof(double));

	for (i = 0; i < N; i++)
	{
		DD[i] = 1.0 / Diag[i];
	}

	// {r0} = {b} - [A] {x_ini}
	for (neib = 0; neib < NeibPETot; neib++)
	{
		for (k = export_index[neib]; k < export_index[neib + 1]; k++)
		{
			kk = export_item[k];
			SendBuf[k] = PHI[kk];
		}
	}

	for (neib = 0; neib < NeibPETot; neib++)
	{
		is = export_index[neib];
		len_s = export_index[neib + 1] - export_index[neib];
		MPI_Isend(&SendBuf[is], len_s, MPI_DOUBLE, NeibPE[neib], 0, MPI_COMM_WORLD, &RequestSend[neib]);
	}

	for (neib = 0; neib < NeibPETot; neib++)
	{
		ir = import_index[neib];
		len_r = import_index[neib + 1] - import_index[neib];
		MPI_Irecv(&RecvBuf[ir], len_r, MPI_DOUBLE, NeibPE[neib], 0, MPI_COMM_WORLD, &RequestRecv[neib]);
	}

	MPI_Waitall(NeibPETot, RequestRecv, StatRecv);

	for (neib = 0; neib < NeibPETot; neib++)
	{
		for (k = import_index[neib]; k < import_index[neib + 1]; k++)
		{
			kk = import_item[k];
			PHI[kk] = RecvBuf[k];
		}
	}

	MPI_Waitall(NeibPETot, RequestSend, StatSend);

	for (i = 0; i < N; i++)
	{
		R[i] = Diag[i] * PHI[i];
		for (j = Index[i]; j < Index[i + 1]; j++)
		{
			R[i] += AMat[j] * PHI[Item[j]];
		}
	}

	BNorm20 = 0.0;
	for (i = 0; i < N; i++)
	{
		BNorm20 += Rhs[i] * Rhs[i];
		R[i] = Rhs[i] - R[i];
	}
	ierr = MPI_Allreduce(&BNorm20, &BNorm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	StartTime = MPI_Wtime();

	for (iter = 1; iter <= IterMax; iter++)
	{
		// {z} = [Minv] {r}
		for (i = 0; i < N; i++)
		{
			Z[i] = DD[i] * R[i];
		}

		// printf("1\n");

		// RHO = {r} {z}
		Rho0 = 0.0;
		for (i = 0; i < N; i++)
		{
			Rho0 += R[i] * Z[i];
		}
		ierr = MPI_Allreduce(&Rho0, &Rho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// printf("2\n");

		// {p} = {z} or BETA = RHO / RHO1
		if (iter == 1)
		{
			for (i = 0; i < N; i++)
			{
				P[i] = Z[i];
			}
		}
		else
		{
			Beta = Rho / Rho1;
			for (i = 0; i < N; i++)
			{
				P[i] = Z[i] + Beta * P[i];
			}
		}

		// printf("3\n");

		// {q} = [A] {p}
		for (neib = 0; neib < NeibPETot; neib++)
		{
			for (k = export_index[neib]; k < export_index[neib + 1]; k++)
			{
				kk = export_item[k];
				SendBuf[k] = P[kk];
			}
		}

		for (neib = 0; neib < NeibPETot; neib++)
		{
			is = export_index[neib];
			len_s = export_index[neib + 1] - export_index[neib];
			MPI_Isend(&SendBuf[is], len_s, MPI_DOUBLE, NeibPE[neib], 0, MPI_COMM_WORLD, &RequestSend[neib]);
		}

		for (neib = 0; neib < NeibPETot; neib++)
		{
			ir = import_index[neib];
			len_r = import_index[neib + 1] - import_index[neib];
			MPI_Irecv(&RecvBuf[ir], len_r, MPI_DOUBLE, NeibPE[neib], 0, MPI_COMM_WORLD, &RequestRecv[neib]);
		}

		MPI_Waitall(NeibPETot, RequestRecv, StatRecv);

		for (neib = 0; neib < NeibPETot; neib++)
		{
			for (k = import_index[neib]; k < import_index[neib + 1]; k++)
			{
				kk = import_item[k];
				P[kk] = RecvBuf[k];
			}
		}

		MPI_Waitall(NeibPETot, RequestSend, StatSend);

		for (i = 0; i < N; i++)
		{
			Q[i] = Diag[i] * P[i];
			for (j = Index[i]; j < Index[i + 1]; j++)
			{
				Q[i] += AMat[j] * P[Item[j]];
			}
		}

		// printf("4\n");

		// ALPHA = RHO / {p} {q}
		C10 = 0.0;
		for (i = 0; i < N; i++)
		{
			C10 += P[i] * Q[i];
		}
		ierr = MPI_Allreduce(&C10, &C1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		Alpha = Rho / C1;

		// printf("5\n");

		// {x} += ALPHA * {p}
		// {r} -= ALPHA * {q}
		for (i = 0; i < N; i++)
		{
			PHI[i] += Alpha * P[i];
			R[i] -= Alpha * Q[i];
		}

		DNorm20 = 0.0;
		for (i = 0; i < N; i++)
		{
			DNorm20 += R[i] * R[i];
		}

		// printf("6\n");

		ierr = MPI_Allreduce(&DNorm20, &DNorm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		Resid = sqrt(DNorm2 / BNorm2);

		// if (MyRank == 0)
		// {
		// 	printf("%8d%s%16.6e\n", iter, " iters, RESID=", Resid);
		// }

		if (Resid < Eps)
		{
			ierr = 0;
			break;
		}

		// printf("7\n");

		Rho1 = Rho;
	}

	EndTime = MPI_Wtime();

	double tick = MPI_Wtick();

	// output
	if (MyRank == 0)
	{
		printf("Elapsed Time: %13.3e\n", (EndTime - StartTime));
		fp = fopen("./out/out-parallel.dat", "w");
		fprintf(fp, "%s\n", "# rank, node id, temperature");
		fclose(fp);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for (i = 0; i < PETot; i++)
	{
		if (MyRank == i)
		{
			fp = fopen("./out/out-parallel.dat", "a");
			for (j = 0; j < N; j++)
			{
				fprintf(fp, "%6d, %7d, %11.3e\n", MyRank, j + 1, PHI[j]);
			}
			fclose(fp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	ierr = MPI_Finalize();

	return ierr;
}