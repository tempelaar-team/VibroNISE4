#include "GlobalsMod.h"
#include "VOverlapsMod.h"

#include <math.h>


static void SaveVOverlaps(const int PriState, const int SecState, const RealType *VOverlaps, const int VTotal);

static RealType VOverlap(const int VShifted, const int VUnshifted, const RealType VShift, const RealType VShiftSquared);
static RealType Laguerre(const int K, const RealType X);
static RealType Function1(const int N, const int M);
static RealType Function2(const int N, const int M, const RealType X);
static long int Factorial(const int LowBound, const int UpBound);


/////////////
// Globals //
/////////////


void BuildVOverlaps(RealType ***VOverlaps, const RealType *HuangRhys, const int Characters, const int VTotal, const int Embarrassing)
{
	int PriCCount, SecCCount, PriVCount, SecVCount, SmartCount;
	RealType VShift, VShiftSquared;
	FILE *FileID;
	
	*VOverlaps = malloc(Characters * Characters * sizeof(RealType *));
	
	if (!Embarrassing)
		FileID = fopen("VOverlaps.dat", "w"),
		fclose(FileID);
	
	for (PriCCount = 0; PriCCount < Characters; PriCCount++)
		for(SecCCount = 0; SecCCount < Characters; SecCCount++) {
			SmartCount = PriCCount * Characters + SecCCount;
			
			(*VOverlaps)[SmartCount] = malloc(VTotal * VTotal * sizeof(RealType));
			
			VShift = (RealType) (sqrt(HuangRhys[PriCCount]) - sqrt(HuangRhys[SecCCount]));
			VShiftSquared = VShift * VShift;
			for (PriVCount = 0; PriVCount < VTotal; PriVCount++)
				for (SecVCount = 0; SecVCount < VTotal; SecVCount++)
					(*VOverlaps)[SmartCount][PriVCount * VTotal + SecVCount] = VOverlap(PriVCount, SecVCount, VShift, VShiftSquared);
			if (!Embarrassing) SaveVOverlaps(PriCCount, SecCCount, (*VOverlaps)[SmartCount], VTotal);
		}
}


void FreeVOverlaps(RealType ***VOverlaps, const int Characters)
{
	int Count;
	
	for (Count = 0; Count < Characters * Characters; Count++) free((*VOverlaps)[Count]);
	free(*VOverlaps);
}


////////////
// Locals //
////////////


static RealType VOverlap(const int VShifted, const int VUnshifted, const RealType VShift, const RealType VShiftSquared)
{
	RealType Output = (RealType) exp(-VShiftSquared / 2);
	
	if (VShifted == VUnshifted)
		Output *= Laguerre(VUnshifted, VShiftSquared);
	else if (VShifted > VUnshifted)
		Output *= pow(VShift, VShifted - VUnshifted) * Function1(VUnshifted, VShifted) * Function2(VUnshifted, VShifted - VUnshifted, VShiftSquared);
	else
		Output *= pow(-VShift, VUnshifted - VShifted) * Function1(VShifted, VUnshifted) * Function2(VShifted, VUnshifted - VShifted, VShiftSquared);
	
	return Output;
}


static RealType Laguerre(const int K, const RealType X)
{
	RealType Output, Dummy1, Dummy2;
	int Count;
	
	if (K == 0)			Output = 1;
	else if (K == 1)	Output = 1 - X;
	else {
		Output = 0;
		Dummy2 = 1;
		Dummy1 = 1 - X;
		for (Count = 2; Count <= K; Count++) {
			Output = (RealType) ((2. * Count - 1 - X) * Dummy1 - (Count - 1) * Dummy2) / Count;
			Dummy2 = Dummy1;
			Dummy1 = Output;
		}		
	}
	
	return Output;
}


static RealType Function1(const int N, const int M)
{
	RealType Output;
	
	if (N >= (M - N))	Output = (RealType) sqrt((RealType) Factorial(N + 1, M)) / Factorial(1, M - N);
	else if (N > 0)		Output = (RealType) sqrt((RealType) Factorial(M - N + 1, M) / (Factorial(1, N) * Factorial(1, M - N)));
	else				Output = (RealType) (1 / sqrt(Factorial(1, M)));
	
	return Output;
}


static RealType Function2(const int N, const int M, const RealType X)
{
	RealType Output = 1;
	int Count;
	
	for (Count = N; Count >= 1; Count--) Output = 1 - (RealType) (N - Count + 1) / (Count * (M + Count)) * X * Output;
	
	return Output;
}


static long int Factorial(const int LowBound, const int UpBound)
{
	long int Output = LowBound;
	int Count;
	
	for (Count = LowBound + 1; Count <= UpBound; Count++) Output *= Count;
	
	return Output;
}


static void SaveVOverlaps(const int PriState, const int SecState, const RealType *VOverlaps, const int VTotal)
{
	int RowCount, ColumnCount;
	FILE *FileID;
	
	FileID = fopen("VOverlaps.dat", "a");
	fprintf(FileID, "%d -> %d\n", PriState, SecState);
	for (RowCount = 0; RowCount < VTotal; RowCount++) {
		for (ColumnCount = 0; ColumnCount < VTotal; ColumnCount++) fprintf(FileID, "%f, ", VOverlaps[RowCount * VTotal + ColumnCount]);
		fprintf(FileID, "\n");
	}
	fprintf(FileID, "\n");
	fclose(FileID);
}