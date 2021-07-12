#include "GlobalsMod.h"
#include "ToolsMod.h"

#include <time.h>
#include <math.h>


#ifdef DoublePrec
extern void dsyev_(char *Jobz, char *UpLo, int *N, double *A, int *LDA, double *W, double *Work, int *LWork, int *Info);
extern void dsyevd_(char *Jobz, char *UpLo, int *N, double *A, int *LDA, double *W, double *Work, int *LWork, int *IWork, int *LIWork, int *Info);
#else
extern void ssyev_(char *Jobz, char *UpLo, int *N, float *A, int *LDA, float *W, float *Work, int *LWork, int *Info);
extern void ssyevd_(char *Jobz, char *UpLo, int *N, float *A, int *LDA, float *W, float *Work, int *LWork, int *IWork, int *LIWork, int *Info);
#endif

RealType *Work;
int *IWork;
int LIWork, LWork;


void InitializeLog(void)
{
	FILE *FileID;
	char FileName[20];
	
	CreateFileName(FileName, "Log"); FileID = fopen(FileName, "w");
	fclose(FileID);
}


void PrintTitle(void)
{
	
#ifdef DoublePrec
	char PrecisionMode[] = "double";
#else
	char PrecisionMode[] = "single";
#endif

	PrintMessage("", "", "");
	PrintMessage("     // // // ////  ////   ///  /  // //  //// ////           ", "", "");
	PrintMessage("    // // // // // // // // // // // // //    //              ", "", "");
	PrintMessage("   // // // ////  ////  // // ///// //  ///  ///              ", "", "");
	PrintMessage("   ///  // // // // // // // // // //    // //                ", "", "");
	PrintMessage("   /   // ////  //  // ///  //  / // ////  ////  Version 4.26 ", "", "");
	PrintMessage("", "", "");
	PrintMessage("                 CRYSTAL EDITION ", "", "");
	PrintMessage("", "", "");
	PrintTime();
	PrintMessage("   Running in ", PrecisionMode, " precision mode.");
	PrintMessage("", "", "");
}


void CreateFileName(char *FileName, const char *Name)
{
	sprintf(FileName, "%s", Name);
	if (Embarrassing) sprintf(FileName, "%s_%d", FileName, Embarrassing);
	sprintf(FileName, "%s%s", FileName, ".dat");
}


void PrintTime(void)
{
	time_t TimeRawFormat;
	struct tm *GetTime;
	char Time[25], FileName[20];
	time (&TimeRawFormat);
	FILE *FileID;
	
	CreateFileName(FileName, "Log"); FileID = fopen(FileName, "a");
	GetTime = localtime(&TimeRawFormat);
	strftime(Time, 25, "%Y:%m:%d %H:%M:%S", GetTime);
	printf("%s - ", Time);
	fprintf(FileID, "%s - ", Time);
	fclose(FileID);
}


void PrintMessage(const char *Left, const char *Middle, const char *Right)
{
	char FileName[20];
	FILE *FileID;
	
	CreateFileName(FileName, "Log"); FileID = fopen(FileName, "a");
	printf("%s%s%s\n", Left, Middle, Right);
	fprintf(FileID, "%s%s%s\n", Left, Middle, Right);
	fclose(FileID);
}


long int SetTimer(void)
{
	time_t TimeRawFormat;
	struct tm *GetTime;
	char ContrChar[25];
	time (&TimeRawFormat);
	int ContrInt;
	long int Timer;
	
	GetTime = localtime(&TimeRawFormat);
	
	strftime(ContrChar, 25, "%S", GetTime);
	sscanf(ContrChar, "%d", &ContrInt);
	Timer = ContrInt;
	
	strftime(ContrChar, 25, "%M", GetTime);
	sscanf(ContrChar, "%d", &ContrInt);
	Timer += ContrInt * 60;
	
	strftime(ContrChar, 25, "%H", GetTime);
	sscanf(ContrChar, "%d", &ContrInt);
	Timer += ContrInt * 3600;
	
	strftime(ContrChar, 25, "%d", GetTime);
	sscanf(ContrChar, "%d", &ContrInt);
	Timer += ContrInt * 3600 * 24;
	
	return Timer;
}


int Minus1Pow(const int Input)
{
	return Input % 2 == 0 ? 1 : -1;
}


int IntPow2(const int Input)
{
	return Input * Input;
}


int IntPowN(const int Input, const int N)
{
	int Count, ReturnValue = Input;

	for (Count = 1; Count < N; Count++) ReturnValue *= Input;

	return ReturnValue;
}


RealType RealPow2(const RealType Input)
{
	return Input * Input;
}


int IntMin(const int PriInput, const int SecInput)
{
	return (PriInput < SecInput) ? PriInput : SecInput;
}


int IntMax(const int PriInput, const int SecInput)
{
	return (PriInput > SecInput) ? PriInput : SecInput;
}


RealType RArrayMax(const RealType *Array, const int Dimension)
{
	const RealType *ArrayElement = Array;
	RealType MaxValue = (RealType) (-1. / 0.);
	int Count;

	for (Count = 0; Count < Dimension; Count ++) {
		if (*ArrayElement > MaxValue) MaxValue = *ArrayElement;
		ArrayElement++;
	}

	return MaxValue;
}


int FourDimId_(const int PriId, const int SecId, const int TerId, const int QuaId, const int SecDimension, const int TerDimension, const int QuaDimension)
{
	return ((PriId * SecDimension + SecId) * TerDimension + TerId) * QuaDimension + QuaId;
}


RealType Lorentzian(const RealType Value, const RealType Width)
{
	return Width / (RealPow2(Value) + RealPow2(Width));
}


RealType Gaussian(const RealType Value, const RealType StandardDeviation)
{
	return (RealType) exp(-.5 * RealPow2(Value) / RealPow2(StandardDeviation)) / StandardDeviation;
}


int SumIArray(const int *Array, const int Dimension)
{
	int Count, Value = 0;
	const int *ArrayPos = Array;

	for (Count = 0; Count < Dimension; Count++) Value += (*ArrayPos++);

	return Value;
}


RealType SumRArray(const RealType *Array, const int Dimension)
{
	int Count;
	RealType Value = 0;
	const RealType *ArrayPos = Array;

	for (Count = 0; Count < Dimension; Count++) Value += (*ArrayPos++);

	return Value;
}


void AllocateMaskArray(MaskRArrayStruct *MaskArray, const int Dimension, const int InitArray, const int InitMask)
{
	if (InitArray)	(*MaskArray).Array = calloc((size_t) Dimension, sizeof(RealType));
	else			(*MaskArray).Array = malloc(Dimension * sizeof(RealType));
	if (InitMask)	(*MaskArray).Mask = calloc((size_t) Dimension, sizeof(int));
	else			(*MaskArray).Mask = malloc(Dimension * sizeof(int));
}


void AllocateMaskCArray(MaskCArrayStruct *MaskArray, const int Dimension, const int InitArray, const int InitMask)
{
	if (InitArray)	(*MaskArray).Array = calloc((size_t) Dimension, sizeof(ComplexType));
	else			(*MaskArray).Array = malloc(Dimension * sizeof(ComplexType));
	if (InitMask)	(*MaskArray).Mask = calloc((size_t) Dimension, sizeof(int));
	else			(*MaskArray).Mask = malloc(Dimension * sizeof(int));
}


void DismissMaskCArray(MaskCArrayStruct *MaskArray)
{
	(*MaskArray).Array = NULL; (*MaskArray).Mask = NULL;
}


void FreeMaskArray(MaskRArrayStruct *MaskArray)
{
	free((*MaskArray).Array); free((*MaskArray).Mask);
}


void FreeMaskCArray(MaskCArrayStruct *MaskArray)
{
	if ((*MaskArray).Mask == NULL) return;
	free((*MaskArray).Array); free((*MaskArray).Mask);
}


int CompareIArray(const int *PriArray, const int *SecArray, const int Dimension)
{
	const int *PriArrayElement = PriArray, *SecArrayElement = SecArray;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) if ((*PriArrayElement++) != (*SecArrayElement++)) return false;
	
	return true;
}


ComplexType Trace(const ComplexType *Matrix, const int Dimension)
{
	int Count;
	ComplexType Value = 0;

	for (Count = 0; Count < Dimension; Count++) Value += Matrix[Count * (Dimension + 1)];

	return Value;
}


void FillRArray(RealType *Array, const RealType Value, const int Dimension)
{
	RealType *ArrayElement = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) = Value;
}


void FillCArray(ComplexType *Array, const ComplexType Value, const int Dimension)
{
	ComplexType *ArrayElement = Array;
	int Count;

	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) = Value;
}


void RArrayToCArray(ComplexType *CArray, const RealType *RArray, const int Dimension)
{
	ComplexType *CArrayElement = CArray;
	const RealType *RArrayElement = RArray;
	int Count;

	for (Count = 0; Count < Dimension; Count++) (*CArrayElement++) = (*RArrayElement++);
}


void MirrorIMatrix(int *Matrix, const int Dimension)
{
	int *MatrixElement;
	int PriCount, SecCount;
	
	for (PriCount = 0; PriCount < Dimension; PriCount++) {
		MatrixElement = &Matrix[PriCount * Dimension];
		for (SecCount = 0; SecCount < PriCount; SecCount++) Matrix[SecCount * Dimension + PriCount] = (*MatrixElement++);
	}
}


void MirrorRMatrix(RealType *Matrix, const int Dimension)
{
	RealType *MatrixElement;
	int PriCount, SecCount;
	
	for (PriCount = 0; PriCount < Dimension; PriCount++) {
		MatrixElement = &Matrix[PriCount * Dimension];
		for (SecCount = 0; SecCount < PriCount; SecCount++) Matrix[SecCount * Dimension + PriCount] = (*MatrixElement++);
	}
}


void TransposeRMatrix(RealType *Transpose, const RealType *Matrix, const int Dimension)
{
	RealType *TransposeElement;
	const RealType *MatrixElement;
	int PriCount, SecCount;

	for (PriCount = 0; PriCount < Dimension; PriCount++) {
		TransposeElement = &Transpose[PriCount];
		MatrixElement = &Matrix[PriCount * Dimension];
		for (SecCount = 0; SecCount < Dimension; SecCount++) *TransposeElement = (*MatrixElement++), TransposeElement += Dimension;
	}
}


void NormalizeRArray(RealType *Array, const int Dimension)
{
	RealType *ArrayElement;
	RealType Maximum = 0;
	int Count;

	ArrayElement = Array;
	for (Count = 0; Count < Dimension; Count++) {
		if (fabs(*ArrayElement) > Maximum) Maximum = (RealType) fabs(*ArrayElement);
		ArrayElement++;
	}

	ArrayElement = Array;
	if (Maximum > 0) for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) /= Maximum;
}


void RealPlusRArray(RealType *Array, const RealType Addition, int Dimension)
{
	RealType *ArrayElement = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) += Addition;
}


void RealTimesRArray(RealType *Array, const RealType Factor, int Dimension)
{
	RealType *ArrayElement = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) *= Factor;
}


void RealTimesCArray(ComplexType *Array, const RealType Factor, int Dimension)
{
	ComplexType *ArrayElement = Array;
	int Count;
	
	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) *= Factor;
}


void CArrayTimesCArray(ComplexType *Array, const ComplexType *FactorArray, const int Dimension)
{
	const ComplexType *FactorArrayElement = &FactorArray[0];
	ComplexType *ArrayElement = &Array[0];
	int Count;

	for (Count = 0; Count < Dimension; Count++) (*ArrayElement++) *= (*FactorArrayElement++);
}


ComplexType CInnerProduct(const ComplexType *Bra, const ComplexType *Ket, const int Dimension)
{
	const ComplexType *BraEntry = Bra, *KetEntry = Ket;
	ComplexType Value = 0;
	int Count;

	for (Count = 0; Count < Dimension; Count++) Value += conj(*BraEntry++) * (*KetEntry++);

	return Value;
}


void RMatrixTimesCMatrix(ComplexType *Product, const RealType *RMatrix, const ComplexType *CMatrix, const int Dimension)
{
	ComplexType *ProductElement = Product;
	const ComplexType *CMatrixElement;
	const RealType *RMatrixElement;
	int PriCount, SecCount, SharedCount;

	for (PriCount = 0; PriCount < Dimension; PriCount++)
		for (SecCount = 0; SecCount < Dimension; SecCount++) {
			RMatrixElement = &RMatrix[PriCount * Dimension]; CMatrixElement = &CMatrix[SecCount];
			*ProductElement = 0;
			for (SharedCount = 0; SharedCount < Dimension; SharedCount++) *ProductElement += (*RMatrixElement++) * *CMatrixElement, CMatrixElement += Dimension;
			ProductElement++;
		}
}


void RMatrixTimesCMatrixCutOff(ComplexType *Product, const RealType *RMatrix, const ComplexType *CMatrix, const int Dimension, const int CutOff)
{
	ComplexType *ProductElement = Product;
	const ComplexType *CMatrixElement;
	const RealType *RMatrixElement;
	int PriCount, SecCount, SharedCount;

	for (PriCount = 0; PriCount < Dimension; PriCount++)
		for (SecCount = 0; SecCount < Dimension; SecCount++) {
			RMatrixElement = &RMatrix[PriCount * Dimension]; CMatrixElement = &CMatrix[SecCount];
			*ProductElement = 0;
			for (SharedCount = 0; SharedCount < CutOff; SharedCount++) *ProductElement += (*RMatrixElement++) * *CMatrixElement, CMatrixElement += Dimension;
			ProductElement++;
		}
}


void CMatrixTimesRMatrix(ComplexType *Product, const ComplexType *CMatrix, const RealType *RMatrix, const int Dimension)
{
	ComplexType *ProductElement = Product;
	const ComplexType *CMatrixElement;
	const RealType *RMatrixElement;
	int PriCount, SecCount, SharedCount;

	for (PriCount = 0; PriCount < Dimension; PriCount++)
		for (SecCount = 0; SecCount < Dimension; SecCount++) {
			CMatrixElement = &CMatrix[PriCount * Dimension]; RMatrixElement = &RMatrix[SecCount];
			*ProductElement = 0;
			for (SharedCount = 0; SharedCount < Dimension; SharedCount++) *ProductElement += (*CMatrixElement++) * *RMatrixElement, RMatrixElement += Dimension;
			ProductElement++;
		}
}


void CMatrixTimesRMatrixCutOff(ComplexType *Product, const ComplexType *CMatrix, const RealType *RMatrix, const int Dimension, const int CutOff)
{
	ComplexType *ProductElement = Product;
	const ComplexType *CMatrixElement;
	const RealType *RMatrixElement;
	int PriCount, SecCount, SharedCount;

	for (PriCount = 0; PriCount < Dimension; PriCount++)
		for (SecCount = 0; SecCount < Dimension; SecCount++) {
			*ProductElement = 0;
			if (PriCount < CutOff) {
				CMatrixElement = &CMatrix[PriCount * Dimension]; RMatrixElement = &RMatrix[SecCount];
				for (SharedCount = 0; SharedCount < CutOff; SharedCount++) *ProductElement += (*CMatrixElement++) * *RMatrixElement, RMatrixElement += Dimension;
			}
			ProductElement++;
		}
}


void Conjugate(ComplexType *Conjugate, const ComplexType *Array, const int Dimension)
{
	ComplexType *ConjugateElement = Conjugate;
	int PriCount, SecCount;

	for (PriCount = 0; PriCount < Dimension; PriCount++)
		for (SecCount = 0; SecCount < Dimension; SecCount++)
			(*ConjugateElement++) = (ComplexType) conj(Array[SecCount * Dimension + PriCount]);
}


void IntToString(char *Output, const int Input)
{
	sprintf(Output, "%d", Input);
}


void RealToString(char *Output, const RealType Input)
{
	sprintf(Output, "%5.2E", Input);
}


void CreateGrid(RealType *Grid, const RealType Min, const RealType Max, const int Steps)
{
	int Count;
	
	for (Count = 0; Count < Steps; Count++) Grid[Count] = Min + (Max - Min) * (RealType) Count / (Steps - 1);
}


void PrintIArray(const int *Array, const int RowTotal, const int ColumnTotal)
{
	int RowCount, ColumnCount;

	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) printf("%d, ", Array[RowCount * ColumnTotal + ColumnCount]);
		printf("\n");
	}
}


void PrintRArray(const RealType *Array, const int RowTotal, const int ColumnTotal)
{
	int RowCount, ColumnCount;

	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) printf("%f, ", Array[RowCount * ColumnTotal + ColumnCount]);
		printf("\n");
	}
}


void PrintCArray(const ComplexType *Array, const int RowTotal, const int ColumnTotal)
{
	int RowCount, ColumnCount;

	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) printf("%f + %fi, ", creal(Array[RowCount * ColumnTotal + ColumnCount]), cimag(Array[RowCount * ColumnTotal + ColumnCount]));
		printf("\n");
	}
}


void SaveRArray(const char *FileName, const RealType *Array, const int RowTotal, const int ColumnTotal, int RowTrunc, int ColumnTrunc, const RealType *RowIdVector, const int Index)
{
	int RowCount, ColumnCount;
	char ModFileName[20];
	FILE *FileID;

	if (RowTrunc == 0) RowTrunc = 10000;
	if (ColumnTrunc == 0) ColumnTrunc = 20;

	CreateFileName(ModFileName, FileName); FileID = fopen(ModFileName, "w");
	for (RowCount = 0; RowCount < IntMin(RowTrunc, RowTotal); RowCount++) {
		if (Index) fprintf(FileID, "%d, ", RowCount);
		if (RowIdVector != NULL) fprintf(FileID, "%f, ", RowIdVector[RowCount]);
		for (ColumnCount = 0; ColumnCount < IntMin(ColumnTrunc, ColumnTotal); ColumnCount++) fprintf(FileID, "%f, ", Array[RowCount * ColumnTotal + ColumnCount]);
		fprintf(FileID, "\n");
	}
	fclose(FileID);
}


void SaveCArray(const char *FileName, const ComplexType *Array, const int RowTotal, const int ColumnTotal, const RealType *RowIdVector)
{
	int RowCount, ColumnCount;
	char ModFileName[20];
	FILE *FileID;

	if (RowTotal * ColumnTotal > 100000) return;

	CreateFileName(ModFileName, FileName); FileID = fopen(ModFileName, "w");
	for (RowCount = 0; RowCount < RowTotal; RowCount++) {
		if (RowIdVector != NULL) fprintf(FileID, "%f, ", RowIdVector[RowCount]);
		for (ColumnCount = 0; ColumnCount < ColumnTotal; ColumnCount++) fprintf(FileID, "%f, %f, ", creal(Array[RowCount * ColumnTotal + ColumnCount]), cimag(Array[RowCount * ColumnTotal + ColumnCount]));
		fprintf(FileID, "\n");
	}
	fclose(FileID);
}


void InitializeDiag(RealType *Matrix, RealType *EigenValues, int Dimension, const int DAndC)
// You might wanna check out what's going on at the bottom of this subroutine...
{
	int Info;

	PrintTime(); PrintMessage("Diagonalizing. ", "", "");
	
	LWork = -1;
	Work = calloc(1, sizeof(RealType));
	if (DAndC) {
		LIWork = -1;
		IWork = calloc(1, sizeof(int));
#ifdef DoublePrec
		dsyevd_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, IWork, &LIWork, &Info);			
#else
		ssyevd_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, IWork, &LIWork, &Info);
#endif
		LIWork = IWork[0];
		free(IWork);
		IWork = calloc((size_t) LIWork, sizeof(int));
	}
	else {
#ifdef DoublePrec
		dsyev_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, &Info);
#else
		ssyev_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, &Info);
#endif
	}
	LWork = (int) Work[0];
	free(Work);
	if (LWork < 1 + 6 * Dimension + 2 * IntPow2(Dimension)) LWork = 1 + 6 * Dimension + 2 * IntPow2(Dimension);
	Work = calloc((size_t) LWork, sizeof(RealType));
}


void DiagonalizeMatrix(RealType *Matrix, RealType *EigenValues, int Dimension, const int DAndC)
{
	int Info;
	
	if (Dimension == 0) return;
	
	if (DAndC) {
#ifdef DoublePrec
		dsyevd_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, IWork, &LIWork, &Info);
#else
		ssyevd_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, IWork, &LIWork, &Info);
#endif
	}
	else {
#ifdef DoublePrec
		dsyev_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, &Info);
#else
		ssyev_("V", "U", &Dimension, Matrix, &Dimension, EigenValues, Work, &LWork, &Info);
#endif
	}
}


void CloseDiag(const int DAndC)
{
	free(Work);
	if (DAndC) free(IWork);
}
