#ifndef ToolsMod
#define ToolsMod


// Rename to MaskRArrayStruct

typedef struct {

	RealType *Array;
	int *Mask;

} MaskRArrayStruct;


typedef struct {

	ComplexType *Array;
	int *Mask;

} MaskCArrayStruct;


int Embarrassing;


void InitializeLog(void);
void PrintTitle(void);
void CreateFileName(char *FileName, const char *Name);
void PrintTime(void);
void PrintMessage(const char *Left, const char *Middle, const char *Right);
void AllocateMaskArray(MaskRArrayStruct *MaskArray, const int Dimension, const int InitArray, const int InitMask);
void AllocateMaskCArray(MaskCArrayStruct *MaskArray, const int Dimension, const int InitArray, const int InitMask);
void DismissMaskCArray(MaskCArrayStruct *MaskArray);
void FreeMaskArray(MaskRArrayStruct *MaskArray);
void FreeMaskCArray(MaskCArrayStruct *MaskArray);
void FillRArray(RealType *Array, const RealType Value, const int Dimension);
void FillCArray(ComplexType *Array, const ComplexType Value, const int Dimension);
void RArrayToCArray(ComplexType *CArray, const RealType *RArray, const int Dimension);
void MirrorIMatrix(int *Matrix, const int Dimension);
void MirrorRMatrix(RealType *Matrix, const int Dimension);
void TransposeRMatrix(RealType *Transpose, const RealType *Matrix, const int Dimension);
void NormalizeRArray(RealType *Array, const int Dimension);
void RealPlusRArray(RealType *Array, const RealType Addition, int Dimension);
void RealTimesRArray(RealType *Array, const RealType Factor, int Dimension);
void RealTimesCArray(ComplexType *Array, const RealType Factor, int Dimension);
void CArrayTimesCArray(ComplexType *Array, const ComplexType *FactorArray, const int Dimension);
void RMatrixTimesCMatrix(ComplexType *Product, const RealType *RMatrix, const ComplexType *CMatrix, const int Dimension);
void RMatrixTimesCMatrixCutOff(ComplexType *Product, const RealType *RMatrix, const ComplexType *CMatrix, const int Dimension, const int CutOff);
void CMatrixTimesRMatrix(ComplexType *Product, const ComplexType *CMatrix, const RealType *RMatrix, const int Dimension);
void CMatrixTimesRMatrixCutOff(ComplexType *Product, const ComplexType *CMatrix, const RealType *RMatrix, const int Dimension, const int CutOff);
void Conjugate(ComplexType *Conjugate, const ComplexType *Array, const int Dimension);
void IntToString(char *Output, const int Input);
void RealToString(char *Output, const RealType Input);
void CreateGrid(RealType *Grid, const RealType Min, const RealType Max, const int Steps);
void PrintIArray(const int *Array, const int RowTotal, const int ColumnTotal);
void PrintRArray(const RealType *Array, const int RowTotal, const int ColumnTotal);
void PrintCArray(const ComplexType *Array, const int RowTotal, const int ColumnTotal);
void SaveRArray(const char *FileName, const RealType *Array, const int RowTotal, const int ColumnTotal, int RowTrunc, int ColumnTrunc, const RealType *RowIdVector, const int Index);
void SaveCArray(const char *FileName, const ComplexType *Array, const int RowTotal, const int ColumnTotal, const RealType *RowIdVector);
void InitializeDiag(RealType *Matrix, RealType *EigenValues, int Dimension, const int DAndC);
void DiagonalizeMatrix(RealType *Matrix, RealType *EigenValues, int Dimension, const int DAndC);
void CloseDiag(const int DAndC);

long int SetTimer(void);
int Minus1Pow(const int Input);
int IntPow2(const int Input);
int IntPowN(const int Input, const int N);
RealType RealPow2(const RealType Input);
int IntMin(const int PriInput, const int SecInput);
int IntMax(const int PriInput, const int SecInput);
RealType RArrayMax(const RealType *Array, const int Dimension);
int FourDimId_(const int PriId, const int SecId, const int TerId, const int QuaId, const int SecDimension, const int TerDimension, const int QuaDimension);
RealType Lorentzian(const RealType Value, const RealType Width);
RealType Gaussian(const RealType Value, const RealType StandardDeviation);
int SumIArray(const int *Array, const int Dimension);
RealType SumRArray(const RealType *Array, const int Dimension);
int CompareIArray(const int *PriArray, const int *SecArray, const int Dimension);
ComplexType Trace(const ComplexType *Matrix, const int Dimension);
ComplexType CInnerProduct(const ComplexType *Bra, const ComplexType *Ket, const int Dimension);

#endif