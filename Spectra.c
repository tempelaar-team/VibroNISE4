#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "ToolsMod.h"
#include "ParmsHandleMod.h"

#include <math.h>
#include <fftw3.h>


static void Generate1DSpec(const ParmStruct *Parms);
static void FourierTransform1D(ComplexType *Output, const ComplexType *Input, const int DimensionIn, const int DimensionOut);
static void Generate2DSpec(ComplexType *TotalSpectrum, ComplexType *PriSubSpectrum, ComplexType *SecSubSpectrum, const RealType Rescaling, const char *IdChar, const ParmStruct *Parms, const int Rephasing);
static void FourierTransform2D(ComplexType *Output, const ComplexType *Input, const int DimensionIn, const int DimensionOut);
static void ReadResponse(ComplexType *Response, const char *FileName, const ParmStruct *Parms);
static void SaveEnergyGrid(const ParmStruct *Parms);
static void SaveSpectrum(const ComplexType *Spectrum, const char *FileName, const ParmStruct *Parms);


int main(int ArgsTotal, char *Args[])
{
	ParmStruct *Parms;
	ComplexType *Spectrum2D, *SpectrumRR, *SpectrumNR, *SpectrumGB, *SpectrumSE, *SpectrumEA;

	PrintTitle();
	PrintMessage("              // Spectra Generator //", "", "");
	PrintMessage("", "", "");
	
	Parms = calloc(1, sizeof(ParmStruct));

	InitializeParms(Parms);
	ReadParms(Args, Parms);
	DerivedParms(Parms);

	Parms->SpectrumTotal = 2 * *Parms->TimeTotal;

	if (Parms->Mode == ModeResponse1D) Generate1DSpec(Parms);
	else if (Parms->Mode == ModeResponse2D) {

		Spectrum2D = (ComplexType *) calloc((size_t) IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumRR = (ComplexType *) calloc((size_t) IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumNR = (ComplexType *) calloc((size_t) IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumGB = (ComplexType *) calloc((size_t) IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumSE = (ComplexType *) calloc((size_t) IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));
		SpectrumEA = (ComplexType *) calloc((size_t) IntPow2(Parms->SpectrumTotal), sizeof(ComplexType));

		SaveEnergyGrid(Parms);

		Generate2DSpec(Spectrum2D, SpectrumRR, SpectrumGB, 1.0, "GBRR", Parms, true);
		Generate2DSpec(Spectrum2D, SpectrumNR, SpectrumGB, 1.0, "GBNR", Parms, false);
		if (!Parms->NoSE) {
			Generate2DSpec(Spectrum2D, SpectrumRR, SpectrumSE, 1.0, "SERR", Parms, true);
			Generate2DSpec(Spectrum2D, SpectrumNR, SpectrumSE, 1.0, "SENR", Parms, false);
		}
		Generate2DSpec(Spectrum2D, SpectrumRR, SpectrumEA, Parms->RescaleEA, "EARR", Parms, true);
		Generate2DSpec(Spectrum2D, SpectrumNR, SpectrumEA, Parms->RescaleEA, "EANR", Parms, false);

		SaveSpectrum(Spectrum2D, "Spectrum2D.dat", Parms);
		SaveSpectrum(SpectrumRR, "SpectrumRR.dat", Parms);
		SaveSpectrum(SpectrumNR, "SpectrumNR.dat", Parms);
		SaveSpectrum(SpectrumGB, "SpectrumGB.dat", Parms);
		if (!Parms->NoSE) SaveSpectrum(SpectrumSE, "SpectrumSE.dat", Parms);
		SaveSpectrum(SpectrumEA, "SpectrumEA.dat", Parms);

		free(Spectrum2D); free(SpectrumRR); free(SpectrumNR); free(SpectrumGB); free(SpectrumSE); free(SpectrumEA);

	}
	
	PrintMessage("Finished!", "", "");
	
	free(Parms);
	
	return 0;
}


static void Generate1DSpec(const ParmStruct *Parms)
{
	double Dummy, ReadReal, ReadImag;
	RealType Energy;
	int Count;
	FILE *FileID;
	
	ComplexType *ResponseX = malloc(*Parms->TimeTotal * sizeof(ComplexType)), *SpectrumX = malloc(Parms->SpectrumTotal * sizeof(ComplexType));
	ComplexType *ResponseY = malloc(*Parms->TimeTotal * sizeof(ComplexType)), *SpectrumY = malloc(Parms->SpectrumTotal * sizeof(ComplexType));
	ComplexType *ResponseZ = malloc(*Parms->TimeTotal * sizeof(ComplexType)), *SpectrumZ = malloc(Parms->SpectrumTotal * sizeof(ComplexType));
	RealType *Spectrum = malloc(4 * Parms->SpectrumTotal * sizeof(RealType));

	FileID = fopen("Response1D.dat", "r");
	for (Count = 0; Count < *Parms->TimeTotal; Count++) {
		fscanf(FileID, "%lf, %lf, %lf,", &Dummy, &ReadReal, &ReadImag);
		ResponseX[Count] = ReadReal + I * ReadImag;
		if (Parms->LifeTime[0]) ResponseX[Count] *= exp(-Parms->TimeStep * Count / (2 * Parms->LifeTime[0]));
		fscanf(FileID, "%lf, %lf,", &ReadReal, &ReadImag);
		ResponseY[Count] = ReadReal + I * ReadImag;
		if (Parms->LifeTime[1]) ResponseY[Count] *= exp(-Parms->TimeStep * Count / (2 * Parms->LifeTime[1]));
		fscanf(FileID, "%lf, %lf,", &ReadReal, &ReadImag);
		ResponseZ[Count] = ReadReal + I * ReadImag;
		if (Parms->LifeTime[2]) ResponseZ[Count] *= exp(-Parms->TimeStep * Count / (2 * Parms->LifeTime[2]));
	}
	fclose(FileID);

	FourierTransform1D(SpectrumX, ResponseX, *Parms->TimeTotal, Parms->SpectrumTotal);
	FourierTransform1D(SpectrumY, ResponseY, *Parms->TimeTotal, Parms->SpectrumTotal);
	FourierTransform1D(SpectrumZ, ResponseZ, *Parms->TimeTotal, Parms->SpectrumTotal);

	for (Count = 0; Count < Parms->SpectrumTotal; Count++) {
		Spectrum[Parms->SpectrumTotal + Count] = (RealType) creal(SpectrumX[Count]);
		Spectrum[2 * Parms->SpectrumTotal + Count] = (RealType) creal(SpectrumY[Count]);
		Spectrum[3 * Parms->SpectrumTotal + Count] = (RealType) creal(SpectrumZ[Count]);
		Spectrum[Count] = (RealType) (creal(SpectrumX[Count]) + creal(SpectrumY[Count]) + creal(SpectrumZ[Count]));
	}

	if (Parms->Normalize) NormalizeRArray(Spectrum, 4 * Parms->SpectrumTotal);

	FileID = fopen("Spectrum1D.dat", "w");
	for (Count = 0; Count < Parms->SpectrumTotal; Count++) {

		Energy = 2 * Pi * (Count - Parms->SpectrumTotal / 2) / (Parms->EnergyToFreq * Parms->TimeStep * Parms->SpectrumTotal);
		Energy += Parms->SEnergyMean;

		if (Energy >= Parms->EnergyBounds[0] && Energy <= Parms->EnergyBounds[1])
			fprintf(FileID, "%f, %f, %f, %f, %f,\n", Energy, Spectrum[Count], Spectrum[Parms->SpectrumTotal + Count], Spectrum[2 * Parms->SpectrumTotal + Count], Spectrum[3 * Parms->SpectrumTotal + Count]);
	}
	fclose(FileID);
	
	free(ResponseX); free(ResponseY); free(ResponseZ); free(SpectrumX); free(SpectrumY); free(SpectrumZ); free(Spectrum);
}


static void FourierTransform1D(ComplexType *Output, const ComplexType *Input, const int DimensionIn, const int DimensionOut)
{
	fftw_complex *In, *Out;
	fftw_plan Plan;
	int Count;
	
	In = fftw_malloc(sizeof(fftw_complex) * DimensionOut);
	Out = fftw_malloc(sizeof(fftw_complex) * DimensionOut);
	Plan = fftw_plan_dft_1d(DimensionOut, In, Out, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	for (Count = 0; Count < DimensionOut; Count++) In[Count] = 0;
	for (Count = 0; Count < DimensionIn; Count++) {
		In[Count] = Input[Count];
		if (Count > 0) In[DimensionOut - Count] = conj(In[Count]);
	}
	
	fftw_execute(Plan);
	
	for (Count = 0; Count < DimensionOut; Count++)
		Output[Count] = (ComplexType) Out[(Count < (int) ceil(.5 * DimensionOut)) * (Count + (int) floor(.5 * DimensionOut)) + (Count > (int) floor(.5 * DimensionOut)) * (Count - (int) ceil(.5 * DimensionOut))];

	fftw_destroy_plan(Plan);
	fftw_free(In); fftw_free(Out);
	fftw_cleanup();
}


static void Generate2DSpec(ComplexType *TotalSpectrum, ComplexType *PriSubSpectrum, ComplexType *SecSubSpectrum, const RealType Rescaling, const char *IdChar, const ParmStruct *Parms, const int Rephasing)
{
	ComplexType *Response = malloc(IntPow2(*Parms->TimeTotal) * sizeof(ComplexType));
	ComplexType *SpectrumRaw = malloc(IntPow2(Parms->SpectrumTotal) * sizeof(ComplexType)), *Spectrum = malloc(IntPow2(Parms->SpectrumTotal) * sizeof(ComplexType));
	ComplexType SpectralValue;
	RealType PriEnergy, SecEnergy;
	int PriCount, SecCount;
	char FileName[100];

	sprintf(FileName, "%s%s%s", "Response", IdChar, ".dat");
	ReadResponse(Response, FileName, Parms);

	FourierTransform2D(SpectrumRaw, Response, *Parms->TimeTotal, Parms->SpectrumTotal);

	for (PriCount = 0; PriCount < Parms->SpectrumTotal; PriCount++) {

		PriEnergy = (PriCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->SEnergyMean;
		if (PriEnergy < Parms->EnergyBounds[0] || PriEnergy > Parms->EnergyBounds[1]) continue;

		for (SecCount = 0; SecCount < Parms->SpectrumTotal; SecCount++) {
			SecEnergy = (SecCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->SEnergyMean;
			if (SecEnergy < Parms->EnergyBounds[0] || SecEnergy > Parms->EnergyBounds[1]) continue;

			SpectralValue = Rescaling * SpectrumRaw[(Rephasing * Parms->SpectrumTotal + Minus1Pow(Rephasing) * PriCount) * Parms->SpectrumTotal + SecCount];
			Spectrum[PriCount * Parms->SpectrumTotal + SecCount] = SpectralValue;
			TotalSpectrum[PriCount * Parms->SpectrumTotal + SecCount] += SpectralValue;
			PriSubSpectrum[PriCount * Parms->SpectrumTotal + SecCount] += SpectralValue;
			SecSubSpectrum[PriCount * Parms->SpectrumTotal + SecCount] += SpectralValue;
		}
	}

	sprintf(FileName, "%s%s%s", "Spectrum", IdChar, ".dat");
	SaveSpectrum(Spectrum, FileName, Parms);

	free(Response); free(SpectrumRaw); free(Spectrum);
}


static void FourierTransform2D(ComplexType *Output, const ComplexType *Input, const int DimensionIn, const int DimensionOut)
{
	fftw_complex *In, *Out;
	fftw_plan Plan;
	int PriCount, SecCount;

	In = fftw_malloc(sizeof(fftw_complex) * IntPow2(DimensionOut));
	Out = fftw_malloc(sizeof(fftw_complex) * IntPow2(DimensionOut));
	Plan = fftw_plan_dft_2d(DimensionOut, DimensionOut, In, Out, FFTW_BACKWARD, FFTW_ESTIMATE);

	for (PriCount = 0; PriCount < IntPow2(DimensionOut); PriCount++) In[PriCount] = 0;
	for (PriCount = 0; PriCount < DimensionIn; PriCount++)
		for (SecCount = 0; SecCount < DimensionIn; SecCount++) In[PriCount * DimensionOut + SecCount] = Input[PriCount * DimensionIn + SecCount];

	fftw_execute(Plan);

	for (PriCount = 0; PriCount < DimensionOut; PriCount++)
		for (SecCount = 0; SecCount < DimensionOut; SecCount++)
			Output[PriCount * DimensionOut + SecCount] = (ComplexType) Out[((PriCount < (int) ceil(.5 * DimensionOut)) * (PriCount + (int) floor(.5 * DimensionOut)) + (PriCount > (int) floor(.5 * DimensionOut)) * (PriCount - (int) ceil(.5 * DimensionOut))) * DimensionOut + (SecCount < (int) ceil(.5 * DimensionOut)) * (SecCount + (int) floor(.5 * DimensionOut)) + (SecCount > (int) floor(.5 * DimensionOut)) * (SecCount - (int) ceil(.5 * DimensionOut))];

	fftw_destroy_plan(Plan);
	fftw_free(In); fftw_free(Out);
	fftw_cleanup();
}


static void ReadResponse(ComplexType *Response, const char *FileName, const ParmStruct *Parms)
{
	double ReadPriTime, ReadSecTime, ReadReal, ReadImag;
	int PriCount, SecCount, NewFormat = false;
	FILE *FileID;

	FileID = fopen(FileName, "r");
	for (SecCount = 0; SecCount < 2; SecCount++) {
		fscanf(FileID, "%lf %lf %lf %lf\n", &ReadPriTime, &ReadSecTime, &ReadImag, &ReadReal);
		if (ReadPriTime != 0 || ReadSecTime != SecCount * Parms->TimeStep) {
			NewFormat = true;
			break;
		}
	}
	fclose(FileID);

	FileID = fopen(FileName, "r");
	for (PriCount = 0; PriCount < *Parms->TimeTotal; PriCount++)
		for (SecCount = 0; SecCount < *Parms->TimeTotal; SecCount++) {
			if (NewFormat)	fscanf(FileID, "%lf %lf", &ReadImag, &ReadReal);
			else 			fscanf(FileID, "%lf %lf %lf %lf", &ReadPriTime, &ReadSecTime, &ReadImag, &ReadReal);
			Response[PriCount * *Parms->TimeTotal + SecCount] = (-ReadReal - I * ReadImag);
			if (*Parms->LifeTime) Response[PriCount * *Parms->TimeTotal + SecCount] *= exp(-Parms->TimeStep * (PriCount + SecCount) / (2 * *Parms->LifeTime));
		}
	fclose(FileID);
}


static void SaveEnergyGrid(const ParmStruct *Parms)
{
	FILE *FileID;
	RealType Energy;
	int Count;

	FileID = fopen("EnergyGrid.dat", "w");
	for (Count = 0; Count < Parms->SpectrumTotal; Count++) {
		Energy = (Count - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->SEnergyMean;
		if (Energy < Parms->EnergyBounds[0] || Energy > Parms->EnergyBounds[1]) continue;
		fprintf(FileID, "%f ", Energy);
	}
	fclose(FileID);
}


static void SaveSpectrum(const ComplexType *Spectrum, const char *FileName, const ParmStruct *Parms)
{
	FILE *FileID;
	RealType Energy;
	int PriCount, SecCount;

	FileID = fopen(FileName, "w");
	for (PriCount = 0; PriCount < Parms->SpectrumTotal; PriCount++) {
		Energy = (PriCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->SEnergyMean;
		if (Energy < Parms->EnergyBounds[0] || Energy > Parms->EnergyBounds[1]) continue;
		for (SecCount = 0; SecCount < Parms->SpectrumTotal; SecCount++) {
			Energy = (SecCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->SEnergyMean;
			if (Energy < Parms->EnergyBounds[0] || Energy > Parms->EnergyBounds[1]) continue;
			fprintf(FileID, "%f ", -creal(Spectrum[PriCount * Parms->SpectrumTotal + SecCount]));
		}
		fprintf(FileID, "\n");
	}
	for (PriCount = 0; PriCount < Parms->SpectrumTotal; PriCount++) {
		Energy = (PriCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->SEnergyMean;
		if (Energy < Parms->EnergyBounds[0] || Energy > Parms->EnergyBounds[1]) continue;
		for (SecCount = 0; SecCount < Parms->SpectrumTotal; SecCount++) {
			Energy = (SecCount - Parms->SpectrumTotal / 2) / (Parms->TimeStep * Parms->SpectrumTotal * LightSpeed) + Parms->SEnergyMean;
			if (Energy < Parms->EnergyBounds[0] || Energy > Parms->EnergyBounds[1]) continue;
			fprintf(FileID, "%f ", -cimag(Spectrum[PriCount * Parms->SpectrumTotal + SecCount]));
		}
		fprintf(FileID, "\n");
	}
	fclose(FileID);
}