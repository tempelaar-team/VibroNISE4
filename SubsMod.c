#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "RandomMod.h"
#include "ToolsMod.h"
#include "BasisMod.h"
#include "SubsMod.h"

#include <string.h>
#include <math.h>


static void SaveCrystal(const CrystalStruct *Crystal, const ParmStruct *Parms);
static void CoefIncrement(ComplexType *CoefIncrement, const ComplexType *Coef, const MaskRArrayStruct Hamiltonian, const int BasisTotal);
static void CoefIncrement_Redfield(ComplexType *CoefIncrement, const ComplexType *Coef, const RealType *KetEigenValues, const RealType *BraEigenValues, const MaskCArrayStruct RedfieldTensor, const int KetBasisTotal, const int BraBasisTotal, const int KetCutOff, const int BraCutOff, const RealType TimeStepAsInvEnergy);

static RealType Distance(const RealType ADistance, const RealType BDistance, const RealType TwoCosABAngle);
static int VSum(const BasisStruct BasisState);
static RealType DiagonalEnergy(const BasisStruct *BasisState, const CouplingsStruct *Couplings, const ParmStruct *Parms);
static RealType LineStrength(const int GStateId, const int PolId, const MaskRArrayStruct *GSDipoles,  const RealType *PermDipoles, const RealType *SEigenVector, const ParmStruct *Parms);
static RealType TransitionDipole(const int UnitId, const int GStateId, const MaskRArrayStruct *GSDipoles, const RealType *SEigenVector, const ParmStruct *Parms);
static int ClubCountCycle(int *PriClubCount, int *SecClubCount, const int *BasisBound);
static int ClubCountInc(int *ClubCount, const int *BasisBound, const int Limit);


/////////////
// Globals //
/////////////


void InitProgress(ProgressStruct *Progress, const PolStruct *Pols, const ParmStruct *Parms)
{
	(*Progress).Timer = SetTimer();
	(*Progress).Count = 0;
	(*Progress).Total = 100;
	(*Progress).RefCount = 0;
	if (Parms->Mode == ModeDynamics)
		(*Progress).RefTotal = Parms->SampleTotal * *Parms->TimeTotal;
	else if (Parms->Mode == ModeResponse1D)
		(*Progress).RefTotal = *Parms->TimeTotal * SumIArray(Pols->Active, 3);
	else if (Parms->Mode == ModeResponse2D) {
		if (Parms->DoubleCrossPol)	(*Progress).RefTotal = *Parms->TimeTotal * SumIArray(Pols->Active, 3) * 4;
		else						(*Progress).RefTotal = *Parms->TimeTotal * (1 + 7 * (SumIArray(Pols->Active, 3) == 2) + 20 * (SumIArray(Pols->Active, 3) == 3));
	}
}


void UpdateProgress(ProgressStruct *Progress)
{
	char DisplayNumber[3];

	if ((RealType) ((*Progress).Count + 1) / (RealType) (*Progress).Total < (RealType) ++(*Progress).RefCount / (RealType) (*Progress).RefTotal) {
		if ((*Progress).Count == 0 && (abs((int) (SetTimer() - (*Progress).Timer)) < 1)) (*Progress).Total /= 10;
		else {
			(*Progress).Count++;
			IntToString(DisplayNumber, 100 * (*Progress).Count / (*Progress).Total);
			PrintMessage("Progress ", DisplayNumber, "%.");
		}
	}
}


void BuildCrystal(CrystalStruct **Crystal, const ParmStruct *Parms)
{
	int ACount, BCount, SubCount, SiteCount = 0;
	
	*Crystal = malloc(Parms->SiteTotal * sizeof(CrystalStruct));
	
	for (ACount = 0; ACount < Parms->ACellTotal; ACount++)
		for (BCount = 0; BCount < Parms->BCellTotal; BCount++)
			for (SubCount = 0; SubCount < Parms->SubTotal; SubCount++) {
				(*Crystal)[SiteCount].TwoTimesPosA = 2 * ACount + (SubCount == 1);
				(*Crystal)[SiteCount].TwoTimesPosB = 2 * BCount + (SubCount == 1);
				(*Crystal)[SiteCount].CellA = ACount; (*Crystal)[SiteCount].CellB = BCount;
				(*Crystal)[SiteCount].Sub = SubCount;
				SiteCount++;
			}

	if (!Embarrassing) SaveCrystal(*Crystal, Parms);
}


void BuildDCrystal(DCrystalStruct **DCrystal, RealType **Distances, const CrystalStruct *Crystal, const ParmStruct *Parms)
{
	RealType MinDist, Dist, HalfA = (RealType) .5 * Parms->AConstant, HalfB = (RealType) .5 * Parms->BConstant, TwoCosABAngle = (RealType) (2.0 * cos(Parms->ABAngle));
	int PriCount, SecCount, TransCount, SmartCount, TwoTimesDPA, TwoTimesDPB;
	
	int Translations[18] = { 0, 0, Parms->ACellTotal, 0, Parms->ACellTotal, Parms->BCellTotal, 0, Parms->BCellTotal, -1 * Parms->ACellTotal, Parms->BCellTotal, -1 * Parms->ACellTotal, 0, -1 * Parms->ACellTotal, -1 * Parms->BCellTotal, 0, -1 * Parms->BCellTotal, Parms->ACellTotal, -1 * Parms->BCellTotal };
	
	*DCrystal = calloc((size_t) IntPow2(Parms->SiteTotal), sizeof(DCrystalStruct));
	*Distances = calloc((size_t) 2 * IntPow2(Parms->SiteTotal), sizeof(RealType));
	
	for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++) {
		for (SecCount = 0; SecCount < Parms->SiteTotal; SecCount++) {
			SmartCount = PriCount * Parms->SiteTotal + SecCount;
			MinDist = (RealType) (1.0 / 0.0);
			for (TransCount = 0; TransCount < 9; TransCount++) {
				TwoTimesDPA = Crystal[SecCount].TwoTimesPosA + 2 * Translations[2 * TransCount] - Crystal[PriCount].TwoTimesPosA;
				TwoTimesDPB = Crystal[SecCount].TwoTimesPosB + 2 * Translations[2 * TransCount + 1] - Crystal[PriCount].TwoTimesPosB;
				Dist = Distance(TwoTimesDPA * HalfA, TwoTimesDPB * HalfB, TwoCosABAngle);
				if (Dist < MinDist) {
					MinDist = Dist;
					(*DCrystal)[SmartCount].TwoTimesDPosA = TwoTimesDPA;
					(*DCrystal)[SmartCount].TwoTimesDPosB = TwoTimesDPB;
					(*DCrystal)[SmartCount].DCellA = Crystal[SecCount].CellA + Translations[2 * TransCount] - Crystal[PriCount].CellA;
					(*DCrystal)[SmartCount].DCellB = Crystal[SecCount].CellB + Translations[2 * TransCount + 1] - Crystal[PriCount].CellB;
				}
				if (!Parms->Periodic) break;
			}
			(*Distances)[SmartCount] = MinDist;
			if ((*DCrystal)[SmartCount].TwoTimesDPosA >= 0 && (*DCrystal)[SmartCount].TwoTimesDPosA <= Parms->ACellTotal && (*DCrystal)[SmartCount].TwoTimesDPosB >= 0 && (*DCrystal)[SmartCount].TwoTimesDPosB <= Parms->BCellTotal)
				(*Distances)[IntPow2(Parms->SiteTotal) + SmartCount] = true;
		}
	}
}


void BuildTransDistances(int **TransDistances, const CrystalStruct *Crystal, const ParmStruct *Parms)
{
	int PriCount, SecCount, ScanCount, TwoTimesDistanceA, TwoTimesDistanceB;
	int TwoTimesACellTotal = 2 * Parms->ACellTotal, TwoTimesBCellTotal = 2 * Parms->BCellTotal;
	const CrystalStruct *PriCrystalElement, *SecCrystalElement, *ScanCrystalElement;
	
	*TransDistances = calloc((size_t) IntPow2(Parms->SiteTotal), sizeof(int));
	
	PriCrystalElement = Crystal;
	for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++) {
		
		SecCrystalElement = Crystal;
		for (SecCount = 0; SecCount < Parms->SiteTotal; SecCount++) {
			TwoTimesDistanceA = (*SecCrystalElement).TwoTimesPosA - (*PriCrystalElement).TwoTimesPosA;
			TwoTimesDistanceB = (*SecCrystalElement++).TwoTimesPosB - (*PriCrystalElement).TwoTimesPosB;
			if (TwoTimesDistanceA < 0) TwoTimesDistanceA += TwoTimesACellTotal;
			if (TwoTimesDistanceB < 0) TwoTimesDistanceB += TwoTimesBCellTotal;
			
			ScanCrystalElement = Crystal;
			for (ScanCount = 0; ScanCount < Parms->SiteTotal; ScanCount++) {
				if ((*ScanCrystalElement).TwoTimesPosA == TwoTimesDistanceA && (*ScanCrystalElement).TwoTimesPosB == TwoTimesDistanceB)
					(*TransDistances)[PriCount * Parms->SiteTotal + SecCount] = ScanCount;
				ScanCrystalElement++;
			}
		}
		PriCrystalElement++;
	}
}


void BuildPermDipoles(RealType **PermDipoles, PolStruct *Pols, const CrystalStruct *Crystal, const ParmStruct *Parms)
{
	int SiteCount, SubCount, PolCount, STCount;
	double ReadDouble;
	FILE *FileID;
	
	*PermDipoles = calloc((size_t) 6 * Parms->SiteTotal, sizeof(RealType));

	if (strncmp(Parms->DipolesFile, "|", 1) == 0) {
		Pols->Active[0] = true;
		for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++)
			(*PermDipoles)[SiteCount] = 1;
	}
	else if (Parms->DDCouplingsFormat == DDCouplingsFormatOld) {
		FileID = fopen(Parms->DipolesFile, "r");
		for (STCount = 0; STCount < 2; STCount++) {
			for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++)
				for (PolCount = 0; PolCount < 3; PolCount++) {
					fscanf(FileID, "%lf ", &ReadDouble);
					(*PermDipoles)[(PolCount * 2 + STCount) * Parms->SiteTotal + SiteCount] = (RealType) ReadDouble;
					if (fabs(ReadDouble) > 0) Pols->Active[PolCount] = true;
				}
			if (feof(FileID)) break;
		}
		fclose(FileID);
	}
	else {
		FileID = fopen(Parms->DipolesFile, "r");
		for (STCount = 0; STCount < 2; STCount++) {
			for (SubCount = 0; SubCount < Parms->SubTotal; SubCount++)
				for (PolCount = 0; PolCount < 3; PolCount++) {
					fscanf(FileID, "%lf", &ReadDouble);
					if (ReadDouble) {
						Pols->Active[PolCount] = true;
						for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++)
							if (Crystal[SiteCount].Sub == SubCount)
								(*PermDipoles)[(PolCount * 2 + STCount) * Parms->SiteTotal + SiteCount] = (RealType) ReadDouble;
					}
				}
			if (feof(FileID)) break;
		}
		fclose(FileID);
	}
}


void BuildDDCouplings(MaskRArrayStruct *DDCouplings, const CrystalStruct *Crystal, const DCrystalStruct *DCrystal, const ParmStruct *Parms)
{
	RealType *List;
	int *ListId, *ListIdCompare;
	double ReadDouble;
	int PriCount, SecCount, SmartCount, ListCount = 0, ListTotal, Dummy;
	FILE *FileID;

	if (Parms->DDCouplingsFormat == DDCouplingsFormatNick && strncmp(Parms->DDCouplingsFile, "|", 1) != 0) {

		AllocateMaskArray(DDCouplings, IntPow2(Parms->SiteTotal), true, true);
		List = malloc(10000 * sizeof(RealType));
		ListId = malloc(4 * 10000 * sizeof(int));
		ListIdCompare = malloc(4 * sizeof(int));

		FileID = fopen(Parms->DDCouplingsFile, "r");
		while (!feof(FileID) && ListCount < 10000) {
			SmartCount = 4 * ListCount;
			fscanf(FileID, "%d %d %d %d %d %lf\n", &ListId[SmartCount], &ListId[SmartCount + 1],
				   &ListId[SmartCount + 2], &ListId[SmartCount + 3], &Dummy, &ReadDouble);
			ListId[SmartCount]--;
			ListId[SmartCount + 1]--; // Convert from Nick's convention to mine.
			List[ListCount] = (RealType) ReadDouble;
			ListCount++;
		}
		ListTotal = ListCount;
		fclose(FileID);

		for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++)
			for (SecCount = 0; SecCount < PriCount; SecCount++) {

				SmartCount = PriCount * Parms->SiteTotal + SecCount;
				ListIdCompare[0] = Crystal[PriCount].Sub;
				ListIdCompare[1] = Crystal[SecCount].Sub;
				ListIdCompare[2] = DCrystal[SmartCount].DCellA;
				ListIdCompare[3] = DCrystal[SmartCount].DCellB;

				for (ListCount = 0; ListCount < ListTotal; ListCount++)
					if (CompareIArray(&ListId[4 * ListCount], ListIdCompare, 4)) {
						(*DDCouplings).Array[SmartCount] = List[ListCount];
						(*DDCouplings).Mask[SmartCount] = true;
						break;
					}
			}

		MirrorRMatrix((*DDCouplings).Array, Parms->SiteTotal);
		MirrorIMatrix((*DDCouplings).Mask, Parms->SiteTotal);

		free(List);
		free(ListId);
		free(ListIdCompare);
	}
	else if (Parms->DDCouplingsFormat == DDCouplingsFormatOld && strncmp(Parms->DDCouplingsFile, "|", 1) != 0) {

		AllocateMaskArray(DDCouplings, IntPow2(Parms->SiteTotal), false, false);

		FileID = fopen(Parms->DDCouplingsFile, "r");
		for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++) {
			for (SecCount = 0; SecCount < Parms->SiteTotal; SecCount++) {
				fscanf(FileID, "%lf ", &ReadDouble);
				(*DDCouplings).Array[PriCount * Parms->SiteTotal + SecCount] = (RealType) ReadDouble;
				(*DDCouplings).Mask[PriCount * Parms->SiteTotal + SecCount] = true;
			}
			(*DDCouplings).Array[PriCount * (Parms->SiteTotal + 1)] -= Parms->SEnergyMean;
		}
		fclose(FileID);

	}
	else BuildCouplings(DDCouplings, Parms->DDCouplingsFile, Parms->DDCouplingNN, Crystal, DCrystal, Parms);

	if (Parms->OptDielectric != 1)
		for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++)
			for (SecCount = 0; SecCount < Parms->SiteTotal; SecCount++)
				if (PriCount != SecCount) (*DDCouplings).Array[PriCount * Parms->SiteTotal + SecCount] /= Parms->OptDielectric;
}


void BuildCouplings(MaskRArrayStruct *Couplings, const char *FileName, const RealType CouplingNN, const CrystalStruct *Crystal, const DCrystalStruct *DCrystal, const ParmStruct *Parms)
{
	RealType *List;
	int *ListId, *ListIdCompare;
	double ReadDouble;
	int PriCount, SecCount, SmartCount, ListCount = 0, ListInversion = true, ListTotal;
	FILE *FileID;

	AllocateMaskArray(Couplings, IntPow2(Parms->SiteTotal), true, true);

	if (strncmp(FileName, "|", 1) == 0) {
		for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++)
			for (SecCount = 0; SecCount < PriCount; SecCount++) {
				if (Crystal[PriCount].Sub != Crystal[SecCount].Sub) continue;
				SmartCount = PriCount * Parms->SiteTotal + SecCount;
				if ((DCrystal[SmartCount].DCellA == 0 && abs(DCrystal[SmartCount].DCellB) == 1) || (DCrystal[SmartCount].DCellB == 0 && abs(DCrystal[SmartCount].DCellA) == 1)) {
					(*Couplings).Array[SmartCount] = CouplingNN;
					(*Couplings).Mask[SmartCount] = true;
				}
			}
		MirrorRMatrix((*Couplings).Array, Parms->SiteTotal);
		MirrorIMatrix((*Couplings).Mask, Parms->SiteTotal);
	}
	else {

		List = malloc(10000 * sizeof(RealType));
		ListId = malloc(3 * 10000 * sizeof(int));
		ListIdCompare = malloc(3 * sizeof(int));

		FileID = fopen(FileName, "r");
		while (!feof(FileID) && ListCount < 10000) {
			fscanf(FileID, "%lf", &ReadDouble);
			ListId[3 * ListCount] = (int) ReadDouble;
			fscanf(FileID, "%lf", &ReadDouble);
			if (ListCount == 0 && ReadDouble < 0) ListInversion = false;
			ListId[3 * ListCount + 1] = (int) (2 * ReadDouble); // Doubling for conversion to integer.
			fscanf(FileID, "%lf", &ReadDouble);
			ListId[3 * ListCount + 2] = (int) (2 * ReadDouble); // Idem.
			fscanf(FileID, "%lf\n", &ReadDouble);
			List[ListCount] = (RealType) ReadDouble;
			ListCount++;
		}
		ListTotal = ListCount;
		fclose(FileID);

		for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++)
			for (SecCount = 0; SecCount < Parms->SiteTotal; SecCount++) {

				SmartCount = PriCount * Parms->SiteTotal + SecCount;
				ListIdCompare[0] = Crystal[PriCount].Sub;
				ListIdCompare[1] = DCrystal[SmartCount].TwoTimesDPosA;
				ListIdCompare[2] = DCrystal[SmartCount].TwoTimesDPosB;
				for (ListCount = 0; ListCount < ListTotal; ListCount++)
					if (CompareIArray(&ListId[3 * ListCount], ListIdCompare, 3)) {
						(*Couplings).Array[SmartCount] = List[ListCount];
						(*Couplings).Mask[SmartCount] = true;
						SmartCount = -1;
						break;
					}

				if (SmartCount == -1 || !ListInversion) continue;

				ListIdCompare[1] *= -1;
				ListIdCompare[2] *= -1;
				for (ListCount = 0; ListCount < ListTotal; ListCount++)
					if (CompareIArray(&ListId[3 * ListCount], ListIdCompare, 3)) {
						(*Couplings).Array[SmartCount] = List[ListCount];
						(*Couplings).Mask[SmartCount] = true;
						break;
					}
			}

		free(List); free(ListId); free(ListIdCompare);
	}
}


void BuildEHEnergies(RealType **EHEnergies, const DCrystalStruct *DCrystal, const RealType *Distances, ParmStruct *Parms)
{
	RealType *List;
	int *ListId, *ListIdCompare;
	double ReadDouble;
	int PriCount, SecCount, SmartCount, ListCount = 0, ListTotal;
	FILE *FileID;

	*EHEnergies = calloc((size_t) IntPow2(Parms->SiteTotal), sizeof(RealType));
	
	List = malloc(10000 * sizeof(RealType));
	ListId = malloc(2 * 10000 * sizeof(int));
	ListIdCompare = malloc(2 * sizeof(int));
	
	if (strncmp(Parms->EHEnergiesFile, "|", 1) != 0) {
		
		FileID = fopen(Parms->EHEnergiesFile, "r");
		while (!feof(FileID) && ListCount < 10000) {
			fscanf(FileID, "%lf", &ReadDouble);
			ListId[2 * ListCount] = (int) (2 * ReadDouble);		// Doubling for conversion to integer.
			fscanf(FileID, "%lf", &ReadDouble);
			ListId[2 * ListCount + 1] = (int) (2 * ReadDouble); // Idem.
			fscanf(FileID, "%lf\n", &ReadDouble);
			List[ListCount] = (RealType) ReadDouble;
			ListCount++;
		}
		fclose(FileID);
		
	}
	ListTotal = ListCount;

	for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++)
		for (SecCount = 0; SecCount < Parms->SiteTotal; SecCount++) {
			
			if (PriCount == SecCount || (Parms->EHRadius >= 0 && Distances[PriCount * Parms->SiteTotal + SecCount] > Parms->EHRadius)) continue;

			SmartCount = PriCount * Parms->SiteTotal + SecCount;
			
			ListIdCompare[0] = DCrystal[SmartCount].TwoTimesDPosA;
			ListIdCompare[1] = DCrystal[SmartCount].TwoTimesDPosB;
			for (ListCount = 0; ListCount < ListTotal; ListCount++)
				if (CompareIArray(&ListId[2 * ListCount], ListIdCompare, 2)) {
					(*EHEnergies)[SmartCount] = List[ListCount];
					SmartCount = -1;
					break;
				}
			if (SmartCount == -1) continue;
			
			ListIdCompare[0] *= -1;
			ListIdCompare[1] *= -1;
			for (ListCount = 0; ListCount < ListTotal; ListCount++)
				if (CompareIArray(&ListId[2 * ListCount], ListIdCompare, 2)) {
					(*EHEnergies)[SmartCount] = List[ListCount];
					SmartCount = -1;
					break;
				}
			if (SmartCount == -1) continue;

			(*EHEnergies)[SmartCount] = Parms->EHSeparated - Parms->EHScaling / Distances[SmartCount];
		}

	free(List); free(ListId); free(ListIdCompare);
}


void BuildK0Vectors(int **K0Vectors, const int *TransDistances, const CrystalStruct *Crystal, const BasisStruct *K0Basis, const BasisStruct *SBasis, const ParmStruct *Parms)
{
	BasisStruct K0BasisState, SBasisState;
	int K0BasisCount, SBasisCount, ClubCount, Particles, PreFactorId;

	*K0Vectors = calloc((size_t) *Parms->K0BasisBound * *Parms->SBasisBound, sizeof(int));

	for (ClubCount = 1; ClubCount < ClubTotal - 1; ClubCount++) {
		if (Parms->K0BasisBound[ClubCount] == *Parms->K0BasisBound) break;
		Particles = CountParticles(&K0Basis[Parms->K0BasisBound[ClubCount]]);

		if (Particles > 2) continue;

		for (K0BasisCount = Parms->K0BasisBound[ClubCount]; K0BasisCount < Parms->K0BasisBound[ClubCount + 1]; K0BasisCount++) {
			memcpy(&K0BasisState, &K0Basis[K0BasisCount], sizeof(BasisStruct));

			for (SBasisCount = Parms->SBasisBound[ClubCount]; SBasisCount < Parms->SBasisBound[ClubCount + 1]; SBasisCount++) {
				memcpy(&SBasisState, &SBasis[SBasisCount], sizeof(BasisStruct));

				PreFactorId = 1;

				if (SBasisState.P[0].C == Electron || SBasisState.P[0].C == Ground) SwitchFirstTwoParticles(&SBasisState);
				else if (SBasisState.P[0].C == Triplet) {
					if (Crystal[SBasisState.P[0].S].Sub != Crystal[SBasisState.P[1].S].Sub) {
						if (Crystal[SBasisState.P[0].S].Sub == 1) SwitchFirstTwoParticles(&SBasisState);
					}
					else if (TransDistances[SBasisState.P[0].S * Parms->SiteTotal + SBasisState.P[1].S] > TransDistances[SBasisState.P[1].S * Parms->SiteTotal + SBasisState.P[0].S]) SwitchFirstTwoParticles(&SBasisState);
					else if (TransDistances[SBasisState.P[0].S * Parms->SiteTotal + SBasisState.P[1].S] == TransDistances[SBasisState.P[1].S * Parms->SiteTotal + SBasisState.P[0].S]) PreFactorId = 2;
				}

				if (K0BasisState.P[0].S == Crystal[SBasisState.P[0].S].Sub && IdentifyFirstTwoVibrations(&K0BasisState, &SBasisState)) {
					if (Particles == 1)
						(*K0Vectors)[K0BasisCount * *Parms->SBasisBound + SBasisCount] = PreFactorId;
					else if (TransDistances[SBasisState.P[0].S * Parms->SiteTotal + SBasisState.P[1].S] == K0BasisState.P[1].S)
						(*K0Vectors)[K0BasisCount * *Parms->SBasisBound + SBasisCount] = PreFactorId;
				}
			}
		}
	}
}


void BuildEnergyGrid(RealType **EnergyGrid, const ParmStruct *Parms)
{
	int Count;
	
	*EnergyGrid = malloc(Parms->SpectrumTotal * sizeof(RealType));
	
	for (Count = 0; Count < Parms->SpectrumTotal; Count++) (*EnergyGrid)[Count] = Parms->EnergyBounds[0] + (Parms->EnergyBounds[1] - Parms->EnergyBounds[0]) * Count / (Parms->SpectrumTotal - 1) - Parms->SEnergyMean;
}


void BuildDipoles(const char *FileName, MaskRArrayStruct *Dipoles, const BasisStruct *LoBasis, const BasisStruct *HiBasis, const int LoBasisTotal, const int HiBasisTotal, RealType * const *VOverlaps, const ParmStruct *Parms)
{
	BasisStruct HiBasisState, LoBasisState;
	RealType EffectiveVOverlap;
	int LoBasisCount, HiBasisCount, SiteCount;

	PrintMessage("Constructing ", FileName, ".");

	AllocateMaskArray(Dipoles, HiBasisTotal * LoBasisTotal, true, false);
	memset((*Dipoles).Mask, -1, HiBasisTotal * LoBasisTotal * sizeof(int));

	for (HiBasisCount = 0; HiBasisCount < HiBasisTotal; HiBasisCount++)
		for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++) {

			memcpy(&HiBasisState, &HiBasis[HiBasisCount], sizeof(BasisStruct));
			if (AnnihilateS(&HiBasisState, SiteCount)) {

				memcpy(&LoBasisState, &HiBasisState, sizeof(BasisStruct));
				for (LoBasisCount = 0; LoBasisCount < LoBasisTotal; LoBasisCount++)
					if (IdentifyBasisState(&EffectiveVOverlap, &LoBasisState, &LoBasis[LoBasisCount], VOverlaps, Parms->VTotal)) {
						(*Dipoles).Array[LoBasisCount * HiBasisTotal + HiBasisCount] = EffectiveVOverlap;
						(*Dipoles).Mask[LoBasisCount * HiBasisTotal + HiBasisCount] = SiteCount;
					}
			}

			memcpy(&HiBasisState, &HiBasis[HiBasisCount], sizeof(BasisStruct));
			if (AnnihilateT_n(&HiBasisState, SiteCount)) {
				memcpy(&LoBasisState, &HiBasisState, sizeof(BasisStruct));
				for (LoBasisCount = 0; LoBasisCount < LoBasisTotal; LoBasisCount++)
					if (IdentifyBasisState(&EffectiveVOverlap, &LoBasisState, &LoBasis[LoBasisCount], VOverlaps, Parms->VTotal)) {
						(*Dipoles).Array[LoBasisCount * HiBasisTotal + HiBasisCount] = Parms->TT_nDipole * EffectiveVOverlap;
						(*Dipoles).Mask[LoBasisCount * HiBasisTotal + HiBasisCount] = Parms->SiteTotal + SiteCount;
					}
			}
		}

	SaveOperator(FileName, (*Dipoles).Array, LoBasis, HiBasis, LoBasisTotal, HiBasisTotal);
}


void BuildGEigenValues(RealType **GEigenValues, const BasisStruct *GBasis, const ParmStruct *Parms)
{
	RealType *GEigenEnergiesElement;
	int Count;

	*GEigenValues = malloc(*Parms->GBasisBound * sizeof(RealType));

	GEigenEnergiesElement = *GEigenValues;

	for (Count = 0; Count < *Parms->GBasisBound; Count++) (*GEigenEnergiesElement++) = Parms->VEnergy * VSum(GBasis[Count]);
}


void BuildGPropagator(ComplexType **GPropagator, const BasisStruct *GBasis, const ParmStruct *Parms)
{
	RealType TimeStepAsInvEnergy = Parms->EnergyToFreq * Parms->TimeStep;
	ComplexType *GPropagatorElement;
	int Count;

	*GPropagator = malloc(*Parms->GBasisBound * sizeof(ComplexType));

	GPropagatorElement = *GPropagator;
	for (Count = 0; Count < *Parms->GBasisBound; Count++) (*GPropagatorElement++) = cexp(Parms->VEnergy * VSum(GBasis[Count]) * TimeStepAsInvEnergy * -I);
}


void BuildPropagator(ComplexType **Propagator, const RealType *EigenValues, const RealType TimeStepAsInvEnergy, const int BasisTotal)
{
	ComplexType *PropagatorElement;
	int Count;

	*Propagator = malloc(BasisTotal * sizeof(ComplexType));

	PropagatorElement = *Propagator;
	for (Count = 0; Count < BasisTotal; Count++) (*PropagatorElement++) = cexp(EigenValues[Count] * TimeStepAsInvEnergy * -I);
}


void BuildHamiltonian(const char *FileName, MaskRArrayStruct *Hamiltonian, RealType **ToDiag, RealType **EigenValues, int *ToDiagTotal, const CouplingsStruct *Couplings, RealType *const *VOverlaps, const BasisStruct *Basis, const int *BasisBound, const ParmStruct *Parms)
{
	BasisStruct PriBasisState, SecBasisState;
	RealType EffectiveVOverlap;
	int PriBasisCount, SecBasisCount, PriSiteId, SecSiteCount, PriClubCount, SecClubCount, PriPCount, SecPCount;

	PrintTime(); PrintMessage("Constructing ", FileName, ".");

	AllocateMaskArray(Hamiltonian, IntPow2(*BasisBound), true, true);
	if (EigenValues != NULL) *EigenValues = malloc(*BasisBound * sizeof(RealType));

	for (PriBasisCount = 0; PriBasisCount < *BasisBound; PriBasisCount++) {
		(*Hamiltonian).Array[PriBasisCount * (*BasisBound + 1)] = DiagonalEnergy(&Basis[PriBasisCount], Couplings, Parms);
		(*Hamiltonian).Mask[PriBasisCount * (*BasisBound + 1)] = true;
	}

	PriClubCount = 0; SecClubCount = 0;
	while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {
		if (!Present(&Basis[BasisBound[PriClubCount]], Singlet) || !Present(&Basis[BasisBound[SecClubCount]], Singlet)) continue;

		for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++)
			for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++)
				if (Basis[PriBasisCount].P[PriPCount].C == Singlet) {

					memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
					PriBasisState.P[PriPCount].C = Ground;
					PriSiteId = PriBasisState.P[PriPCount].S;

					for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
						if (!(*Couplings).DD.Mask[PriSiteId * Parms->SiteTotal + SecSiteCount] || PriSiteId == SecSiteCount) continue;

						memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
						if (!CreateS(&SecBasisState, SecSiteCount)) continue;

						for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < IntMin(PriBasisCount, BasisBound[SecClubCount + 1]); SecBasisCount++)
							if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
								(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).DD.Array[PriSiteId * Parms->SiteTotal + SecSiteCount] * EffectiveVOverlap;
								(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
							}
					}
				}
	}

	PriClubCount = 0; SecClubCount = 0;
	while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {
		if (!Present(&Basis[BasisBound[PriClubCount]], Triplet_n) || !Present(&Basis[BasisBound[SecClubCount]], Triplet_n)) continue;

		for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++)
			for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++)
				if (Basis[PriBasisCount].P[PriPCount].C == Triplet_n) {

					memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
					PriBasisState.P[PriPCount].C = Triplet;
					PriSiteId = PriBasisState.P[PriPCount].S;

					for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
						if (!(*Couplings).TripletDD.Mask[PriSiteId * Parms->SiteTotal + SecSiteCount] || PriSiteId == SecSiteCount) continue;

						memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
						if (!DipoleCreateT_n(&SecBasisState, SecSiteCount)) continue;

						for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < IntMin(PriBasisCount, BasisBound[SecClubCount + 1]); SecBasisCount++)
							if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
								(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).TripletDD.Array[PriSiteId * Parms->SiteTotal + SecSiteCount] * EffectiveVOverlap;
								(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
							}
					}
				}
	}


	PriClubCount = 0; SecClubCount = 0;
	while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {

		if (Parms->NoCTCTCoupling && !Present(&Basis[BasisBound[SecClubCount]], Singlet)) continue;
		if (!Present(&Basis[BasisBound[PriClubCount]], Electron) && !Present(&Basis[BasisBound[SecClubCount]], Electron)) continue;

		for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++)
			for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++)
				if (Basis[PriBasisCount].P[PriPCount].C == Singlet || Basis[PriBasisCount].P[PriPCount].C == Electron) {

					memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
					if (PriBasisState.P[PriPCount].C == Electron)	PriBasisState.P[PriPCount].C = Ground;
					else											PriBasisState.P[PriPCount].C = Hole;
					PriSiteId = PriBasisState.P[PriPCount].S;

					for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
						if (!(*Couplings).LL.Mask[PriSiteId * Parms->SiteTotal + SecSiteCount] || PriSiteId == SecSiteCount) continue;

						memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
						if (!CreateE(&SecBasisState, SecSiteCount)) continue;

						for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
							if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
								(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).LL.Array[PriSiteId * Parms->SiteTotal + SecSiteCount] * EffectiveVOverlap;
								(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
							}
					}
				}
	}

	PriClubCount = 0; SecClubCount = 0;
	while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {

		if (Parms->NoCTCTCoupling && !Present(&Basis[BasisBound[SecClubCount]], Singlet)) continue;
		if (!Present(&Basis[BasisBound[PriClubCount]], Hole) && !Present(&Basis[BasisBound[SecClubCount]], Hole)) continue;

		for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++)
			for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++)
				if (Basis[PriBasisCount].P[PriPCount].C == Singlet || Basis[PriBasisCount].P[PriPCount].C == Hole) {

					memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
					if (PriBasisState.P[PriPCount].C == Hole)	PriBasisState.P[PriPCount].C = Ground;
					else										PriBasisState.P[PriPCount].C = Electron;
					PriSiteId = PriBasisState.P[PriPCount].S;

					for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
						if (!(*Couplings).HH.Mask[PriSiteId * Parms->SiteTotal + SecSiteCount] || PriSiteId == SecSiteCount) continue;

						memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
						if (!CreateH(&SecBasisState, SecSiteCount)) continue;

						for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
							if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
								(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).HH.Array[PriSiteId * Parms->SiteTotal + SecSiteCount] * EffectiveVOverlap;
								(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
							}
					}
				}
	}

	PriClubCount = 0; SecClubCount = 0;
	while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {

		if (Present(&Basis[BasisBound[PriClubCount]], Hole))
			if (!Present(&Basis[BasisBound[SecClubCount]], Triplet)) continue;
		else if (Present(&Basis[BasisBound[PriClubCount]], Triplet))
			if (!Present(&Basis[BasisBound[SecClubCount]], Hole)) continue;
		else continue;

		for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++) {
			memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
			for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++) {

				PriSiteId = PriBasisState.P[PriPCount].S;

				if (PriBasisState.P[PriPCount].C == Hole)
					for (SecPCount = PriPCount + 1; SecPCount < PTotal && PriBasisState.P[SecPCount].S != -1; SecPCount++)
						if (PriBasisState.P[SecPCount].C == Electron && (*Couplings).HL.Mask[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId]) {

							memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
							SecBasisState.P[PriPCount].C = Triplet;
							SecBasisState.P[SecPCount].C = Triplet;
							for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
								if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
									(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += SqrtHalfOfThree * (*Couplings).HL.Array[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId] * EffectiveVOverlap;
									(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
								}
						}

				if (PriBasisState.P[PriPCount].C == Electron)
					for (SecPCount = PriPCount + 1; SecPCount < PTotal && PriBasisState.P[SecPCount].S != -1; SecPCount++)
						if (PriBasisState.P[SecPCount].C == Hole && (*Couplings).HL.Mask[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S]) {

							memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
							SecBasisState.P[PriPCount].C = Triplet;
							SecBasisState.P[SecPCount].C = Triplet;
							for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
								if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
									(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += SqrtHalfOfThree * (*Couplings).HL.Array[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S] * EffectiveVOverlap;
									(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
								}
						}

				if (PriBasisState.P[PriPCount].C == Triplet) {
					for (SecPCount = PriPCount + 1; SecPCount < PTotal && PriBasisState.P[SecPCount].S != -1; SecPCount++)
						if (PriBasisState.P[SecPCount].C == Triplet) {

							if ((*Couplings).HL.Mask[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S]) {
								memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
								SecBasisState.P[PriPCount].C = Electron;
								SecBasisState.P[SecPCount].C = Hole;
								for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
									if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
										(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += SqrtHalfOfThree * (*Couplings).HL.Array[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S] * EffectiveVOverlap;
										(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
									}
							}

							if ((*Couplings).HL.Mask[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId]) {
								memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
								SecBasisState.P[PriPCount].C = Hole;
								SecBasisState.P[SecPCount].C = Electron;
								for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
									if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
										(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += SqrtHalfOfThree * (*Couplings).HL.Array[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId] * EffectiveVOverlap;
										(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
									}
							}
						}
				}
			}
		}
	}

	PriClubCount = 0; SecClubCount = 0;
	while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {

		if (Present(&Basis[BasisBound[PriClubCount]], Triplet) && Present(&Basis[BasisBound[SecClubCount]], Triplet))
			for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++) {

				for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++) {

					memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
					PriSiteId = PriBasisState.P[PriPCount].S;

					if (PriBasisState.P[PriPCount].C == Triplet) {
						PriBasisState.P[PriPCount].C = Ground;

						for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
							if (!(*Couplings).LLHH.Mask[PriSiteId * Parms->SiteTotal + SecSiteCount] || PriSiteId == SecSiteCount) continue;

							memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
							if (!CreateT(&SecBasisState, SecSiteCount)) continue;

							for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < IntMin(BasisBound[SecClubCount + 1], PriBasisCount); SecBasisCount++)
								if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
									(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).LLHH.Array[PriSiteId * Parms->SiteTotal + SecSiteCount] * EffectiveVOverlap;
									(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
								}
						}
					}
				}
			}

		if (Present(&Basis[BasisBound[PriClubCount]], Triplet_n) && Present(&Basis[BasisBound[SecClubCount]], Triplet_n))
			for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++) {

				for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++) {

					memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
					PriSiteId = PriBasisState.P[PriPCount].S;

					if (PriBasisState.P[PriPCount].C == Triplet_n) {
						PriBasisState.P[PriPCount].C = Ground;

						for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
							if (!(*Couplings).LnLnHH.Mask[PriSiteId * Parms->SiteTotal + SecSiteCount] || PriSiteId == SecSiteCount) continue;

							memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
							if (!CreateT_n(&SecBasisState, SecSiteCount)) continue;
							for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < IntMin(BasisBound[SecClubCount + 1], PriBasisCount); SecBasisCount++)
								if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
									(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).LnLnHH.Array[PriSiteId * Parms->SiteTotal + SecSiteCount]; // * EffectiveVOverlap;
									(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
								}
						}
					}
				}
			}

		if (Present(&Basis[BasisBound[PriClubCount]], Singlet) && Present(&Basis[BasisBound[SecClubCount]], Triplet))
			for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++) {
				memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
				for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++) {

					PriSiteId = PriBasisState.P[PriPCount].S;

					if (PriBasisState.P[PriPCount].C == Singlet) {
						PriBasisState.P[PriPCount].C = Triplet;

						for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
							if (!(*Couplings).LLHL.Mask[SecSiteCount * Parms->SiteTotal + PriSiteId] || PriSiteId == SecSiteCount) continue;

							memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
							if (!CreateT(&SecBasisState, SecSiteCount)) continue;

							for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
								if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
									(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).LLHL.Array[SecSiteCount * Parms->SiteTotal + PriSiteId] * EffectiveVOverlap;
									(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
								}
						}
					}
				}
			}

		if (!Present(&Basis[BasisBound[PriClubCount]], Triplet) || !Present(&Basis[BasisBound[SecClubCount]], Singlet)) continue;

		for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++) {
			memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
			for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++) {

				PriSiteId = PriBasisState.P[PriPCount].S;

				if (PriBasisState.P[PriPCount].C == Triplet) {
					for (SecPCount = PriPCount + 1; SecPCount < PTotal && PriBasisState.P[SecPCount].S != -1; SecPCount++)
						if (PriBasisState.P[SecPCount].C == Triplet) {

							if ((*Couplings).LLHL.Mask[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S]) {
								memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
								SecBasisState.P[PriPCount].C = Ground;
								SecBasisState.P[SecPCount].C = Singlet;
								for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
									if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
										(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).LLHL.Array[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S] * EffectiveVOverlap;
										(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
									}
							}

							if ((*Couplings).LLHL.Mask[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId]) {
								memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
								SecBasisState.P[PriPCount].C = Singlet;
								SecBasisState.P[SecPCount].C = Ground;
								for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
									if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal)) {
										(*Hamiltonian).Array[PriBasisCount * *BasisBound + SecBasisCount] += (*Couplings).LLHL.Array[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId] * EffectiveVOverlap;
										(*Hamiltonian).Mask[PriBasisCount * *BasisBound + SecBasisCount] = true;
									}
							}
						}
				}
			}
		}
	}

	if (ToDiag != NULL) *ToDiag = (*Hamiltonian).Array, *ToDiagTotal = *BasisBound;

	SaveOperator(FileName, (*Hamiltonian).Array, Basis, Basis, *BasisBound, *BasisBound);
}


void DrawOffDiagDisorder(MaskRArrayStruct *OffDiagDisorder, const MaskRArrayStruct Couplings, const RealType Sigma, const ParmStruct *Parms, const int Symmetric)
{
    int PriCount, SecCount, SmartCount;

    for (PriCount = 0; PriCount < Parms->SiteTotal; PriCount++)
        for (SecCount = 0; SecCount < (Symmetric ? PriCount : Parms->SiteTotal); SecCount++) {
            SmartCount = PriCount * Parms->SiteTotal + SecCount;
            if (Couplings.Mask[SmartCount]) (*OffDiagDisorder).Array[SmartCount] = (RealType) RandomGaussian(0, Sigma * fabs(Couplings.Array[SmartCount]));
            //if (Couplings.Mask[SmartCount]) (*OffDiagDisorder).Array[SmartCount] = (RealType) RandomGaussian(0, Sigma);
        }

    if (Symmetric) MirrorRMatrix((*OffDiagDisorder).Array, Parms->SiteTotal);
}


void AddOffDiagDisorder(MaskRArrayStruct Hamiltonian, const CouplingsStruct *OffDiagDisorder, RealType *const *VOverlaps, const BasisStruct *Basis, const int *BasisBound, const ParmStruct *Parms)
{
    BasisStruct PriBasisState, SecBasisState;
    RealType EffectiveVOverlap;
    int PriBasisCount, SecBasisCount, PriSiteId, SecSiteCount, PriClubCount, SecClubCount, PriPCount, SecPCount;

    PriClubCount = 0; SecClubCount = 0;
    while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {

        if (Parms->NoCTCTCoupling && !Present(&Basis[BasisBound[SecClubCount]], Singlet)) continue;
        if (!Present(&Basis[BasisBound[PriClubCount]], Electron) && !Present(&Basis[BasisBound[SecClubCount]], Electron)) continue;

        for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++)
            for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++)
                if (Basis[PriBasisCount].P[PriPCount].C == Singlet || Basis[PriBasisCount].P[PriPCount].C == Electron) {

                    memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
                    if (PriBasisState.P[PriPCount].C == Electron)	PriBasisState.P[PriPCount].C = Ground;
                    else											PriBasisState.P[PriPCount].C = Hole;
                    PriSiteId = PriBasisState.P[PriPCount].S;

                    for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
                        if (!(*OffDiagDisorder).LL.Mask[PriSiteId * Parms->SiteTotal + SecSiteCount] || PriSiteId == SecSiteCount) continue;

                        memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
                        if (!CreateE(&SecBasisState, SecSiteCount)) continue;

                        for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
                            if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal))
                                Hamiltonian.Array[PriBasisCount * *BasisBound + SecBasisCount] += (*OffDiagDisorder).LL.Array[PriSiteId * Parms->SiteTotal + SecSiteCount] * EffectiveVOverlap;
                    }
                }
    }

    PriClubCount = 0; SecClubCount = 0;
    while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {

        if (Parms->NoCTCTCoupling && !Present(&Basis[BasisBound[SecClubCount]], Singlet)) continue;
        if (!Present(&Basis[BasisBound[PriClubCount]], Hole) && !Present(&Basis[BasisBound[SecClubCount]], Hole)) continue;

        for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++)
            for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++)
                if (Basis[PriBasisCount].P[PriPCount].C == Singlet || Basis[PriBasisCount].P[PriPCount].C == Hole) {

                    memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
                    if (PriBasisState.P[PriPCount].C == Hole)	PriBasisState.P[PriPCount].C = Ground;
                    else										PriBasisState.P[PriPCount].C = Electron;
                    PriSiteId = PriBasisState.P[PriPCount].S;

                    for (SecSiteCount = 0; SecSiteCount < Parms->SiteTotal; SecSiteCount++) {
                        if (!(*OffDiagDisorder).HH.Mask[PriSiteId * Parms->SiteTotal + SecSiteCount] || PriSiteId == SecSiteCount) continue;

                        memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
                        if (!CreateH(&SecBasisState, SecSiteCount)) continue;

                        for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
                            if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal))
                                Hamiltonian.Array[PriBasisCount * *BasisBound + SecBasisCount] += (*OffDiagDisorder).HH.Array[PriSiteId * Parms->SiteTotal + SecSiteCount] * EffectiveVOverlap;
                    }
                }
    }

    PriClubCount = 0; SecClubCount = 0;
    while (ClubCountCycle(&PriClubCount, &SecClubCount, BasisBound)) {

        if (Present(&Basis[BasisBound[PriClubCount]], Hole))
            if (!Present(&Basis[BasisBound[SecClubCount]], Triplet)) continue;
            else if (Present(&Basis[BasisBound[PriClubCount]], Triplet))
                if (!Present(&Basis[BasisBound[SecClubCount]], Hole)) continue;
                else continue;

        for (PriBasisCount = BasisBound[PriClubCount]; PriBasisCount < BasisBound[PriClubCount + 1]; PriBasisCount++) {
            memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
            for (PriPCount = 0; PriPCount < PTotal && Basis[PriBasisCount].P[PriPCount].S != -1; PriPCount++) {

                PriSiteId = PriBasisState.P[PriPCount].S;

                if (PriBasisState.P[PriPCount].C == Hole)
                    for (SecPCount = PriPCount + 1; SecPCount < PTotal && PriBasisState.P[SecPCount].S != -1; SecPCount++)
                        if (PriBasisState.P[SecPCount].C == Electron && (*OffDiagDisorder).HL.Mask[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId]) {

                            memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
                            SecBasisState.P[PriPCount].C = Triplet;
                            SecBasisState.P[SecPCount].C = Triplet;
                            for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
                                if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal))
                                    Hamiltonian.Array[PriBasisCount * *BasisBound + SecBasisCount] += SqrtHalfOfThree * (*OffDiagDisorder).HL.Array[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId] * EffectiveVOverlap;
                        }

                if (PriBasisState.P[PriPCount].C == Electron)
                    for (SecPCount = PriPCount + 1; SecPCount < PTotal && PriBasisState.P[SecPCount].S != -1; SecPCount++)
                        if (PriBasisState.P[SecPCount].C == Hole && (*OffDiagDisorder).HL.Mask[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S]) {

                            memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
                            SecBasisState.P[PriPCount].C = Triplet;
                            SecBasisState.P[SecPCount].C = Triplet;
                            for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
                                if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal))
                                    Hamiltonian.Array[PriBasisCount * *BasisBound + SecBasisCount] += SqrtHalfOfThree * (*OffDiagDisorder).HL.Array[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S] * EffectiveVOverlap;
                        }

                if (PriBasisState.P[PriPCount].C == Triplet) {
                    for (SecPCount = PriPCount + 1; SecPCount < PTotal && PriBasisState.P[SecPCount].S != -1; SecPCount++)
                        if (PriBasisState.P[SecPCount].C == Triplet) {

                            if ((*OffDiagDisorder).HL.Mask[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S]) {
                                memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
                                SecBasisState.P[PriPCount].C = Electron;
                                SecBasisState.P[SecPCount].C = Hole;
                                for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
                                    if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal))
                                        Hamiltonian.Array[PriBasisCount * *BasisBound + SecBasisCount] += SqrtHalfOfThree * (*OffDiagDisorder).HL.Array[PriSiteId * Parms->SiteTotal + PriBasisState.P[SecPCount].S] * EffectiveVOverlap;
                            }

                            if ((*OffDiagDisorder).HL.Mask[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId]) {
                                memcpy(&SecBasisState, &PriBasisState, sizeof(BasisStruct));
                                SecBasisState.P[PriPCount].C = Hole;
                                SecBasisState.P[SecPCount].C = Electron;
                                for (SecBasisCount = BasisBound[SecClubCount]; SecBasisCount < BasisBound[SecClubCount + 1]; SecBasisCount++)
                                    if (IdentifyBasisState(&EffectiveVOverlap, &SecBasisState, &Basis[SecBasisCount], VOverlaps, Parms->VTotal))
                                        Hamiltonian.Array[PriBasisCount * *BasisBound + SecBasisCount] += SqrtHalfOfThree * (*OffDiagDisorder).HL.Array[PriBasisState.P[SecPCount].S * Parms->SiteTotal + PriSiteId] * EffectiveVOverlap;
                            }
                        }
                }
            }
        }
    }
}


void BuildK0Hamiltonian(RealType **K0Hamiltonian, RealType **ToDiag, RealType **SEigenValues, int *ToDiagTotal, MaskRArrayStruct *SHamiltonian, const int *K0Vectors, const ParmStruct *Parms)
{
	RealType *K0HamiltonianElement;
	const int *K0VectorsPriElement, *K0VectorsSecElement;
	RealType InvCellTotal = (RealType) 1.0 / Parms->CellTotal;
	RealType K0PreFactors[] = { 1.0, (RealType) sqrt(2.0), 2.0 };
	int PriClubCount, SecClubCount, PriK0Count, SecK0Count, PriSCount, SecSCount;

	PrintTime(); PrintMessage("Constructing k = 0 Hamiltonian.", "", "");

	*K0Hamiltonian = calloc((size_t) IntPow2(*Parms->K0BasisBound), sizeof(RealType));
	*SEigenValues = realloc(*SEigenValues, *Parms->K0BasisBound * sizeof(RealType));

	MirrorRMatrix((*SHamiltonian).Array, *Parms->SBasisBound);
	MirrorIMatrix((*SHamiltonian).Mask, *Parms->SBasisBound);

	for (PriClubCount = 1; PriClubCount < ClubTotal - 1; PriClubCount++) {
		if (Parms->K0BasisBound[PriClubCount] == *Parms->K0BasisBound) break;

		for (PriK0Count = Parms->K0BasisBound[PriClubCount]; PriK0Count < Parms->K0BasisBound[PriClubCount + 1]; PriK0Count++) {
			K0VectorsPriElement = &K0Vectors[PriK0Count * *Parms->SBasisBound + Parms->SBasisBound[PriClubCount]];
			for (PriSCount = Parms->SBasisBound[PriClubCount]; PriSCount < Parms->SBasisBound[PriClubCount + 1]; PriSCount++) {

				if (*K0VectorsPriElement)
					for (SecClubCount = 1; SecClubCount < ClubTotal - 1; SecClubCount++) {
						if (Parms->K0BasisBound[SecClubCount] == *Parms->K0BasisBound) break;

						K0HamiltonianElement = &(*K0Hamiltonian)[PriK0Count * *Parms->K0BasisBound + Parms->K0BasisBound[SecClubCount]];

						for (SecK0Count = Parms->K0BasisBound[SecClubCount]; SecK0Count < Parms->K0BasisBound[SecClubCount + 1]; SecK0Count++) {
							K0VectorsSecElement = &K0Vectors[SecK0Count * *Parms->SBasisBound + Parms->SBasisBound[SecClubCount]];
							for (SecSCount = Parms->SBasisBound[SecClubCount]; SecSCount < Parms->SBasisBound[SecClubCount + 1]; SecSCount++) {
								if (*K0VectorsSecElement && (*SHamiltonian).Mask[PriSCount * *Parms->SBasisBound + SecSCount])
									*K0HamiltonianElement += (*SHamiltonian).Array[PriSCount * *Parms->SBasisBound + SecSCount] * K0PreFactors[*K0VectorsPriElement + *K0VectorsSecElement - 2];
								K0VectorsSecElement++;
							}
							K0HamiltonianElement++;
						}
					}
				K0VectorsPriElement++;
			}
		}
	}

	RealTimesRArray(*K0Hamiltonian, InvCellTotal, IntPow2(*Parms->K0BasisBound));

	*ToDiag = *K0Hamiltonian;
	*ToDiagTotal = *Parms->K0BasisBound;
}


int DetermineStateCutOff(const RealType *EigenValues, const int BasisTotal, const ParmStruct *Parms)
{
	int Count, StateCutOff = BasisTotal;

	if (Parms->EnergyBounds[1])
		for (Count = 0; Count < BasisTotal; Count++)
			if (EigenValues[Count] > Parms->EnergyBounds[1] - Parms->SEnergyMean) {
				StateCutOff = Count;
				break;
			}

	return StateCutOff;
}


void BuildBathFT(ComplexType **BathFT, const RealType *EigenValues, const BasisStruct *Basis, const int BasisTotal, const int CutOff, const ParmStruct *Parms)
{
	ComplexType *BathFTElement;
	const RealType *EigenValuesPriElement, *EigenValuesSecElement;
	double ReadTime, ReadBathReal, ReadBathImag, Dummy;
	RealType TimeAsInvEnergy, TimeStep = 0;
	ComplexType BathCorrelation;
	int PriCount, SecCount;
	FILE *FileID;

	*BathFT = calloc((size_t) IntPow2(CutOff), sizeof(ComplexType));

	if (strncmp(Parms->BathFile, "|", 1) == 0) return;

	BathFTElement = *BathFT;

	FileID = fopen(Parms->BathFile, "r");
	while (!feof(FileID)) {
		fscanf(FileID, "%lf %lf %lf\n", &ReadTime, &ReadBathReal, &ReadBathImag);
		TimeAsInvEnergy = (RealType) (ReadTime * Parms->EnergyToFreq);
		if (TimeStep == 0 && ReadTime > 0) TimeStep = (RealType) ReadTime;
		BathCorrelation = (ReadBathReal + I * ReadBathImag);

		BathFTElement = *BathFT; EigenValuesPriElement = EigenValues;
		for (PriCount = 0; PriCount < CutOff; PriCount++) {
			EigenValuesSecElement = EigenValues;
			for (SecCount = 0; SecCount < CutOff; SecCount++)
				(*BathFTElement++) += BathCorrelation * cexp(I * ((*EigenValuesSecElement++) - *EigenValuesPriElement) * TimeAsInvEnergy);
			EigenValuesPriElement++;
		}
	}
	fclose(FileID);

	RealTimesCArray(*BathFT, Parms->BathReorganization * TimeStep, IntPow2(CutOff));
}


void BuildReperator(RealType **Reperator, const RealType *EigenVectors, const BasisStruct *Basis, const int BasisTotal, const int CutOff, const ParmStruct *Parms)
{
	BasisStruct PriBasisState;
	const RealType *PriEigenVectorsElement, *SecEigenVectorsElement, *TerEigenVectorsElement, *QuaEigenVectorsElement;
	const RealType *PriEigenVectorsInit = EigenVectors, *SecEigenVectorsInit;
	RealType *ReperatorElement;
	RealType PriSecEigenVectors, PriSecTerEigenVectors;
	int PriCount, SecCount, TerCount, QuaCount, PriBasisCount, SecBasisCount;

	if (strncmp(Parms->BathFile, "|", 1) == 0) {
		*Reperator = malloc(sizeof(RealType));
		return;
	}
	else *Reperator = calloc((size_t) IntPowN(CutOff, 4), sizeof(RealType));

	for (PriBasisCount = 0; PriBasisCount < BasisTotal; PriBasisCount++) {
		SecEigenVectorsInit = EigenVectors;
		memcpy(&PriBasisState, &Basis[PriBasisCount], sizeof(BasisStruct));
		for (SecBasisCount = 0; SecBasisCount < BasisTotal; SecBasisCount++) {
			if (IdentifyPurelyElectronic(&PriBasisState, &Basis[SecBasisCount])) {
				ReperatorElement = *Reperator;
				PriEigenVectorsElement = PriEigenVectorsInit;
				for (PriCount = 0; PriCount < CutOff; PriCount++) {
					SecEigenVectorsElement = PriEigenVectorsInit;
					for (SecCount = 0; SecCount < CutOff; SecCount++) {
						PriSecEigenVectors = *PriEigenVectorsElement * *SecEigenVectorsElement;
						TerEigenVectorsElement = SecEigenVectorsInit;
						for (TerCount = 0; TerCount < CutOff; TerCount++) {
							PriSecTerEigenVectors = PriSecEigenVectors * *TerEigenVectorsElement;
							QuaEigenVectorsElement = SecEigenVectorsInit;
							for (QuaCount = 0; QuaCount < CutOff; QuaCount++)
								(*ReperatorElement++) += PriSecTerEigenVectors * *QuaEigenVectorsElement, QuaEigenVectorsElement += BasisTotal;
							TerEigenVectorsElement += BasisTotal;
						}
						SecEigenVectorsElement += BasisTotal;
					}
					PriEigenVectorsElement += BasisTotal;
				}
			}
			SecEigenVectorsInit++;
		}
		PriEigenVectorsInit++;
	}
}


void BuildRedfieldTensor(MaskCArrayStruct *RedfieldTensor, const ComplexType *KetBathFT, const ComplexType *BraBathFT, const RealType *KetReperator, const RealType *BraReperator, const RealType *KetEigenValues, const RealType *BraEigenValues, const int KetBasisTotal, const int BraBasisTotal, const int KetCutOff, const int BraCutOff, const ParmStruct *Parms)
{
	RealType PriDiff, SecBraEigenValue;
	int PriKetCount, PriBraCount, SecKetCount, SecBraCount, TerCount, TensorId;

	if (strncmp(Parms->BathFile, "|", 1) == 0) {
		DismissMaskCArray(RedfieldTensor);
		return;
	}
	else AllocateMaskCArray(RedfieldTensor, IntPow2(KetCutOff) * IntPow2(BraCutOff), true, true);

	for (PriKetCount = 0; PriKetCount < KetCutOff; PriKetCount++)
		for (PriBraCount = 0; PriBraCount < BraCutOff; PriBraCount++) {
			PriDiff = KetEigenValues[PriKetCount] - BraEigenValues[PriBraCount];
			for (SecKetCount = 0; SecKetCount < KetCutOff; SecKetCount++) {
				SecBraEigenValue = KetEigenValues[SecKetCount] - PriDiff;
				for (SecBraCount = 0; SecBraCount < BraCutOff; SecBraCount++) {
					if (fabs(BraEigenValues[SecBraCount] - SecBraEigenValue) > .0001 * fabs(SecBraEigenValue)) continue; // This is unstable, should be modified

					TensorId = FourDimId_(PriKetCount, PriBraCount, SecKetCount, SecBraCount, BraCutOff, KetCutOff, BraCutOff);

					if (KetReperator == BraReperator) {
						(*RedfieldTensor).Array[TensorId] += (KetBathFT[PriKetCount * KetCutOff + SecKetCount] + conj(BraBathFT[PriBraCount * BraCutOff + SecBraCount])) * KetReperator[FourDimId_(SecBraCount, PriBraCount, PriKetCount, SecKetCount, BraCutOff, KetCutOff, KetCutOff)];
						(*RedfieldTensor).Mask[TensorId] = true;
					}
					
					if (PriBraCount == SecBraCount && KetReperator != NULL) {
						for (TerCount = 0; TerCount < KetCutOff; TerCount++)
							(*RedfieldTensor).Array[TensorId] -= KetBathFT[TerCount * KetCutOff + SecKetCount] * KetReperator[FourDimId_(PriKetCount, TerCount, TerCount, SecKetCount, KetCutOff, KetCutOff, KetCutOff)];
						(*RedfieldTensor).Mask[TensorId] = true;
					}

					if (PriKetCount == SecKetCount && BraReperator != NULL) {
						for (TerCount = 0; TerCount < BraCutOff; TerCount++)
							(*RedfieldTensor).Array[TensorId] -= conj(BraBathFT[TerCount * BraCutOff + SecBraCount]) * BraReperator[FourDimId_(SecBraCount, TerCount, TerCount, PriBraCount, BraCutOff, BraCutOff, BraCutOff)];
						(*RedfieldTensor).Mask[TensorId] = true;
					}
				}
			}
		}
}


void BuildAdiabaticDipoles(RealType **ADipoles, const MaskRArrayStruct Dipoles, const RealType *PermDipoles, const RealType *LoEigenVectors, const RealType *HiEigenVectors, const RealType *LoEigenValues, const RealType *HiEigenValues, const int LoBasisTotal, const int HiBasisTotal, const PolStruct *Pols, const ParmStruct *Parms)
{
	int LoStateCount, LoSiteCount, HiStateCount, HiSiteCount, TransitionSite;
	RealType *PriADipolesElement, *SecADipolesElement, *TerADipolesElement;
	const RealType *LoEigenVectorsElement, *HiEigenVectorsElement, *DipolesArrayElement;
	const RealType Unity = 1.0;

	(*ADipoles) = calloc((size_t) 3 * LoBasisTotal * HiBasisTotal, sizeof(RealType));

	PriADipolesElement = *ADipoles;
	SecADipolesElement = &(*ADipoles)[LoBasisTotal * HiBasisTotal];
	TerADipolesElement = &(*ADipoles)[2 * LoBasisTotal * HiBasisTotal];

	for (LoStateCount = 0; LoStateCount < LoBasisTotal; LoStateCount++)
		for (HiStateCount = 0; HiStateCount < HiBasisTotal; HiStateCount++) {

			if (*Parms->LaserBounds && LoEigenValues != NULL && HiEigenValues != NULL)
				if (Parms->SEnergyMean + HiEigenValues[HiStateCount] - LoEigenValues[LoStateCount] < Parms->LaserBounds[0] || Parms->SEnergyMean + HiEigenValues[HiStateCount] - LoEigenValues[LoStateCount] > Parms->LaserBounds[1]) {
					PriADipolesElement++; SecADipolesElement++; TerADipolesElement++;
					continue;
				}

			*PriADipolesElement = 0; *SecADipolesElement = 0; *TerADipolesElement = 0;
			if (LoEigenVectors != NULL) LoEigenVectorsElement = &LoEigenVectors[LoStateCount * LoBasisTotal];
			else LoEigenVectorsElement = &Unity;
			DipolesArrayElement = Dipoles.Array;
			for (LoSiteCount = 0; LoSiteCount < LoBasisTotal; LoSiteCount++) {
				if (LoEigenVectors == NULL) {
					LoSiteCount = LoStateCount;
					DipolesArrayElement = &Dipoles.Array[LoStateCount * HiBasisTotal];
				}
				HiEigenVectorsElement = &HiEigenVectors[HiStateCount * HiBasisTotal];
				for (HiSiteCount = 0; HiSiteCount < HiBasisTotal; HiSiteCount++) {
					if (Dipoles.Mask[LoSiteCount * HiBasisTotal + HiSiteCount] != -1) {
						TransitionSite = Dipoles.Mask[LoSiteCount * HiBasisTotal + HiSiteCount];
						if (Pols->Active[0])
							*PriADipolesElement += PermDipoles[TransitionSite] * *LoEigenVectorsElement * *DipolesArrayElement * *HiEigenVectorsElement;
						if (Pols->Active[1])
							*SecADipolesElement += PermDipoles[2 * Parms->SiteTotal + TransitionSite] * *LoEigenVectorsElement * *DipolesArrayElement * *HiEigenVectorsElement;
						if (Pols->Active[2])
							*TerADipolesElement += PermDipoles[4 * Parms->SiteTotal + TransitionSite] * *LoEigenVectorsElement * *DipolesArrayElement * *HiEigenVectorsElement;
					}
					DipolesArrayElement++; HiEigenVectorsElement++;
				}

				if (LoEigenVectors == NULL) break;
				LoEigenVectorsElement++;
			}
			PriADipolesElement++; SecADipolesElement++; TerADipolesElement++;
		}
}


void K0ToSEigenVectors(RealType **SEigenVectors, const RealType *K0EigenVectors, const int *K0Vectors, const ParmStruct *Parms)
{
	const RealType *K0EigenVectorsElement;
	RealType *SEigenVectorsElement;
	int SBasisCount, K0BasisCount, StateCount, ClubCount;

	PrintTime(); PrintMessage("Transforming back from k-space.", "", "");

	*SEigenVectors = realloc(*SEigenVectors, *Parms->K0BasisBound * *Parms->SBasisBound * sizeof(RealType));

	for (ClubCount = 1; ClubCount < ClubTotal - 1 && Parms->SBasisBound[ClubCount] != *Parms->SBasisBound; ClubCount++) {
		for (StateCount = 0; StateCount < *Parms->K0BasisBound; StateCount++) {
			SEigenVectorsElement = &(*SEigenVectors)[StateCount * *Parms->SBasisBound + Parms->SBasisBound[ClubCount]];
			for (SBasisCount = Parms->SBasisBound[ClubCount]; SBasisCount < Parms->SBasisBound[ClubCount + 1]; SBasisCount++) {
				*SEigenVectorsElement = 0;
				K0EigenVectorsElement = &K0EigenVectors[StateCount * *Parms->K0BasisBound];
				for (K0BasisCount = 0; K0BasisCount < *Parms->K0BasisBound; K0BasisCount++)
					*SEigenVectorsElement += sqrt(K0Vectors[K0BasisCount * *Parms->SBasisBound + SBasisCount]) * *K0EigenVectorsElement++;
				(*SEigenVectorsElement++) /= (RealType) sqrt(Parms->CellTotal);
			}
		}
	}
}


void SetPolarizations(PolStruct *Pols, const int DoubleCrossPol)
{
	if (DoubleCrossPol) {
		Pols->Weight = .0833333;
		if (Pols->Active[0] && Pols->Active[1]) {
			if (Pols->Id < 5)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 1; Pols->Pulses[2] = 0; Pols->Pulses[3] = 1; Pols->Id = 5;  return; }
			if (Pols->Id < 6)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 0; Pols->Pulses[2] = 1; Pols->Pulses[3] = 0; Pols->Id = 6;  return; }
			if (Pols->Id < 7)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 1; Pols->Pulses[2] = 1; Pols->Pulses[3] = 0; Pols->Id = 7;  Pols->Weight *= -1; return; }
			if (Pols->Id < 8)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 0; Pols->Pulses[2] = 0; Pols->Pulses[3] = 1; Pols->Id = 8;  Pols->Weight *= -1; return; }
		}
		if (Pols->Active[0] && Pols->Active[2]) {
			if (Pols->Id < 11)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 2; Pols->Pulses[2] = 0; Pols->Pulses[3] = 2; Pols->Id = 11; return; }
			if (Pols->Id < 12)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 0; Pols->Pulses[2] = 2; Pols->Pulses[3] = 0; Pols->Id = 12; return; }
			if (Pols->Id < 13)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 2; Pols->Pulses[2] = 2; Pols->Pulses[3] = 0; Pols->Id = 13; Pols->Weight *= -1; return; }
			if (Pols->Id < 14)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 0; Pols->Pulses[2] = 0; Pols->Pulses[3] = 2; Pols->Id = 14; Pols->Weight *= -1; return; }
		}
		if (Pols->Active[1] && Pols->Active[2]) {
			if (Pols->Id < 17)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 2; Pols->Pulses[2] = 1; Pols->Pulses[3] = 2; Pols->Id = 17; return; }
			if (Pols->Id < 18)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 1; Pols->Pulses[2] = 2; Pols->Pulses[3] = 1; Pols->Id = 18; return; }
			if (Pols->Id < 19)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 2; Pols->Pulses[2] = 2; Pols->Pulses[3] = 1; Pols->Id = 19; Pols->Weight *= -1; return; }
			if (Pols->Id < 20)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 1; Pols->Pulses[2] = 1; Pols->Pulses[3] = 2; Pols->Id = 20; Pols->Weight *= -1; return; }
		}
	}
	else {
		Pols->Weight = .2;
		if (Pols->Id < 0 && Pols->Active[0])	{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 0; Pols->Pulses[2] = 0; Pols->Pulses[3] = 0; Pols->Id = 0;  return; }
		if (Pols->Id < 1 && Pols->Active[1])	{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 1; Pols->Pulses[2] = 1; Pols->Pulses[3] = 1; Pols->Id = 1;  return; }
		if (Pols->Id < 2 && Pols->Active[2])	{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 2; Pols->Pulses[2] = 2; Pols->Pulses[3] = 2; Pols->Id = 2;  return; }

		Pols->Weight = .0666667;
		if (Pols->Active[1] && Pols->Active[0]) {
			if (Pols->Id < 3)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 0; Pols->Pulses[2] = 1; Pols->Pulses[3] = 1; Pols->Id = 3;  return; }
			if (Pols->Id < 4)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 1; Pols->Pulses[2] = 0; Pols->Pulses[3] = 0; Pols->Id = 4;  return; }
			if (Pols->Id < 5)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 1; Pols->Pulses[2] = 0; Pols->Pulses[3] = 1; Pols->Id = 5;  return; }
			if (Pols->Id < 6)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 0; Pols->Pulses[2] = 1; Pols->Pulses[3] = 0; Pols->Id = 6;  return; }
			if (Pols->Id < 7)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 1; Pols->Pulses[2] = 1; Pols->Pulses[3] = 0; Pols->Id = 7;  return; }
			if (Pols->Id < 8)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 0; Pols->Pulses[2] = 0; Pols->Pulses[3] = 1; Pols->Id = 8;  return; }
		}
		if (Pols->Active[2] && Pols->Active[0]) {
			if (Pols->Id < 9)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 0; Pols->Pulses[2] = 2; Pols->Pulses[3] = 2; Pols->Id = 9;  return; }
			if (Pols->Id < 10)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 2; Pols->Pulses[2] = 0; Pols->Pulses[3] = 0; Pols->Id = 10; return; }
			if (Pols->Id < 11)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 2; Pols->Pulses[2] = 0; Pols->Pulses[3] = 2; Pols->Id = 11; return; }
			if (Pols->Id < 12)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 0; Pols->Pulses[2] = 2; Pols->Pulses[3] = 0; Pols->Id = 12; return; }
			if (Pols->Id < 13)					{ Pols->Pulses[0] = 0; Pols->Pulses[1] = 2; Pols->Pulses[2] = 2; Pols->Pulses[3] = 0; Pols->Id = 13; return; }
			if (Pols->Id < 14)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 0; Pols->Pulses[2] = 0; Pols->Pulses[3] = 2; Pols->Id = 14; return; }
		}
		if (Pols->Active[2] && Pols->Active[1]) {
			if (Pols->Id < 15)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 1; Pols->Pulses[2] = 2; Pols->Pulses[3] = 2; Pols->Id = 15; return; }
			if (Pols->Id < 16)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 2; Pols->Pulses[2] = 1; Pols->Pulses[3] = 1; Pols->Id = 16; return; }
			if (Pols->Id < 17)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 2; Pols->Pulses[2] = 1; Pols->Pulses[3] = 2; Pols->Id = 17; return; }
			if (Pols->Id < 18)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 1; Pols->Pulses[2] = 2; Pols->Pulses[3] = 1; Pols->Id = 18; return; }
			if (Pols->Id < 19)					{ Pols->Pulses[0] = 1; Pols->Pulses[1] = 2; Pols->Pulses[2] = 2; Pols->Pulses[3] = 1; Pols->Id = 19; return; }
			if (Pols->Id < 20)					{ Pols->Pulses[0] = 2; Pols->Pulses[1] = 1; Pols->Pulses[2] = 1; Pols->Pulses[3] = 2; Pols->Id = 20; return; }
		}
	}
	Pols->Id = 21;
}


void CreatePermDipolesSeq(RealType *PermDipolesSeq, PolStruct *Pols, const RealType *PermDipoles, const ParmStruct *Parms)
{
	int SeqCount, SiteCount;

	if (*Parms->Polarizations != -1) {
		for (SeqCount = 0; SeqCount < 4; SeqCount++)
			for (SiteCount = 0; SiteCount < 2 * Parms->SiteTotal; SiteCount++)
				PermDipolesSeq[SeqCount * 2 * Parms->SiteTotal + SiteCount]
						= Parms->Polarizations[SeqCount] * PermDipoles[SiteCount] + Parms->Polarizations[4 + SeqCount] * PermDipoles[2 * Parms->SiteTotal + SiteCount];
		Pols->Id = 21;
	}
	else
		for (SeqCount = 0; SeqCount < 4; SeqCount++)
			memcpy(&PermDipolesSeq[SeqCount * 2 * Parms->SiteTotal], &PermDipoles[Pols->Pulses[SeqCount] * 2 * Parms->SiteTotal], 2 * Parms->SiteTotal * sizeof(RealType));
}


void DirectAbsorption(RealType **Absorption, RealType **LineStrengths, const RealType *EnergyGrid, const RealType *PermDipoles, const RealType *SEigenVectors, const RealType *SEigenValues, const MaskRArrayStruct *GSDipoles, const PolStruct *Pols, const ParmStruct *Parms, const int StateTotal)
{
	RealType *LineStrengthsPol = malloc(StateTotal * sizeof(RealType));
	int SpectrumCount, StateCount, PolCount;

	*Absorption = calloc((size_t) Parms->SpectrumTotal * 4, sizeof(RealType));
	*LineStrengths = calloc((size_t) StateTotal, sizeof(RealType));

	for (PolCount = 0; PolCount < 3; PolCount++) {
		if (!Pols->Active[PolCount]) continue;

		for (StateCount = 0; StateCount < StateTotal; StateCount++) {
			LineStrengthsPol[StateCount] = LineStrength(0, PolCount, GSDipoles, PermDipoles, &SEigenVectors[StateCount * *Parms->SBasisBound], Parms);
			(*LineStrengths)[StateCount] += LineStrengthsPol[StateCount];
		}

		for (SpectrumCount = 0; SpectrumCount < Parms->SpectrumTotal; SpectrumCount++) {
			for (StateCount = 0; StateCount < StateTotal; StateCount++) {
				if (Parms->AbsMode == AbsModeNick)
					(*Absorption)[SpectrumCount * 4 + PolCount + 1] += (Parms->SEnergyMean + SEigenValues[StateCount]) * LineStrengthsPol[StateCount] * Lorentzian(EnergyGrid[SpectrumCount] - SEigenValues[StateCount], Parms->HBroad[PolCount]);
				else if (Parms->Gaussian || Parms->AbsMode == AbsModeTim)
					(*Absorption)[SpectrumCount * 4 + PolCount + 1] += LineStrengthsPol[StateCount] * Gaussian(EnergyGrid[SpectrumCount] - SEigenValues[StateCount], Parms->HBroad[PolCount]);
				else
					(*Absorption)[SpectrumCount * 4 + PolCount + 1] += LineStrengthsPol[StateCount] * Lorentzian(EnergyGrid[SpectrumCount] - SEigenValues[StateCount], Parms->HBroad[PolCount]);
			}

			if (Parms->AbsMode == AbsModeTim)
				(*Absorption)[SpectrumCount * 4 + PolCount + 1] /= RealPow2((Parms->SEnergyMean + EnergyGrid[SpectrumCount]) / Parms->SEnergyMean);

			(*Absorption)[SpectrumCount * 4] += (*Absorption)[SpectrumCount * 4 + PolCount + 1];
		}
	}

	free(LineStrengthsPol);
}


void DirectFluorescence(RealType **Fluorescence, RealType **FluStrengths, const RealType *EnergyGrid, const RealType *PermDipoles, const RealType *SEigenVectors, const RealType *SEigenValues, const MaskRArrayStruct *GSDipoles, const BasisStruct *GBasis, const PolStruct *Pols, const ParmStruct *Parms, const int StateTotal)
{
    RealType *LineStrengthsPol = malloc(Parms->VTotal * sizeof(RealType));
    RealType BoltzmannFactor, PartitionFunction = 0;
    int SpectrumCount, GCount, SCount, PolCount, VCount;

    *Fluorescence = calloc((size_t) Parms->SpectrumTotal * 4, sizeof(RealType));
    *FluStrengths = calloc((size_t) StateTotal, sizeof(RealType));

    if (Parms->Temperature != 0)
        for (SCount = 0; SCount < StateTotal; SCount++)
            PartitionFunction += (RealType) exp((double) -(SEigenValues[SCount] - *SEigenValues) / Parms->Temperature);

    for (SCount = 0; SCount < StateTotal; SCount++) {
        FillRArray(LineStrengthsPol, 0, Parms->VTotal);
        BoltzmannFactor = (Parms->Temperature == 0) ? (SCount == 0) : (RealType) exp((double) -(SEigenValues[SCount] - *SEigenValues) / Parms->Temperature) / PartitionFunction;

        if (BoltzmannFactor < .00001) break;

        for (PolCount = 0; PolCount < 3; PolCount++) {
            if (!Pols->Active[PolCount]) continue;

            for (GCount = 0; GCount < *Parms->GBasisBound; GCount++) {
                VCount = IntMax(GBasis[GCount].P[0].V[Ground], 0) + IntMax(GBasis[GCount].P[1].V[Ground], 0);
                LineStrengthsPol[VCount] += LineStrength(GCount, PolCount, GSDipoles, PermDipoles, &SEigenVectors[SCount * *Parms->SBasisBound], Parms);
            }

            //if (PolCount == 0)
                (*FluStrengths)[SCount] += BoltzmannFactor * SumRArray(LineStrengthsPol, Parms->VTotal);

            for (SpectrumCount = 0; SpectrumCount < Parms->SpectrumTotal; SpectrumCount++)
                for (VCount = 0; VCount < Parms->VTotal; VCount++)
                    (*Fluorescence)[SpectrumCount * 4 + PolCount + 1] += BoltzmannFactor * LineStrengthsPol[VCount] * Gaussian(EnergyGrid[SpectrumCount] - SEigenValues[SCount] + Parms->StokesShift + VCount * Parms->VEnergy, Parms->HBroad[PolCount]);
        }
    }

    for (SpectrumCount = 0; SpectrumCount < Parms->SpectrumTotal; SpectrumCount++)
        (*Fluorescence)[SpectrumCount * 4] += SumRArray(&(*Fluorescence)[SpectrumCount * 4 + 1], 3);

    free(LineStrengthsPol);
}


void DirectCPopulations(RealType **CPopulations, const RealType *SEigenVectors, const BasisStruct *SBasis, const ParmStruct *Parms, const int StateTotal)
{
	int StateCount, SBasisCount;

	*CPopulations = calloc((size_t) StateTotal * 3, sizeof(RealType));

	for (StateCount = 0; StateCount < StateTotal; StateCount++)
		for (SBasisCount = 0; SBasisCount < *Parms->SBasisBound; SBasisCount++) {
			if (Present(&SBasis[SBasisCount], Singlet))			(*CPopulations)[3 * StateCount] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
			else if (Present(&SBasis[SBasisCount], Electron))	(*CPopulations)[3 * StateCount + 1] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
			else if (Present(&SBasis[SBasisCount], Triplet))	(*CPopulations)[3 * StateCount + 2] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
		}
}


void DirectCPopulationsExtended(const RealType *SEigenVectors, const RealType *LineStrengths, const RealType *SEigenValues, const BasisStruct *SBasis, const ParmStruct *Parms, const int StateTotal)
{
	int CCount, StateCount, SBasisCount, LineCount = 0;
	RealType LineStrengthTrunc = (RealType) .00001 * RArrayMax(LineStrengths, StateTotal);
	char FileName[20];
	FILE *FileID;

	RealType *CPopulations = calloc((size_t) StateTotal * 9, sizeof(RealType));

	for (StateCount = 0; StateCount < StateTotal; StateCount++)
		for (SBasisCount = 0; SBasisCount < *Parms->SBasisBound; SBasisCount++) {
			if (Present(&SBasis[SBasisCount], Singlet)) {
				if (CountVibrations(&SBasis[SBasisCount]) == 0)
					CPopulations[9 * StateCount] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
				else if (CountVibrations(&SBasis[SBasisCount]) == 1)
					CPopulations[9 * StateCount + 1] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
				else
					CPopulations[9 * StateCount + 2] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
			}
			else if (Present(&SBasis[SBasisCount], Electron)) {
				if (CountVibrations(&SBasis[SBasisCount]) == 0)
					CPopulations[9 * StateCount + 3] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
				else if (CountVibrations(&SBasis[SBasisCount]) == 1)
					CPopulations[9 * StateCount + 4] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
				else
					CPopulations[9 * StateCount + 5] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
			}
			else if (Present(&SBasis[SBasisCount], Triplet)) {
				if (CountVibrations(&SBasis[SBasisCount]) == 0)
					CPopulations[9 * StateCount + 6] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
				else if (CountVibrations(&SBasis[SBasisCount]) == 1)
					CPopulations[9 * StateCount + 7] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
				else
					CPopulations[9 * StateCount + 8] += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
			}
		}

	CreateFileName(FileName, "DirectCPops"); FileID = fopen(FileName, "w");
	for (StateCount = 0; StateCount < StateTotal; StateCount++) {
		if (LineCount > 100000) {
			fprintf(FileID, "File truncated...");
			break;
		}
		if (LineStrengths[StateCount] > LineStrengthTrunc) {
			fprintf(FileID, "%d, %f, ", StateCount, SEigenValues[StateCount]);
			for (CCount = 0; CCount < 9; CCount++) fprintf(FileID, "%f, ", CPopulations[StateCount * 9 + CCount]);
			fprintf(FileID, "\n");
			LineCount++;
		}
	}
	fclose(FileID);
}


void DirectTDensity(const RealType *SEigenVectors, const RealType *LineStrengths, const int *TransDistances, const DCrystalStruct *DCrystal, const BasisStruct *SBasis, const ParmStruct *Parms, const int StateTotal)
// Only works for bare TT basis states
{
	RealType *TDensity = calloc((size_t) StateTotal * Parms->SiteTotal, sizeof(RealType));
	RealType LineStrengthTrunc = (RealType) .00001 * RArrayMax(LineStrengths, StateTotal);
	int StateCount, SiteCount, SBasisCount, PrintCount = 0, PrintTotal = 150;
	FILE *FileID;

	for (StateCount = 0; StateCount < StateTotal; StateCount++)
		if (LineStrengths[StateCount] > LineStrengthTrunc)
			for (SBasisCount = 0; SBasisCount < *Parms->SBasisBound; SBasisCount++)
				if (SBasis[SBasisCount].P[0].C == Triplet)
					TDensity[StateCount * Parms->SiteTotal + TransDistances[SBasis[SBasisCount].P[0].S * Parms->SiteTotal + SBasis[SBasisCount].P[1].S]]
							+= RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);

	FileID = fopen("DirectTDensity.dat", "w");
	fprintf(FileID, "              ");
	for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++) fprintf(FileID, "%5.1f   ", .5 * DCrystal[SiteCount].TwoTimesDPosA);
	fprintf(FileID, "\n              ");
	for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++) fprintf(FileID, "%5.1f   ", .5 * DCrystal[SiteCount].TwoTimesDPosB);
	fprintf(FileID, "\n");
	for (StateCount = 0; StateCount < StateTotal; StateCount++) {
		if (LineStrengths[StateCount] <= LineStrengthTrunc) continue;
		fprintf(FileID, "%3d  %7.4f  ", StateCount, LineStrengths[StateCount]);
		for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++) fprintf(FileID, "  %6.4f", TDensity[StateCount * Parms->SiteTotal + SiteCount]);
		fprintf(FileID, "\n");
		if (++PrintCount > PrintTotal) break;
	}
	fclose(FileID);

	free(TDensity);
}


void DirectCTDensity(const RealType *SEigenVectors, const RealType *LineStrengths, const int *TransDistances, const DCrystalStruct *DCrystal, const BasisStruct *SBasis, const ParmStruct *Parms, const int StateTotal)
// Only works for bare CT basis states
{
    RealType *CTDensity = calloc((size_t) StateTotal * Parms->SiteTotal, sizeof(RealType));
    RealType LineStrengthTrunc = (RealType) .00001 * RArrayMax(LineStrengths, StateTotal);
    int StateCount, SiteCount, SBasisCount, PrintCount = 0, PrintTotal = 150;
    FILE *FileID;

    for (StateCount = 0; StateCount < StateTotal; StateCount++)
        if (LineStrengths[StateCount] > LineStrengthTrunc)
            for (SBasisCount = 0; SBasisCount < *Parms->SBasisBound; SBasisCount++) {
                if (SBasis[SBasisCount].P[0].C == Electron)
                    CTDensity[StateCount * Parms->SiteTotal + TransDistances[SBasis[SBasisCount].P[0].S * Parms->SiteTotal + SBasis[SBasisCount].P[1].S]]
                            += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
                if (SBasis[SBasisCount].P[0].C == Hole)
                    CTDensity[StateCount * Parms->SiteTotal + TransDistances[SBasis[SBasisCount].P[1].S * Parms->SiteTotal + SBasis[SBasisCount].P[0].S]]
                            += RealPow2(SEigenVectors[StateCount * *Parms->SBasisBound + SBasisCount]);
            }

    FileID = fopen("DirectCTDensity.dat", "w");
    fprintf(FileID, "              ");
    for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++) fprintf(FileID, "%5.1f   ", .5 * DCrystal[SiteCount].TwoTimesDPosA);
    fprintf(FileID, "\n              ");
    for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++) fprintf(FileID, "%5.1f   ", .5 * DCrystal[SiteCount].TwoTimesDPosB);
    fprintf(FileID, "\n");
    for (StateCount = 0; StateCount < StateTotal; StateCount++) {
        if (LineStrengths[StateCount] <= LineStrengthTrunc) continue;
        fprintf(FileID, "%3d  %7.4f  ", StateCount, LineStrengths[StateCount]);
        for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++) fprintf(FileID, "  %6.4f", CTDensity[StateCount * Parms->SiteTotal + SiteCount]);
        fprintf(FileID, "\n");
        if (++PrintCount > PrintTotal) break;
    }
    fclose(FileID);

    free(CTDensity);
}


void DirectLineStrengths(RealType **LineStrengths, const RealType *SEigenValues, const RealType *GSAdiabDipoles, const ParmStruct *Parms, const int StateTotal)
{
	int PolCount, StateCount;

	*LineStrengths = calloc((size_t) StateTotal, sizeof(RealType));

	for (StateCount = 0; StateCount < StateTotal; StateCount++)
		for (PolCount = 0; PolCount < 3; PolCount++)
			(*LineStrengths)[StateCount] += RealPow2(GSAdiabDipoles[PolCount * *Parms->GBasisBound * StateTotal + StateCount]);
}


void DirectFTMapAnalysis(PolStruct *Pols, const RealType *LineStrengths, const RealType *CPopulations, const RealType *GSAdiabDipoles, const RealType *SDAdiabDipoles, const RealType *SEigenValues, const RealType *DEigenValues, const BasisStruct *GBasis, const ParmStruct *Parms)
{
	int GCount, PriSCount, SecSCount, DCount, IdS1 = 0, IdTTMax = 0;
	int GSBasisTotal = *Parms->GBasisBound * *Parms->SBasisBound, SDBasisTotal = *Parms->SBasisBound * *Parms->DBasisBound;
	RealType TTMax = 0, SDLineStrength, GB = 0, EA_1 = 0, EA_2 = 0;
	char FileName[20];
	FILE *FileID;

	PrintTime(); PrintMessage("FT amplitude map analysis in progress...", "", "");

	for (PriSCount = 0; PriSCount < *Parms->SBasisBound; PriSCount++)
		if (CPopulations[3 * PriSCount + 2] < .5) {
			IdS1 = PriSCount;
			break;
		}

	for (PriSCount = 0; PriSCount < IdS1; PriSCount++)
		if (LineStrengths[PriSCount] > TTMax) {
			TTMax = LineStrengths[PriSCount];
			IdTTMax = PriSCount;
		}

	CreateFileName(FileName, "TripletTransitions"); FileID = fopen(FileName, "w");
	for (DCount = 0; DCount < *Parms->DBasisBound; DCount++) {
		SDLineStrength = 0;
		for (Pols->Id = 0; Pols->Id < 3; Pols->Id++)
			SDLineStrength += RealPow2(SDAdiabDipoles[Pols->Id * SDBasisTotal + IdTTMax * *Parms->DBasisBound + DCount]);
		if (SDLineStrength != 0) fprintf(FileID, "%d %f %f\n", DCount, DEigenValues[DCount] - SEigenValues[IdTTMax], SDLineStrength);
	}
	fclose(FileID);

	Pols->Id = -1;
	while (true) {
		SetPolarizations(Pols, Parms->DoubleCrossPol);

		if (Pols->Id > 20) break;

		for (GCount = 0; GCount < Parms->GBasisBound[3]; GCount++)
			if (VSum(GBasis[GCount]) == 1)
				GB += Pols->Weight
					  * GSAdiabDipoles[Pols->Pulses[0] * GSBasisTotal + IdS1]
					  * GSAdiabDipoles[Pols->Pulses[1] * GSBasisTotal + GCount * *Parms->SBasisBound + IdS1]
					  * GSAdiabDipoles[Pols->Pulses[2] * GSBasisTotal + GCount * *Parms->SBasisBound + IdS1]
					  * GSAdiabDipoles[Pols->Pulses[3] * GSBasisTotal + IdS1];

		for (PriSCount = 0; PriSCount < IdS1; PriSCount++)
			for (DCount = 0; DCount < *Parms->DBasisBound; DCount++)
				EA_1 -= Pols->Weight
						* GSAdiabDipoles[Pols->Pulses[0] * GSBasisTotal + PriSCount]
						* SDAdiabDipoles[Pols->Pulses[1] * SDBasisTotal + PriSCount * *Parms->DBasisBound + DCount]
						* SDAdiabDipoles[Pols->Pulses[2] * SDBasisTotal + IdS1 * *Parms->DBasisBound + DCount]
						* GSAdiabDipoles[Pols->Pulses[3] * GSBasisTotal + IdS1];

		for (PriSCount = 0; PriSCount < IdS1; PriSCount++)
			for (SecSCount = IdS1 + 1; SecSCount < *Parms->SBasisBound; SecSCount++)
				for (DCount = 0; DCount < *Parms->DBasisBound; DCount++)
					EA_2 -= Pols->Weight
						  * GSAdiabDipoles[Pols->Pulses[0] * GSBasisTotal + PriSCount]
						  * SDAdiabDipoles[Pols->Pulses[1] * SDBasisTotal + PriSCount * *Parms->DBasisBound + DCount]
						  * SDAdiabDipoles[Pols->Pulses[2] * SDBasisTotal + SecSCount * *Parms->DBasisBound + DCount]
						  * GSAdiabDipoles[Pols->Pulses[3] * GSBasisTotal + SecSCount];
	}

	CreateFileName(FileName, "FTMapAnalysis"); FileID = fopen(FileName, "w");
	fprintf(FileID, "%f\n%f\n%f\n", GB, EA_1, EA_2);
	fclose(FileID);
}


void Excite(ComplexType *HiCoef, const ComplexType *LoCoef, const RealType *LoHiDipoles, const int HiTotal, const int LoTotal)
{
	ComplexType *HiCoefElement = HiCoef;
	int HiCount, LoCount;

	for (HiCount = 0; HiCount < HiTotal; HiCount++) {
		*HiCoefElement = 0;
		for (LoCount = 0; LoCount < LoTotal; LoCount++) *HiCoefElement += LoHiDipoles[LoCount * HiTotal + HiCount] * LoCoef[LoCount];
		HiCoefElement++;
	}
}


void DeExcite(ComplexType *LoCoef, const ComplexType *HiCoef, const RealType *LoHiDipoles, const int LoTotal, const int HiTotal)
{
	ComplexType *LoCoefElement = LoCoef;
	const ComplexType *HiCoefElement;
	const RealType *LoHiDipolesElement = LoHiDipoles;
	int LoCount, HiCount;

	for (LoCount = 0; LoCount < LoTotal; LoCount++) {
		*LoCoefElement = 0;
		HiCoefElement = HiCoef;
		for (HiCount = 0; HiCount < HiTotal; HiCount++) *LoCoefElement += (*LoHiDipolesElement++) * (*HiCoefElement++);
		LoCoefElement++;
	}
}


ComplexType DeExciteVacuum(const ComplexType *Coef, const RealType *LoHiDipoles, const int Total)
{
	ComplexType ReturnValue = 0;
	const ComplexType *CoefElement = Coef;
	const RealType *LoHiDipolesElement = LoHiDipoles;
	int Count;

	for (Count = 0; Count < Total; Count++) ReturnValue += (*LoHiDipolesElement++) * (*CoefElement++);

	return ReturnValue;
}


void ExciteVacuumDiabatic(ComplexType *Coef, const RealType *PermDipoles, const MaskRArrayStruct *LoHiDipoles, const int Total)
{
	ComplexType *CoefElement = Coef;
	int Count;

	for (Count = 0; Count < Total; Count++) {
		if ((*LoHiDipoles).Mask[Count] != -1) 	(*CoefElement++) = (*LoHiDipoles).Array[Count] * PermDipoles[(*LoHiDipoles).Mask[Count]];
		else									(*CoefElement++) = 0;
	}
}


ComplexType DeExciteVacuumDiabatic(const ComplexType *Coef, const RealType *PermDipoles, const MaskRArrayStruct *LoHiDipoles, const int Total)
{
	ComplexType ReturnValue = 0;
	const ComplexType *CoefElement = Coef;
	const int *LoHiDipolesMaskElement = (*LoHiDipoles).Mask;
	const RealType *LoHiDipolesArrayElement = (*LoHiDipoles).Array;
	int Count;

	for (Count = 0; Count < Total; Count++) {
		if (*LoHiDipolesMaskElement != -1) ReturnValue += (*LoHiDipolesArrayElement) * PermDipoles[*LoHiDipolesMaskElement] * (*CoefElement);
		LoHiDipolesMaskElement++; LoHiDipolesArrayElement++; CoefElement++;
	}

	return ReturnValue;
}


void ExciteDiabatic(ComplexType *HiCoef, const ComplexType *LoCoef, const RealType *PermDipoles, const MaskRArrayStruct *LoHiDipoles, const int HiTotal, const int LoTotal)
{
	ComplexType *HiCoefElement = HiCoef;
	int HiCount, LoCount;

	for (HiCount = 0; HiCount < HiTotal; HiCount++) {
		*HiCoefElement = 0;
		for (LoCount = 0; LoCount < LoTotal; LoCount++)
			if ((*LoHiDipoles).Mask[LoCount * HiTotal + HiCount] != -1)
				*HiCoefElement += (*LoHiDipoles).Array[LoCount * HiTotal + HiCount] * PermDipoles[(*LoHiDipoles).Mask[LoCount * HiTotal + HiCount]] * LoCoef[LoCount];
		HiCoefElement++;
	}
}


void DeExciteDiabatic(ComplexType *LoCoef, const ComplexType *HiCoef, const RealType *PermDipoles, const MaskRArrayStruct *LoHiDipoles, const int LoTotal, const int HiTotal)
{
	ComplexType *LoCoefElement = LoCoef;
	const ComplexType *HiCoefElement;
	const int *LoHiDipolesMaskElement = (*LoHiDipoles).Mask;
	const RealType *LoHiDipolesArrayElement = (*LoHiDipoles).Array;
	int LoCount, HiCount;

	for (LoCount = 0; LoCount < LoTotal; LoCount++) {
		*LoCoefElement = 0;
		HiCoefElement = HiCoef;
		for (HiCount = 0; HiCount < HiTotal; HiCount++) {
			if (*LoHiDipolesMaskElement != -1) *LoCoefElement += (*LoHiDipolesArrayElement) * PermDipoles[*LoHiDipolesMaskElement] * (*HiCoefElement);
			LoHiDipolesMaskElement++; LoHiDipolesArrayElement++; HiCoefElement++;
		}
		LoCoefElement++;
	}
}


void Propagate(ComplexType *Coef, ComplexType *CoefWork, const MaskRArrayStruct Hamiltonian, const int BasisTotal)
{
	ComplexType *CoefPos, *PriCoefWorkPos, *SecCoefWorkPos, *TerCoefWorkPos, *QuaCoefWorkPos;
	int Iteration, Count, TwoTimesFactor;

	CoefIncrement(CoefWork, Coef, Hamiltonian, BasisTotal); // k1

	for (Iteration = 0; Iteration < 3; Iteration++) {

		TwoTimesFactor = 1 + (Iteration == 2);

		CoefPos = Coef; PriCoefWorkPos = &CoefWork[(2 * Iteration + 1) * BasisTotal]; SecCoefWorkPos = &CoefWork[2 * Iteration * BasisTotal];
		for (Count = 0; Count < BasisTotal; Count++) (*PriCoefWorkPos++) = (*CoefPos++) + (ComplexType) (.5 * TwoTimesFactor) * (*SecCoefWorkPos++); // y + c_I * k_I, where I = Iteration + 1

		CoefIncrement(&CoefWork[(2 * Iteration + 2) * BasisTotal], &CoefWork[(2 * Iteration + 1) * BasisTotal],
					  Hamiltonian, BasisTotal); // k_(I+1)
	}

	CoefPos = Coef; PriCoefWorkPos = CoefWork; SecCoefWorkPos = &CoefWork[2 * BasisTotal]; TerCoefWorkPos = &CoefWork[4 * BasisTotal]; QuaCoefWorkPos = &CoefWork[6 * BasisTotal];
	for (Count = 0; Count < BasisTotal; Count++) (*CoefPos++) += 0.16666666667 * ((*PriCoefWorkPos++) + 2 * (*SecCoefWorkPos++) + 2 * (*TerCoefWorkPos++) + (*QuaCoefWorkPos++));
}


void Propagate_Redfield(ComplexType *Coef, ComplexType *CoefWork, const RealType *KetEigenValues, const RealType *BraEigenValues, const MaskCArrayStruct RedfieldTensor, const int KetBasisTotal, const int BraBasisTotal, const int KetCutOff, const int BraCutOff, const ParmStruct *Parms)
{
	ComplexType *CoefPos, *PriCoefWorkPos, *SecCoefWorkPos, *TerCoefWorkPos, *QuaCoefWorkPos;
	int Iteration, Count, Dimension = KetBasisTotal * BraBasisTotal, TwoTimesFactor;

	CoefIncrement_Redfield(CoefWork, Coef, KetEigenValues, BraEigenValues, RedfieldTensor, KetBasisTotal, BraBasisTotal, KetCutOff, BraCutOff, Parms->EnergyToFreq * Parms->TimeStep); // k1

	for (Iteration = 0; Iteration < 3; Iteration++) {

		TwoTimesFactor = 1 + (Iteration == 2);

		CoefPos = Coef; PriCoefWorkPos = &CoefWork[(2 * Iteration + 1) * Dimension]; SecCoefWorkPos = &CoefWork[2 * Iteration * Dimension];
		for (Count = 0; Count < Dimension; Count++) (*PriCoefWorkPos++) = (*CoefPos++) + (ComplexType) .5 * (TwoTimesFactor * (*SecCoefWorkPos++)); // y + c_I * k_I, where I = Iteration + 1

		CoefIncrement_Redfield(&CoefWork[(2 * Iteration + 2) * Dimension], &CoefWork[(2 * Iteration + 1) * Dimension], KetEigenValues, BraEigenValues, RedfieldTensor, KetBasisTotal, BraBasisTotal, KetCutOff, BraCutOff, Parms->EnergyToFreq * Parms->TimeStep); // k_(I+1)
	}

	CoefPos = Coef; PriCoefWorkPos = CoefWork; SecCoefWorkPos = &CoefWork[2 * Dimension]; TerCoefWorkPos = &CoefWork[4 * Dimension]; QuaCoefWorkPos = &CoefWork[6 * Dimension];
	for (Count = 0; Count < Dimension; Count++) (*CoefPos++) += 0.16666666667 * ((*PriCoefWorkPos++) + 2 * (*SecCoefWorkPos++) + 2 * (*TerCoefWorkPos++) + (*QuaCoefWorkPos++));
}


void ExtractCPopulations(RealType  *CPopulations, const RealType *Populations, const BasisStruct *SBasis, const ParmStruct *Parms)
{
	int TimeCount, SCount;

	for (TimeCount = 0; TimeCount < *Parms->TimeTotal; TimeCount++)
		for (SCount = 0; SCount < *Parms->SBasisBound; SCount++) {
			if (Present(&SBasis[SCount], Singlet))	CPopulations[3 * TimeCount] += Populations[*Parms->SBasisBound * TimeCount + SCount];
			if (Present(&SBasis[SCount], Electron))	CPopulations[3 * TimeCount + 1] += Populations[*Parms->SBasisBound * TimeCount + SCount];
			if (Present(&SBasis[SCount], Triplet))	CPopulations[3 * TimeCount + 2] += Populations[*Parms->SBasisBound * TimeCount + SCount];
		}
}


void ExtractSPopulations(RealType  *SPopulations, const RealType *Populations, const BasisStruct *SBasis, const ParmStruct *Parms)
{
    int TimeCount, SCount, PCount;

    for (SCount = 0; SCount < *Parms->SBasisBound; SCount++)
        for (PCount = 0; PCount < PTotal && SBasis[SCount].P[PCount].S != -1; PCount++)
            if (SBasis[SCount].P[PCount].C == Singlet)
                for (TimeCount = 0; TimeCount < *Parms->TimeTotal; TimeCount++)
                    SPopulations[Parms->SiteTotal * TimeCount + SBasis[SCount].P[PCount].S] += Populations[*Parms->SBasisBound * TimeCount + SCount];
}


void SaveLineStrengths(const RealType *LineStrengths, const RealType *FluStrengths, const RealType *SEigenValues, const int StateTotal)
{
	int StateCount, LineCount = 0;
	RealType LineStrengthTrunc = (RealType) .00001 * RArrayMax(LineStrengths, StateTotal);
    RealType FluStrengthTrunc = (FluStrengths == NULL ? 0 : (RealType) .00001 * RArrayMax(FluStrengths, StateTotal));
	char FileName[20];
	FILE *FileID;

	CreateFileName(FileName, "DirectSticks"); FileID = fopen(FileName, "w");
	for (StateCount = 0; StateCount < StateTotal; StateCount++) {
		if (LineCount > 100000) {
			fprintf(FileID, "File truncated...");
			break;
		}
		if (LineStrengths[StateCount] > LineStrengthTrunc || (FluStrengths == NULL ? false : FluStrengths[StateCount] > FluStrengthTrunc)) {
            if (FluStrengths == NULL)
			    fprintf(FileID, "%d, %f, %f,\n", StateCount, SEigenValues[StateCount], LineStrengths[StateCount]);
            else
                fprintf(FileID, "%d, %f, %f, %f,\n", StateCount, SEigenValues[StateCount], LineStrengths[StateCount], FluStrengths[StateCount]);
			LineCount++;
		}
	}
	fclose(FileID);
}


void SaveCPopulations(const RealType *CPopulations, const RealType *LineStrengths, const RealType *SEigenValues, const int StateTotal)
{
	int StateCount, LineCount = 0;
	RealType LineStrengthTrunc = (RealType) .00001 * RArrayMax(LineStrengths, StateTotal);
	char FileName[20];
	FILE *FileID;

	CreateFileName(FileName, "DirectCPops"); FileID = fopen(FileName, "w");
	for (StateCount = 0; StateCount < StateTotal; StateCount++) {
		if (LineCount > 100000) {
			fprintf(FileID, "File truncated...");
			break;
		}
		//if (LineStrengths[StateCount] > LineStrengthTrunc) {
			fprintf(FileID, "%d, %f, %f, %f, %f,\n", StateCount, SEigenValues[StateCount], CPopulations[StateCount * 3], CPopulations[StateCount * 3 + 1], CPopulations[StateCount * 3 + 2]);
			LineCount++;
		//}
	}
	fclose(FileID);
}


void SaveResponse2D(const char *FileName, const ComplexType *Response, const ParmStruct *Parms, const int MemSkip, const int Minus)
{
	int PriCount, SecCount;
	char ModFileName[20];
	FILE *FileID;

	CreateFileName(ModFileName, FileName); FileID = fopen(ModFileName, "w");
	for (PriCount = 0; PriCount < *Parms->TimeTotal; PriCount++) {
		for (SecCount = 0; SecCount < *Parms->TimeTotal; SecCount++)
			fprintf(FileID, "%lf %lf ", Minus1Pow(1 + Minus) * cimag(Response[(PriCount * *Parms->TimeTotal + SecCount) * MemSkip]), Minus1Pow(1 + Minus) * creal(Response[(PriCount * *Parms->TimeTotal + SecCount) * MemSkip]));
		fprintf(FileID, "\n");
	}
	fclose(FileID);
}


////////////
// Locals //
////////////


static void SaveCrystal(const CrystalStruct *Crystal, const ParmStruct *Parms)
{
	int SiteCount;
	FILE *FileID;

	FileID = fopen("Crystal.dat", "w");
	fprintf(FileID, "Id:     Pos A:  Pos B:  Cell A: Cell B: Sub:\n");
	for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++)
		fprintf(FileID, "%6d  %6.4f  %6.4f  %6d  %6d  %6d\n", SiteCount, .5 * Crystal[SiteCount].TwoTimesPosA, .5 * Crystal[SiteCount].TwoTimesPosB, Crystal[SiteCount].CellA, Crystal[SiteCount].CellB, Crystal[SiteCount].Sub);
	fclose(FileID);
}


static RealType Distance(const RealType ADistance, const RealType BDistance, const RealType TwoCosABAngle)
{
	return (RealType) sqrt(RealPow2(ADistance) + RealPow2(BDistance) + ADistance * BDistance * TwoCosABAngle);
}


static RealType DiagonalEnergy(const BasisStruct *BasisState, const CouplingsStruct *Couplings, const ParmStruct *Parms)
// Only works for a single electron-hole pair and for a single triplet-triplet pair!
{
	RealType Value = 0;
	int PCount;
	int TripletPresent = false, Triplet_nPresent = false, ElectronPresent = -1, HolePresent = -1;

	for (PCount = 0; PCount < PTotal - 1 && (*BasisState).P[PCount].S != -1; PCount++) {
		if ((*BasisState).P[PCount].C == Triplet_n)			Triplet_nPresent = true;
		else if ((*BasisState).P[PCount].C == Triplet)		TripletPresent = true;
		else if ((*BasisState).P[PCount].C == Electron)		ElectronPresent = (*BasisState).P[PCount].S;
		else if ((*BasisState).P[PCount].C == Hole)			HolePresent = (*BasisState).P[PCount].S;
		else if ((*BasisState).P[PCount].C == Singlet)		Value += (*Couplings).DD.Array[(*BasisState).P[PCount].S * (Parms->SiteTotal + 1)];
	}

	if (Triplet_nPresent)		Value += Parms->TT_nEnergyMean;
	else if (TripletPresent)	Value += Parms->TTEnergyMean;

	if (ElectronPresent != -1)	Value += (*Couplings).EHEnergies[ElectronPresent * Parms->SiteTotal + HolePresent];

	return Value + Parms->VEnergy * VSum(*BasisState);
}


static int VSum(const BasisStruct BasisState)
{
	int PCount, Value = 0;

	for (PCount = 0; PCount < PTotal; PCount++) {
		if (BasisState.P[PCount].S == -1) break;
		Value += BasisState.P[PCount].V[BasisState.P[PCount].C];
	}

	return Value;
}


static RealType LineStrength(const int GStateId, const int PolId, const MaskRArrayStruct *GSDipoles,  const RealType *PermDipoles, const RealType *SEigenVector, const ParmStruct *Parms)
{
	RealType Contribution = 0;
	int SiteCount;

	for (SiteCount = 0; SiteCount < Parms->SiteTotal; SiteCount++)
		Contribution += PermDipoles[PolId * 2 * Parms->SiteTotal + SiteCount] * TransitionDipole(SiteCount, GStateId, GSDipoles, SEigenVector, Parms);

	return RealPow2(Contribution);
}


static RealType TransitionDipole(const int UnitId, const int GStateId, const MaskRArrayStruct *GSDipoles, const RealType *SEigenVector, const ParmStruct *Parms)
{
	RealType Output = 0;
	int SCount;

	for (SCount = 0; SCount < *Parms->SBasisBound; SCount++)
		if ((*GSDipoles).Mask[GStateId * *Parms->SBasisBound + SCount] == UnitId)
			Output += (*GSDipoles).Array[GStateId * *Parms->SBasisBound + SCount] * SEigenVector[SCount];

	return Output;
}


static void CoefIncrement(ComplexType *CoefIncrement, const ComplexType *Coef, const MaskRArrayStruct Hamiltonian, const int BasisTotal)
{
	ComplexType *CoefIncrementElement = CoefIncrement;
	const int *HamiltonianMaskElement = Hamiltonian.Mask;
	int PriCount, SecCount;

	for (PriCount = 0; PriCount < BasisTotal; PriCount++) {
		*CoefIncrementElement = 0;
		for (SecCount = 0; SecCount < BasisTotal; SecCount++)
			if (*HamiltonianMaskElement++)
				*CoefIncrementElement -= Hamiltonian.Array[PriCount * BasisTotal + SecCount] * Coef[SecCount] * I;
		CoefIncrementElement++;
	}
}


static void CoefIncrement_Redfield(ComplexType *CoefIncrement, const ComplexType *Coef, const RealType *KetEigenValues, const RealType *BraEigenValues, const MaskCArrayStruct RedfieldTensor, const int KetBasisTotal, const int BraBasisTotal, const int KetCutOff, const int BraCutOff, const RealType TimeStepAsInvEnergy)
{
	ComplexType *CoefIncrementElement = CoefIncrement;
	const ComplexType *RedfieldTensorArrayElement = RedfieldTensor.Array, *PriCoefElement = Coef, *SecCoefElement;
	const int *RedfieldTensorMaskElement = RedfieldTensor.Mask;
	int PriKetCount, PriBraCount, SecKetCount, SecBraCount, PriPreCalc;

	for (PriKetCount = 0; PriKetCount < KetBasisTotal; PriKetCount++)
		for (PriBraCount = 0; PriBraCount < BraBasisTotal; PriBraCount++) (*CoefIncrementElement++) = 0;

	CoefIncrementElement = CoefIncrement;
	for (PriKetCount = 0; PriKetCount < KetCutOff; PriKetCount++) {
		for (PriBraCount = 0; PriBraCount < BraCutOff; PriBraCount++) {

			if (KetEigenValues == BraEigenValues && PriKetCount == PriBraCount) *CoefIncrementElement = 0;
			else if (BraEigenValues == NULL)	*CoefIncrementElement = KetEigenValues[PriKetCount] * TimeStepAsInvEnergy * *PriCoefElement * -I;
			else								*CoefIncrementElement = (BraEigenValues[PriBraCount] - KetEigenValues[PriKetCount]) * TimeStepAsInvEnergy * *PriCoefElement * I;

			if (RedfieldTensorArrayElement != NULL) {
				SecCoefElement = Coef;
				for (SecKetCount = 0; SecKetCount < KetCutOff; SecKetCount++) {
					for (SecBraCount = 0; SecBraCount < BraCutOff; SecBraCount++) {
						if (*RedfieldTensorMaskElement++) *CoefIncrementElement += TimeStepAsInvEnergy * *RedfieldTensorArrayElement * *SecCoefElement;
						RedfieldTensorArrayElement++; SecCoefElement++;
					}
					SecCoefElement += BraBasisTotal - BraCutOff;
				}
			}

			PriCoefElement++; CoefIncrementElement++;
		}
		PriCoefElement += BraBasisTotal - BraCutOff; CoefIncrementElement += BraBasisTotal - BraCutOff;
	}
}


static int ClubCountCycle(int *PriClubCount, int *SecClubCount, const int *BasisBound)
{
	if (*PriClubCount == 0) {
		if (!ClubCountInc(PriClubCount, BasisBound, ClubTotal - 1)) return false;
		*SecClubCount = *PriClubCount;
		return true;
	}

	if (!ClubCountInc(SecClubCount, BasisBound, *PriClubCount + 1)) {
		*SecClubCount = 0;
		if (!ClubCountInc(PriClubCount, BasisBound, ClubTotal - 1)) return false;
		if (!ClubCountInc(SecClubCount, BasisBound, *PriClubCount + 1)) return false;
	}

	return true;
}


static int ClubCountInc(int *ClubCount, const int *BasisBound, const int Limit)
{
	int Count;

	for (Count = *ClubCount + 1; Count < Limit; Count++) {
		if (BasisBound[Count + 1] != BasisBound[Count]) {
			*ClubCount = Count;
			return true;
		}
	}
	return false;
}