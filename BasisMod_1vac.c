#include "GlobalsMod.h"
#include "ToolsMod.h"
#include "BasisMod.h"

#include <string.h>


static void BuildGBasis(BasisStruct *GBasis, int *GBasisBound, const int SiteTotal, const int VTotal);
static void BuildSBasis(BasisStruct *SBasis, int *SBasisBound, const RealType *Distances, const RealType VSlope, const RealType VRadius, const RealType EHRadius, const RealType TTRadius, const int SiteTotal, const int VTotal);
static void BuildDBasis(BasisStruct *DBasis, int *DBasisBound, const RealType *Distances, const RealType EHRadius, const RealType TTRadius, const int SiteTotal, const int VTotal, const int TT_nOnly);
static void BuildK0Basis(BasisStruct *K0Basis, int *K0BasisBound, const RealType *Distances, const RealType VRadius, const RealType EHRadius, const RealType TTRadius, const int SiteTotal, const int VTotal, const int CellSiteTotal);
static void ClearBasisState(BasisStruct *BS);
static void AssignBasisState(BasisStruct *Basis, const BasisStruct BS, int *Id, const int Sort);
static void SortBasisState(BasisStruct *BS);
static void SaveGBasis(const BasisStruct *GBasis, const int *GBasisBound);
static void SaveSBasis(const char *FileName, const BasisStruct *SBasis, const int *SBasisBound);
static void SaveDBasis(const BasisStruct *SBasis, const int *SBasisBound);
static void BasisStateToString(char *Output, const BasisStruct BS);
static void BasisStateToStringAdvanced(char *Output, const BasisStruct BS);
static void ParticleToString(char *Output, int *CharCount, int *Printed, const char *Character, const int S, const int V);

static RealType FindVOverlap(const ParticleStruct Particle, RealType * const *VOverlaps, const enum Character C, const int VId, const int VTotal);
static int SaveBasisState(FILE *FileID, const BasisStruct BS);


/////////////
// Globals //
/////////////


void BuildBases(BasisStruct **GBasis, BasisStruct **SBasis, BasisStruct **DBasis, BasisStruct **K0Basis, int *GBasisBound, int *SBasisBound, int *DBasisBound, int *K0BasisBound, const RealType *Distances, const RealType VSlope, const RealType VRadius, const RealType EHRadius, const RealType TTRadius, const int SiteTotal, const int VTotal, const int CellSiteTotal, const int TT_nOnly)
{
	PrintMessage("Constructing the bases.", "", "");
	
	if (GBasis != NULL) {
		
		BuildGBasis(NULL, GBasisBound, SiteTotal, VTotal);
		*GBasis = malloc(*GBasisBound * sizeof(BasisStruct));
		BuildGBasis(*GBasis, NULL, SiteTotal, VTotal);
		
		if (!Embarrassing) SaveGBasis(*GBasis, GBasisBound);
	}
	
	if (SBasis != NULL) {
		
		BuildSBasis(NULL, SBasisBound, Distances, VSlope, VRadius, EHRadius, TTRadius, SiteTotal, VTotal);
		*SBasis = malloc(*SBasisBound * sizeof(BasisStruct));
		BuildSBasis(*SBasis, NULL, Distances, VSlope, VRadius, EHRadius, TTRadius, SiteTotal, VTotal);
		
		if (!Embarrassing) SaveSBasis("SBasis.dat", *SBasis, SBasisBound);
	}
	
	if (DBasis != NULL) {

		BuildDBasis(NULL, DBasisBound, Distances, EHRadius, TTRadius, SiteTotal, VTotal, TT_nOnly);
		*DBasis = malloc(*DBasisBound * sizeof(BasisStruct));
		BuildDBasis(*DBasis, NULL, Distances, EHRadius, TTRadius, SiteTotal, VTotal, TT_nOnly);
		
		if (!Embarrassing) SaveDBasis(*DBasis, DBasisBound);
	}
	
	if (K0Basis != NULL) {
		
		BuildK0Basis(NULL, K0BasisBound, Distances, VRadius, EHRadius, TTRadius, SiteTotal, VTotal, CellSiteTotal);
		*K0Basis = malloc(*K0BasisBound * sizeof(BasisStruct));
		BuildK0Basis(*K0Basis, NULL, Distances, VRadius, EHRadius, TTRadius, SiteTotal, VTotal, CellSiteTotal);
		
		if (!Embarrassing) SaveSBasis("K0Basis.dat", *K0Basis, K0BasisBound);
	}
}


int Present(const BasisStruct *BS, const int C)
{
	int Count;

	for (Count = 0; Count < PTotal && (*BS).P[Count].S != -1; Count++)
		if ((*BS).P[Count].C == C)
			return true;

	return false;
}


int CreateS(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < PTotal && (*BS).P[Count].S >= UnitId; Count++)
		if ((*BS).P[Count].S == UnitId) {
			if ((*BS).P[Count].C == Ground) {
				(*BS).P[Count].C = Singlet;
				return true;
			}
			else return false;
		}
	
	(*BS).P[PTotal - 1].S = UnitId;
	(*BS).P[PTotal - 1].C = Singlet;
	(*BS).P[PTotal - 1].V[Ground] = 0;
	SortBasisState(BS);
	return true;
}


int AnnihilateS(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < PTotal && (*BS).P[Count].S >= UnitId; Count++)
		if ((*BS).P[Count].S == UnitId && (*BS).P[Count].C == Singlet) {
			(*BS).P[Count].C = Ground;
			return true;
		}
	
	return false;
}


int CreateE(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < PTotal && (*BS).P[Count].S >= UnitId; Count++)
		if ((*BS).P[Count].S == UnitId) {
			if ((*BS).P[Count].C == Ground) {
				(*BS).P[Count].C = Electron;
				return true;
			}
			else if ((*BS).P[Count].C == Hole) {
				(*BS).P[Count].C = Singlet;
				return true;
			}
			else return false;
		}
	
	(*BS).P[PTotal - 1].S = UnitId;
	(*BS).P[PTotal - 1].C = Electron;
	(*BS).P[PTotal - 1].V[Ground] = 0;
	SortBasisState(BS);
	return true;
}


int CreateH(BasisStruct *BS, const int UnitId)
{
	int Count;
	
	for (Count = 0; Count < PTotal && (*BS).P[Count].S >= UnitId; Count++)
		if ((*BS).P[Count].S == UnitId) {
			if ((*BS).P[Count].C == Ground) {
				(*BS).P[Count].C = Hole;
				return true;
			}
			else if ((*BS).P[Count].C == Electron) {
				(*BS).P[Count].C = Singlet;
				return true;
			}
			else return false;
		}
	
	(*BS).P[PTotal - 1].S = UnitId;
	(*BS).P[PTotal - 1].C = Hole;
	(*BS).P[PTotal - 1].V[Ground] = 0;
	SortBasisState(BS);
	return true;
}


int CreateT(BasisStruct *BS, const int UnitId)
{
	int Count;

	for (Count = 0; Count < PTotal && (*BS).P[Count].S >= UnitId; Count++)
		if ((*BS).P[Count].S == UnitId) {
			if ((*BS).P[Count].C == Ground) {
				(*BS).P[Count].C = Triplet;
				return true;
			}
			else return false;
		}

	(*BS).P[PTotal - 1].S = UnitId;
	(*BS).P[PTotal - 1].C = Triplet;
	(*BS).P[PTotal - 1].V[Ground] = 0;
	SortBasisState(BS);
	return true;
}


int DipoleCreateT_n(BasisStruct *BS, const int UnitId)
{
	int Count;

	for (Count = 0; Count < PTotal && (*BS).P[Count].S >= UnitId; Count++)
		if ((*BS).P[Count].S == UnitId) {
			if ((*BS).P[Count].C == Triplet) {
				(*BS).P[Count].C = Triplet_n;
				return true;
			}
			else return false;
		}

	return false;
}


int CreateT_n(BasisStruct *BS, const int UnitId)
{
	int Count;

	for (Count = 0; Count < PTotal && (*BS).P[Count].S >= UnitId; Count++)
		if ((*BS).P[Count].S == UnitId) {
			if ((*BS).P[Count].C == Ground) {
				(*BS).P[Count].C = Triplet_n;
				return true;
			}
			else return false;
		}

	(*BS).P[PTotal - 1].S = UnitId;
	(*BS).P[PTotal - 1].C = Triplet_n;
	(*BS).P[PTotal - 1].V[Ground] = 0;
	SortBasisState(BS);
	return true;
}


int AnnihilateT_n(BasisStruct *BS, const int UnitId)
{
	int Count;

	for (Count = 0; Count < PTotal && (*BS).P[Count].S >= UnitId; Count++)
		if ((*BS).P[Count].S == UnitId && (*BS).P[Count].C == Triplet_n) {
			(*BS).P[Count].C = Triplet;
			return true;
		}

	return false;
}


int CountParticles(const BasisStruct *BS)
{
	int Value;

	for (Value = 0; Value< PTotal; Value++)
		if ((*BS).P[Value].S == -1) break;

	return Value;
}


int CountVibrations(const BasisStruct *BS)
{
	int PCount, Value = 0;

	for (PCount = 0; PCount< PTotal; PCount++) {
		if ((*BS).P[PCount].S == -1) break;
		Value += (*BS).P[PCount].V[(*BS).P[PCount].C];
	}

	return Value;
}


int IdentifyBasisState(RealType *EffectiveVOverlap, const BasisStruct *BS, const BasisStruct *RefBS, RealType * const *VOverlaps, const int VTotal)
{
	int Count, RefCount, Offset = 0, VOverlapNonUnity = false;
	
	for (Count = 0; Count < PTotal; Count++) {
		RefCount = Count - Offset;
		if ((*BS).P[Count].S == -1 && (*RefBS).P[RefCount].S == -1) break;
		else if ((*BS).P[Count].S != (*RefBS).P[RefCount].S) {
			if ((*BS).P[Count].C == Ground && (*BS).P[Count].V[Ground] == -1) Offset++, VOverlapNonUnity = true;
			else return false;
		}
		else if ((*BS).P[Count].C != (*RefBS).P[RefCount].C) return false;
		else if ((*BS).P[Count].V[(*BS).P[Count].C] != (*RefBS).P[RefCount].V[(*RefBS).P[RefCount].C]) {
			if ((*BS).P[Count].V[(*BS).P[Count].C] == -1) VOverlapNonUnity = true;
			else return false;
		}
	}
	*EffectiveVOverlap = 1;
	
	if (VOverlapNonUnity) {
		Offset = 0;
		for (Count = 0; Count < PTotal; Count++) {
			RefCount = Count - Offset;
			if ((*BS).P[Count].S == -1 && (*RefBS).P[RefCount].S == -1) break;
			else if ((*BS).P[Count].S != (*RefBS).P[RefCount].S) {
				if ((*BS).P[Count].C == Ground && (*BS).P[Count].V[Ground] == -1) {
					*EffectiveVOverlap *= FindVOverlap((*BS).P[Count], VOverlaps, Ground, 0, VTotal);
					Offset++;
				}
				else return false;
			}
			else if ((*BS).P[Count].C != (*RefBS).P[RefCount].C) return false;
			else if ((*BS).P[Count].V[(*BS).P[Count].C] != (*RefBS).P[RefCount].V[(*RefBS).P[RefCount].C]) {
				if ((*BS).P[Count].V[(*BS).P[Count].C] == -1)
					*EffectiveVOverlap *= FindVOverlap((*BS).P[Count], VOverlaps, (*RefBS).P[RefCount].C, (*RefBS).P[RefCount].V[(*RefBS).P[RefCount].C], VTotal);
				else return false;
			}
		}
	}
	
	return true;
}


int IdentifyFirstTwoVibrations(const BasisStruct *BS, const BasisStruct *RefBS)
{
	int CCount;
	
	for (CCount = 0; CCount < CTotal; CCount++) {
		if ((*BS).P[0].V[CCount] != (*RefBS).P[0].V[CCount]) return false;
		if ((*BS).P[1].V[CCount] != (*RefBS).P[1].V[CCount]) return false;
	}
	
	return true;
}


int IdentifyPurelyElectronic(const BasisStruct *BS, const BasisStruct *RefBS)
{
	int Count = 0, RefCount = 0;

	while (Count < PTotal && RefCount < PTotal) {
		if ((*BS).P[Count].S == -1 && (*RefBS).P[RefCount].S == -1) break;
		else if ((*BS).P[Count].C == Ground) Count++;
		else if ((*RefBS).P[RefCount].C == Ground) RefCount++;
		else if ((*BS).P[Count].S != (*RefBS).P[RefCount].S || (*BS).P[Count].C != (*RefBS).P[RefCount].C) return false;
		else Count++, RefCount++;
	}

	return true;
}


void SwitchFirstTwoParticles(BasisStruct *BS)
{
	ParticleStruct *WorkParticle = malloc(sizeof(ParticleStruct));
	
	memcpy(WorkParticle, &(*BS).P[0], sizeof(ParticleStruct));
	memcpy(&(*BS).P[0], &(*BS).P[1], sizeof(ParticleStruct));
	memcpy(&(*BS).P[1], WorkParticle, sizeof(ParticleStruct));
	
	free(WorkParticle);
}


void PrintBasisState(const BasisStruct BS)
{
	char PrintLine[100];
	
	BasisStateToString(PrintLine, BS);
	printf("%s\n", PrintLine);
}


void PrintBasisStateAdvanced(const BasisStruct BS)
{
	char PrintLine[100];

	BasisStateToStringAdvanced(PrintLine, BS);
	printf("%s\n", PrintLine);
}


void SaveOperator(const char *FileName, const RealType *Operator, const BasisStruct *PriBasis, const BasisStruct *SecBasis, const int PriDimension, const int SecDimension)
{
	char ModFileName[20];
	int PriCount, SecCount, PrintCount;
	FILE *FileID;
	
	if (Embarrassing) return;
	
	sprintf(ModFileName, "%s%s", FileName, ".dat");
	
	FileID = fopen(ModFileName, "w");
	
	if (IntMax(PriDimension, SecDimension) > 70) fprintf(FileID, "Operator too large for storage!");
	else if (IntMin(PriDimension, SecDimension) == 0) fprintf(FileID, "Operator empty!");
	else {
		
		fprintf(FileID, "                ");
		for (PriCount = 0; PriCount < PriDimension; PriCount++) {
			if (PriBasis == NULL) fprintf(FileID, "%2d               ", PriCount);
			else {
				fprintf(FileID, "|");
				PrintCount = SaveBasisState(FileID, PriBasis[PriCount]);
				fprintf(FileID, ">");
				for (; PrintCount < 15; PrintCount++) fprintf(FileID, " ");
			}
		}
		fprintf(FileID, "\n");
		for (SecCount = 0; SecCount < SecDimension; SecCount++) {
			if (SecBasis == NULL) fprintf(FileID, "%2d               ", SecCount);
			else {
				fprintf(FileID, "|");
				PrintCount = SaveBasisState(FileID, SecBasis[SecCount]);
				fprintf(FileID, ">");
				for (; PrintCount < 15; PrintCount++) fprintf(FileID, " ");
			}
			for (PriCount = 0; PriCount < PriDimension; PriCount++) fprintf(FileID, "%10.3f       ", Operator[PriCount * SecDimension + SecCount]);
			fprintf(FileID, "\n");
		}
		
	}
	
	fclose(FileID);
}


////////////
// Locals //
////////////


static void BuildGBasis(BasisStruct *GBasis, int *GBasisBound, const int SiteTotal, const int VTotal)
{
	BasisStruct BS;
	int Id = 0;

	if (GBasisBound != NULL) GBasisBound[1] = Id;

	ClearBasisState(&BS);
	
	AssignBasisState(GBasis, BS, &Id, false);

	if (GBasisBound != NULL) GBasisBound[2] = Id;

	BS.P[0].C = Ground;
	for (BS.P[0].S = 0; BS.P[0].S < SiteTotal; BS.P[0].S++)
		for (BS.P[0].V[Ground] = 1; BS.P[0].V[Ground] < VTotal; BS.P[0].V[Ground]++)
			AssignBasisState(GBasis, BS, &Id, false);
	
	if (GBasisBound != NULL) GBasisBound[3] = Id;
	
	BS.P[1].C = Ground;
	for (BS.P[0].S = 1; BS.P[0].S < SiteTotal; BS.P[0].S++)
		for (BS.P[1].S = 0; BS.P[1].S < BS.P[0].S; BS.P[1].S++)
			for (BS.P[0].V[Ground] = 1; BS.P[0].V[Ground] < VTotal - 1; BS.P[0].V[Ground]++)
				for (BS.P[1].V[Ground] = 1; BS.P[1].V[Ground] < VTotal - BS.P[0].V[Ground]; BS.P[1].V[Ground]++)
					AssignBasisState(GBasis, BS, &Id, false);
	
	if (GBasisBound != NULL) GBasisBound[4] = Id;
	if (GBasisBound != NULL) GBasisBound[5] = Id;
	if (GBasisBound != NULL) GBasisBound[6] = Id;
	if (GBasisBound != NULL) GBasisBound[7] = Id;
	if (GBasisBound != NULL) GBasisBound[0] = Id;
}


static void BuildSBasis(BasisStruct *SBasis, int *SBasisBound, const RealType *Distances, const RealType VSlope, const RealType VRadius, const RealType EHRadius, const RealType TTRadius, const int SiteTotal, const int VTotal)
{
	BasisStruct BS;
	int Id = 0;

	if (SBasisBound != NULL) SBasisBound[1] = Id;

	ClearBasisState(&BS);
	
	BS.P[0].C = Singlet;
	for (BS.P[0].S = 0; BS.P[0].S < SiteTotal; BS.P[0].S++)
		for (BS.P[0].V[Singlet] = 0; BS.P[0].V[Singlet] < VTotal; BS.P[0].V[Singlet]++)
			AssignBasisState(SBasis, BS, &Id, false);
	
	if (SBasisBound != NULL) SBasisBound[2] = Id;

	BS.P[1].C = Ground;
	for (BS.P[0].S = 0; BS.P[0].S < SiteTotal; BS.P[0].S++)
		for (BS.P[1].S = 0; BS.P[1].S < SiteTotal; BS.P[1].S++) {
			if (BS.P[0].S == BS.P[1].S || (VRadius >= 0 && Distances[BS.P[0].S * SiteTotal + BS.P[1].S] > VRadius)) continue;
			for (BS.P[0].V[Singlet] = 0; BS.P[0].V[Singlet] < VTotal - 1; BS.P[0].V[Singlet]++)
				for (BS.P[1].V[Ground] = 1; BS.P[1].V[Ground] < VTotal - BS.P[0].V[Singlet]; BS.P[1].V[Ground]++)
					if (!VSlope || BS.P[0].V[Singlet] + BS.P[1].V[Ground] < VTotal - (int) (Distances[BS.P[0].S * SiteTotal + BS.P[1].S] / VSlope))
						AssignBasisState(SBasis, BS, &Id, true);
		}
	
	if (SBasisBound != NULL) SBasisBound[3] = Id;

	ClearBasisState(&BS);

	BS.P[0].C = Electron;
	BS.P[1].C = Hole;
	for (BS.P[0].S = 0; BS.P[0].S < SiteTotal; BS.P[0].S++)
		for (BS.P[1].S = 0; BS.P[1].S < SiteTotal; BS.P[1].S++) {
			if (BS.P[0].S == BS.P[1].S || (EHRadius >= 0 && Distances[BS.P[0].S * SiteTotal + BS.P[1].S] > EHRadius)) continue;
			for (BS.P[0].V[Electron] = 0; BS.P[0].V[Electron] < VTotal; BS.P[0].V[Electron]++)
				for (BS.P[1].V[Hole] = 0; BS.P[1].V[Hole] < VTotal - BS.P[0].V[Electron]; BS.P[1].V[Hole]++)
					AssignBasisState(SBasis, BS, &Id, true);
		}
	
	if (SBasisBound != NULL) SBasisBound[4] = Id;

	ClearBasisState(&BS);

	BS.P[0].C = Triplet;
	BS.P[1].C = Triplet;
	for (BS.P[0].S = 1; BS.P[0].S < SiteTotal; BS.P[0].S++)
		for (BS.P[1].S = 0; BS.P[1].S < BS.P[0].S; BS.P[1].S++) {
			if (TTRadius >= 0 && Distances[BS.P[0].S * SiteTotal + BS.P[1].S] > TTRadius) continue;
			for (BS.P[0].V[Triplet] = 0; BS.P[0].V[Triplet] < VTotal; BS.P[0].V[Triplet]++)
				for (BS.P[1].V[Triplet] = 0; BS.P[1].V[Triplet] < VTotal - BS.P[0].V[Triplet]; BS.P[1].V[Triplet]++)
					AssignBasisState(SBasis, BS, &Id, false);
		}

	if (SBasisBound != NULL) SBasisBound[5] = Id;
	if (SBasisBound != NULL) SBasisBound[6] = Id;
	if (SBasisBound != NULL) SBasisBound[7] = Id;
	if (SBasisBound != NULL) SBasisBound[0] = Id;
}


static void BuildK0Basis(BasisStruct *K0Basis, int *K0BasisBound, const RealType *Distances, const RealType VRadius, const RealType EHRadius, const RealType TTRadius, const int SiteTotal, const int VTotal, const int CellSiteTotal)
{
	BasisStruct BS;
	int Id = 0;

	if (K0BasisBound != NULL) K0BasisBound[1] = Id;

	ClearBasisState(&BS);
	
	BS.P[0].C = Singlet;
	for (BS.P[0].S = 0; BS.P[0].S < CellSiteTotal; BS.P[0].S++)
		for (BS.P[0].V[Singlet] = 0; BS.P[0].V[Singlet] < VTotal; BS.P[0].V[Singlet]++)
			AssignBasisState(K0Basis, BS, &Id, false);
	
	if (K0BasisBound != NULL) K0BasisBound[2] = Id;
	
	BS.P[0].C = Singlet;
	BS.P[1].C = Ground;
	for (BS.P[0].S = 0; BS.P[0].S < CellSiteTotal; BS.P[0].S++)
		for (BS.P[1].S = 1; BS.P[1].S < SiteTotal; BS.P[1].S++) {
			if (VRadius >= 0 && Distances[BS.P[1].S] > VRadius) continue;
			for (BS.P[0].V[Singlet] = 0; BS.P[0].V[Singlet] < VTotal; BS.P[0].V[Singlet]++)
				for (BS.P[1].V[Ground] = 1; BS.P[1].V[Ground] < VTotal - BS.P[0].V[Singlet]; BS.P[1].V[Ground]++)
					AssignBasisState(K0Basis, BS, &Id, false);
		}
	
	if (K0BasisBound != NULL) K0BasisBound[3] = Id;

	ClearBasisState(&BS);

	BS.P[0].C = Hole;
	BS.P[1].C = Electron;
	for (BS.P[0].S = 0; BS.P[0].S < CellSiteTotal; BS.P[0].S++)
		for (BS.P[1].S = 1; BS.P[1].S < SiteTotal; BS.P[1].S++) {
			if (EHRadius >= 0 && Distances[BS.P[1].S] > EHRadius) continue;
			for (BS.P[0].V[Hole] = 0; BS.P[0].V[Hole] < VTotal; BS.P[0].V[Hole]++)
				for (BS.P[1].V[Electron] = 0; BS.P[1].V[Electron] < VTotal - BS.P[0].V[Hole]; BS.P[1].V[Electron]++)
					AssignBasisState(K0Basis, BS, &Id, false);
		}
	
	if (K0BasisBound != NULL) K0BasisBound[4] = Id;

	ClearBasisState(&BS);

	BS.P[0].C = Triplet;
	BS.P[1].C = Triplet;
	BS.P[0].S = 0;
	if (CellSiteTotal == 2)
		for (BS.P[1].S = 1; BS.P[1].S < SiteTotal; BS.P[1].S += 2) {
			if (TTRadius >= 0 && Distances[BS.P[1].S] > TTRadius) continue;
			for (BS.P[0].V[Triplet] = 0; BS.P[0].V[Triplet] < VTotal; BS.P[0].V[Triplet]++)
				for (BS.P[1].V[Triplet] = 0; BS.P[1].V[Triplet] < VTotal - BS.P[0].V[Triplet]; BS.P[1].V[Triplet]++)
					AssignBasisState(K0Basis, BS, &Id, false);
		}

	for (BS.P[0].S = 0; BS.P[0].S < CellSiteTotal; BS.P[0].S++)
		for (BS.P[1].S = CellSiteTotal; BS.P[1].S < SiteTotal; BS.P[1].S += CellSiteTotal) {
			if ((TTRadius >= 0 && Distances[BS.P[1].S] > TTRadius) || !Distances[IntPow2(SiteTotal) + BS.P[1].S]) continue;
			for (BS.P[0].V[Triplet] = 0; BS.P[0].V[Triplet] < VTotal; BS.P[0].V[Triplet]++)
				for (BS.P[1].V[Triplet] = 0; BS.P[1].V[Triplet] < VTotal - BS.P[0].V[Triplet]; BS.P[1].V[Triplet]++)
					AssignBasisState(K0Basis, BS, &Id, false);
		}

	if (K0BasisBound != NULL) K0BasisBound[5] = Id;
	if (K0BasisBound != NULL) K0BasisBound[6] = Id;
	if (K0BasisBound != NULL) K0BasisBound[7] = Id;
	if (K0BasisBound != NULL) K0BasisBound[0] = Id;
}


static void BuildDBasis(BasisStruct *DBasis, int *DBasisBound, const RealType *Distances, const RealType EHRadius, const RealType TTRadius, const int SiteTotal, const int VTotal, const int TT_nOnly)
{
	BasisStruct BS;
	int Id = 0;

	if (DBasisBound != NULL) DBasisBound[1] = Id;

	ClearBasisState(&BS);

	if (!TT_nOnly) {
		BS.P[0].C = Singlet;
		BS.P[1].C = Singlet;
		for (BS.P[0].S = 1; BS.P[0].S < SiteTotal; BS.P[0].S++)
			for (BS.P[1].S = 0; BS.P[1].S < BS.P[0].S; BS.P[1].S++)
				for (BS.P[0].V[Singlet] = 0; BS.P[0].V[Singlet] < VTotal; BS.P[0].V[Singlet]++)
					for (BS.P[1].V[Singlet] = 0; BS.P[1].V[Singlet] < VTotal - BS.P[0].V[Singlet]; BS.P[1].V[Singlet]++)
						AssignBasisState(DBasis, BS, &Id, false);
	}

	if (DBasisBound != NULL) DBasisBound[2] = Id;

	if (!TT_nOnly) {
		BS.P[0].C = Singlet;
		BS.P[1].C = Electron;
		BS.P[2].C = Hole;
		for (BS.P[0].S = 0; BS.P[0].S < SiteTotal; BS.P[0].S++)
			for (BS.P[1].S = 0; BS.P[1].S < SiteTotal; BS.P[1].S++)
				for (BS.P[2].S = 0; BS.P[2].S < SiteTotal; BS.P[2].S++) {
					if (BS.P[0].S == BS.P[1].S || BS.P[0].S == BS.P[2].S || BS.P[1].S == BS.P[2].S || (EHRadius >= 0 && Distances[BS.P[1].S * SiteTotal + BS.P[2].S] > EHRadius)) continue;
					for (BS.P[0].V[Singlet] = 0; BS.P[0].V[Singlet] < VTotal; BS.P[0].V[Singlet]++)
						for (BS.P[1].V[Electron] = 0; BS.P[1].V[Electron] < VTotal - BS.P[0].V[Singlet]; BS.P[1].V[Electron]++)
							for (BS.P[1].V[Hole] = 0; BS.P[1].V[Hole] < VTotal - BS.P[0].V[Singlet] - BS.P[1].V[Electron]; BS.P[1].V[Hole]++)
								AssignBasisState(DBasis, BS, &Id, true);
				}
	}

	if (DBasisBound != NULL) DBasisBound[3] = Id;

	if (!TT_nOnly) {
		BS.P[0].C = Singlet;
		BS.P[1].C = Triplet;
		BS.P[2].C = Triplet;
		for (BS.P[0].S = 0; BS.P[0].S < SiteTotal; BS.P[0].S++)
			for (BS.P[1].S = 1; BS.P[1].S < SiteTotal; BS.P[1].S++)
				for (BS.P[2].S = 0; BS.P[2].S < BS.P[1].S; BS.P[2].S++) {
					if (BS.P[0].S == BS.P[1].S || BS.P[0].S == BS.P[2].S || (TTRadius >= 0 && Distances[BS.P[1].S * SiteTotal + BS.P[2].S] > TTRadius)) continue;
					for (BS.P[0].V[Singlet] = 0; BS.P[0].V[Singlet] < VTotal; BS.P[0].V[Singlet]++)
						for (BS.P[1].V[Triplet] = 0; BS.P[1].V[Triplet] < VTotal - BS.P[0].V[Singlet]; BS.P[1].V[Triplet]++)
							for (BS.P[1].V[Triplet] = 0; BS.P[1].V[Triplet] < VTotal - BS.P[0].V[Singlet] - BS.P[1].V[Triplet]; BS.P[1].V[Triplet]++)
								AssignBasisState(DBasis, BS, &Id, true);
				}
	}

	if (DBasisBound != NULL) DBasisBound[4] = Id;

	ClearBasisState(&BS);

	BS.P[0].C = Triplet;
	BS.P[1].C = Triplet_n;
	for (BS.P[0].S = 0; BS.P[0].S < SiteTotal; BS.P[0].S++)
		for (BS.P[1].S = 0; BS.P[1].S < SiteTotal; BS.P[1].S++) {
			if (BS.P[0].S == BS.P[1].S || (TTRadius >= 0 && Distances[BS.P[0].S * SiteTotal + BS.P[1].S] > TTRadius)) continue;
			for (BS.P[0].V[Triplet] = 0; BS.P[0].V[Triplet] < VTotal; BS.P[0].V[Triplet]++)
				for (BS.P[1].V[Triplet_n] = 0; BS.P[1].V[Triplet_n] < VTotal - BS.P[0].V[Triplet]; BS.P[1].V[Triplet_n]++)
					AssignBasisState(DBasis, BS, &Id, true);
		}

	if (DBasisBound != NULL) DBasisBound[5] = Id;
	if (DBasisBound != NULL) DBasisBound[6] = Id;
	if (DBasisBound != NULL) DBasisBound[7] = Id;
	if (DBasisBound != NULL) DBasisBound[0] = Id;
}


static void ClearBasisState(BasisStruct *BS)
{
	memset(BS, -1, sizeof(BasisStruct));
}


static void AssignBasisState(BasisStruct *Basis, const BasisStruct BS, int *Id, const int Sort)
{
	int PCount;
	
	for (PCount = 0; PCount < PTotal; PCount++)
		if (BS.P[PCount].S == 0) return;

	if (Basis != NULL) {
		memcpy(&Basis[*Id], &BS, sizeof(BasisStruct));
		if (Sort) SortBasisState(&Basis[*Id]);
	}
	(*Id)++;
}


static void SortBasisState(BasisStruct *BS)
{
	ParticleStruct *WorkParticle = malloc(sizeof(ParticleStruct));
	int PriCount, SecCount;
	
	for (PriCount = 1; PriCount < PTotal; PriCount++) {
		memcpy(WorkParticle, &(*BS).P[PriCount], sizeof(ParticleStruct));
		
		for (SecCount = PriCount - 1; SecCount >= 0 && (*BS).P[SecCount].S < (*WorkParticle).S; SecCount--)
			memcpy(&(*BS).P[SecCount + 1], &(*BS).P[SecCount], sizeof(ParticleStruct));
		
		memcpy(&(*BS).P[SecCount + 1], WorkParticle, sizeof(ParticleStruct));
	}
	
	free(WorkParticle);
}


static RealType FindVOverlap(const ParticleStruct Particle, RealType * const *VOverlaps, const enum Character C, const int VId, const int VTotal)
{
	int CCount;
	
	for (CCount = 0; CCount < CTotal; CCount++) {
		if (CCount == C) continue;
		if (Particle.V[CCount] != -1) return VOverlaps[C * CTotal + CCount][VId * VTotal + Particle.V[CCount]];
	}
	return 0;
}


static void SaveGBasis(const BasisStruct *GBasis, const int *GBasisBound)
{
	int Count;
	FILE *FileID;

	FileID = fopen("GBasis.dat", "w");

	fprintf(FileID, "Total: %d\n", *GBasisBound);

	if (*GBasisBound > 5000) {
		fprintf(FileID, "\nBasis too large for storage!\n");
		fclose(FileID);
		return;
	}

	fprintf(FileID, "\n   Vacuum state\n");
	for (Count = 0; Count < *GBasisBound; Count++) {
		if (Count == GBasisBound[3])		fprintf(FileID, "\n   Two vibrations: %d\n", GBasisBound[4] - GBasisBound[3]);
		else if (Count == GBasisBound[2])	fprintf(FileID, "\n   Single vibration: %d\n", GBasisBound[3] - GBasisBound[2]);
		fprintf(FileID, "%8d: |", Count);
		SaveBasisState(FileID, GBasis[Count]);
		fprintf(FileID, ">\n");
	}
	
	fclose(FileID);
}


static void SaveSBasis(const char *FileName, const BasisStruct *SBasis, const int *SBasisBound)
{
	int Count;
	FILE *FileID;

	FileID = fopen(FileName, "w");

	fprintf(FileID, "Total: %d\n", *SBasisBound);

	if (*SBasisBound > 5000) {
		fprintf(FileID, "\nBasis too large for storage!\n");
		fclose(FileID);
		return;
	}

	for (Count = 0; Count < *SBasisBound; Count++) {
		if (Count == SBasisBound[4])	fprintf(FileID, "\n   Triplet-triplet pairs: %d\n", SBasisBound[5] - SBasisBound[4]);
		else if (Count == SBasisBound[3])	fprintf(FileID, "\n   Charge-transfer states: %d\n", SBasisBound[4] - SBasisBound[3]);
		else if (Count == SBasisBound[2])	fprintf(FileID, "\n   Vibronic + vibration: %d\n", SBasisBound[3] - SBasisBound[2]);
		else if (Count == SBasisBound[1])	fprintf(FileID, "\n   Purely vibronic: %d\n", SBasisBound[2] - SBasisBound[1]);
		fprintf(FileID, "%8d: |", Count);
		SaveBasisState(FileID, SBasis[Count]);
		fprintf(FileID, ">\n");
	}
	
	fclose(FileID);
}


static void SaveDBasis(const BasisStruct *DBasis, const int *DBasisBound)
{
	int Count;
	FILE *FileID;

	FileID = fopen("DBasis.dat", "w");

	fprintf(FileID, "Total: %d\n", *DBasisBound);

	if (*DBasisBound > 5000) {
		fprintf(FileID, "\nBasis too large for storage!\n");
		fclose(FileID);
		return;
	}

	for (Count = 0; Count < *DBasisBound; Count++) {
		if (Count == DBasisBound[4])		fprintf(FileID, "\n   Higher-lying triplet-triplet pairs: %d\n", DBasisBound[5] - DBasisBound[4]);
		else if (Count == DBasisBound[3])	fprintf(FileID, "\n   Singlets and triplet-triplet pairs: %d\n", DBasisBound[4] - DBasisBound[3]);
		else if (Count == DBasisBound[2])	fprintf(FileID, "\n   Singlets and charge transfer states: %d\n", DBasisBound[3] - DBasisBound[2]);
		else if (Count == DBasisBound[1])	fprintf(FileID, "\n   Double singlets: %d\n", DBasisBound[2] - DBasisBound[1]);
		fprintf(FileID, "%8d: |", Count);
		SaveBasisState(FileID, DBasis[Count]);
		fprintf(FileID, ">\n");
	}

	fclose(FileID);
}


static int SaveBasisState(FILE *FileID, const BasisStruct BS)
{
	char PrintLine[100];
	
	BasisStateToString(PrintLine, BS);
	fprintf(FileID, "%s", PrintLine);
	
	return (int) strlen(PrintLine);
}


static void BasisStateToString(char *Output, const BasisStruct BS)
{
	int PCount, CharCount = 2, Printed = false;

	for (PCount = 0; PCount < PTotal; PCount++) {
		if (BS.P[PCount].C == Ground)			ParticleToString(Output, &CharCount, &Printed, "g", BS.P[PCount].S, BS.P[PCount].V[Ground]);
		else if (BS.P[PCount].C == Singlet)		ParticleToString(Output, &CharCount, &Printed, "s", BS.P[PCount].S, BS.P[PCount].V[Singlet]);
		else if (BS.P[PCount].C == Electron)	ParticleToString(Output, &CharCount, &Printed, "e", BS.P[PCount].S, BS.P[PCount].V[Electron]);
		else if (BS.P[PCount].C == Hole)		ParticleToString(Output, &CharCount, &Printed, "h", BS.P[PCount].S, BS.P[PCount].V[Hole]);
		else if (BS.P[PCount].C == Triplet)		ParticleToString(Output, &CharCount, &Printed, "t", BS.P[PCount].S, BS.P[PCount].V[Triplet]);
		else if (BS.P[PCount].C == Triplet_n)	ParticleToString(Output, &CharCount, &Printed, "tn", BS.P[PCount].S, BS.P[PCount].V[Triplet_n]), CharCount++;
	}
	if (!Printed) sprintf(Output, "VAC"), CharCount += 3;
}


static void BasisStateToStringAdvanced(char *Output, const BasisStruct BS)
{
	int PCount, VId = 0, CharCount = 4, Printed = false;
	char Char[4];

	for (PCount = 0; PCount < PTotal; PCount++) {
		if (BS.P[PCount].C == Ground) 			sprintf(Char, "g/");
		else if (BS.P[PCount].C == Singlet)		sprintf(Char, "s/");
		else if (BS.P[PCount].C == Electron)	sprintf(Char, "e/");
		else if (BS.P[PCount].C == Hole)		sprintf(Char, "h/");
		else if (BS.P[PCount].C == Triplet)		sprintf(Char, "t/");
		else if (BS.P[PCount].C == Triplet_n)	sprintf(Char, "tn/"), CharCount++;
		else break;

		if (BS.P[PCount].V[Ground] != -1) 			sprintf(Char, "%sg", Char), VId = Ground;
		else if (BS.P[PCount].V[Singlet] != -1)		sprintf(Char, "%ss", Char), VId = Singlet;
		else if (BS.P[PCount].V[Electron] != -1)	sprintf(Char, "%se", Char), VId = Electron;
		else if (BS.P[PCount].V[Hole] != -1)		sprintf(Char, "%sh", Char), VId = Hole;
		else if (BS.P[PCount].V[Triplet] != -1)		sprintf(Char, "%st", Char), VId = Triplet;
		else if (BS.P[PCount].V[Triplet_n] != -1)	sprintf(Char, "%stn", Char), VId = Triplet_n, CharCount++;

		ParticleToString(Output, &CharCount, &Printed, Char, BS.P[PCount].S, BS.P[PCount].V[VId]);
	}
	if (!Printed) sprintf(Output, "VAC"), CharCount += 1;
}


static void ParticleToString(char *Output, int *CharCount, int *Printed, const char *Character, const int S, const int V)
{
	if (*Printed)	sprintf(Output, "%s;%s%d,%d", Output, Character, S, V), (*CharCount) += 5;
	else			sprintf(Output, "%s%d,%d", Character, S, V), (*CharCount) += 4;
	if (S > 9) (*CharCount)++;
	if (S > 99) (*CharCount)++;
	*Printed = true;
}
