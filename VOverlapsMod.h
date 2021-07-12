#ifndef VOverlapsMod
#define VOverlapsMod

void BuildVOverlaps(RealType ***VOverlaps, const RealType *HuangRhys, const int Characters, const int VTotal, const int Embarrassing);
void FreeVOverlaps(RealType ***VOverlaps, const int Characters);

#endif