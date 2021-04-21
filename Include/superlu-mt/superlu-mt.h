#ifndef SUPERLU_MT_H
#define SUPERLU_MT_H

using superlu::SuperMatrix;
using superlu::Stype_t;
using superlu::Dtype_t;
using superlu::Mtype_t;

#ifdef __cplusplus
extern "C" {
#endif

void dCreate_CompCol_Matrix(SuperMatrix*, int, int, int, double*, int*, int*, Stype_t, Dtype_t, Mtype_t);
void sCreate_CompCol_Matrix(SuperMatrix*, int, int, int, float*, int*, int*, Stype_t, Dtype_t, Mtype_t);

#ifdef __cplusplus
}
#endif

#ifdef SUANPAN_SUPERLUMT

extern int SUANPAN_NUM_THREADS;

#ifdef __cplusplus
extern "C" {
#endif

void pdgssv(int, SuperMatrix*, int*, int*, SuperMatrix*, SuperMatrix*, SuperMatrix*, int*);
void psgssv(int, SuperMatrix*, int*, int*, SuperMatrix*, SuperMatrix*, SuperMatrix*, int*);
void Destroy_CompCol_NCP(SuperMatrix*);
void Destroy_SuperNode_SCP(SuperMatrix*);

#ifdef __cplusplus
}
#endif

#endif

#endif
