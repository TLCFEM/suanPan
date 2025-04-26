#ifndef SUPERLU_MT_H
#define SUPERLU_MT_H

// superlu still needs lp64 version of BLAS
// #ifdef SUANSPAN_64BIT_INT
// typedef long long superlu_int_t;
// #else
// typedef int superlu_int_t;
// #endif
typedef int superlu_int_t;
typedef float flops_t;

#ifndef SUANPAN_SUPERLUMT

using superlu::Dtype_t;
using superlu::Mtype_t;
using superlu::Stype_t;
using superlu::superlu_options_t;
using superlu::SuperLUStat_t;
using superlu::SuperMatrix;

#else

extern int SUANPAN_NUM_THREADS;

typedef enum {
    SLU_NC,
    SLU_NCP,
    SLU_NR,
    SLU_SC,
    SLU_SCP,
    SLU_SR,
    SLU_DN,
    SLU_NR_loc
} Stype_t;

typedef enum {
    SLU_S,
    SLU_D,
    SLU_C,
    SLU_Z
} Dtype_t;

typedef enum {
    SLU_GE,
    SLU_TRLU,
    SLU_TRUU,
    SLU_TRL,
    SLU_TRU,
    SLU_SYL,
    SLU_SYU,
    SLU_HEL,
    SLU_HEU
} Mtype_t;

typedef struct {
    Stype_t Stype;
    Dtype_t Dtype;
    Mtype_t Mtype;
    superlu_int_t nrow;
    superlu_int_t ncol;
    void* Store;
} SuperMatrix;

typedef struct {
    superlu_int_t panels;
    float fcops;
    double fctime;
    superlu_int_t skedwaits;
    double skedtime;
    double cs_time;
    double spintime;
    superlu_int_t pruned;
    superlu_int_t unpruned;
} procstat_t;

typedef struct {
    superlu_int_t size;
    superlu_int_t pnum;
    double starttime;
    double fctime;
    float flopcnt;
    superlu_int_t pipewaits;
    double spintime;
} panstat_t;

typedef struct {
    flops_t flops;
    superlu_int_t nzs;
    double fctime;
} stat_relax_t;

typedef struct {
    flops_t flops;
    superlu_int_t nzs;
    double fctime;
} stat_col_t;

typedef struct {
    superlu_int_t ncols;
    flops_t flops;
    superlu_int_t nzs;
    double fctime;
} stat_snode_t;

typedef struct {
    float est;
    float pdiv;
} cp_panel_t;

typedef struct {
    float eft;
    float pmod;
} desc_eft_t;

typedef struct {
    superlu_int_t* panel_histo;
    double* utime;
    flops_t* ops;
    procstat_t* procstat;
    panstat_t* panstat;
    superlu_int_t num_panels;
    float dom_flopcnt;
    float flops_last_P_panels;
    stat_relax_t* stat_relax;
    stat_col_t* stat_col;
    stat_snode_t* stat_snode;
    superlu_int_t* panhows;
    cp_panel_t* cp_panel;
    desc_eft_t* desc_eft;
    superlu_int_t *cp_firstkid, *cp_nextkid;
    superlu_int_t* height;
    float* flops_by_height;
} Gstat_t;

typedef enum {
    NOTRANS,
    TRANS,
    CONJ
} trans_t;

#ifdef __cplusplus
extern "C" {
#endif

superlu_int_t sp_ienv(superlu_int_t);

void StatAlloc(superlu_int_t, superlu_int_t, superlu_int_t, superlu_int_t, Gstat_t*);
void StatInit(superlu_int_t, superlu_int_t, Gstat_t*);
void StatFree(Gstat_t*);

void get_perm_c(superlu_int_t, SuperMatrix*, superlu_int_t*);

void pdgssv(superlu_int_t, SuperMatrix*, superlu_int_t*, superlu_int_t*, SuperMatrix*, SuperMatrix*, SuperMatrix*, superlu_int_t*);
void psgssv(superlu_int_t, SuperMatrix*, superlu_int_t*, superlu_int_t*, SuperMatrix*, SuperMatrix*, SuperMatrix*, superlu_int_t*);
void dgstrs(trans_t, SuperMatrix*, SuperMatrix*, superlu_int_t*, superlu_int_t*, SuperMatrix*, Gstat_t*, superlu_int_t*);
void sgstrs(trans_t, SuperMatrix*, SuperMatrix*, superlu_int_t*, superlu_int_t*, SuperMatrix*, Gstat_t*, superlu_int_t*);

void Destroy_SuperMatrix_Store(SuperMatrix*);
void Destroy_CompCol_NCP(SuperMatrix*);
void Destroy_SuperNode_SCP(SuperMatrix*);

#ifdef __cplusplus
}
#endif

#endif

#ifdef __cplusplus
extern "C" {
#endif
void dCreate_CompCol_Matrix(SuperMatrix*, superlu_int_t, superlu_int_t, superlu_int_t, double*, superlu_int_t*, superlu_int_t*, Stype_t, Dtype_t, Mtype_t);
void sCreate_CompCol_Matrix(SuperMatrix*, superlu_int_t, superlu_int_t, superlu_int_t, float*, superlu_int_t*, superlu_int_t*, Stype_t, Dtype_t, Mtype_t);
void dCreate_Dense_Matrix(SuperMatrix*, superlu_int_t, superlu_int_t, double*, superlu_int_t, Stype_t, Dtype_t, Mtype_t);
void sCreate_Dense_Matrix(SuperMatrix*, superlu_int_t, superlu_int_t, float*, superlu_int_t, Stype_t, Dtype_t, Mtype_t);
#ifdef __cplusplus
}
#endif

#endif
