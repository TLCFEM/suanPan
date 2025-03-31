#ifndef SUPERLU_MT_H
#define SUPERLU_MT_H

#ifndef SUANPAN_SUPERLUMT

using superlu::Dtype_t;
using superlu::GlobalLU_t;
using superlu::mem_usage_t;
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
    int nrow;
    int ncol;
    void* Store;
} SuperMatrix;

typedef struct {
    int panels;
    float fcops;
    double fctime;
    int skedwaits;
    double skedtime;
    double cs_time;
    double spintime;
    int pruned;
    int unpruned;
} procstat_t;

typedef struct {
    int size;
    int pnum;
    double starttime;
    double fctime;
    float flopcnt;
    int pipewaits;
    double spintime;
} panstat_t;

typedef struct {
    float flops;
    int nzs;
    double fctime;
} stat_relax_t;

typedef struct {
    float flops;
    int nzs;
    double fctime;
} stat_col_t;

typedef struct {
    int ncols;
    float flops;
    int nzs;
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
    int* panel_histo;
    double* utime;
    float* ops;
    procstat_t* procstat;
    panstat_t* panstat;
    int num_panels;
    float dom_flopcnt;
    float flops_last_P_panels;
    stat_relax_t* stat_relax;
    stat_col_t* stat_col;
    stat_snode_t* stat_snode;
    int* panhows;
    cp_panel_t* cp_panel;
    desc_eft_t* desc_eft;
    int *cp_firstkid, *cp_nextkid;
    int* height;
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

int sp_ienv(int);
void* superlu_malloc(size_t);
void superlu_free(void*);
void get_perm_c(int, SuperMatrix*, int*);
void pdgssv(int, SuperMatrix*, int*, int*, SuperMatrix*, SuperMatrix*, SuperMatrix*, int*);
void psgssv(int, SuperMatrix*, int*, int*, SuperMatrix*, SuperMatrix*, SuperMatrix*, int*);
void dgstrs(trans_t, SuperMatrix*, SuperMatrix*, int*, int*, SuperMatrix*, Gstat_t*, int*);
void sgstrs(trans_t, SuperMatrix*, SuperMatrix*, int*, int*, SuperMatrix*, Gstat_t*, int*);
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

void dCreate_CompCol_Matrix(SuperMatrix*, int, int, int, double*, int*, int*, Stype_t, Dtype_t, Mtype_t);
void sCreate_CompCol_Matrix(SuperMatrix*, int, int, int, float*, int*, int*, Stype_t, Dtype_t, Mtype_t);
void dCreate_Dense_Matrix(SuperMatrix*, int, int, double*, int, Stype_t, Dtype_t, Mtype_t);
void sCreate_Dense_Matrix(SuperMatrix*, int, int, float*, int, Stype_t, Dtype_t, Mtype_t);

#ifdef __cplusplus
}
#endif

#endif
