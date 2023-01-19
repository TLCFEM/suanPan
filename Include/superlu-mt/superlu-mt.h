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
    NO,
    YES
} yes_no_t;

typedef enum {
    DOFACT,
    SamePattern,
    SamePattern_SameRowPerm,
    FACTORED
} fact_t;

typedef enum {
    NOROWPERM,
    LargeDiag_MC64,
    LargeDiag_HWPM,
    MY_PERMR
} rowperm_t;

typedef enum {
    NATURAL,
    MMD_ATA,
    MMD_AT_PLUS_A,
    COLAMD,
    METIS_AT_PLUS_A,
    PARMETIS,
    ZOLTAN,
    MY_PERMC
} colperm_t;

typedef enum {
    NOTRANS,
    TRANS,
    CONJ
} trans_t;

typedef enum {
    NOREFINE,
    SLU_SINGLE = 1,
    SLU_DOUBLE,
    SLU_EXTRA
} IterRefine_t;

typedef enum {
    SYSTEM,
    USER
} LU_space_t;

typedef enum {
    ONE_NORM,
    TWO_NORM,
    INF_NORM
} norm_t;

typedef enum {
    SILU,
    SMILU_1,
    SMILU_2,
    SMILU_3
} milu_t;

typedef struct {
    fact_t Fact;
    yes_no_t Equil;
    colperm_t ColPerm;
    trans_t Trans;
    IterRefine_t IterRefine;
    double DiagPivotThresh;
    yes_no_t SymmetricMode;
    yes_no_t PivotGrowth;
    yes_no_t ConditionNumber;
    rowperm_t RowPerm;
    int ILU_DropRule;
    double ILU_DropTol;
    double ILU_FillFactor;
    norm_t ILU_Norm;
    double ILU_FillTol;
    milu_t ILU_MILU;
    double ILU_MILU_Dim;
    yes_no_t ParSymbFact;
    yes_no_t ReplaceTinyPivot;
    yes_no_t SolveInitialized;
    yes_no_t RefineInitialized;
    yes_no_t PrintStat;
    int nnzL, nnzU;
    int num_lookaheads;
    yes_no_t lookahead_etree;
    yes_no_t SymPattern;
} superlu_options_t;

typedef float flops_t;

typedef struct {
    int* panel_histo;
    double* utime;
    flops_t* ops;
    int TinyPivots;
    int RefineSteps;
    int expansions;
} SuperLUStat_t;

typedef struct {
    float for_lu;
    float total_needed;
} mem_usage_t;

typedef struct e_node {
    int size;
    void* mem;
} ExpHeader;

typedef struct {
    int size;
    int used;
    int top1;
    int top2;
    void* array;
} LU_stack_t;

typedef struct {
    int* xsup;
    int* supno;
    int* lsub;
    int* xlsub;
    void* lusup;
    int* xlusup;
    void* ucol;
    int* usub;
    int* xusub;
    int nzlmax;
    int nzumax;
    int nzlumax;
    int n;
    LU_space_t MemModel;
    int num_expansions;
    ExpHeader* expanders;
    LU_stack_t stack;
} GlobalLU_t;

#ifdef __cplusplus
extern "C" {
#endif

void* superlu_malloc(size_t);
void superlu_free(void*);
int sp_ienv(int);
void StatAlloc(int, int, int, int, Gstat_t*);
void StatInit(int, int, Gstat_t*);
void StatFree(Gstat_t*);
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
void dgsisx(superlu_options_t*, SuperMatrix*, int*, int*, int*, char*, double*, double*, SuperMatrix*, SuperMatrix*, void*, int, SuperMatrix*, SuperMatrix*, double*, double*, GlobalLU_t*, mem_usage_t*, SuperLUStat_t*, int*);
void sgsisx(superlu_options_t*, SuperMatrix*, int*, int*, int*, char*, float*, float*, SuperMatrix*, SuperMatrix*, void*, int, SuperMatrix*, SuperMatrix*, float*, float*, GlobalLU_t*, mem_usage_t*, SuperLUStat_t*, int*);
void ilu_set_default_options(superlu_options_t*);

#ifdef __cplusplus
}
#endif

#endif
