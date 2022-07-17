#ifndef TESTSOLVER_H
#define TESTSOLVER_H

#include <suanPan.h>

struct System {
    mat A;

    explicit System(mat&& AA)
        : A(std::forward<mat>(AA)) {}

    [[nodiscard]] vec evaluate(const vec& x) const { return A * x; }
};

struct NonlinearSystem {
    [[nodiscard]] vec evaluate(const vec& x) const { return square(x); }
};

template<sp_d data_t> class DiagonalPreconditioner {
    vec diag_reciprocal;
public:
    explicit DiagonalPreconditioner(const Mat<data_t>& in_mat)
        : diag_reciprocal(1. / in_mat.diag()) {}

    [[nodiscard]] Col<data_t> apply(const Col<data_t>& in) const { return diag_reciprocal % in; }
};

#endif // TESTSOLVER_H
