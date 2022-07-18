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

template<typename Container> class Jacobi {
    vec diag_reciprocal;
public:
    explicit Jacobi(Container& in_mat)
        : diag_reciprocal(1. / in_mat.diag()) {}

    [[nodiscard]] vec apply(const vec& in) const { return diag_reciprocal % in; }
};

#endif // TESTSOLVER_H
