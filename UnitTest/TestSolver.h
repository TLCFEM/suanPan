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
    const vec diag_reciprocal;
public:
    explicit Jacobi(const Container& in_mat)
        : diag_reciprocal(1. / vec(in_mat.diag()).replace(0., 1.)) {}

    [[nodiscard]] vec apply(const vec& in) const { return diag_reciprocal % in; }
};

class Quadratic {
public:
    [[nodiscard]] vec evaluate_residual(const vec& x) const { return square(x) - 1.; }

    [[nodiscard]] mat evaluate_jacobian(const vec& x) const { return diagmat(2. * x); }
};

#endif // TESTSOLVER_H
