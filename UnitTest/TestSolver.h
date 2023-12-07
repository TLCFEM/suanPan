#ifndef TESTSOLVER_H
#define TESTSOLVER_H

#include <suanPan.h>

struct System {
    mat A;

    explicit System(mat&& AA)
        : A(std::move(AA)) {}

    [[nodiscard]] vec evaluate(const vec& x) const { return A * x; }
};

struct NonlinearSystem {
    [[nodiscard]] vec evaluate(const vec& x) const { return square(x); }
};

class Quadratic {
public:
    [[nodiscard]] vec evaluate_residual(const vec& x) const { return square(x) - 1.; }

    [[nodiscard]] mat evaluate_jacobian(const vec& x) const { return diagmat(2. * x); }
};

#endif // TESTSOLVER_H
