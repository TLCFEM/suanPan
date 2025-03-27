#ifndef TESTSOLVER_H
#define TESTSOLVER_H

#include <suanPan.h>

class Quadratic {
public:
    [[nodiscard]] vec evaluate_residual(const vec& x) const { return square(x) - 1.; }

    [[nodiscard]] mat evaluate_jacobian(const vec& x) const { return diagmat(2. * x); }
};

#endif // TESTSOLVER_H
