#include "sat.h"

namespace suanpan {
    /**
     * @brief Test intersection between two convex polygons using the Separating Axis Theorem (SAT).
     * For each edge of both polygons, constructs an axis and compares the projections
     * of the polygons onto the axis perpendicular to that edge. If a separating axis
     * is found, the polygons do not intersect; otherwise, they intersect or touch.
     *
     * @param a A 2xNa matrix of vertices for polygon A, columns are ordered 2D points.
     * @param b A 2xNb matrix of vertices for polygon B, columns are ordered 2D points.
     * @return true if polygons intersect or touch; false if a separating axis is found.
     *
     * @pre Polygons should be convex and vertices ordered (counterclockwise).
     */
    bool sat(const mat& a, const mat& b) {
        const auto project = [](const vec& axis, const mat& nodes) -> rowvec { return rowvec{axis(1), -axis(0)} * nodes; };

        const auto overlap = [&](const mat& target) {
            for(auto I = 0llu, J = 1llu; I < target.n_cols; ++I, ++J) {
                const vec axis = target.col(J % target.n_cols) - target.col(I);
                const auto proj_a = project(axis, a), proj_b = project(axis, b);
                if(std::max(proj_a.min(), proj_b.min()) > std::min(proj_a.max(), proj_b.max())) return false;
            }
            return true;
        };

        return overlap(a) && overlap(b);
    }

    auto intersect(const vec& x1, const vec& x2, const vec& x3, const vec& x4) {
        const vec d12 = x2 - x1;
        const vec d34 = x4 - x3;
        const vec d13 = x3 - x1;

        const auto cross = [](const vec& a, const vec& b) { return a(0) * b(1) - a(1) * b(0); };

        static const mat rotation{{0., -1.}, {1., 0.}};

        const auto numerator = cross(d13, d34);
        const auto denominator = cross(d12, d34);
        const auto p = numerator / denominator;

        const rowvec pnp1 = rotation * d34;
        const rowvec pnp4 = rotation * d13;
        const rowvec pnp3 = pnp1 + pnp4;

        const rowvec pdp1 = rotation * d34;
        const rowvec pdp2 = -pdp1;
        const rowvec pdp4 = rotation * d12;
        const rowvec pdp3 = -pdp4;

        const rowvec ppp1 = (pnp1 - p * pdp1) / denominator;
        const rowvec ppp2 = -p / denominator * pdp2;
        const rowvec ppp3 = (pnp3 - p * pdp3) / denominator;
        const rowvec ppp4 = (pnp4 - p * pdp4) / denominator;

        const vec intersection = x1 + p * d12;

        const mat pip1 = eye(2, 2) * (1. - p) + d12 * ppp1;
        const mat pip2 = eye(2, 2) * p + d12 * ppp2;
        const mat pip3 = d12 * ppp3;
        const mat pip4 = d12 * ppp4;

        return intersection;
    }
} // namespace suanpan
