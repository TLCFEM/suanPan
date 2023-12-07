/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "TabularSpline.h"
#include <Domain/DomainBase.h>
#include <Domain/MetaMat/BandMat.hpp>

void TabularSpline::initialize(const shared_ptr<DomainBase>& D) {
    Tabular::initialize(D);

    if(!is_active()) return;

    dt = diff(time);

    if(!all(dt)) {
        D->disable_amplitude(get_tag());
        suanpan_warning("Repeated data points detected.\n");
        return;
    }

    dy = diff(magnitude);

    const auto np = time.n_elem;
    const auto n = np - 1llu;
    const auto nm = n - 1llu;

    BandMat<double> system(np, 1, 1);

    suanpan::for_each(0llu, np, [&](const uword I) { system.at(I, I) = 2.; });

    system.at(0llu, 1llu) = 1.;
    system.at(n, nm) = 1.;

    vec b(np, fill::none);
    b(0) = dy(0) / dt(0) / dt(0);
    b(n) = -dy(nm) / dt(nm) / dt(nm);

    suanpan::for_each(1llu, n, [&](const uword I) {
        const auto J = I - 1llu, K = I + 1llu;
        const auto denom = time(K) - time(J);
        const auto mu = dt(J) / denom;
        system.at(I, J) = mu;
        system.at(I, I) = 2.;
        system.at(I, K) = 1. - mu;
        b(I) = (dy(I) / dt(I) - dy(J) / dt(J)) / denom;
    });

    system.solve(m, b);
}

double TabularSpline::get_amplitude(const double T) {
    const auto step_time = T - start_time;

    if(step_time <= time.front()) return magnitude.front();

    if(step_time >= time.back()) return magnitude.back();

    const auto I = std::distance(time.cbegin(), std::lower_bound(time.cbegin(), time.cend(), step_time));
    const auto J = I - 1;

    double y = (m(J) * pow(time(I) - step_time, 3.) + m(I) * pow(step_time - time(J), 3.)) / dt(J);
    y += (magnitude(J) / dt(J) - m(J) * dt(J)) * (time(I) - step_time);
    y += (magnitude(I) / dt(J) - m(I) * dt(J)) * (step_time - time(J));

    return y;
}

void TabularSpline::print() {
    suanpan_info("TabularSpline.\n");
}
