/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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

#include "Frequency.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Arnoldi.h>
#include <Solver/Integrator/Integrator.h>

Frequency::Frequency(const unsigned T, const unsigned N, const char TP)
    : Step(T, 0.)
    , eigen_number(N)
    , eigen_type(TP) {}

int Frequency::initialize() {
    configure_storage_scheme();

    factory->set_analysis_type(AnalysisType::EIGEN);

    const auto t_domain = database.lock();

    // integrator
    modifier = make_shared<Integrator>();
    modifier->set_domain(t_domain);

    // solver
    if(nullptr == solver) solver = make_shared<Arnoldi>(0, eigen_number, eigen_type);
    solver->set_integrator(modifier);

    if(SUANPAN_SUCCESS != modifier->initialize()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != solver->initialize()) return SUANPAN_FAIL;

    return t_domain->restart();
}

int Frequency::analyze() {
    auto& G = get_integrator();

    if(SUANPAN_SUCCESS != solver->analyze()) {
        suanpan_warning("Fail to decompose the system, try to increase the number of eigenvalues.\n");
        return SUANPAN_SUCCESS;
    }

    G->record();

    return SUANPAN_SUCCESS;
}

void Frequency::set_eigen_number(const unsigned N) const { access::rw(eigen_number) = N; }

unsigned Frequency::get_eigen_number() const { return eigen_number; }
