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
/**
 * @class EnergyEvolution
 * @brief A EnergyEvolution class.
 *
 * The EnergyEvolution class.
 *
 * @author tlc
 * @date 13/09/2020
 * @version 0.1.0
 * @file EnergyEvolution.h
 * @addtogroup Criterion
 * @{
 */

#ifndef ENERGYEVOLUTION_H
#define ENERGYEVOLUTION_H

#include <Constraint/Criterion/Criterion.h>
#include <Toolbox/container.h>

class Element;

class EnergyEvolution : public Criterion {
    const unsigned iteration, reactive_ratio;

    const unsigned incre_level, final_level;

    const double weight, propagation_weight;

    const double tolerance;

    suanpan::graph<unsigned> map;
    std::vector<unsigned> index;

    unsigned current_level = incre_level;

    vec energy;

    double total_energy = 0.;

    unsigned balanced_iteration = 0;

protected:
    std::function<double(const Element*)> get_energy;

public:
    EnergyEvolution(unsigned,      // tag
                    unsigned,      // step tag
                    unsigned,      // incre level
                    unsigned,      // final level
                    double = 1.,   // centre weight
                    unsigned = 2,  // propagation iteration
                    unsigned = 10, // reactivation ratio
                    double = .5,   // propagation weight
                    double = 1E-5  // tolerance
        );

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
