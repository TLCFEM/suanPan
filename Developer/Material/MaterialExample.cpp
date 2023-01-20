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

#include "MaterialExample.h"
#include <Toolbox/utility.h>

SUANPAN_EXPORT void new_materialexample(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("A valid elastic modulus is required.\n");
        return;
    }

    double yield_stress;
    if(!get_input(command, yield_stress)) {
        suanpan_error("A valid yield stress is required.\n");
        return;
    }

    auto hardening_ratio = 0.;
    if(command.eof())
        suanpan_debug("Zero hardening ratio assumed.\n");
    else if(!get_input(command, hardening_ratio)) {
        suanpan_error("A valid hardening ratio is required.\n");
        return;
    }

    auto beta = 0.;
    if(command.eof())
        suanpan_debug("Isotropic hardening assumed.\n");
    else if(!get_input(command, beta)) {
        suanpan_error("A valid beta is required.\n");
        return;
    }

    auto density = 0.;
    if(command.eof())
        suanpan_debug("Zero density assumed.\n");
    else if(!get_input(command, density)) {
        suanpan_error("A valid density is required.\n");
        return;
    }

    return_obj = make_unique<MaterialExample>(tag, elastic_modulus, yield_stress, hardening_ratio, beta, density);
}

/**
 * \brief The constructor of material model depends on the specific model to be
 * implemented. But for the base `Material` class, it takes only three input arguments:
 *
 * - Unique Material Tag
 * - Material Type
 * - Density
 *
 * Here we are implementing a uniaxial bilinear hardening model, hence `MaterialType::D1`
 * is passed to the base. The material type will be used to validate with associated
 * elements/sections to ensure the consistency of the sizes of data passed between objects.
 *
 * \param T Unique Material Tag
 * \param E Elastic Modulus
 * \param Y Yield Stress
 * \param H Hardening Ratio
 * \param B Beta
 * \param R Density
 */
MaterialExample::MaterialExample(const unsigned T, const double E, const double Y, const double H, const double B, const double R)
    : MaterialExampleData{E, Y, H, B, E * H / (1. - H)}
    , Material(T, MaterialType::D1, R) {}

/**
 * \brief Unless the target material model involves other material models to compute
 * responses, in general, it is not necessary to get additional information from
 * other parts of the model.
 *
 * In general cases, history variables shall be initialised and initial stiffness
 * (and initial damping if appropriate) shall be set to proper value.
 *
 * To enable initial values for history variables, build-in `trial_history` and
 * `current_history` shall be used, developers can initialise them via method
 * `initialise_history()` so that initial values set by `initial` command will not be
 * overwritten.
 *
 * In this example, instead of using build-in history variables, we manage history
 * variables, namely back stress and plastic strain, by ourselves.
 *
 */
int MaterialExample::initialize(const shared_ptr<DomainBase>&) {
    current_back_stress = trial_back_stress = 0.;
    current_plastic_strain = trial_plastic_strain = 0.;

    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

/**
 * \brief The `get_copy()` method should always be implemented with `make_unique`.
 * In case the model defines other memory management, developers may need to further
 * provide a copy ctor to make `make_unique` work.
 *
 * **!!!NEVER DO A MANUAL COPY OF DATA IN THIS METHOD!!!**
 *
 * \return a copy of material model
 */
unique_ptr<Material> MaterialExample::get_copy() { return make_unique<MaterialExample>(*this); }

/**
 * \brief There are two states we are managing at any time point.
 * The current state is the converged state from the last time substep. Since it is
 * converged, all data shall be valid and accurate.
 * The trial state stores the response computed based on converged state and new trial
 * strain. It may be discarded, committed or overwritten with new trial values.
 *
 * **WE ALWAYS COMPUTE TRIAL STATE BASED ON CURRENT STATE AND NEW TRIAL STRAIN**
 * **NEVER COMPUTE RESPONSE BASED ON ANY INFORMATION FROM UNCONVERGED STATE**
 *
 * Developers who are not familiar with classic plasticity theory may consult textbooks
 * for details.
 *
 * \param t_strain trial strain
 * \return error flag
 */
int MaterialExample::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    incre_strain = trial_strain - current_strain;

    if(fabs(incre_strain(0)) <= tolerance) return 0;

    trial_back_stress = current_back_stress;
    trial_plastic_strain = current_plastic_strain;
    trial_stiffness = initial_stiffness;

    trial_stress = current_stress + elastic_modulus * incre_strain;

    const auto shifted_stress = trial_stress(0) - current_back_stress;

    const auto yield_func = fabs(shifted_stress) - yield_stress - (1. - beta) * plastic_modulus * current_plastic_strain;

    if(yield_func > 0.) {
        const auto incre_plastic_strain = yield_func / (elastic_modulus + plastic_modulus);
        trial_stress -= suanpan::sign(shifted_stress) * elastic_modulus * incre_plastic_strain;
        trial_stiffness *= hardening_ratio;
        trial_back_stress += suanpan::sign(shifted_stress) * beta * plastic_modulus * incre_plastic_strain;
        trial_plastic_strain += incre_plastic_strain;
    }

    return 0;
}

/**
 * \brief Operations are required to achieve the following objective.
 *
 * current state -> 0
 * trial state -> 0
 *
 * \return error flag
 */
int MaterialExample::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_stiffness = initial_stiffness;
    current_back_stress = 0.;
    current_plastic_strain = 0.;
    return reset_status();
}

/**
 * \brief Operations are required to achieve the following objective.
 *
 * current state <- trial state
 *
 * \return error flag
 */
int MaterialExample::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_back_stress = trial_back_stress;
    current_plastic_strain = trial_plastic_strain;
    return 0;
}

/**
 * \brief Operations are required to achieve the following objective.
 *
 * current state -> trial state
 *
 * \return error flag
 */
int MaterialExample::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_back_stress = current_back_stress;
    trial_plastic_strain = current_plastic_strain;
    return 0;
}

void MaterialExample::print() {
    suanpan_info("A material example based on uniaxial J2 bilinear mixed hardening model.\n");
    suanpan_info("Strain: {:.5E}\tStress: {:.5E}\n", current_strain.at(0), current_stress.at(0));
}
