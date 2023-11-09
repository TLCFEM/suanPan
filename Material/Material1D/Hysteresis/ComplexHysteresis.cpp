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

#include "ComplexHysteresis.h"

ComplexHysteresis::ComplexHysteresis(const unsigned T, const double E, const double R)
    : Material1D(T, R)
    , elastic_modulus(fabs(E)) {}

int ComplexHysteresis::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(20);

    return SUANPAN_SUCCESS;
}

int ComplexHysteresis::update_trial_status(const vec& n_strain) {
    incre_strain = (trial_strain = n_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_load_status = current_load_status;
    trial_history = current_history;

    const auto& unload_c_strain = trial_history(0);
    const auto& unload_c_stress = trial_history(1);
    const auto& residual_c_strain = trial_history(4);
    const auto& unload_t_strain = trial_history(6);
    const auto& unload_t_stress = trial_history(7);
    const auto& residual_t_strain = trial_history(10);
    const auto& connect_c_stress = trial_history(12);
    const auto& connect_t_stress = trial_history(14);

    // const auto& residual_c_stiffness = trial_history(5);
    // const auto& residual_t_stiffness = trial_history(11);
    // const auto& connect_c_stiffness = trial_history(13);
    // const auto& connect_t_stiffness = trial_history(15);

    auto& reverse_c_strain = trial_history(2);
    auto& reverse_c_stress = trial_history(3);
    auto& reverse_t_strain = trial_history(8);
    auto& reverse_t_stress = trial_history(9);
    auto& inter_strain = trial_history(16);
    auto& inter_stress = trial_history(17);
    auto& reload_c_stiffness = trial_history(18);
    auto& reload_t_stiffness = trial_history(19);

    const auto tension_sign = incre_strain(0) > 0.;

    podarray<double> response(2);

    if(Status::CBACKBONE == trial_load_status) {
        if(!tension_sign) response = compute_compression_backbone(trial_strain(0));
        else {
            update_compression_unload(current_strain(0));
            response = compute_compression_unload(trial_strain(0));
        }
    }
    else if(Status::TBACKBONE == trial_load_status) {
        if(tension_sign) response = compute_tension_backbone(trial_strain(0));
        else {
            update_tension_unload(current_strain(0));
            response = compute_tension_unload(trial_strain(0));
        }
    }
    else if(Status::CUNLOAD == trial_load_status) {
        if(tension_sign) response = compute_compression_unload(trial_strain(0));
        else if(current_strain(0) > residual_t_strain) {
            reverse_t_strain = current_strain(0);
            reverse_t_stress = current_stress(0);
            response = compute_tension_subunload(trial_strain(0));
        }
        else if(current_strain(0) > residual_c_strain)
            if(trial_strain(0) < residual_c_strain) response = compute_tension_unload(trial_strain(0));
            else {
                trial_load_status = Status::CTRANS;
                inter_strain = current_strain(0);
                inter_stress = current_stress(0);
                response(1) = (connect_c_stress - inter_stress) / (residual_c_strain - inter_strain);
                response(0) = current_stress(0) + incre_strain(0) * response(1);
            }
        else {
            reload_c_stiffness = (unload_c_stress - current_stress(0)) / (unload_c_strain - current_strain(0));
            response = compute_compression_reload(trial_strain(0));
        }
    }
    else if(Status::TUNLOAD == trial_load_status) {
        if(!tension_sign) response = compute_tension_unload(trial_strain(0));
        else if(current_strain(0) < residual_c_strain) {
            reverse_c_strain = current_strain(0);
            reverse_c_stress = current_stress(0);
            response = compute_compression_subunload(trial_strain(0));
        }
        else if(current_strain(0) < residual_t_strain)
            if(trial_strain(0) > residual_t_strain) response = compute_compression_unload(trial_strain(0));
            else {
                trial_load_status = Status::TTRANS;
                inter_strain = current_strain(0);
                inter_stress = current_stress(0);
                response(1) = (connect_t_stress - inter_stress) / (residual_t_strain - inter_strain);
                response(0) = current_stress(0) + incre_strain(0) * response(1);
            }
        else {
            reload_t_stiffness = (unload_t_stress - current_stress(0)) / (unload_t_strain - current_strain(0));
            response = compute_tension_reload(trial_strain(0));
        }
    }
    else if(Status::CSUBUNLOAD == trial_load_status) {
        if(tension_sign) response = compute_compression_subunload(trial_strain(0));
        else {
            reload_c_stiffness = (unload_c_stress - current_stress(0)) / (unload_c_strain - current_strain(0));
            response = compute_compression_reload(trial_strain(0));
        }
    }
    else if(Status::TSUBUNLOAD == trial_load_status) {
        if(!tension_sign) response = compute_tension_subunload(trial_strain(0));
        else {
            reload_t_stiffness = (unload_t_stress - current_stress(0)) / (unload_t_strain - current_strain(0));
            response = compute_tension_reload(trial_strain(0));
        }
    }
    else if(Status::CRELOAD == trial_load_status) {
        if(!tension_sign) response = compute_compression_reload(trial_strain(0));
        else {
            reverse_c_strain = current_strain(0);
            reverse_c_stress = current_stress(0);
            response = compute_compression_subunload(trial_strain(0));
        }
    }
    else if(Status::TRELOAD == trial_load_status) {
        if(tension_sign) response = compute_tension_reload(trial_strain(0));
        else {
            reverse_t_strain = current_strain(0);
            reverse_t_stress = current_stress(0);
            response = compute_tension_subunload(trial_strain(0));
        }
    }
    else if(Status::CTRANS == trial_load_status) {
        if(tension_sign)
            if(trial_strain(0) > residual_t_strain) response = compute_compression_unload(trial_strain(0));
            else {
                trial_load_status = Status::TTRANS;
                inter_strain = current_strain(0);
                inter_stress = current_stress(0);
                response(1) = (connect_t_stress - inter_stress) / (residual_t_strain - inter_strain);
                response(0) = current_stress(0) + incre_strain(0) * response(1);
            }
        else if(trial_strain(0) < residual_c_strain) response = compute_tension_unload(trial_strain(0));
        else {
            response(1) = (connect_c_stress - inter_stress) / (residual_c_strain - inter_strain);
            response(0) = current_stress(0) + incre_strain(0) * response(1);
        }
    }
    else if(Status::TTRANS == trial_load_status) {
        if(tension_sign)
            if(trial_strain(0) > residual_t_strain) response = compute_compression_unload(trial_strain(0));
            else {
                response(1) = (connect_t_stress - inter_stress) / (residual_t_strain - inter_strain);
                response(0) = current_stress(0) + incre_strain(0) * response(1);
            }
        else if(trial_strain(0) < residual_c_strain) response = compute_tension_unload(trial_strain(0));
        else {
            trial_load_status = Status::CTRANS;
            inter_strain = current_strain(0);
            inter_stress = current_stress(0);
            response(1) = (connect_c_stress - inter_stress) / (residual_c_strain - inter_strain);
            response(0) = current_stress(0) + incre_strain(0) * response(1);
        }
    }
    else {
        if(tension_sign) {
            trial_load_status = Status::TBACKBONE;
            response = compute_tension_backbone(trial_strain(0));
        }
        else {
            trial_load_status = Status::CBACKBONE;
            response = compute_compression_backbone(trial_strain(0));
        }
    }

    trial_stress = response(0);
    trial_stiffness = response(1);

    return SUANPAN_SUCCESS;
}

int ComplexHysteresis::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history.zeros();
    current_load_status = Status::NONE;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int ComplexHysteresis::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_load_status = trial_load_status;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int ComplexHysteresis::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_load_status = current_load_status;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}
