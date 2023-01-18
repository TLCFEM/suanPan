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

#include "SimpleHysteresis.h"

podarray<double> SimpleHysteresis::compute_compression_inner(const double t_strain) const {
    podarray<double> response(2);

    auto& reverse_c_strain = trial_history(2);  // unloading point strain compression side
    auto& reverse_c_stress = trial_history(3);  // unloading point stress compression side
    auto& reverse_t_strain = trial_history(4);  // unloading point strain tension side
    auto& reverse_t_stress = trial_history(5);  // unloading point stress tension side
    auto& residual_c_strain = trial_history(6); // residual strain in compression unloading path
    auto& residual_t_strain = trial_history(7); // residual strain in compression unloading path

    if(t_strain < residual_c_strain) {
        response(1) = reverse_c_stress / (reverse_c_strain - residual_c_strain);
        response(0) = response(1) * (t_strain - residual_c_strain);
    }
    else if(const auto middle_strain = (1. - middle_point) * residual_t_strain + middle_point * reverse_t_strain; t_strain < middle_strain) {
        response(1) = middle_point * reverse_t_stress / (middle_strain - residual_c_strain);
        response(0) = response(1) * (t_strain - residual_c_strain);
    }
    else {
        response(1) = reverse_t_stress / (reverse_t_strain - residual_t_strain);
        response(0) = response(1) * (t_strain - residual_t_strain);
    }

    return response;
}

podarray<double> SimpleHysteresis::compute_tension_inner(const double t_strain) const {
    podarray<double> response(2);

    auto& reverse_c_strain = trial_history(2);  // unloading point strain compression side
    auto& reverse_c_stress = trial_history(3);  // unloading point stress compression side
    auto& reverse_t_strain = trial_history(4);  // unloading point strain tension side
    auto& reverse_t_stress = trial_history(5);  // unloading point stress tension side
    auto& residual_c_strain = trial_history(6); // residual strain in compression unloading path
    auto& residual_t_strain = trial_history(7); // residual strain in compression unloading path

    if(t_strain > residual_t_strain) {
        response(1) = reverse_t_stress / (reverse_t_strain - residual_t_strain);
        response(0) = response(1) * (t_strain - residual_t_strain);
    }
    else if(const auto middle_strain = (1. - middle_point) * residual_c_strain + middle_point * reverse_c_strain; t_strain > middle_strain) {
        response(1) = middle_point * reverse_c_stress / (middle_strain - residual_t_strain);
        response(0) = response(1) * (t_strain - residual_t_strain);
    }
    else {
        response(1) = reverse_c_stress / (reverse_c_strain - residual_c_strain);
        response(0) = response(1) * (t_strain - residual_c_strain);
    }

    return response;
}

SimpleHysteresis::SimpleHysteresis(const unsigned T, const double M, const double R)
    : DataSimpleHysteresis{std::min(std::max(datum::eps, M), 1.)}
    , Material1D(T, R) {}

int SimpleHysteresis::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    if(current_history.is_empty() || !any(current_history)) {
        current_history.zeros(8);
        auto point = compute_compression_initial_reverse();
        current_history(2) = point(0);
        current_history(3) = point(1);
        point = compute_tension_initial_reverse();
        current_history(4) = point(0);
        current_history(5) = point(1);
        initial_history = current_history;
    }
    else if(current_history.size() != 8) {
        current_history.resize(8);
        initial_history = current_history;
    }

    trial_history = current_history;
    auto& max_c_strain = trial_history(0);      // maximum compression strain recorded
    auto& max_t_strain = trial_history(1);      // maximum tension strain recorded
    auto& reverse_c_strain = trial_history(2);  // unloading point strain compression side
    auto& reverse_c_stress = trial_history(3);  // unloading point stress compression side
    auto& reverse_t_strain = trial_history(4);  // unloading point strain tension side
    auto& reverse_t_stress = trial_history(5);  // unloading point stress tension side
    auto& residual_c_strain = trial_history(6); // residual strain in compression unloading path
    auto& residual_t_strain = trial_history(7); // residual strain in compression unloading path

    podarray<double> response;

    if(Status::NONE == current_flag) {
        if(incre_strain(0) > 0.) {
            trial_flag = Status::TBACKBONE;
            response = compute_tension_backbone(max_t_strain = trial_strain(0));
        }
        else {
            trial_flag = Status::CBACKBONE;
            response = compute_compression_backbone(max_c_strain = trial_strain(0));
        }
    }
    else if(Status::CBACKBONE == current_flag) {
        if(incre_strain(0) > 0.) {
            residual_c_strain = compute_compression_residual(reverse_c_strain = current_strain(0), reverse_c_stress = current_stress(0));
            if(trial_strain(0) > reverse_t_strain) {
                trial_flag = Status::TBACKBONE;
                response = compute_tension_backbone(max_t_strain = trial_strain(0));
            }
            else {
                trial_flag = Status::CINNER;
                response = compute_compression_inner(trial_strain(0));
            }
        }
        else {
            trial_flag = Status::CBACKBONE;
            response = compute_compression_backbone(max_c_strain = trial_strain(0));
        }
    }
    else if(Status::TBACKBONE == current_flag) {
        if(incre_strain(0) < 0.) {
            residual_t_strain = compute_tension_residual(reverse_t_strain = current_strain(0), reverse_t_stress = current_stress(0));
            if(trial_strain(0) < reverse_c_strain) {
                trial_flag = Status::CBACKBONE;
                response = compute_compression_backbone(max_c_strain = trial_strain(0));
            }
            else {
                trial_flag = Status::TINNER;
                response = compute_tension_inner(trial_strain(0));
            }
        }
        else {
            trial_flag = Status::TBACKBONE;
            response = compute_tension_backbone(max_t_strain = trial_strain(0));
        }
    }
    else if(Status::CINNER == current_flag) {
        if(trial_strain(0) > reverse_t_strain) {
            trial_flag = Status::TBACKBONE;
            response = compute_tension_backbone(max_t_strain = trial_strain(0));
        }
        else if(trial_strain(0) < reverse_c_strain) {
            trial_flag = Status::CBACKBONE;
            response = compute_compression_backbone(max_c_strain = trial_strain(0));
        }
        else {
            trial_flag = Status::CINNER;
            response = compute_compression_inner(trial_strain(0));
        }
    }
    else if(Status::TINNER == current_flag) {
        if(trial_strain(0) > reverse_t_strain) {
            trial_flag = Status::TBACKBONE;
            response = compute_tension_backbone(max_t_strain = trial_strain(0));
        }
        else if(trial_strain(0) < reverse_c_strain) {
            trial_flag = Status::CBACKBONE;
            response = compute_compression_backbone(max_c_strain = trial_strain(0));
        }
        else {
            trial_flag = Status::TINNER;
            response = compute_tension_inner(trial_strain(0));
        }
    }

    trial_stress = response(0);
    trial_stiffness = response(1);

    suanpan_assert([&] { if(!trial_stress.is_finite() || !trial_stiffness.is_finite()) throw invalid_argument("infinite number detected"); });

    return SUANPAN_SUCCESS;
}

int SimpleHysteresis::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_flag = Status::NONE;
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int SimpleHysteresis::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_flag = trial_flag;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int SimpleHysteresis::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_flag = current_flag;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}
