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

#include "ConcreteCM.h"
#include <Toolbox/utility.h>

pod2 ConcreteCM::compute_compression_backbone(const double n_strain) {
    pod2 response;

    suanpan_assert([&] { if(n_strain > 0.) throw std::invalid_argument("argument is not acceptable"); });

    const auto normal_strain = std::max(datum::eps, n_strain / c_strain);

    const auto tmp_a = pow(normal_strain, c_n);
    const auto tmp_b = c_n == 1. ? 1. + (c_m - 1. + log(normal_strain)) * normal_strain : 1. + (c_m - c_n / (c_n - 1.)) * normal_strain + tmp_a / (c_n - 1.);
    response[0] = c_stress * c_m * normal_strain / tmp_b;
    response[1] = initial_stiffness(0) * (1. - tmp_a) / tmp_b / tmp_b;

    return response;
}

pod2 ConcreteCM::compute_tension_backbone(const double n_strain) {
    pod2 response;

    suanpan_assert([&] { if(n_strain < 0.) throw std::invalid_argument("argument is not acceptable"); });

    const auto normal_strain = std::max(datum::eps, n_strain / t_strain);

    const auto tmp_a = pow(normal_strain, t_n);
    const auto tmp_b = t_n == 1. ? 1. + (t_m - 1. + log(normal_strain)) * normal_strain : 1. + (t_m - t_n / (t_n - 1.)) * normal_strain + tmp_a / (t_n - 1.);
    response[0] = t_stress * t_m * normal_strain / tmp_b;
    response[1] = initial_stiffness(0) * (1. - tmp_a) / tmp_b / tmp_b;

    return response;
}

pod2 ConcreteCM::compute_compression_unload(const double n_strain) {
    const auto& unload_c_strain = trial_history(0);
    const auto& unload_c_stress = trial_history(1);
    const auto& residual_c_strain = trial_history(4);
    const auto& residual_c_stiffness = trial_history(5);
    const auto& unload_t_strain = trial_history(6);
    const auto& unload_t_stress = trial_history(7);
    const auto& reload_t_stiffness = trial_history(19);

    pod2 response;

    if(n_strain > unload_t_strain) {
        trial_load_status = Status::TBACKBONE;
        response = compute_tension_backbone(n_strain);
    }
    else if(n_strain > residual_c_strain) {
        trial_load_status = Status::CUNLOAD;
        response = compute_transition(n_strain, residual_c_strain, 0., residual_c_stiffness, unload_t_strain, unload_t_stress, reload_t_stiffness);
    }
    else {
        trial_load_status = Status::CUNLOAD;
        response = compute_transition(n_strain, unload_c_strain, unload_c_stress, initial_stiffness(0), residual_c_strain, 0., residual_c_stiffness);
    }

    return response;
}

pod2 ConcreteCM::compute_tension_unload(const double n_strain) {
    const auto& unload_c_strain = trial_history(0);
    const auto& unload_c_stress = trial_history(1);
    const auto& unload_t_strain = trial_history(6);
    const auto& unload_t_stress = trial_history(7);
    const auto& residual_t_strain = trial_history(10);
    const auto& residual_t_stiffness = trial_history(11);
    const auto& reload_c_stiffness = trial_history(18);

    pod2 response;

    if(n_strain < unload_c_strain) {
        trial_load_status = Status::CBACKBONE;
        response = compute_compression_backbone(n_strain);
    }
    else if(n_strain < residual_t_strain) {
        trial_load_status = Status::TUNLOAD;
        response = compute_transition(n_strain, residual_t_strain, 0., residual_t_stiffness, unload_c_strain, unload_c_stress, reload_c_stiffness);
    }
    else {
        trial_load_status = Status::TUNLOAD;
        response = compute_transition(n_strain, unload_t_strain, unload_t_stress, initial_stiffness(0), residual_t_strain, 0., residual_t_stiffness);
    }

    return response;
}

pod2 ConcreteCM::compute_compression_reload(const double n_strain) {
    const auto& unload_c_strain = trial_history(0);
    const auto& reload_c_stiffness = trial_history(18);

    pod2 response;

    if(n_strain < unload_c_strain) {
        trial_load_status = Status::CBACKBONE;
        response = compute_compression_backbone(n_strain);
    }
    else {
        trial_load_status = Status::CRELOAD;
        response[1] = reload_c_stiffness;
        response[0] = current_stress(0) + (n_strain - current_strain(0)) * response[1];
    }

    return response;
}

pod2 ConcreteCM::compute_tension_reload(const double n_strain) {
    const auto& unload_t_strain = trial_history(6);
    const auto& reload_t_stiffness = trial_history(19);

    pod2 response;

    if(n_strain > unload_t_strain) {
        trial_load_status = Status::TBACKBONE;
        response = compute_tension_backbone(n_strain);
    }
    else {
        trial_load_status = Status::TRELOAD;
        response[1] = reload_t_stiffness;
        response[0] = current_stress(0) + (n_strain - current_strain(0)) * response[1];
    }

    return response;
}

pod2 ConcreteCM::compute_compression_subunload(const double n_strain) {
    const auto& reverse_c_strain = trial_history(2);
    const auto& reverse_c_stress = trial_history(3);
    const auto& residual_c_strain = trial_history(4);

    pod2 response;

    if(n_strain > residual_c_strain) response = compute_compression_unload(n_strain);
    else {
        trial_load_status = Status::CSUBUNLOAD;
        response[1] = reverse_c_stress / (reverse_c_strain - residual_c_strain);
        response[0] = current_stress(0) + (n_strain - current_strain(0)) * response[1];
    }

    return response;
}

pod2 ConcreteCM::compute_tension_subunload(const double n_strain) {
    const auto& reverse_t_strain = trial_history(8);
    const auto& reverse_t_stress = trial_history(9);
    const auto& residual_t_strain = trial_history(10);

    pod2 response;

    if(n_strain < residual_t_strain) response = compute_tension_unload(n_strain);
    else {
        trial_load_status = Status::TSUBUNLOAD;
        response[1] = reverse_t_stress / (reverse_t_strain - residual_t_strain);
        response[0] = current_stress(0) + (n_strain - current_strain(0)) * response[1];
    }

    return response;
}

pod2 ConcreteCM::compute_transition(const double EM, const double EA, const double SA, const double KA, const double EB, const double SB, const double KB) const {
    pod2 response;

    if(fabs(EM - EA) <= 1E-15) {
        response[0] = SA;
        response[1] = KA;
    }
    else if(fabs(EM - EB) <= 1E-15) {
        response[0] = SB;
        response[1] = KB;
    }
    else if(linear_trans) {
        response[1] = (SB - SA) / (EB - EA);
        response[0] = SA + response[1] * (EM - EA);
    }
    else {
        const auto i_strain = EM - EA;
        const auto d_strain = EB - EA;
        const auto secant = (SB - SA) / d_strain;
        const auto tmp_a = secant - KA;
        const auto ratio = (KB - secant) / tmp_a;
        suanpan_assert([&] { if(ratio <= 0.) throw std::invalid_argument("argument is not acceptable"); });
        const auto tmp_b = tmp_a * pow(i_strain / d_strain, ratio);

        response[0] = SA + i_strain * (KA + tmp_b);
        response[1] = KA + (ratio + 1.) * tmp_b;
    }

    return response;
}

void ConcreteCM::update_compression_unload(const double r_strain) {
    auto& unload_c_strain = trial_history(0);
    auto& unload_c_stress = trial_history(1);
    auto& residual_c_strain = trial_history(4);
    auto& residual_c_stiffness = trial_history(5);
    auto& residual_t_stiffness = trial_history(11);
    auto& reload_c_stiffness = trial_history(18);
    const auto& unload_t_strain = trial_history(6);
    const auto& unload_t_stress = trial_history(7);
    const auto& residual_t_strain = trial_history(10);

    // inputs shall be on the backbone
    auto response = compute_compression_backbone(r_strain);

    // update unload point
    unload_c_strain = r_strain;    // -
    unload_c_stress = response[0]; // -

    const auto normal_strain = std::max(datum::eps, unload_c_strain / c_strain);                                                                                                        // +
    const auto secant_stiffness = std::max((std::max(datum::eps, unload_c_stress / c_strain) + .57 * initial_stiffness(0)) / (normal_strain + .57), unload_c_stress / unload_c_strain); // +

    // update residual strain
    residual_c_strain = std::min(-datum::eps, unload_c_strain - unload_c_stress / secant_stiffness);              // -
    residual_c_stiffness = std::min(.8 * secant_stiffness, .1 * initial_stiffness(0) * exp(-2. * normal_strain)); // +

    reload_c_stiffness = secant_stiffness;

    if(unload_t_strain == 0.) update_tension_unload(1.2 * t_strain);

    residual_c_stiffness = std::min(residual_c_stiffness, .8 * unload_t_stress / (unload_t_strain - residual_c_strain));
    residual_t_stiffness = std::min(residual_t_stiffness, .8 * unload_c_stress / (unload_c_strain - residual_t_strain));

    update_connect();
}

void ConcreteCM::update_tension_unload(const double r_strain) {
    auto& residual_c_stiffness = trial_history(5);
    auto& unload_t_strain = trial_history(6);
    auto& unload_t_stress = trial_history(7);
    auto& residual_t_strain = trial_history(10);
    auto& residual_t_stiffness = trial_history(11);
    auto& reload_t_stiffness = trial_history(19);
    const auto& unload_c_strain = trial_history(0);
    const auto& unload_c_stress = trial_history(1);
    const auto& residual_c_strain = trial_history(4);

    // inputs shall be on the backbone
    auto response = compute_tension_backbone(r_strain);

    // update unload point
    unload_t_strain = r_strain;
    unload_t_stress = response[0];

    const auto normal_strain = std::max(datum::eps, unload_t_strain / t_strain);                                                                                                        // +
    const auto secant_stiffness = std::max((std::max(datum::eps, unload_t_stress / t_strain) + .67 * initial_stiffness(0)) / (normal_strain + .67), unload_t_stress / unload_t_strain); // +

    residual_t_strain = std::max(datum::eps, unload_t_strain - unload_t_stress / secant_stiffness);                // +
    residual_t_stiffness = std::min(.8 * secant_stiffness, initial_stiffness(0) / (1. + pow(normal_strain, 1.1))); // +

    reload_t_stiffness = secant_stiffness;

    if(unload_c_strain == 0.) update_compression_unload(1.2 * c_strain);

    residual_t_stiffness = std::min(residual_t_stiffness, .8 * unload_c_stress / (unload_c_strain - residual_t_strain));
    residual_c_stiffness = std::min(residual_c_stiffness, .8 * unload_t_stress / (unload_t_strain - residual_c_strain));

    update_connect();
}

void ConcreteCM::update_connect() {
    auto& connect_c_stress = trial_history(12);
    auto& connect_c_stiffness = trial_history(13);
    auto& connect_t_stress = trial_history(14);
    auto& connect_t_stiffness = trial_history(15);
    const auto& unload_c_strain = trial_history(0);
    const auto& unload_c_stress = trial_history(1);
    const auto& residual_c_strain = trial_history(4);
    const auto& residual_c_stiffness = trial_history(5);
    const auto& unload_t_strain = trial_history(6);
    const auto& unload_t_stress = trial_history(7);
    const auto& residual_t_strain = trial_history(10);
    const auto& residual_t_stiffness = trial_history(11);
    const auto& reload_c_stiffness = trial_history(18);
    const auto& reload_t_stiffness = trial_history(19);

    auto response = compute_transition(residual_t_strain, residual_c_strain, 0., residual_c_stiffness, unload_t_strain, unload_t_stress, reload_t_stiffness);
    connect_t_stress = response[0];
    connect_t_stiffness = response[1];

    response = compute_transition(residual_c_strain, residual_t_strain, 0., residual_t_stiffness, unload_c_strain, unload_c_stress, reload_c_stiffness);
    connect_c_stress = response[0];
    connect_c_stiffness = response[1];
}

ConcreteCM::ConcreteCM(const unsigned T, const double E, const double SC, const double ST, const double NCC, const double NTT, const double EC, const double ET, const bool LT, const double R)
    : ComplexHysteresis(T, E, R)
    , c_stress(-perturb(fabs(SC)))
    , c_strain(-perturb(fabs(EC)))
    , t_stress(perturb(fabs(ST)))
    , t_strain(perturb(fabs(ET)))
    , c_m(elastic_modulus * c_strain / c_stress)
    , c_n(std::max(perturb(1.), NCC))
    , t_m(elastic_modulus * t_strain / t_stress)
    , t_n(std::max(perturb(1.), NTT))
    , linear_trans(LT) {}

unique_ptr<Material> ConcreteCM::get_copy() { return make_unique<ConcreteCM>(*this); }

double ConcreteCM::get_parameter(const ParameterType P) const {
    if(ParameterType::ELASTICMODULUS == P) return elastic_modulus;
    if(ParameterType::PEAKSTRAIN == P) return c_strain;
    if(ParameterType::CRACKSTRAIN == P) return t_strain;
    return 0.;
}

void ConcreteCM::print() {
    suanpan_info("A concrete model based on Chang & Mander's concrete model.\n");
    ComplexHysteresis::print();
}
