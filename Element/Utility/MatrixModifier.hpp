/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @fn MatrixModifier
 * @brief The MatrixModifier class.
 *
 * @author tlc
 * @date 27/10/2017
 * @version 0.1.0
 * @file MatrixModifier.hpp
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef MATRIXMODIFIER_H
#define MATRIXMODIFIER_H

#include <Element/Element.h>
#include <suanPan.h>

namespace suanpan {
    namespace mass {
        struct lumped_simple {
            template<typename T> static void apply(Mat<T>&);
        };

        struct lumped_scale {
            template<typename T> static void apply(Mat<T>&, unsigned);
        };
    } // namespace mass
    namespace damping {
        struct rayleigh {
            template<typename T> static void apply(const shared_ptr<Element>&, T, T, T, T);
        };

        struct elemental {
            template<typename T> static void apply(const shared_ptr<Element>&, T);
        };
    } // namespace damping
}     // namespace suanpan

template<typename T> void suanpan::mass::lumped_simple::apply(Mat<T>& mass) { mass = diagmat(sum(mass)); }

template<typename T> void suanpan::mass::lumped_scale::apply(Mat<T>& mass, const unsigned dim) {
    Col<T> diag_mass(mass.n_rows, fill::zeros);

    for(unsigned I = 0; I < dim; ++I) {
        auto total_mass = 0.;
        auto true_mass = 0.;
        for(auto J = I; J < diag_mass.n_elem; J += dim) {
            true_mass += mass(J, J);
            total_mass += sum(mass.row(J));
        }
        if(fabs(true_mass) > 1E-14) {
            auto factor = total_mass / true_mass;
            for(auto J = I; J < diag_mass.n_elem; J += dim) diag_mass(J) = mass(J, J) * factor;
        }
    }

    mass = diagmat(diag_mass);
}

template<typename T> void suanpan::damping::rayleigh::apply(const shared_ptr<Element>& element_obj, const T alpha, const T beta, const T zeta, const T eta) {
    auto& ele_damping = access::rw(element_obj->get_trial_damping());

    if(auto& ele_force = access::rw(element_obj->get_trial_damping_force()); element_obj->if_update_damping()) {
        mat damping(element_obj->get_total_number(), element_obj->get_total_number(), fill::zeros);

        if(0. != alpha && !element_obj->get_current_mass().is_empty()) damping += alpha * element_obj->get_current_mass();
        if(0. != beta && !element_obj->get_current_stiffness().is_empty()) {
            damping += beta * element_obj->get_current_stiffness();
            if(element_obj->is_nlgeom() && !element_obj->get_current_geometry().is_empty()) damping += beta * element_obj->get_current_geometry();
        }
        if(0. != zeta && !element_obj->get_initial_stiffness().is_empty()) damping += zeta * element_obj->get_initial_stiffness();
        if(0. != eta && !element_obj->get_trial_stiffness().is_empty()) {
            damping += eta * element_obj->get_trial_stiffness();
            if(element_obj->is_nlgeom() && !element_obj->get_trial_geometry().is_empty()) damping += eta * element_obj->get_trial_geometry();
        }

        ele_damping = damping;

        ele_force = damping * get_trial_velocity(element_obj.get());
    }
    else if(!ele_damping.is_empty()) ele_force = ele_damping * get_trial_velocity(element_obj.get());
}

template<typename T> void suanpan::damping::elemental::apply(const shared_ptr<Element>& element_obj, const T damping_ratio) {
    const auto& t_stiffness = element_obj->get_current_stiffness();
    const auto& t_mass = element_obj->get_current_mass();

    if(t_stiffness.is_empty() || t_mass.is_empty()) return;

    cx_vec eig_val;
    cx_mat eig_vec;

    auto num_mode = std::min(damping_ratio.n_elem, rank(t_mass));

    if(!eig_pair(eig_val, eig_vec, t_stiffness, t_mass)) return;

    const mat theta = t_mass * abs(eig_vec.head_cols(num_mode));
    const mat damping = theta * diagmat(2. * sqrt(abs(eig_val.head(num_mode)) * damping_ratio)) * theta.t();

    access::rw(element_obj->get_trial_damping()) = damping;

    access::rw(element_obj->get_trial_damping_force()) = damping * get_trial_velocity(element_obj.get());
}

#endif

//! @}
