/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "Contact2D.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Group/Group.h>
#include <Domain/Node.h>

const mat Contact2D::rotation{{0., -1.}, {1., 0.}};

void Contact2D::update_position() {
    for(auto& [node, position, axis, norm] : master) {
        position.zeros(2);
        auto& t_coor = node.lock()->get_coordinate();
        for(uword K = 0; K < std::min(2llu, t_coor.n_elem); ++K) position(K) += t_coor(K);
        auto& t_disp = node.lock()->get_trial_displacement();
        for(uword K = 0; K < std::min(2llu, t_disp.n_elem); ++K) position(K) += t_disp(K);
    }

    for(decltype(master.size()) I = 0, J = 1; I < master.size() - 1; ++I, ++J) master[I].norm = rotation * normalise(master[I].axis = master[J].position - master[I].position);

    for(auto& [node, position] : slave) {
        position.zeros(2);
        auto& t_coor = node.lock()->get_coordinate();
        for(uword K = 0; K < std::min(2llu, t_coor.n_elem); ++K) position(K) += t_coor(K);
        auto& t_disp = node.lock()->get_trial_displacement();
        for(uword K = 0; K < std::min(2llu, t_disp.n_elem); ++K) position(K) += t_disp(K);
    }
}

Contact2D::Contact2D(const unsigned T, const unsigned M, const unsigned S, const double P)
    : Element(T, c_dof, {M, S})
    , master_tag(M)
    , slave_tag(S)
    , alpha(P) {}

int Contact2D::initialize(const shared_ptr<DomainBase>& D) {
    const auto& m_pool = D->get<Group>(master_tag)->get_pool();
    master.reserve(m_pool.n_elem);
    for(auto& I : m_pool) master.emplace_back(MasterNode{D->get<Node>(I), {}, {}, {}});
    master.shrink_to_fit();

    const auto& s_pool = D->get<Group>(slave_tag)->get_pool();
    slave.reserve(s_pool.n_elem);
    for(auto& I : s_pool) slave.emplace_back(SlaveNode{D->get<Node>(I), {}});
    slave.shrink_to_fit();

    return SUANPAN_SUCCESS;
}

int Contact2D::update_status() {
    update_position();

    trial_stiffness.zeros(get_total_number(), get_total_number());
    trial_resistance.zeros(get_total_number());

    for(decltype(master.size()) I = 0; I < master.size() - 1; ++I)
        for(decltype(slave.size()) J = 0; J < slave.size(); ++J) {
            const vec s = slave[J].position - master[I].position; // 1-3
            const auto tn = dot(s, master[I].axis);               // numerator, 1-2-3
            const auto td = dot(master[I].axis, master[I].axis);  // denominator, 1-2
            const auto t = tn / td;                               // numerator/denominator ==> ranges from zero to one, 1-2-3
            if(const auto u = dot(s, master[I].norm); u < datum::eps && t <= 1. && t >= 0.) {
                auto K = c_dof * I;
                const auto span_a = span(K, K + c_dof - 1);
                K += c_dof;
                const auto span_b = span(K, K + c_dof - 1);
                K = c_dof * (master.size() + J);
                const auto span_c = span(K, K + c_dof - 1);

                const vec reaction = u * alpha * master[I].norm;
                trial_resistance(span_a) += (t - 1.) * reaction;
                trial_resistance(span_b) += -t * reaction;
                trial_resistance(span_c) += reaction;

                const mat dn = rotation / sqrt(td) - master[I].norm * master[I].axis.t() / td;

                const rowvec pupc = master[I].norm.t();
                const rowvec pupb = s.t() * dn;
                const rowvec pupa = -pupb - pupc;

                const rowvec plpc = u * master[I].axis.t() / td;
                const rowvec plpb = u / td * s.t() - 2. * t * plpc;
                const rowvec plpa = -plpb - plpc;

                trial_stiffness(span_a, span_a) += alpha * (master[I].norm * (plpa + (t - 1.) * pupa) - (t - 1.) * u * dn);
                trial_stiffness(span_a, span_b) += alpha * (master[I].norm * (plpb + (t - 1.) * pupb) + (t - 1.) * u * dn);
                trial_stiffness(span_a, span_c) += alpha * master[I].norm * (plpc + (t - 1.) * pupc);
                trial_stiffness(span_b, span_a) += alpha * (master[I].norm * (-plpa - t * pupa) + t * u * dn);
                trial_stiffness(span_b, span_b) += alpha * (master[I].norm * (-plpb - t * pupb) - t * u * dn);
                trial_stiffness(span_b, span_c) += alpha * master[I].norm * (-plpc - t * pupc);
                trial_stiffness(span_c, span_a) += alpha * (master[I].norm * pupa - u * dn);
                trial_stiffness(span_c, span_b) += alpha * (master[I].norm * pupb + u * dn);
                trial_stiffness(span_c, span_c) += alpha * master[I].norm * pupc;
            }
        }

    return SUANPAN_SUCCESS;
}

int Contact2D::clear_status() {
    current_stiffness.reset();
    trial_stiffness.reset();
    current_resistance.reset();
    trial_resistance.reset();
    return SUANPAN_SUCCESS;
}

int Contact2D::commit_status() {
    current_stiffness = trial_stiffness;
    current_resistance = trial_resistance;
    return SUANPAN_SUCCESS;
}

int Contact2D::reset_status() {
    trial_stiffness = current_stiffness;
    trial_resistance = current_resistance;
    return SUANPAN_SUCCESS;
}
