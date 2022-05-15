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

#include "Contact3D.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Group.h>
#include <Domain/Node.h>
#include <Toolbox/tensorToolbox.h>

void Contact3D::update_position() {
    for(auto& [node, span, position] : slave) {
        position.zeros(3);
        auto& t_coor = node.lock()->get_coordinate();
        for(auto K = 0llu; K < std::min(3llu, t_coor.n_elem); ++K) position(K) += t_coor(K);
        auto& t_disp = node.lock()->get_trial_displacement();
        for(auto K = 0llu; K < std::min(3llu, t_disp.n_elem); ++K) position(K) += t_disp(K);
    }

    for(auto& [m_node, f_norm, f_area] : master) {
        for(auto& [node, span, position, e_norm] : m_node) {
            position.zeros(3);
            auto& t_coor = node.lock()->get_coordinate();
            for(auto K = 0llu; K < std::min(3llu, t_coor.n_elem); ++K) position(K) += t_coor(K);
            auto& t_disp = node.lock()->get_trial_displacement();
            for(auto K = 0llu; K < std::min(3llu, t_disp.n_elem); ++K) position(K) += t_disp(K);
        }

        constexpr unsigned i = 0, j = 1, k = 2;

        const vec edge_i = m_node[k].position - m_node[j].position;
        const vec edge_j = m_node[i].position - m_node[k].position;
        const vec edge_k = m_node[j].position - m_node[i].position;

        f_norm = cross(edge_j, -edge_i);
        f_area = norm(f_norm);
        f_norm /= f_area;

        m_node[i].outer_norm = cross(edge_i, f_norm);
        m_node[j].outer_norm = cross(edge_j, f_norm);
        m_node[k].outer_norm = cross(edge_k, f_norm);
    }
}

void Contact3D::check_contact(const MasterFacet& m, const SlaveNode& s) {
    constexpr unsigned i = 0, j = 1, k = 2;

    const vec s_i = s.position - m.node[i].position;
    const vec s_j = s.position - m.node[j].position;
    const vec s_k = s.position - m.node[k].position;

    const double pen = dot(m.facet_outer_norm, s_i);

    vec h(3);

    h(i) = dot(s_j, m.node[i].outer_norm);
    h(j) = dot(s_k, m.node[j].outer_norm);
    h(k) = dot(s_i, m.node[k].outer_norm);

    // check if the slave node penetrates the master facet
    if(pen > datum::eps) return;

    // check if penetration point is located inside of the facet
    if(h(i) > 0. || h(j) > 0. || h(k) > 0.) return;

    const vec n = h / m.facet_area;

    vec resistance = pen * m.facet_outer_norm;

    trial_resistance(s.local_span) += resistance;
    trial_resistance(m.node[i].local_span) += n(i) * resistance;
    trial_resistance(m.node[j].local_span) += n(j) * resistance;
    trial_resistance(m.node[k].local_span) += n(k) * resistance;

    const mat skew_fn = transform::skew_symm(m.facet_outer_norm);

    const mat skew_edge_i = transform::skew_symm(vec(m.node[k].position - m.node[j].position));
    const mat skew_edge_j = transform::skew_symm(vec(m.node[i].position - m.node[k].position));
    const mat skew_edge_k = transform::skew_symm(vec(m.node[j].position - m.node[i].position));

    const rowvec dadi = m.facet_outer_norm.t() * skew_edge_i;
    const rowvec dadj = m.facet_outer_norm.t() * skew_edge_j;
    const rowvec dadk = m.facet_outer_norm.t() * skew_edge_k;

    const mat dfna = eye(3, 3) - m.facet_outer_norm * m.facet_outer_norm.t();
    const mat dfndi = dfna * skew_edge_i / m.facet_area;
    const mat dfndj = dfna * skew_edge_j / m.facet_area;
    const mat dfndk = dfna * skew_edge_k / m.facet_area;

    const mat dra = m.facet_outer_norm * s_i.t() + pen * eye(3, 3);
    const mat drds = m.facet_outer_norm * m.facet_outer_norm.t();
    const mat drdi = dra * dfndi - drds;
    const mat drdj = dra * dfndj;
    const mat drdk = dra * dfndk;

    const rowvec si_fn = s_i.t() * skew_fn;
    const rowvec sj_fn = s_j.t() * skew_fn;
    const rowvec sk_fn = s_k.t() * skew_fn;

    const rowvec si_skew = s_i.t() * skew_edge_k;
    const rowvec sj_skew = s_j.t() * skew_edge_i;
    const rowvec sk_skew = s_k.t() * skew_edge_j;

    h /= m.facet_area;
    resistance /= m.facet_area;

    trial_stiffness(s.local_span, s.local_span) += drds;
    trial_stiffness(s.local_span, m.node[i].local_span) += drdi;
    trial_stiffness(s.local_span, m.node[j].local_span) += drdj;
    trial_stiffness(s.local_span, m.node[k].local_span) += drdk;

    trial_stiffness(m.node[i].local_span, s.local_span) += n(i) * drds + resistance * m.node[i].outer_norm.t();
    trial_stiffness(m.node[i].local_span, m.node[i].local_span) += n(i) * drdi + resistance * (sj_skew * dfndi - h(i) * dadi);
    trial_stiffness(m.node[i].local_span, m.node[j].local_span) += n(i) * drdj + resistance * (sj_skew * dfndj + sj_fn - m.node[i].outer_norm.t() - h(i) * dadj);
    trial_stiffness(m.node[i].local_span, m.node[k].local_span) += n(i) * drdk + resistance * (sj_skew * dfndk - sj_fn - h(i) * dadk);

    trial_stiffness(m.node[j].local_span, s.local_span) += n(j) * drds + resistance * m.node[j].outer_norm.t();
    trial_stiffness(m.node[j].local_span, m.node[i].local_span) += n(j) * drdi + resistance * (sk_skew * dfndi - sk_fn - h(j) * dadi);
    trial_stiffness(m.node[j].local_span, m.node[j].local_span) += n(j) * drdj + resistance * (sk_skew * dfndj - h(j) * dadj);
    trial_stiffness(m.node[j].local_span, m.node[k].local_span) += n(j) * drdk + resistance * (sk_skew * dfndk + sk_fn - m.node[j].outer_norm.t() - h(j) * dadk);

    trial_stiffness(m.node[k].local_span, s.local_span) += n(k) * drds + resistance * m.node[k].outer_norm.t();
    trial_stiffness(m.node[k].local_span, m.node[i].local_span) += n(k) * drdi + resistance * (si_skew * dfndi + si_fn - m.node[k].outer_norm.t() - h(k) * dadi);
    trial_stiffness(m.node[k].local_span, m.node[j].local_span) += n(k) * drdj + resistance * (si_skew * dfndj - si_fn - h(k) * dadj);
    trial_stiffness(m.node[k].local_span, m.node[k].local_span) += n(k) * drdk + resistance * (si_skew * dfndk - h(k) * dadk);
}

Contact3D::Contact3D(const unsigned T, const unsigned M, const unsigned S, const double P)
    : Element(T, c_dof, {M, S})
    , master_tag(M)
    , slave_tag(S)
    , alpha(P) {}

int Contact3D::initialize(const shared_ptr<DomainBase>& D) {
    const auto& m_pool = D->get<Group>(master_tag)->get_pool();

    if(0 != m_pool.n_elem % 3) {
        suanpan_error("Contact3D %u master group has wrong number of nodes.", get_tag());
        return SUANPAN_FAIL;
    }

    master.reserve(m_pool.n_elem / 3);
    for(auto I = 0llu, J = 1llu, K = 2llu; K < m_pool.n_elem; I += 3llu, J += 3llu, K += 3llu) master.emplace_back(MasterFacet{{MasterNode{D->get<Node>(m_pool(I)), {}, {}, {}}, MasterNode{D->get<Node>(m_pool(J)), {}, {}, {}}, MasterNode{D->get<Node>(m_pool(K)), {}, {}, {}}}, {}, {}});

    master.shrink_to_fit();

    const auto& s_pool = D->get<Group>(slave_tag)->get_pool();
    slave.reserve(s_pool.n_elem);
    for(auto& I : s_pool) slave.emplace_back(SlaveNode{D->get<Node>(I), {}, {}});
    slave.shrink_to_fit();

    auto counter = 0llu;

    for(auto& [m_node, f_norm, f_area] : master)
        for(auto& [node, span, pos, e_norm] : m_node) {
            span.set_size(3);
            span.imbue([&] { return counter++; });
        }

    for(auto& [node, span, pos] : slave) {
        span.set_size(3);
        span.imbue([&] { return counter++; });
    }

    return SUANPAN_SUCCESS;
}

int Contact3D::update_status() {
    update_position();

    trial_stiffness.zeros(get_total_number(), get_total_number());
    trial_resistance.zeros(get_total_number());

    for(auto& I : master) for(auto& J : slave) check_contact(I, J);

    trial_stiffness *= alpha;
    trial_resistance *= alpha;

    return SUANPAN_SUCCESS;
}

int Contact3D::clear_status() {
    current_stiffness.reset();
    trial_stiffness.reset();
    current_resistance.reset();
    trial_resistance.reset();
    return SUANPAN_SUCCESS;
}

int Contact3D::commit_status() {
    current_stiffness = trial_stiffness;
    current_resistance = trial_resistance;
    return SUANPAN_SUCCESS;
}

int Contact3D::reset_status() {
    trial_stiffness = current_stiffness;
    trial_resistance = current_resistance;
    return SUANPAN_SUCCESS;
}
