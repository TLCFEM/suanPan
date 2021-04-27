/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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
	auto counter = 0llu;

	for(auto& [node, pos, idx] : slave) {
		pos.zeros(3);
		auto& t_coor = node.lock()->get_coordinate();
		for(uword K = 0; K < std::min(3llu, t_coor.n_elem); ++K) pos(K) += t_coor(K);
		auto& t_disp = node.lock()->get_trial_displacement();
		for(uword K = 0; K < std::min(3llu, t_disp.n_elem); ++K) pos(K) += t_disp(K);
		idx.set_size(3);
		idx.imbue([&]() { return counter++; });
	}

	for(auto& [m_node, edge_outer_norm, facet_outer_norm] : master) {
		for(auto& [node, pos, idx] : m_node) {
			pos.zeros(3);
			auto& t_coor = node.lock()->get_coordinate();
			for(auto K = 0llu; K < std::min(3llu, t_coor.n_elem); ++K) pos(K) += t_coor(K);
			auto& t_disp = node.lock()->get_trial_displacement();
			for(auto K = 0llu; K < std::min(3llu, t_disp.n_elem); ++K) pos(K) += t_disp(K);
			idx.set_size(3);
			idx.imbue([&]() { return counter++; });
		}

		const vec edge_i = m_node[1].position - m_node[0].position;
		const vec edge_j = m_node[2].position - m_node[1].position;
		const vec edge_k = m_node[0].position - m_node[2].position;

		facet_outer_norm = normalise(cross(edge_i, edge_j));

		edge_outer_norm[0] = normalise(cross(edge_i, facet_outer_norm));
		edge_outer_norm[1] = normalise(cross(edge_j, facet_outer_norm));
		edge_outer_norm[2] = normalise(cross(edge_k, facet_outer_norm));
	}
}

void Contact3D::check_contact(const MasterFacet& m, const SlaveNode& s) {
	const vec s_i = s.position - m.node[0].position;

	const double pen = dot(m.facet_outer_norm, s_i);

	if(pen > datum::eps) return;

	const vec s_j = s.position - m.node[1].position;
	const vec s_k = s.position - m.node[2].position;

	vec h(3);

	h(0) = -dot(s_i, m.edge_outer_norm[0]);
	h(1) = -dot(s_j, m.edge_outer_norm[1]);
	h(2) = -dot(s_k, m.edge_outer_norm[2]);

	if(h(0) < datum::eps || h(1) < datum::eps || h(2) < datum::eps) return;

	mat edge(3, 3);

	edge.col(0) = m.node[1].position - m.node[0].position;
	edge.col(1) = m.node[2].position - m.node[1].position;
	edge.col(2) = m.node[0].position - m.node[2].position;

	vec l(3);

	const vec unit_i = edge.col(0) / (l(0) = norm(edge.col(0)));
	const vec unit_j = edge.col(1) / (l(1) = norm(edge.col(1)));
	const vec unit_k = edge.col(2) / (l(2) = norm(edge.col(2)));

	const auto area = norm(cross(edge.col(2), edge.col(0))); // double the area

	const vec n = (h % l / area).eval()(uvec{1, 2, 0});

	const vec resistance = alpha * pen * m.facet_outer_norm;

	trial_resistance.head(3) -= resistance;
	for(auto I = 0; I < 3; ++I) trial_resistance(m.node[I].local_span) += n(I) * resistance;

	const rowvec dlidj = unit_i.t();
	const rowvec dlidi = -dlidj;

	const rowvec dljdk = unit_j.t();
	const rowvec dljdj = -dljdk;

	const rowvec dlkdi = unit_k.t();
	const rowvec dlkdk = -dlkdi;

	const mat skew_i = transform::skew_symm(edge.col(0).eval());
	const mat skew_j = transform::skew_symm(edge.col(1).eval());
	const mat skew_k = transform::skew_symm(edge.col(2).eval());

	const rowvec dadi = m.facet_outer_norm.t() * skew_j;
	const rowvec dadj = m.facet_outer_norm.t() * skew_k;
	const rowvec dadk = m.facet_outer_norm.t() * skew_i;

	const mat dneidj = (eye(3, 3) - unit_i * unit_i.t()) / l(0);
	const mat dneidi = -dneidj;

	const mat dnejdk = (eye(3, 3) - unit_j * unit_j.t()) / l(1);
	const mat dnejdj = -dnejdk;

	const mat dnekdi = (eye(3, 3) - unit_k * unit_k.t()) / l(2);
	const mat dnekdk = -dnekdi;

	const mat skew_nei = transform::skew_symm(unit_i);
	const mat skew_nej = transform::skew_symm(unit_j);
	const mat skew_nek = transform::skew_symm(unit_k);
	const mat skew_nf = transform::skew_symm(m.facet_outer_norm);

	const mat dnfd = (ones(3, 3) - m.facet_outer_norm.t() * m.facet_outer_norm) / area;

	const mat dnfdi = dnfd * skew_j;
	const mat dnfdj = dnfd * skew_k;
	const mat dnfdk = dnfd * skew_i;

	const mat dnoidi = skew_nei * dnfdi - skew_nf * dneidi;
	const mat dnoidj = skew_nei * dnfdj - skew_nf * dneidj;
	const mat dnoidk = skew_nei * dnfdk;

	const mat dnojdi = skew_nej * dnfdi;
	const mat dnojdj = skew_nej * dnfdj - skew_nf * dnejdj;
	const mat dnojdk = skew_nej * dnfdk - skew_nf * dnejdk;

	const mat dnokdi = skew_nek * dnfdi - skew_nf * dnekdi;
	const mat dnokdj = skew_nek * dnfdj;
	const mat dnokdk = skew_nek * dnfdk - skew_nf * dnekdk;

	const rowvec dhids = -m.edge_outer_norm[0].t();
	const rowvec dhjds = -m.edge_outer_norm[1].t();
	const rowvec dhkds = -m.edge_outer_norm[2].t();

	const rowvec dhidi = -s_i.t() * dnoidi + m.edge_outer_norm[0].t();
	const rowvec dhidj = -s_i.t() * dnoidj;
	const rowvec dhidk = -s_i.t() * dnoidk;

	const rowvec dhjdi = -s_j.t() * dnojdi;
	const rowvec dhjdj = -s_j.t() * dnojdj + m.edge_outer_norm[1].t();
	const rowvec dhjdk = -s_j.t() * dnojdk;

	const rowvec dhkdi = -s_k.t() * dnokdi;
	const rowvec dhkdj = -s_k.t() * dnokdj;
	const rowvec dhkdk = -s_k.t() * dnokdk + m.edge_outer_norm[2].t();

	const rowvec dnids = l(1) / area * dhjds;
	const rowvec dnjds = l(2) / area * dhkds;
	const rowvec dnkds = l(0) / area * dhids;

	const rowvec dnidi = (dhjdi - h(1) / area * dadi) * l(1) / area;
	const rowvec dnidj = (l(1) * dhjdj + h(1) * dljdj - h(1) * l(1) / area * dadj) / area;
	const rowvec dnidk = (l(1) * dhjdk + h(1) * dljdk - h(1) * l(1) / area * dadk) / area;

	const rowvec dnjdi = (l(2) * dhkdi + h(2) * dlkdi - h(2) * l(2) / area * dadi) / area;
	const rowvec dnjdj = (dhkdj - h(2) / area * dadj) * l(2) / area;
	const rowvec dnjdk = (l(2) * dhkdk + h(2) * dlkdk - h(2) * l(2) / area * dadk) / area;

	const rowvec dnkdi = (l(0) * dhidi + h(0) * dlidi - h(0) * l(0) / area * dadi) / area;
	const rowvec dnkdj = (l(0) * dhidj + h(0) * dlidj - h(0) * l(0) / area * dadj) / area;
	const rowvec dnkdk = (dhidk - h(0) / area * dadk) * l(0) / area;

	const rowvec dpds = m.facet_outer_norm.t();

	const rowvec dpdi = -m.facet_outer_norm.t() + s_i.t() * dnfdi;
	const rowvec dpdj = s_i.t() * dnfdj;
	const rowvec dpdk = s_i.t() * dnfdk;

	const mat drds = alpha * m.facet_outer_norm * dpds;
	const mat drdi = alpha * (pen * dnfdi + m.facet_outer_norm * dpdi);
	const mat drdj = alpha * (pen * dnfdj + m.facet_outer_norm * dpdj);
	const mat drdk = alpha * (pen * dnfdk + m.facet_outer_norm * dpdk);

	trial_stiffness(s.local_span, s.local_span) -= drds;
	trial_stiffness(s.local_span, m.node[0].local_span) -= drdi;
	trial_stiffness(s.local_span, m.node[1].local_span) -= drdj;
	trial_stiffness(s.local_span, m.node[2].local_span) -= drdk;

	trial_stiffness(m.node[0].local_span, s.local_span) += n(0) * drds + resistance * dnids;
	trial_stiffness(m.node[0].local_span, m.node[0].local_span) += n(0) * drdi + resistance * dnidi;
	trial_stiffness(m.node[0].local_span, m.node[1].local_span) += n(0) * drdj + resistance * dnidj;
	trial_stiffness(m.node[0].local_span, m.node[2].local_span) += n(0) * drdk + resistance * dnidk;

	trial_stiffness(m.node[1].local_span, s.local_span) += n(1) * drds + resistance * dnjds;
	trial_stiffness(m.node[1].local_span, m.node[0].local_span) += n(1) * drdi + resistance * dnjdi;
	trial_stiffness(m.node[1].local_span, m.node[1].local_span) += n(1) * drdj + resistance * dnjdj;
	trial_stiffness(m.node[1].local_span, m.node[2].local_span) += n(1) * drdk + resistance * dnjdk;

	trial_stiffness(m.node[2].local_span, s.local_span) += n(2) * drds + resistance * dnkds;
	trial_stiffness(m.node[2].local_span, m.node[0].local_span) += n(2) * drdi + resistance * dnkdi;
	trial_stiffness(m.node[2].local_span, m.node[1].local_span) += n(2) * drdj + resistance * dnkdj;
	trial_stiffness(m.node[2].local_span, m.node[2].local_span) += n(2) * drdk + resistance * dnkdk;
}

Contact3D::Contact3D(const unsigned T, const unsigned M, const unsigned S, const double P)
	: Element(T, c_dof, {M, S})
	, master_tag(M)
	, slave_tag(S)
	, alpha(P) {}

void Contact3D::initialize(const shared_ptr<DomainBase>& D) {
	if(!D->find<Group>(master_tag) || !D->find<Group>(slave_tag)) {
		suanpan_error("Contact3D: cannot find the given groups for element %u.\n", get_tag());
		D->disable_element(get_tag());
		return;
	}

	const auto& m_pool = D->get<Group>(master_tag)->get_pool();

	if(0 != m_pool.n_elem % 3) {
		suanpan_error("Contact3D: master group has wrong number of nodes, now disable element %u and return.", get_tag());
		D->disable_element(get_tag());
		return;
	}

	master.reserve(m_pool.n_elem / 3);
	for(auto I = 0llu, J = 1llu, K = 2llu; K < m_pool.n_elem; I += 3llu, J += 3llu, K += 3llu) {
		if(!D->find<Node>(m_pool(I)) || !D->find<Node>(m_pool(J)) || !D->find<Node>(m_pool(K))) {
			suanpan_error("Contact3D: invalid Node %llu detected, now disable element %u and return.", I, get_tag());
			D->disable_element(get_tag());
			return;
		}
		master.emplace_back(MasterFacet{{ContactNode{D->get<Node>(m_pool(I)), {}, {}}, ContactNode{D->get<Node>(m_pool(J)), {}, {}}, ContactNode{D->get<Node>(m_pool(K)), {}, {}}}, {}, {}});
	}

	master.shrink_to_fit();

	const auto& s_pool = D->get<Group>(slave_tag)->get_pool();
	slave.reserve(s_pool.n_elem);
	for(auto& I : s_pool) {
		if(!D->find<Node>(I)) {
			suanpan_error("Contact3D: invalid Node %llu detected, now disable element %u and return.", I, get_tag());
			D->disable_element(get_tag());
			return;
		}
		slave.emplace_back(SlaveNode{D->get<Node>(I), {}, {}});
	}
	slave.shrink_to_fit();
}

int Contact3D::update_status() {
	update_position();

	trial_stiffness.zeros(get_total_number(), get_total_number());
	trial_resistance.zeros(get_total_number());

	for(auto& I : master) for(auto& J : slave) check_contact(I, J);

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
