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

#include "EnergyEvolution.h"
#include <Domain/DomainBase.h>
#include <Element/Element.h>
#include <Step/Step.h>

EnergyEvolution::EnergyEvolution(const unsigned T, const unsigned ST, const unsigned IL, const unsigned FL, const double WT, const unsigned IT, const unsigned RR, const double PW, const double TL)
    : Criterion(T, ST)
    , iteration(IT)
    , reactive_ratio(RR)
    , incre_level(IL)
    , final_level(FL)
    , weight(WT)
    , propagation_weight(PW)
    , tolerance(TL)
    , get_energy(nullptr) {}

int EnergyEvolution::initialize(const shared_ptr<DomainBase>& D) {
    const auto node_pool = D->get_node_connectivity();
    auto element_pool = D->get_element_connectivity();

    map.clear();
    map.resize(element_pool.size(), std::vector<uword>{});

    for(decltype(element_pool.size()) element = 0; element < element_pool.size(); ++element) for(const auto connected_node : element_pool[element]) for(const auto connected_element : node_pool[connected_node]) if(element != connected_element) map[element].emplace_back(connected_element);

    suanpan_for_each(map.begin(), map.end(), [](std::vector<uword>& element) { suanpan::unique(element); });

    energy.zeros(map.size());

    return SUANPAN_SUCCESS;
}

int EnergyEvolution::process(const shared_ptr<DomainBase>& D) {
    if(D->get_current_step()->get_time_left() > 1E-10) return SUANPAN_SUCCESS;

    if(20 == balanced_iteration) return SUANPAN_EXIT;

    vec current_energy(map.size());

    Col<int> state(map.size());

    suanpan_for(0llu, current_energy.n_elem, [&](const uword I) {
        if(D->find<Element>(I)) {
            if(const auto& t_element = D->get<Element>(I); !t_element->is_active()) {
                current_energy(I) = 0.;
                state(I) = 0; // disabled
            }
            else {
                current_energy(I) = get_energy(t_element.get()) / t_element->get_characteristic_length();
                state(I) = 1; // enabled
            }
        }
        else {
            current_energy(I) = 0.;
            state(I) = -1; // not defined
        }
    });

    const auto sum_energy = accu(current_energy);

    if(fabs(sum_energy - total_energy) < tolerance * total_energy) return SUANPAN_EXIT;

    total_energy = sum_energy;

    unsigned counter = 0;
    while(counter++ != iteration) {
        const auto previous_energy = current_energy;
        current_energy *= weight;

        suanpan_for(0llu, current_energy.n_elem, [&](const uword I) {
            for(const auto J : map[I]) current_energy(I) += previous_energy(J);
            current_energy(I) /= static_cast<double>(map[I].size()) + weight;
        });
    }

    suanpan_info("current rejection ratio: %.4E.\n", current_level / 100.);

    energy = (1. - propagation_weight) * energy + propagation_weight * current_energy;

    const uvec disabled_list = find(state == 0);
    uvec sorted_list = sort_index(energy(disabled_list));

    sorted_list = disabled_list(sorted_list.tail(reactive_ratio * disabled_list.n_elem / 100));
    suanpan_for_each(sorted_list.cbegin(), sorted_list.cend(), [&](const uword I) {
        D->enable_element(static_cast<unsigned>(I));
        state(I) = 1;
    });

    const uvec enabled_list = find(state == 1);
    sorted_list = sort_index(energy(enabled_list));
    sorted_list = enabled_list(sorted_list.head(enabled_list.n_elem - std::min((100 - current_level) * current_energy.n_elem / 100, enabled_list.n_elem)));
    suanpan_for_each(sorted_list.cbegin(), sorted_list.cend(), [&](const uword I) { D->disable_element(static_cast<unsigned>(I)); });

    if(D->is_updated() && current_level == final_level) return SUANPAN_EXIT;

    if(current_level == final_level) ++balanced_iteration;
    else current_level += std::min(incre_level, final_level - current_level);

    return SUANPAN_SUCCESS;
}
