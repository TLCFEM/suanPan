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

#include "sort_color.h"

unsigned sort_color(vector<vector<unsigned>>& node_register, const unsigned element_size) {
#ifdef SUANPAN_DEBUG
	wall_clock T;
	T.tic();
#endif

	vector<vector<unsigned>> element_register(element_size + 1llu);
	for(const auto& I : node_register) for(const auto& J : I) for(const auto& K : I) element_register[J].emplace_back(K);

	suanpan_for_each(element_register.begin(), element_register.end(), [](vector<unsigned>& element) {
		std::sort(element.begin(), element.end());
		element.erase(std::unique(element.begin(), element.end()), element.end());
	});

	multimap<unsigned, unsigned, std::greater<>> degree;

	for(unsigned I = 0; I != element_register.size(); ++I) if(!element_register[I].empty()) degree.insert({static_cast<unsigned>(element_register[I].size()) - 1, I});

	node_register.clear();

	if(degree.empty()) return 0;

	while(true) {
		node_register.emplace_back();
		auto& c_color = node_register.back();
		c_color.emplace_back(degree.begin()->second);
		for(auto I = degree.erase(degree.begin()); I != degree.end();) {
			auto flag = false;
			for(const auto& K : c_color) {
				if(K == I->second) continue;
				if(const auto& c_link = element_register[K]; c_link.end() != find(c_link.begin(), c_link.end(), I->second)) {
					flag = true;
					break;
				}
			}
			if(flag) ++I;
			else {
				c_color.emplace_back(I->second);
				I = degree.erase(I);
			}
		}

		if(degree.empty()) break;
	}

	node_register.shrink_to_fit();

#ifdef SUANPAN_DEBUG
	suanpan_debug("Coloring algorithm takes %.5E seconds.\n", T.toc());
#endif

	return static_cast<unsigned>(node_register.size());
}
