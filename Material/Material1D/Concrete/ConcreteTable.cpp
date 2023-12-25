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

#include "ConcreteTable.h"

podarray<double> ConcreteTable::compute_compression_initial_reverse() const {
    podarray<double> response(2);

    for(uword I = 0; I < c_table.n_rows; ++I)
        if(c_table(I, 0) != 0.) {
            response(0) = c_table(I, 0);
            response(1) = c_table(I, 1);
            break;
        }

    return response;
}

podarray<double> ConcreteTable::compute_tension_initial_reverse() const {
    podarray<double> response(2);

    for(uword I = 0; I < t_table.n_rows; ++I)
        if(t_table(I, 0) != 0.) {
            response(0) = t_table(I, 0);
            response(1) = t_table(I, 1);
            break;
        }

    return response;
}

podarray<double> ConcreteTable::compute_compression_backbone(const double n_strain) const {
    podarray<double> response(2);

    for(uword I = 0; I < c_table.n_rows; ++I)
        if(c_table(I, 0) < n_strain) {
            if(0 == I) {
                response(1) = c_table(I, 1) / c_table(I, 0);
                response(0) = n_strain * response(1);
            }
            else {
                response(1) = (c_table(I, 1) - c_table(I - 1, 1)) / (c_table(I, 0) - c_table(I - 1, 0));
                response(0) = c_table(I - 1, 1) + (n_strain - c_table(I - 1, 0)) * response(1);
            }
            return response;
        }

    response(1) = 1E-10;
    response(0) = c_table(c_table.n_rows - 1, 0);

    return response;
}

podarray<double> ConcreteTable::compute_tension_backbone(const double n_strain) const {
    podarray<double> response(2);

    for(uword I = 0; I < t_table.n_rows; ++I)
        if(t_table(I, 0) > n_strain) {
            if(0 == I) {
                response(1) = t_table(I, 1) / t_table(I, 0);
                response(0) = n_strain * response(1);
            }
            else {
                response(1) = (t_table(I, 1) - t_table(I - 1, 1)) / (t_table(I, 0) - t_table(I - 1, 0));
                response(0) = t_table(I - 1, 1) + (n_strain - t_table(I - 1, 0)) * response(1);
            }
            return response;
        }

    response(1) = 1E-10;
    response(0) = t_table(t_table.n_rows - 1, 0);

    return response;
}

double ConcreteTable::compute_compression_residual(const double reverse_c_strain, const double reverse_c_stress) const { return std::min(0., reverse_c_strain - reverse_c_stress * (reverse_c_strain / c_strain + .57) / (reverse_c_stress / c_strain + .57 * initial_stiffness(0))); }

double ConcreteTable::compute_tension_residual(const double reverse_t_strain, const double reverse_t_stress) const { return std::max(0., reverse_t_strain - reverse_t_stress * (reverse_t_strain / t_strain + .67) / (reverse_t_stress / t_strain + .67 * initial_stiffness(0))); }

ConcreteTable::ConcreteTable(const unsigned T, mat&& CT, mat&& TT, const double MP, const double R)
    : SimpleHysteresis(T, MP, R)
    , c_table(std::move(CT))
    , t_table(std::move(TT))
    , c_strain(c_table.col(0)(index_min(c_table.col(1))))
    , t_strain(t_table.col(0)(index_max(t_table.col(1)))) {}

int ConcreteTable::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = .5 * ((c_table(1, 1) - c_table(0, 1)) / (c_table(1, 0) - c_table(0, 0)) + (t_table(1, 1) - t_table(0, 1)) / (t_table(1, 0) - t_table(0, 0)));

    return SUANPAN_SUCCESS;
}

double ConcreteTable::get_parameter(const ParameterType P) const {
    if(ParameterType::ELASTICMODULUS == P) return initial_stiffness(0);
    if(ParameterType::PEAKSTRAIN == P) return compute_compression_initial_reverse()(0);
    if(ParameterType::CRACKSTRAIN == P) return compute_tension_initial_reverse()(0);
    return 0.;
}

unique_ptr<Material> ConcreteTable::get_copy() { return make_unique<ConcreteTable>(*this); }
