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

#include "Converger.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

/**
 * \brief the complete constructor.
 * \param T `unique_tag`
 * \param E `tolerance`
 * \param M `maximum_iteration`
 * \param P `print_flag`
 */
Converger::Converger(const unsigned T, const double E, const unsigned M, const bool P)
    : CopiableTag(T)
    , tolerance(E)
    , max_iteration(M)
    , print_flag(P) {}

int Converger::initialize() {
    if(nullptr == database.lock()) {
        suanpan_error("A valid domain is required.\n");
        return SUANPAN_FAIL;
    }

    return SUANPAN_SUCCESS;
}

/**
 * \brief method to set `tolerance`.
 * \param T `tolerance`
 */
void Converger::set_tolerance(const double T) { tolerance = T; }

/**
 * \brief method to return `tolerance`.
 * \return `tolerance`
 */
double Converger::get_tolerance() const { return tolerance; }

void Converger::set_max_iteration(const unsigned M) { max_iteration = M; }

unsigned Converger::get_max_iteration() const { return max_iteration; }

/**
 * \brief method to set `DomainBase`.
 * \param D `DomainBase`
 */
void Converger::set_domain(const std::weak_ptr<DomainBase>& D) {
    if(database.lock() != D.lock()) database = D;
}

/**
 * \brief method to return `DomainBase`.
 * \return `DomainBase`
 */
const std::weak_ptr<DomainBase>& Converger::get_domain() const { return database; }

/**
 * \brief method to set `error`.
 * \param E `error`
 */
void Converger::set_error(const double E) { error = E; }

/**
 * \brief method to return `error`.
 * \return `error`
 */
double Converger::get_error() const { return error; }

/**
 * \brief method to set `conv_flag`.
 * \param C `conv_flag`
 */
void Converger::set_conv_flag(const bool C) { conv_flag = C; }

/**
 * \brief method to return `conv_flag`.
 * \return `conv_flag`
 */
bool Converger::get_conv_flag() const { return conv_flag; }

vec Converger::get_residual() const {
    const auto D = get_domain().lock();
    auto& W = D->get_factory();

    vec residual = W->get_trial_load() - W->get_sushi();
    if(!W->get_reference_load().is_empty() && !W->get_trial_load_factor().is_empty()) residual += W->get_reference_load() * W->get_trial_load_factor();

    for(const auto& t_dof : D->get_restrained_dof()) residual(t_dof) = 0.;

    return residual;
}

/**
 * \brief method to return `print_flag`.
 * \return `print_flag`
 */
bool Converger::is_print() const { return print_flag; }
