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

#include "SectionShell.h"
#include <Domain/DomainBase.h>

SectionShell::SectionShell(const unsigned T, const unsigned MT, vec&& E)
    : SectionShellData{MT, std::forward<vec>(E)}
    , Tag(T) {}

void SectionShell::set_initialized(const bool T) const { access::rw(initialized) = T; }

void SectionShell::set_symmetric(const bool T) const { access::rw(symmetric) = T; }

bool SectionShell::is_initialized() const { return initialized; }

bool SectionShell::is_symmetric() const { return symmetric; }

void SectionShell::set_eccentricity(const vec& E) const { access::rw(eccentricity) = E; }

const vec& SectionShell::get_eccentricity() const { return eccentricity; }

const vec& SectionShell::get_trial_membrane_strain() const { return trial_membrane_strain; }

const vec& SectionShell::get_trial_membrane_strain_rate() const { return trial_membrane_strain_rate; }

const vec& SectionShell::get_trial_plate_strain() const { return trial_plate_strain; }

const vec& SectionShell::get_trial_plate_strain_rate() const { return trial_plate_strain_rate; }

const vec& SectionShell::get_trial_membrane_stress() const { return trial_membrane_stress; }

const vec& SectionShell::get_trial_plate_stress() const { return trial_plate_stress; }

const mat& SectionShell::get_trial_membrane_stiffness() const { return trial_membrane_stiffness; }

const mat& SectionShell::get_trial_plate_stiffness() const { return trial_plate_stiffness; }

const vec& SectionShell::get_current_membrane_strain() const { return current_membrane_strain; }

const vec& SectionShell::get_current_membrane_strain_rate() const { return current_membrane_strain_rate; }

const vec& SectionShell::get_current_plate_strain() const { return current_plate_strain; }

const vec& SectionShell::get_current_plate_strain_rate() const { return current_plate_strain_rate; }

const vec& SectionShell::get_current_membrane_stress() const { return current_membrane_stress; }

const vec& SectionShell::get_current_plate_stress() const { return current_plate_stress; }

const mat& SectionShell::get_current_membrane_stiffness() const { return current_membrane_stiffness; }

const mat& SectionShell::get_current_plate_stiffness() const { return current_plate_stiffness; }

const mat& SectionShell::get_initial_membrane_stiffness() const { return initial_membrane_stiffness; }

const mat& SectionShell::get_initial_plate_stiffness() const { return initial_plate_stiffness; }

int SectionShell::update_incre_status(const double ME, const double PE) {
    const vec m_strain{ME}, p_strain{PE};
    return update_incre_status(m_strain, p_strain);
}

int SectionShell::update_incre_status(const double ME, const double PE, const double MER, const double PER) {
    const vec m_strain{ME}, p_strain{PE}, m_strain_rate{MER}, p_strain_rate{PER};
    return update_incre_status(m_strain, p_strain, m_strain_rate, p_strain_rate);
}

int SectionShell::update_trial_status(const double ME, const double PE) {
    const vec m_strain{ME}, p_strain{PE};
    return update_trial_status(m_strain, p_strain);
}

int SectionShell::update_trial_status(const double ME, const double PE, const double MER, const double PER) {
    const vec m_strain{ME}, p_strain{PE}, m_strain_rate{MER}, p_strain_rate{PER};
    return update_trial_status(m_strain, p_strain, m_strain_rate, p_strain_rate);
}

int SectionShell::update_incre_status(const vec& m_strain, const vec& p_strain) { return update_trial_status(current_membrane_strain + m_strain, current_plate_strain + p_strain); }

int SectionShell::update_incre_status(const vec& m_strain, const vec& p_strain, const vec& m_strain_rate, const vec& p_strain_rate) { return update_trial_status(current_membrane_strain + m_strain, current_plate_strain + p_strain, current_membrane_strain_rate + m_strain_rate, current_plate_strain_rate + p_strain_rate); }

int SectionShell::update_trial_status(const vec&, const vec&) { throw invalid_argument("hidden method update_trial_status() called"); }

int SectionShell::update_trial_status(const vec& m_strain, const vec& p_strain, const vec& m_strain_rate, const vec& p_strain_rate) {
    trial_membrane_strain_rate = m_strain_rate;
    trial_plate_strain_rate = p_strain_rate;
    return update_trial_status(m_strain, p_strain);
}

vector<vec> SectionShell::record(OutputType) { return {}; }

unique_ptr<SectionShell> suanpan::make_copy(const shared_ptr<SectionShell>& S) { return S->get_copy(); }

unique_ptr<SectionShell> suanpan::make_copy(const unique_ptr<SectionShell>& S) { return S->get_copy(); }
