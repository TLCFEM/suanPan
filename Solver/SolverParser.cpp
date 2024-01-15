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

#include "SolverParser.h"
#include <Domain/DomainBase.h>
#include <Solver/Solver>
#include <Step/Step.h>
#include <Toolbox/utility.h>

int create_new_integrator(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string integrator_type;
    if(!get_input(command, integrator_type)) {
        suanpan_error("A valid integrator type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    if(if_contain(suanpan::to_upper(cref(integrator_type)), suanpan::to_upper(string("Newmark")))) {
        auto alpha = .25, beta = .5;
        if(!command.eof()) {
            if(!get_input(command, alpha)) {
                suanpan_error("A valid alpha is required.\n");
                return SUANPAN_SUCCESS;
            }
            if(!get_input(command, beta)) {
                suanpan_error("A valid beta is required.\n");
                return SUANPAN_SUCCESS;
            }
        }

        if(is_equal(integrator_type, "Newmark")) { if(domain->insert(make_shared<Newmark>(tag, alpha, beta))) code = 1; }
        else if(is_equal(integrator_type, "RayleighNewmark")) {
            vec p(4, fill::zeros);
            auto idx = 0llu;
            while(!command.eof() && idx < p.n_elem)
                if(!get_input(command, p(idx++))) {
                    suanpan_error("A valid parameter for Rayleigh damping is required.\n");
                    return SUANPAN_SUCCESS;
                }

            if(suanpan::approx_equal(sum(p), 0.))
                suanpan_warning("It seems all parameters are zeros, either the inputs are wrong or consider use plain Newmark.\n");

            if(domain->insert(make_shared<RayleighNewmark>(tag, alpha, beta, p(0), p(1), p(2), p(3)))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeNewmark")) {
            vector<double> damping_coef, frequency;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid damping coefficient is required.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid frequency is required.\n");
                    return SUANPAN_SUCCESS;
                }
                frequency.emplace_back(t_para);
            }

            if(domain->insert(make_shared<LeeNewmark>(tag, damping_coef, frequency, alpha, beta))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeNewmarkIterative")) {
            vector<double> damping_coef, frequency;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid damping coefficient is required.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid frequency is required.\n");
                    return SUANPAN_SUCCESS;
                }
                frequency.emplace_back(t_para);
            }

            if(domain->insert(make_shared<LeeNewmarkIterative>(tag, damping_coef, frequency, alpha, beta))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeElementalNewmark")) {
            vector<double> damping_coef, frequency;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid damping coefficient is required.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid frequency is required.\n");
                    return SUANPAN_SUCCESS;
                }
                frequency.emplace_back(t_para);
            }

            if(domain->insert(make_shared<LeeElementalNewmark>(tag, damping_coef, frequency, alpha, beta))) code = 1;
        }
        else if(integrator_type.size() >= 14 && is_equal(integrator_type.substr(0, 14), "LeeNewmarkFull")) {
            vector<LeeNewmarkFull::Mode> modes;

            auto omega = 0., zeta = 0., para_a = .0, para_b = .0;

            auto get_basic_input = [&] {
                if(!get_input(command, zeta)) {
                    suanpan_error("A valid zeta_p is required.\n");
                    return SUANPAN_FAIL;
                }
                if(!get_input(command, omega)) {
                    suanpan_error("A valid omega_p is required.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            auto get_first = [&] {
                if(!get_input(command, para_a)) {
                    suanpan_error("A valid parameter is required.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            auto get_second = [&] {
                if(!get_input(command, para_b)) {
                    suanpan_error("A valid parameter is required.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            while(!command.eof()) {
                string type;
                if(!get_input(command, type)) {
                    suanpan_error("A valid type is required.\n");
                    return SUANPAN_SUCCESS;
                }
                if(is_equal("-type0", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T0, vec{}, zeta, omega});
                }
                else if(is_equal("-type1", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T1, vec{static_cast<double>(static_cast<unsigned>(para_a))}, zeta, omega});
                }
                else if(is_equal("-type2", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T2, vec{static_cast<double>(static_cast<unsigned>(para_a)), static_cast<double>(static_cast<unsigned>(para_b))}, zeta, omega});
                }
                else if(is_equal("-type3", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T3, vec{para_a}, zeta, omega});
                }
                else if(is_equal("-type4", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return SUANPAN_SUCCESS;
                    double para_c, para_d, para_e;
                    if(!get_input(command, para_c) || !get_input(command, para_d) || !get_input(command, para_e)) {
                        suanpan_error("A valid parameter is required.\n");
                        return SUANPAN_SUCCESS;
                    }
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T4, vec{static_cast<double>(static_cast<unsigned>(para_a)), static_cast<double>(static_cast<unsigned>(para_b)), static_cast<double>(static_cast<unsigned>(para_c)), static_cast<double>(static_cast<unsigned>(para_d)), para_e}, zeta, omega});
                }
                else {
                    suanpan_error("A valid type is required.\n");
                    return SUANPAN_SUCCESS;
                }
            }

            if(is_equal(integrator_type, "LeeNewmarkFullTrial")) { if(domain->insert(make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::TRIAL))) code = 1; }
            else if(is_equal(integrator_type, "LeeNewmarkFullCurrent") || is_equal(integrator_type, "LeeNewmarkFull")) { if(domain->insert(make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::CURRENT))) code = 1; }
            else if(is_equal(integrator_type, "LeeNewmarkFullInitial")) { if(domain->insert(make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::INITIAL))) code = 1; }
        }
        else if(is_equal(integrator_type, "WilsonPenzienNewmark")) {
            vector<double> damping_coef;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid damping coefficient is required.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
            }

            if(domain->insert(make_shared<WilsonPenzienNewmark>(tag, damping_coef, alpha, beta))) code = 1;
        }
        else if(is_equal(integrator_type, "NonviscousNewmark")) {
            vector<double> m_r, s_r, m_i, s_i;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid damping coefficient is required.\n");
                    return SUANPAN_SUCCESS;
                }
                m_r.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid damping coefficient is required.\n");
                    return SUANPAN_SUCCESS;
                }
                m_i.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid damping coefficient is required.\n");
                    return SUANPAN_SUCCESS;
                }
                s_r.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("A valid damping coefficient is required.\n");
                    return SUANPAN_SUCCESS;
                }
                s_i.emplace_back(t_para);
            }

            auto m_imag = vec{m_i}, s_imag = vec{s_i};
            if(accu(m_imag) + accu(s_imag) > 1E-10) {
                suanpan_error("Parameters should be conjugate pairs.\n");
                return SUANPAN_SUCCESS;
            }

            if(domain->insert(make_shared<NonviscousNewmark>(tag, alpha, beta, cx_vec{vec{m_r}, m_imag}, cx_vec{vec{s_r}, s_imag}))) code = 1;
        }
    }
    else if(is_equal(integrator_type, "GeneralizedAlpha") || is_equal(integrator_type, "GeneralisedAlpha")) {
        vector<double> pool;
        pool.reserve(2);

        double para;
        while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

        if(pool.empty() && domain->insert(make_shared<GeneralizedAlpha>(tag, .5))) code = 1; // NOLINT(bugprone-branch-clone)
        else if(1 == pool.size() && domain->insert(make_shared<GeneralizedAlpha>(tag, std::min(std::max(0., pool[0]), 1.)))) code = 1;
        else if(2 == pool.size() && domain->insert(make_shared<GeneralizedAlpha>(tag, pool[0], pool[1]))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSU0")) {
        vec pool(3);

        for(auto& I : pool)
            if(!get_input(command, I)) {
                suanpan_error("A valid damping radius is required.\n");
                return SUANPAN_SUCCESS;
            }

        if(domain->insert(make_shared<GSSSSU0>(tag, std::move(pool)))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSV0")) {
        vec pool(3);

        for(auto& I : pool)
            if(!get_input(command, I)) {
                suanpan_error("A valid damping radius is required.\n");
                return SUANPAN_SUCCESS;
            }

        if(domain->insert(make_shared<GSSSSV0>(tag, std::move(pool)))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSOptimal")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<GSSSSOptimal>(tag, std::max(0., std::min(radius, 1.))))) code = 1;
    }
    else if(is_equal(integrator_type, "OALTS")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<OALTS>(tag, std::max(0., std::min(radius, 1.))))) code = 1;
    }
    else if(is_equal(integrator_type, "BatheTwoStep")) {
        auto radius = 0.;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }
        radius = std::max(0., std::min(radius, 1.));

        auto gamma = .5;
        if(!get_optional_input(command, gamma)) {
            suanpan_error("A valid gamma is required.\n");
            return SUANPAN_SUCCESS;
        }
        if(gamma <= 0. || gamma >= 1.) gamma = .5;

        if(domain->insert(make_shared<BatheTwoStep>(tag, radius, gamma))) code = 1;
    }
    else if(is_equal(integrator_type, "Tchamwa")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<Tchamwa>(tag, std::max(0., std::min(radius, 1.))))) code = 1;
    }
    else if(is_equal(integrator_type, "BatheExplicit")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<BatheExplicit>(tag, std::max(0., std::min(radius, 1.))))) code = 1;
    }
    else if(is_equal(integrator_type, "GeneralizedAlphaExplicit") || is_equal(integrator_type, "GeneralisedAlphaExplicit")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<GeneralizedAlphaExplicit>(tag, std::max(0., std::min(radius, 1.))))) code = 1;
    }

    if(1 == code) {
        if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_integrator_tag(tag);
        domain->set_current_integrator_tag(tag);
    }
    else
        suanpan_error("Fail to create new integrator via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_solver(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string solver_type;
    if(!get_input(command, solver_type)) {
        suanpan_error("A valid solver type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    if(is_equal(solver_type, "Newton")) { if(domain->insert(make_shared<Newton>(tag))) code = 1; }
    else if(is_equal(solver_type, "modifiedNewton") || is_equal(solver_type, "mNewton")) { if(domain->insert(make_shared<Newton>(tag, true))) code = 1; }
    else if(is_equal(solver_type, "AICN")) {
        auto length = 1.;
        if(!command.eof() && !get_input(command, length)) {
            suanpan_error("A valid length is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<AICN>(tag, length))) code = 1;
    }
    else if(is_equal(solver_type, "BFGS")) { if(domain->insert(make_shared<BFGS>(tag))) code = 1; }
    else if(is_equal(solver_type, "LBFGS")) {
        auto max_history = 20;
        if(!command.eof() && !get_input(command, max_history)) {
            suanpan_error("A valid maximum step number is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<BFGS>(tag, max_history))) code = 1;
    }
    else if(is_equal(solver_type, "FEAST") || is_equal(solver_type, "QuadraticFEAST")) {
        unsigned eigen_number;
        if(!get_input(command, eigen_number)) {
            suanpan_error("A valid number of frequencies is required.\n");
            return SUANPAN_SUCCESS;
        }

        double centre;
        if(!get_input(command, centre)) {
            suanpan_error("A valid radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        auto radius = centre;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<FEAST>(tag, eigen_number, centre, radius, is_equal(solver_type, "QuadraticFEAST")))) code = 1;
    }
    else if(is_equal(solver_type, "DisplacementControl") || is_equal(solver_type, "MPDC")) { if(domain->insert(make_shared<MPDC>(tag))) code = 1; }
    else if(is_equal(solver_type, "Ramm")) { if(domain->insert(make_shared<Ramm>(tag))) code = 1; }
    else
        suanpan_error("Cannot identify the solver type.\n");

    if(1 == code) {
        if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_solver_tag(tag);
        domain->set_current_solver_tag(tag);
    }
    else
        suanpan_error("Fail to create new solver via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
