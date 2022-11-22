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

#include "SolverParser.h"
#include <Domain/DomainBase.h>
#include <Solver/Solver>
#include <Step/Step.h>
#include <Toolbox/utility.h>

int create_new_integrator(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string integrator_type;
    if(!get_input(command, integrator_type)) {
        suanpan_error("create_new_integrator() needs a valid integrator type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_integrator() needs a valid tag.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    if(if_contain(suanpan::to_upper(cref(integrator_type)), suanpan::to_upper(string("Newmark")))) {
        auto alpha = .25, beta = .5;
        if(!command.eof()) {
            if(!get_input(command, alpha)) {
                suanpan_error("create_new_integrator() needs a valid alpha.\n");
                return SUANPAN_SUCCESS;
            }
            if(!get_input(command, beta)) {
                suanpan_error("create_new_integrator() needs a valid beta.\n");
                return SUANPAN_SUCCESS;
            }
        }

        if(is_equal(integrator_type, "Newmark")) { if(domain->insert(make_shared<Newmark>(tag, alpha, beta))) code = 1; }
        else if(is_equal(integrator_type, "RayleighNewmark")) {
            vec p(4, fill::zeros);
            auto idx = 0llu;
            while(!command.eof() && idx < p.n_elem)
                if(!get_input(command, p(idx++))) {
                    suanpan_error("create_new_integrator() needs a valid parameter for Rayleigh damping.\n");
                    return SUANPAN_SUCCESS;
                }

            if(domain->insert(make_shared<RayleighNewmark>(tag, alpha, beta, p(0), p(1), p(2), p(3)))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeNewmark")) {
            vector<double> damping_coef, frequency;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("create_new_integrator() needs a valid damping coefficient.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("create_new_integrator() needs a valid frequency.\n");
                    return SUANPAN_SUCCESS;
                }
                frequency.emplace_back(t_para);
            }

            if(domain->insert(make_shared<LeeNewmark>(tag, damping_coef, frequency, alpha, beta))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeElementalNewmark")) {
            vector<double> damping_coef, frequency;

            while(!command.eof()) {
                double t_para;
                if(!get_input(command, t_para)) {
                    suanpan_error("create_new_integrator() needs a valid damping coefficient.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
                if(!get_input(command, t_para)) {
                    suanpan_error("create_new_integrator() needs a valid frequency.\n");
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
                    suanpan_error("create_new_integrator() needs a valid zeta_p.\n");
                    return SUANPAN_FAIL;
                }
                if(!get_input(command, omega)) {
                    suanpan_error("create_new_integrator() needs a valid omega_p.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            auto get_first = [&] {
                if(!get_input(command, para_a)) {
                    suanpan_error("create_new_integrator() needs a valid parameter.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            auto get_second = [&] {
                if(!get_input(command, para_b)) {
                    suanpan_error("create_new_integrator() needs a valid parameter.\n");
                    return SUANPAN_FAIL;
                }
                return SUANPAN_SUCCESS;
            };

            while(!command.eof()) {
                string type;
                if(!get_input(command, type)) {
                    suanpan_error("create_new_integrator() needs a valid type.\n");
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
                        suanpan_error("create_new_integrator() needs a valid parameter.\n");
                        return SUANPAN_SUCCESS;
                    }
                    modes.emplace_back(LeeNewmarkFull::Mode{LeeNewmarkFull::Type::T4, vec{static_cast<double>(static_cast<unsigned>(para_a)), static_cast<double>(static_cast<unsigned>(para_b)), static_cast<double>(static_cast<unsigned>(para_c)), static_cast<double>(static_cast<unsigned>(para_d)), para_e}, zeta, omega});
                }
                else {
                    suanpan_error("create_new_integrator() needs a valid type.\n");
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
                    suanpan_error("create_new_integrator() needs a valid damping coefficient.\n");
                    return SUANPAN_SUCCESS;
                }
                damping_coef.emplace_back(t_para);
            }

            if(domain->insert(make_shared<WilsonPenzienNewmark>(tag, damping_coef, alpha, beta))) code = 1;
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
                suanpan_error("create_new_integrator() needs a valid damping radius.\n");
                return SUANPAN_SUCCESS;
            }

        if(domain->insert(make_shared<GSSSSU0>(tag, std::move(pool)))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSV0")) {
        vec pool(3);

        for(auto& I : pool)
            if(!get_input(command, I)) {
                suanpan_error("create_new_integrator() needs a valid damping radius.\n");
                return SUANPAN_SUCCESS;
            }

        if(domain->insert(make_shared<GSSSSV0>(tag, std::move(pool)))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSOptimal")) {
        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("create_new_integrator() needs a valid damping radius.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<GSSSSOptimal>(tag, radius))) code = 1;
    }
    else if(is_equal(integrator_type, "BatheTwoStep")) {
        double radius = 0.;
        if(!get_optional_input(command, radius)) {
            suanpan_error("create_new_integrator() needs a valid damping radius.\n");
            return SUANPAN_SUCCESS;
        }
        radius = std::max(0., std::min(radius, 1.));

        double gamma = .5;
        if(!get_optional_input(command, gamma)) {
            suanpan_error("create_new_integrator() needs a valid gamma.\n");
            return SUANPAN_SUCCESS;
        }
        if(suanpan::approx_equal(gamma, 0.) || suanpan::approx_equal(gamma, 1.) || suanpan::approx_equal(gamma, 2. / (1. - radius))) {
            suanpan_warning("BatheTwoStep resets gamma to 0.5.\n");
            gamma = .5;
        }

        if(domain->insert(make_shared<BatheTwoStep>(tag, radius, gamma))) code = 1;
    }

    if(1 == code) {
        if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_integrator_tag(tag);
        domain->set_current_integrator_tag(tag);
    }
    else suanpan_error("create_new_integrator() fails to create the new integrator.\n");

    return SUANPAN_SUCCESS;
}

int create_new_solver(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string solver_type;
    if(!get_input(command, solver_type)) {
        suanpan_error("create_new_solver() requires solver type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_solver() requires a tag.\n");
        return SUANPAN_SUCCESS;
    }

    auto code = 0;
    if(is_equal(solver_type, "Newton")) { if(domain->insert(make_shared<Newton>(tag))) code = 1; }
    else if(is_equal(solver_type, "modifiedNewton") || is_equal(solver_type, "mNewton")) { if(domain->insert(make_shared<Newton>(tag, true))) code = 1; }
    else if(is_equal(solver_type, "BFGS")) { if(domain->insert(make_shared<BFGS>(tag))) code = 1; }
    else if(is_equal(solver_type, "LBFGS")) {
        auto max_history = 20;
        if(!command.eof() && !get_input(command, max_history)) {
            suanpan_error("create_new_solver() requires a valid maximum step for LBFGS algorithm.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<BFGS>(tag, max_history))) code = 1;
    }
    else if(is_equal(solver_type, "Ramm")) {
        auto arc_length = .1;
        string fixed_arc_length = "False";

        if(!command.eof() && !get_input(command, arc_length)) {
            suanpan_error("create_new_solver() requires a valid arc length.\n");
            return SUANPAN_SUCCESS;
        }
        if(!command.eof() && !get_input(command, fixed_arc_length)) {
            suanpan_error("create_new_solver() requires a valid arc length switch.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<Ramm>(tag, arc_length, is_true(fixed_arc_length)))) code = 1;
    }
    else if(is_equal(solver_type, "FEAST")) {
        unsigned eigen_number;
        if(!get_input(command, eigen_number)) {
            suanpan_error("create_new_solver() requires a valid number of frequencies.\n");
            return SUANPAN_SUCCESS;
        }

        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("create_new_solver() requires a valid radius.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<FEAST>(tag, eigen_number, radius, false))) code = 1;
    }
    else if(is_equal(solver_type, "QuadraticFEAST")) {
        unsigned eigen_number;
        if(!get_input(command, eigen_number)) {
            suanpan_error("create_new_solver() requires a valid number of frequencies.\n");
            return SUANPAN_SUCCESS;
        }

        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("create_new_solver() requires a valid radius.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(make_shared<FEAST>(tag, eigen_number, radius, true))) code = 1;
    }
    else if(is_equal(solver_type, "DisplacementControl") || is_equal(solver_type, "MPDC")) { if(domain->insert(make_shared<MPDC>(tag))) code = 1; }
    else suanpan_error("create_new_solver() cannot identify solver type.\n");

    if(1 == code) {
        if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_solver_tag(tag);
        domain->set_current_solver_tag(tag);
    }
    else suanpan_error("create_new_solver() cannot create the new solver.\n");

    return SUANPAN_SUCCESS;
}
