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

#include "SolverParser.h"

#include <Domain/DomainBase.h>
#include <Solver/Solver>
#include <Step/Step.h>
#include <Toolbox/utility.h>

int create_new_integrator(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string integrator_type;
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
    if(if_contain(suanpan::to_upper(cref(integrator_type)), suanpan::to_upper(std::string("Newmark")))) {
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

        if(is_equal(integrator_type, "Newmark")) {
            if(domain->insert(std::make_shared<Newmark>(tag, alpha, beta))) code = 1;
        }
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

            if(domain->insert(std::make_shared<RayleighNewmark>(tag, alpha, beta, p(0), p(1), p(2), p(3)))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeNewmark")) {
            const auto [damping_coef, frequency] = get_remaining<double, double>(command);

            if(domain->insert(std::make_shared<LeeNewmark>(tag, damping_coef, frequency, alpha, beta))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeNewmarkIterative")) {
            using LeeMode = LeeNewmarkIterative::Mode;
            using LeeType = LeeNewmarkIterative::Type;

            std::vector<LeeMode> modes;

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

            auto regularise = [](const double x) { return static_cast<double>(static_cast<unsigned>(x)); };

            while(!command.eof()) {
                std::string type;
                if(!get_input(command, type)) {
                    suanpan_error("A valid type is required.\n");
                    return SUANPAN_SUCCESS;
                }
                if(is_equal("-type0", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeMode{LeeType::T0, vec{}, zeta, omega});
                }
                else if(is_equal("-type1", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeMode{LeeType::T1, vec{regularise(para_a)}, zeta, omega});
                }
                else if(is_equal("-type2", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeMode{LeeType::T2, vec{regularise(para_a), regularise(para_b)}, zeta, omega});
                }
                else if(is_equal("-type3", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeMode{LeeType::T3, vec{para_a}, zeta, omega});
                }
                else if(is_equal("-type4", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return SUANPAN_SUCCESS;
                    double para_c, para_d, para_e;
                    if(!get_input(command, para_c) || !get_input(command, para_d) || !get_input(command, para_e)) {
                        suanpan_error("A valid parameter is required.\n");
                        return SUANPAN_SUCCESS;
                    }
                    modes.emplace_back(LeeMode{LeeType::T4, vec{regularise(para_a), regularise(para_b), regularise(para_c), regularise(para_d), para_e}, zeta, omega});
                }
                else {
                    suanpan_error("A valid type is required.\n");
                    return SUANPAN_SUCCESS;
                }
            }

            if(domain->insert(std::make_shared<LeeNewmarkIterative>(tag, std::move(modes), alpha, beta))) code = 1;
        }
        else if(is_equal(integrator_type, "LeeElementalNewmark")) {
            const auto [damping_coef, frequency] = get_remaining<double, double>(command);

            if(domain->insert(std::make_shared<LeeElementalNewmark>(tag, damping_coef, frequency, alpha, beta))) code = 1;
        }
        else if(integrator_type.size() >= 14 && is_equal(integrator_type.substr(0, 14), "LeeNewmarkFull")) {
            using LeeMode = LeeNewmarkFull::Mode;
            using LeeType = LeeNewmarkFull::Type;

            std::vector<LeeMode> modes;

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

            auto regularise = [](const double x) { return static_cast<double>(static_cast<unsigned>(x)); };

            while(!command.eof()) {
                std::string type;
                if(!get_input(command, type)) {
                    suanpan_error("A valid type is required.\n");
                    return SUANPAN_SUCCESS;
                }
                if(is_equal("-type0", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeMode{LeeType::T0, vec{}, zeta, omega});
                }
                else if(is_equal("-type1", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeMode{LeeType::T1, vec{regularise(para_a)}, zeta, omega});
                }
                else if(is_equal("-type2", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeMode{LeeType::T2, vec{regularise(para_a), regularise(para_b)}, zeta, omega});
                }
                else if(is_equal("-type3", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return SUANPAN_SUCCESS;
                    modes.emplace_back(LeeMode{LeeType::T3, vec{para_a}, zeta, omega});
                }
                else if(is_equal("-type4", type)) {
                    if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return SUANPAN_SUCCESS;
                    double para_c, para_d, para_e;
                    if(!get_input(command, para_c) || !get_input(command, para_d) || !get_input(command, para_e)) {
                        suanpan_error("A valid parameter is required.\n");
                        return SUANPAN_SUCCESS;
                    }
                    modes.emplace_back(LeeMode{LeeType::T4, vec{regularise(para_a), regularise(para_b), regularise(para_c), regularise(para_d), para_e}, zeta, omega});
                }
                else {
                    suanpan_error("A valid type is required.\n");
                    return SUANPAN_SUCCESS;
                }
            }

            if(is_equal(integrator_type, "LeeNewmarkFullTrial")) {
                if(domain->insert(std::make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::TRIAL))) code = 1;
            }
            else if(is_equal(integrator_type, "LeeNewmarkFullCurrent") || is_equal(integrator_type, "LeeNewmarkFull")) {
                if(domain->insert(std::make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::CURRENT))) code = 1;
            }
            else if(is_equal(integrator_type, "LeeNewmarkFullInitial")) {
                if(domain->insert(std::make_shared<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, LeeNewmarkBase::StiffnessType::INITIAL))) code = 1;
            }
        }
        else if(is_equal(integrator_type, "WilsonPenzienNewmark")) {
            if(domain->insert(std::make_shared<WilsonPenzienNewmark>(tag, get_remaining<double>(command), alpha, beta))) code = 1;
        }
        else if(is_equal(integrator_type, "NonviscousNewmark")) {
            const auto [m_r, m_i, s_r, s_i] = get_remaining<double, double, double, double>(command);

            auto m_imag = vec{m_i}, s_imag = vec{s_i};
            if(accu(m_imag) + accu(s_imag) > 1E-12) {
                suanpan_error("Parameters should be conjugate pairs.\n");
                return SUANPAN_SUCCESS;
            }

            auto m = cx_vec{vec{m_r}, m_imag}, s = cx_vec{vec{s_r}, s_imag};

            if(std::abs(accu(m % exp(-1E8 * s))) > 1E-12) {
                suanpan_error("The provided kernel does not converge to zero.\n");
                return SUANPAN_SUCCESS;
            }

            if(domain->insert(std::make_shared<NonviscousNewmark>(tag, alpha, beta, std::move(m), std::move(s)))) code = 1;
        }
    }
    else if(is_equal(integrator_type, "GeneralizedAlpha") || is_equal(integrator_type, "GeneralisedAlpha")) {
        const auto pool = get_remaining<double>(command);

        if(pool.empty() && domain->insert(std::make_shared<GeneralizedAlpha>(tag, .5))) code = 1; // NOLINT(bugprone-branch-clone)
        else if(1 == pool.size() && domain->insert(std::make_shared<GeneralizedAlpha>(tag, suanpan::clamp_unit(pool[0])))) code = 1;
        else if(2 == pool.size() && domain->insert(std::make_shared<GeneralizedAlpha>(tag, pool[0], pool[1]))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSU0")) {
        vec pool(3);

        if(!get_input(command, pool)) {
            suanpan_error("A valid spectral radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<GSSSSU0>(tag, std::move(pool)))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSV0")) {
        vec pool(3);

        if(!get_input(command, pool)) {
            suanpan_error("A valid spectral radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<GSSSSV0>(tag, std::move(pool)))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSSSOptimal")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid spectral radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<GSSSSOptimal>(tag, suanpan::clamp_unit(radius)))) code = 1;
    }
    else if(is_equal(integrator_type, "OALTS")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<OALTS>(tag, suanpan::clamp_unit(radius)))) code = 1;
    }
    else if(is_equal(integrator_type, "BatheTwoStep")) {
        auto radius = 0.;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        auto gamma = .5;
        if(!get_optional_input(command, gamma)) {
            suanpan_error("A valid gamma is required.\n");
            return SUANPAN_SUCCESS;
        }
        if(gamma <= 0. || gamma >= 1.) gamma = .5;

        if(domain->insert(std::make_shared<BatheTwoStep>(tag, suanpan::clamp_unit(radius), gamma))) code = 1;
    }
    else if(is_equal(integrator_type, "Tchamwa")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<Tchamwa>(tag, suanpan::clamp_unit(radius)))) code = 1;
    }
    else if(is_equal(integrator_type, "BatheExplicit")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<BatheExplicit>(tag, suanpan::clamp_unit(radius)))) code = 1;
    }
    else if(is_equal(integrator_type, "ICL")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<ICL>(tag, suanpan::clamp(radius, .5, 1.)))) code = 1;
    }
    else if(is_equal(integrator_type, "GSSE")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<GSSE>(tag, suanpan::clamp_unit(radius)))) code = 1;
    }
    else if(is_equal(integrator_type, "WAT2")) {
        auto para = 1. / 3.;
        if(!get_optional_input(command, para)) {
            suanpan_error("A valid parameter is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<WAT2>(tag, suanpan::clamp_unit(para)))) code = 1;
    }
    else if(is_equal(integrator_type, "GeneralizedAlphaExplicit") || is_equal(integrator_type, "GeneralisedAlphaExplicit")) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<GeneralizedAlphaExplicit>(tag, suanpan::clamp_unit(radius)))) code = 1;
    }

    if(1 == code) {
        if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_integrator_tag(tag);
        domain->set_current_integrator_tag(tag);
    }
    else
        suanpan_error("Fail to create new integrator via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_solver(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string solver_type;
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
    if(is_equal(solver_type, "Newton")) {
        if(domain->insert(std::make_shared<Newton>(tag))) code = 1;
    }
    else if(is_equal(solver_type, "modifiedNewton") || is_equal(solver_type, "mNewton")) {
        if(domain->insert(std::make_shared<Newton>(tag, true))) code = 1;
    }
    else if(is_equal(solver_type, "AICN")) {
        auto length = 1.;
        if(!command.eof() && !get_input(command, length)) {
            suanpan_error("A valid length is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<AICN>(tag, length))) code = 1;
    }
    else if(is_equal(solver_type, "BFGS")) {
        if(domain->insert(std::make_shared<BFGS>(tag))) code = 1;
    }
    else if(is_equal(solver_type, "LBFGS")) {
        auto max_history = 20;
        if(!command.eof() && !get_input(command, max_history)) {
            suanpan_error("A valid maximum step number is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(domain->insert(std::make_shared<BFGS>(tag, max_history))) code = 1;
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

        if(domain->insert(std::make_shared<FEAST>(tag, eigen_number, centre, radius, is_equal(solver_type, "QuadraticFEAST")))) code = 1;
    }
    else if(is_equal(solver_type, "DisplacementControl") || is_equal(solver_type, "MPDC")) {
        if(domain->insert(std::make_shared<MPDC>(tag))) code = 1;
    }
    else if(is_equal(solver_type, "Ramm")) {
        if(domain->insert(std::make_shared<Ramm>(tag))) code = 1;
    }
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
