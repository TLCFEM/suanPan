/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
#include <Toolbox/TokenRegistry.h>

namespace {
    //
    // integrators
    //

    bool get_newmark_parameters(std::istringstream& command, double& alpha, double& beta) {
        alpha = .25;
        beta = .5;

        if(command.eof()) return true;

        if(!get_input(command, alpha)) {
            suanpan_error("A valid alpha is required.\n");
            return false;
        }
        if(!get_input(command, beta)) {
            suanpan_error("A valid beta is required.\n");
            return false;
        }

        return true;
    }

    void new_newmark(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

        return_obj = std::make_unique<Newmark>(tag, alpha, beta);
    }

    void new_rayleigh_newmark(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

        vec p(4, fill::zeros);
        uword idx{0};
        while(!command.eof() && idx < p.n_elem)
            if(!get_input(command, p(idx++))) {
                suanpan_error("A valid parameter for Rayleigh damping is required.\n");
                return;
            }

        if(suanpan::approx_equal(sum(p), 0.)) suanpan_warning("It seems all parameters are zeros, either the inputs are wrong or consider use plain Newmark.\n");

        return_obj = std::make_unique<RayleighNewmark>(tag, alpha, beta, p(0), p(1), p(2), p(3));
    }

    void new_lee_newmark(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

        const auto [damping_coef, frequency] = get_remaining<double, double>(command);

        return_obj = std::make_unique<LeeNewmark>(tag, damping_coef, frequency, alpha, beta);
    }

    void new_lee_newmark_iterative(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

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
                return;
            }
            if(is_equal(type, "-type0")) {
                if(SUANPAN_SUCCESS != get_basic_input()) return;
                modes.emplace_back(LeeMode{.t = LeeType::T0, .p = vec{}, .zeta = zeta, .omega = omega});
            }
            else if(is_equal(type, "-type1")) {
                if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return;
                modes.emplace_back(LeeMode{.t = LeeType::T1, .p = vec{regularise(para_a)}, .zeta = zeta, .omega = omega});
            }
            else if(is_equal(type, "-type2")) {
                if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return;
                modes.emplace_back(LeeMode{.t = LeeType::T2, .p = vec{regularise(para_a), regularise(para_b)}, .zeta = zeta, .omega = omega});
            }
            else if(is_equal(type, "-type3")) {
                if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return;
                modes.emplace_back(LeeMode{.t = LeeType::T3, .p = vec{para_a}, .zeta = zeta, .omega = omega});
            }
            else if(is_equal(type, "-type4")) {
                if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return;
                double para_c, para_d, para_e;
                if(!get_input(command, para_c) || !get_input(command, para_d) || !get_input(command, para_e)) {
                    suanpan_error("A valid parameter is required.\n");
                    return;
                }
                modes.emplace_back(LeeMode{.t = LeeType::T4, .p = vec{regularise(para_a), regularise(para_b), regularise(para_c), regularise(para_d), para_e}, .zeta = zeta, .omega = omega});
            }
            else {
                suanpan_error("A valid type is required.\n");
                return;
            }
        }

        return_obj = std::make_unique<LeeNewmarkIterative>(tag, std::move(modes), alpha, beta);
    }

    void new_lee_elemental_newmark(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

        const auto [damping_coef, frequency] = get_remaining<double, double>(command);

        return_obj = std::make_unique<LeeElementalNewmark>(tag, damping_coef, frequency, alpha, beta);
    }

    template<LeeNewmarkBase::StiffnessType stiffness_type> void new_lee_newmark_full(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

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
                return;
            }
            if(is_equal(type, "-type0")) {
                if(SUANPAN_SUCCESS != get_basic_input()) return;
                modes.emplace_back(LeeMode{.t = LeeType::T0, .p = vec{}, .zeta = zeta, .omega = omega});
            }
            else if(is_equal(type, "-type1")) {
                if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return;
                modes.emplace_back(LeeMode{.t = LeeType::T1, .p = vec{regularise(para_a)}, .zeta = zeta, .omega = omega});
            }
            else if(is_equal(type, "-type2")) {
                if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return;
                modes.emplace_back(LeeMode{.t = LeeType::T2, .p = vec{regularise(para_a), regularise(para_b)}, .zeta = zeta, .omega = omega});
            }
            else if(is_equal(type, "-type3")) {
                if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first()) return;
                modes.emplace_back(LeeMode{.t = LeeType::T3, .p = vec{para_a}, .zeta = zeta, .omega = omega});
            }
            else if(is_equal(type, "-type4")) {
                if(SUANPAN_SUCCESS != get_basic_input() || SUANPAN_SUCCESS != get_first() || SUANPAN_SUCCESS != get_second()) return;
                double para_c, para_d, para_e;
                if(!get_input(command, para_c) || !get_input(command, para_d) || !get_input(command, para_e)) {
                    suanpan_error("A valid parameter is required.\n");
                    return;
                }
                modes.emplace_back(LeeMode{.t = LeeType::T4, .p = vec{regularise(para_a), regularise(para_b), regularise(para_c), regularise(para_d), para_e}, .zeta = zeta, .omega = omega});
            }
            else {
                suanpan_error("A valid type is required.\n");
                return;
            }
        }

        return_obj = std::make_unique<LeeNewmarkFull>(tag, std::move(modes), alpha, beta, stiffness_type);
    }

    void new_wilson_penzien_newmark(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

        return_obj = std::make_unique<WilsonPenzienNewmark>(tag, get_remaining<double>(command), alpha, beta);
    }

    void new_nonviscous_newmark(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

        const auto [m_r, m_i, s_r, s_i] = get_remaining<double, double, double, double>(command);

        const vec m_imag{m_i}, s_imag{s_i};
        if(accu(m_imag) + accu(s_imag) > 1E-12) {
            suanpan_error("Parameters should be conjugate pairs.\n");
            return;
        }

        auto m = cx_vec{vec{m_r}, m_imag}, s = cx_vec{vec{s_r}, s_imag};

        if(std::abs(accu(m % exp(-1E8 * s))) > 1E-12) {
            suanpan_error("The provided kernel does not converge to zero.\n");
            return;
        }

        return_obj = std::make_unique<NonviscousNewmark>(tag, alpha, beta, std::move(m), std::move(s));
    }

    void new_udd_newmark(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

        const auto [m_r, m_i, s_r, s_i] = get_remaining<double, double, double, double>(command);

        const vec m_real{m_r}, s_real{s_r}, m_imag{m_i}, s_imag{s_i};
        if(accu(m_imag) + accu(s_imag) > 1E-12) {
            suanpan_error("Parameters should be conjugate pairs.\n");
            return;
        }

        if(any(m_real > 0.) && any(s_real > 0.)) {
            suanpan_error("UDD requires both weights (m_j) and poles (s_j) to have negative real parts.\n");
            return;
        }

        return_obj = std::make_unique<UDDNewmark>(tag, alpha, beta, cx_vec{m_real, m_imag}, cx_vec{s_real, s_imag});
    }

    void new_uda_newmark(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        double alpha, beta;
        if(!get_newmark_parameters(command, alpha, beta)) return;

        const auto [m_r, m_i, s_r, s_i] = get_remaining<double, double, double, double>(command);

        const vec m_real{m_r}, s_real{s_r}, m_imag{m_i}, s_imag{s_i};
        if(accu(m_imag) + accu(s_imag) > 1E-12) {
            suanpan_error("Parameters should be conjugate pairs.\n");
            return;
        }

        if(any(m_real < 0.) && any(s_real > 0.)) {
            suanpan_error("UDA requires weights (m_j) to have positive real parts and poles (s_j) to have negative real parts.\n");
            return;
        }

        return_obj = std::make_unique<UDANewmark>(tag, alpha, beta, cx_vec{m_real, m_imag}, cx_vec{s_real, s_imag});
    }

    void new_generalized_alpha(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        if(const auto pool = get_remaining<double>(command); pool.empty()) return_obj = std::make_unique<GeneralizedAlpha>(tag, .5);
        else if(1 == pool.size()) return_obj = std::make_unique<GeneralizedAlpha>(tag, suanpan::clamp_unit(pool[0]));
        else if(2 == pool.size()) return_obj = std::make_unique<GeneralizedAlpha>(tag, pool[0], pool[1]);
    }

    template<typename T> void new_gssss(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        vec pool(3);
        if(!get_input(command, pool)) {
            suanpan_error("A valid spectral radius is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag, std::move(pool));
    }

    template<typename T> void new_radius_integrator(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag, suanpan::clamp_unit(radius));
    }

    void new_bathe_two_step(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        auto radius = 0.;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return;
        }

        auto gamma = .5;
        if(!get_optional_input(command, gamma)) {
            suanpan_error("A valid gamma is required.\n");
            return;
        }
        if(gamma <= 0. || gamma >= 1.) gamma = .5;

        return_obj = std::make_unique<BatheTwoStep>(tag, suanpan::clamp_unit(radius), gamma);
    }

    void new_icl(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        auto radius = .5;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid damping radius is required.\n");
            return;
        }

        return_obj = std::make_unique<ICL>(tag, suanpan::clamp(radius, .5, 1.));
    }

    void new_wat2(unique_ptr<Integrator>& return_obj, std::istringstream& command, const unsigned tag) {
        auto para = 1. / 3.;
        if(!get_optional_input(command, para)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<WAT2>(tag, suanpan::clamp_unit(para));
    }

    //
    // solvers
    //

    void new_newton(unique_ptr<Solver>& return_obj, std::istringstream&, const unsigned tag) {
        return_obj = std::make_unique<Newton>(tag);
    }

    void new_modified_newton(unique_ptr<Solver>& return_obj, std::istringstream&, const unsigned tag) {
        return_obj = std::make_unique<Newton>(tag, true);
    }

    void new_aicn(unique_ptr<Solver>& return_obj, std::istringstream& command, const unsigned tag) {
        auto length = 1.;
        if(!command.eof() && !get_input(command, length)) {
            suanpan_error("A valid length is required.\n");
            return;
        }

        return_obj = std::make_unique<AICN>(tag, length);
    }

    void new_bfgs(unique_ptr<Solver>& return_obj, std::istringstream&, const unsigned tag) {
        return_obj = std::make_unique<BFGS>(tag);
    }

    void new_lbfgs(unique_ptr<Solver>& return_obj, std::istringstream& command, const unsigned tag) {
        auto max_history = 20;
        if(!command.eof() && !get_input(command, max_history)) {
            suanpan_error("A valid maximum step number is required.\n");
            return;
        }

        return_obj = std::make_unique<BFGS>(tag, max_history);
    }

    void new_feast(unique_ptr<Solver>& return_obj, std::istringstream& command, const unsigned tag) {
        unsigned eigen_number;
        if(!get_input(command, eigen_number)) {
            suanpan_error("A valid number of frequencies is required.\n");
            return;
        }

        double centre;
        if(!get_input(command, centre)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        auto radius = centre;
        if(!get_optional_input(command, radius)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        return_obj = std::make_unique<FEAST>(tag, eigen_number, centre, radius);
    }

    void new_mpdc(unique_ptr<Solver>& return_obj, std::istringstream&, const unsigned tag) {
        return_obj = std::make_unique<MPDC>(tag);
    }

    void new_ramm(unique_ptr<Solver>& return_obj, std::istringstream&, const unsigned tag) {
        return_obj = std::make_unique<Ramm>(tag);
    }
} // namespace

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

    unique_ptr<Integrator> new_integrator = nullptr;

    static const TokenRegistry<void (*)(unique_ptr<Integrator>&, std::istringstream&, unsigned)> registry{
        {"Newmark", new_newmark},
        {"RayleighNewmark", new_rayleigh_newmark},
        {"LeeNewmark", new_lee_newmark},
        {"LeeNewmarkIterative", new_lee_newmark_iterative},
        {"LeeElementalNewmark", new_lee_elemental_newmark},
        {"LeeNewmarkFullTrial", new_lee_newmark_full<LeeNewmarkBase::StiffnessType::TRIAL>},
        {"LeeNewmarkFullCurrent", new_lee_newmark_full<LeeNewmarkBase::StiffnessType::CURRENT>},
        {"LeeNewmarkFull", new_lee_newmark_full<LeeNewmarkBase::StiffnessType::CURRENT>},
        {"LeeNewmarkFullInitial", new_lee_newmark_full<LeeNewmarkBase::StiffnessType::INITIAL>},
        {"WilsonPenzienNewmark", new_wilson_penzien_newmark},
        {"NonviscousNewmark", new_nonviscous_newmark},
        {"UDDNewmark", new_udd_newmark},
        {"UDANewmark", new_uda_newmark},
        {"GeneralizedAlpha", new_generalized_alpha},
        {"GeneralisedAlpha", new_generalized_alpha},
        {"GSSSSU0", new_gssss<GSSSSU0>},
        {"GSSSSV0", new_gssss<GSSSSV0>},
        {"GSSSSOptimal", new_radius_integrator<GSSSSOptimal>},
        {"OALTS", new_radius_integrator<OALTS>},
        {"BatheTwoStep", new_bathe_two_step},
        {"Tchamwa", new_radius_integrator<Tchamwa>},
        {"BatheExplicit", new_radius_integrator<BatheExplicit>},
        {"ICL", new_icl},
        {"GSSE", new_radius_integrator<GSSE>},
        {"WAT2", new_wat2},
        {"GeneralizedAlphaExplicit", new_radius_integrator<GeneralizedAlphaExplicit>},
        {"GeneralisedAlphaExplicit", new_radius_integrator<GeneralizedAlphaExplicit>},
    };

    if(const auto factory = registry.fetch(integrator_type); nullptr != factory) factory(new_integrator, command, tag);

    if(nullptr == new_integrator || !domain->insert(std::move(new_integrator))) {
        suanpan_error("Fail to create new integrator via \"{}\".\n", command.str());
        return SUANPAN_SUCCESS;
    }

    if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_integrator_tag(tag);
    domain->set_current_integrator_tag(tag);

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

    unique_ptr<Solver> new_solver = nullptr;

    static const TokenRegistry<void (*)(unique_ptr<Solver>&, std::istringstream&, unsigned)> registry{
        {"Newton", new_newton},
        {"modifiedNewton", new_modified_newton},
        {"mNewton", new_modified_newton},
        {"AICN", new_aicn},
        {"BFGS", new_bfgs},
        {"LBFGS", new_lbfgs},
        {"FEAST", new_feast},
        {"DisplacementControl", new_mpdc},
        {"MPDC", new_mpdc},
        {"Ramm", new_ramm},
    };

    if(const auto factory = registry.fetch(solver_type); nullptr != factory) factory(new_solver, command, tag);

    if(nullptr == new_solver || !domain->insert(std::move(new_solver))) {
        suanpan_error("Fail to create new solver via \"{}\".\n", command.str());
        return SUANPAN_SUCCESS;
    }

    if(0 != domain->get_current_step_tag()) domain->get_current_step()->set_solver_tag(tag);
    domain->set_current_solver_tag(tag);

    return SUANPAN_SUCCESS;
}
