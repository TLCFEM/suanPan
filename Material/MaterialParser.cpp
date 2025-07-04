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

// ReSharper disable CppClangTidyCppcoreguidelinesInitVariables
// ReSharper disable IdentifierTypo
// ReSharper disable StringLiteralTypo
#include "MaterialParser.h"

#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Material/Material>
#include <Toolbox/utility.h>

namespace {
    void new_afc01(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double t_yield_stress, t_hardening, t_unloading;
        double c_yield_stress, c_hardening, c_unloading;
        if(!get_input(command, t_yield_stress)) {
            suanpan_error("A valid tension yield stress is required.\n");
            return;
        }
        if(!get_input(command, t_hardening)) {
            suanpan_error("A valid tension hardening modulus is required.\n");
            return;
        }
        if(!get_input(command, t_unloading)) {
            suanpan_error("A valid tension unloading modulus is required.\n");
            return;
        }
        if(!get_input(command, c_yield_stress)) {
            suanpan_error("A valid compression yield stress is required.\n");
            return;
        }
        if(!get_input(command, c_hardening)) {
            suanpan_error("A valid compression hardening modulus is required.\n");
            return;
        }
        if(!get_input(command, c_unloading)) {
            suanpan_error("A valid compression unloading modulus is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<AFC>(tag, elastic_modulus, t_yield_stress, t_hardening, t_unloading, c_yield_stress, c_hardening, c_unloading, 0., density);
    }

    void new_afc02(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double t_yield_stress, t_hardening, t_unloading;
        if(!get_input(command, t_yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }
        if(!get_input(command, t_hardening)) {
            suanpan_error("A valid hardening modulus is required.\n");
            return;
        }
        if(!get_input(command, t_unloading)) {
            suanpan_error("A valid unloading modulus is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<AFC>(tag, elastic_modulus, t_yield_stress, t_hardening, t_unloading, t_yield_stress, t_hardening, t_unloading, 0., density);
    }

    void new_afc03(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double t_yield_stress, t_hardening, t_unloading;
        double c_yield_stress, c_hardening, c_unloading;
        if(!get_input(command, t_yield_stress)) {
            suanpan_error("A valid tension yield stress is required.\n");
            return;
        }
        if(!get_input(command, t_hardening)) {
            suanpan_error("A valid tension hardening modulus is required.\n");
            return;
        }
        if(!get_input(command, t_unloading)) {
            suanpan_error("A valid tension unloading modulus is required.\n");
            return;
        }
        if(!get_input(command, c_yield_stress)) {
            suanpan_error("A valid compression yield stress is required.\n");
            return;
        }
        if(!get_input(command, c_hardening)) {
            suanpan_error("A valid compression hardening modulus is required.\n");
            return;
        }
        if(!get_input(command, c_unloading)) {
            suanpan_error("A valid compression unloading modulus is required.\n");
            return;
        }

        auto degrade = 0.;
        if(command.eof())
            suanpan_debug("Linear degradation assumed.\n");
        else if(!get_input(command, degrade)) {
            suanpan_error("A valid degradation parameter is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<AFC>(tag, elastic_modulus, t_yield_stress, t_hardening, t_unloading, c_yield_stress, c_hardening, c_unloading, degrade, density);
    }

    void new_armstrongfrederick(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec pool{2E5, .2, 4E2, 0., 5E2, 1E1};
        if(!get_optional_input(command, pool)) {
            suanpan_error("Valid inputs are required.\n");
            return;
        }

        const auto all = get_remaining<double>(command);

        auto size = all.size();
        auto density = 0.;
        if(size % 2 == 1) {
            --size;
            density = all.back();
        }

        std::vector<double> ai, bi;
        for(size_t I = 0; I < size;) {
            ai.emplace_back(all.at(I++));
            bi.emplace_back(all.at(I++));
        }

        return_obj = std::make_unique<ArmstrongFrederick>(tag, DataArmstrongFrederick{pool(0), pool(1), pool(2), pool(3), pool(4), pool(5), ai, bi}, density);
    }

    void new_armstrongfrederick1d(unique_ptr<Material>& return_obj, std::istringstream& command, const bool memory = false) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec pa{2E5, 4E2, 0., 5E2, 1E1};
        if(!get_optional_input(command, pa)) {
            suanpan_error("Valid inputs are required.\n");
            return;
        }

        vec pb{.2, 5E2, 1E1};
        if(memory) {
            if(!get_optional_input(command, pb)) {
                suanpan_error("Valid inputs are required.\n");
                return;
            }
        }
        else pb.zeros();

        const auto all = get_remaining<double>(command);

        auto size = all.size();
        auto density = 0.;
        if(size % 2 == 1) {
            --size;
            density = all.back();
        }

        std::vector<double> ai, bi;
        for(size_t I = 0; I < size;) {
            ai.emplace_back(all.at(I++));
            bi.emplace_back(all.at(I++));
        }

        return_obj = std::make_unique<ArmstrongFrederick1D>(tag, DataArmstrongFrederick1D{pa(0), pa(1), pa(2), pa(3), pa(4), pb(0), pb(1), pb(2), ai, bi}, density);
    }

    void new_axisymmetric(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned full_tag;
        if(!get_input(command, full_tag)) {
            suanpan_error("A valid reference material tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Axisymmetric>(tag, full_tag);
    }

    void new_axisymmetricelastic(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<AxisymmetricElastic>(tag, elastic_modulus, poissons_ratio, density);
    }

    void new_bilinear1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto hardening_ratio = 0.;
        if(command.eof())
            suanpan_debug("Zero hardening ratio assumed.\n");
        else if(!get_input(command, hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto beta = 1.;
        if(command.eof())
            suanpan_debug("Isotropic hardening assumed.\n");
        else if(!get_input(command, beta)) {
            suanpan_error("A valid beta is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<Bilinear1D>(tag, elastic_modulus, yield_stress, hardening_ratio, suanpan::clamp_unit(beta), density);
    }

    void new_bilinearcc(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        double beta, m, pt, a, a_slope;
        if(!get_input(command, beta)) {
            suanpan_error("A valid beta is required.\n");
            return;
        }
        if(!get_input(command, m)) {
            suanpan_error("A valid radius ratio is required.\n");
            return;
        }
        if(!get_input(command, pt)) {
            suanpan_error("A valid tensile yield strength is required.\n");
            return;
        }
        if(!get_input(command, a)) {
            suanpan_error("A valid initial size is required.\n");
            return;
        }
        if(!get_input(command, a_slope)) {
            suanpan_error("A valid hardening slope is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearCC>(tag, elastic_modulus, poissons_ratio, beta, m, pt, a, a_slope, density);
    }

    void new_bilineardp(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        double eta_yield, eta_flow, xi, cohesion, cohesion_slope;
        if(!get_input(command, eta_yield)) {
            suanpan_error("A valid eta for yielding criterion is required.\n");
            return;
        }
        if(!get_input(command, eta_flow)) {
            suanpan_error("A valid eta for plasticity flow rule is required.\n");
            return;
        }
        if(!get_input(command, xi)) {
            suanpan_error("A valid xi is required.\n");
            return;
        }
        if(!get_input(command, cohesion)) {
            suanpan_error("A valid cohesion is required.\n");
            return;
        }
        if(!get_input(command, cohesion_slope)) {
            suanpan_error("A valid cohesion is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearDP>(tag, elastic_modulus, poissons_ratio, eta_yield, eta_flow, xi, cohesion, cohesion_slope, density);
    }

    void new_customdp(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        double eta_yield, eta_flow, xi;
        if(!get_input(command, eta_yield)) {
            suanpan_error("A valid eta for yielding criterion is required.\n");
            return;
        }
        if(!get_input(command, eta_flow)) {
            suanpan_error("A valid eta for plasticity flow rule is required.\n");
            return;
        }
        if(!get_input(command, xi)) {
            suanpan_error("A valid xi is required.\n");
            return;
        }

        unsigned expression;
        if(!get_input(command, expression)) {
            suanpan_error("A valid expression tag is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomDP>(tag, elastic_modulus, poissons_ratio, eta_yield, eta_flow, xi, expression, density);
    }

    void new_bilinearelastic1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto hardening_ratio = 0.;
        if(command.eof())
            suanpan_debug("Zero hardening ratio assumed.\n");
        else if(!get_input(command, hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearElastic1D>(tag, elastic_modulus, yield_stress, hardening_ratio, 0., density);
    }

    void new_nle1d01(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        double hardening_ratio;
        if(!get_input(command, hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("A valid radius for transition is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearElastic1D>(tag, elastic_modulus, yield_stress, hardening_ratio, radius, density);
    }

    void new_bilinearorthotropic(unique_ptr<Material>& return_obj, std::istringstream& command, const bool is_hoffman) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec modulus(6);
        if(!get_input(command, modulus)) {
            suanpan_error("A valid modulus is required.\n");
            return;
        }

        vec poissons_ratio(3);
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        vec stress(9);
        if(!get_input(command, stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto hardening = 0.;
        if(command.eof())
            suanpan_debug("Zero hardening ratio assumed.\n");
        else if(!get_input(command, hardening)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        if(is_hoffman) return_obj = std::make_unique<BilinearHoffman>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), hardening, density);
        else return_obj = std::make_unique<BilinearTsaiWu>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), hardening, density);
    }

    void new_bilinearj2(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto hardening_ratio = 0.;
        if(command.eof())
            suanpan_debug("Zero hardening ratio assumed.\n");
        else if(!get_input(command, hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto beta = 1.;
        if(command.eof())
            suanpan_debug("Isotropic hardening assumed.\n");
        else if(!get_input(command, beta)) {
            suanpan_error("A valid beta is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearJ2>(tag, elastic_modulus, poissons_ratio, yield_stress, hardening_ratio, beta, density);
    }

    void new_bilinearmises1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto hardening_ratio = 0.;
        if(command.eof())
            suanpan_debug("Zero hardening ratio assumed.\n");
        else if(!get_input(command, hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearMises1D>(tag, elastic_modulus, yield_stress, hardening_ratio, density);
    }

    void new_bilinearoo(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto pool = get_remaining<double>(command);

        if(3 == pool.size()) {
            pool.insert(pool.end(), pool.begin() + 1, pool.begin() + 3);
            pool.emplace_back(0.);
        }
        else if(4 == pool.size()) pool.insert(pool.end() - 1, pool.begin() + 1, pool.begin() + 3);
        else if(5 == pool.size()) pool.emplace_back(0.);

        if(6 != pool.size()) {
            suanpan_error("3 to 6 parameters are required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearOO>(tag, pool[0], pool[1], pool[2], pool[3], pool[4], pool[5]);
    }

    void new_bilinearpo(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto pool = get_remaining<double>(command);

        if(3 == pool.size()) {
            pool.insert(pool.end(), pool.begin() + 1, pool.begin() + 3);
            pool.emplace_back(0.);
        }
        else if(4 == pool.size()) pool.insert(pool.end() - 1, pool.begin() + 1, pool.begin() + 3);
        else if(5 == pool.size()) pool.emplace_back(0.);

        if(6 != pool.size()) {
            suanpan_error("3 to 6 parameters are required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearPO>(tag, pool[0], pool[1], pool[2], pool[3], pool[4], pool[5]);
    }

    void new_bilinearperic(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        double hardening;
        if(!get_input(command, hardening)) {
            suanpan_error("A valid hardening modulus is required.\n");
            return;
        }

        double mu, epsilon;
        if(!get_input(command, mu)) {
            suanpan_error("A valid mu is required.\n");
            return;
        }
        if(!get_input(command, epsilon)) {
            suanpan_error("A valid epsilon is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearPeric>(tag, elastic_modulus, poissons_ratio, yield_stress, hardening, mu, epsilon, density);
    }

    void new_blatzko(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<BlatzKo>(tag, elastic_modulus, poissons_ratio, density);
    }

    void new_boucwen(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto pool = get_remaining<double>(command);

        if(5 == pool.size()) pool.emplace_back(0.);

        if(6 != pool.size()) {
            suanpan_error("6 or 7 parameters are required.\n");
            return;
        }

        return_obj = std::make_unique<BoucWen>(tag, pool);
    }

    void new_bwbn(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{2E5, 4E2, 1E-2, 5E-1, 1., 1., 0., 1., 0., 1., 0., 0., 0., 0., 0., 1., 0.};

        if(!get_optional_input(command, para_pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<BWBN>(tag, para_pool.head(16), para_pool(16));
    }

    void new_cdp(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{3E4, .2, 3., 30., 5E-4, 5E-2, .2, 2., .5, .65, .2, 1.16, .5, 2400E-12};

        if(!get_optional_input(command, para_pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<CDP>(tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7), para_pool(8), para_pool(9), para_pool(10), para_pool(11), para_pool(12), para_pool(13));
    }

    void new_customcdp(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        Col<unsigned> expressions(2);
        if(!get_input(command, expressions)) {
            suanpan_error("Two valid expression tags are required.\n");
            return;
        }

        vec pool{3E4, .2, 5E-4, 5E-2, .2, 1.16, .5, 2400E-12};
        if(!get_optional_input(command, pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomCDP>(tag, expressions(0), expressions(1), pool(0), pool(1), pool(2), pool(3), pool(4), pool(5), pool(6), pool(7));
    }

    void new_cdpm2(unique_ptr<Material>& return_obj, std::istringstream& command, const CDPM2::DamageType damage_type) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{3E4, .3, 3., 30., .3, .01, .85, .08, .003, 2., 1E-6, 5., 5E-4, 5E-4, 0.};

        if(!get_optional_input(command, para_pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        if(para_pool(4) < 0. || para_pool(4) > 1.) {
            suanpan_error("Initial ratio qh0 must be in [0,1].\n");
            return;
        }
        if(para_pool(5) + para_pool(4) > 1.) {
            suanpan_error("Hardeinng modulus hp must be smaller than 1-qh0.\n");
            return;
        }

        return_obj = std::make_unique<CDPM2>(tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7), para_pool(8), para_pool(9), para_pool(10), para_pool(11), para_pool(12), para_pool(13), damage_type, para_pool(14));
    }

    void new_concrete21(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        const auto para = get_remaining<double>(command);

        if(para.size() == 8) return_obj = std::make_unique<Concrete21>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], 0.);
        else if(para.size() == 9) return_obj = std::make_unique<Concrete21>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8]);
        else
            suanpan_error("Eight or nine double inputs are required.\n");
    }

    void new_concrete22(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        const auto para = get_remaining<double>(command);

        if(para.size() == 10) return_obj = std::make_unique<Concrete22>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8], para[9], 0.);
        else if(para.size() == 11) return_obj = std::make_unique<Concrete22>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8], para[9], para[10]);
        else
            suanpan_error("Ten or eleven double inputs are required.\n");
    }

    void new_concretecm(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double peak_stress;
        if(!get_input(command, peak_stress)) {
            suanpan_error("A valid compression stress is required.\n");
            return;
        }

        double crack_stress;
        if(!get_input(command, crack_stress)) {
            suanpan_error("A valid tension stress is required.\n");
            return;
        }

        double nc, nt;
        if(!get_input(command, nc, nt)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        auto peak_strain = 2E-3;
        if(!command.eof() && !get_input(command, peak_strain)) {
            suanpan_error("A valid tension stress is required.\n");
            return;
        }

        auto crack_strain = 1E-4;
        if(!command.eof() && !get_input(command, crack_strain)) {
            suanpan_error("A valid tension stress is required.\n");
            return;
        }

        std::string linear_trans = "false";
        if(!command.eof() && !get_input(command, linear_trans)) {
            suanpan_error("A valid transition switch is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        if(fabs(peak_strain) <= fabs(peak_stress / elastic_modulus) || fabs(crack_strain) <= fabs(crack_stress / elastic_modulus)) {
            suanpan_error("The secant stiffness at peaks must be smaller than the initial stiffness.\n");
            return;
        }

        return_obj = std::make_unique<ConcreteCM>(tag, elastic_modulus, peak_stress, crack_stress, nc, nt, peak_strain, crack_strain, is_true(linear_trans), density);
    }

    void new_concreteexp(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double f_t, a_t, g_t, f_c, a_c, g_c;
        if(!get_input(command, f_t, a_t, g_t)) {
            suanpan_error("A valid tension parameter is required.\n");
            return;
        }
        if(!get_input(command, f_c, a_c, g_c)) {
            suanpan_error("A valid compression parameter is required.\n");
            return;
        }

        auto middle_point = .2;
        if(!command.eof() && !get_input(command, middle_point)) {
            suanpan_error("A valid middle point is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<ConcreteExp>(tag, elastic_modulus, f_t, a_t, g_t, f_c, a_c, g_c, middle_point, density);
    }

    void new_concretek4(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus, hardening;
        if(!get_input(command, elastic_modulus, hardening)) {
            suanpan_error("A valid modulus is required.\n");
            return;
        }

        vec pool(8);
        if(!get_input(command, pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        auto enable_damage = true, enable_crack_closing = true, objective_damage = false;
        if(!command.eof()) {
            if(!get_input(command, enable_damage, enable_crack_closing, objective_damage)) {
                suanpan_error("A valid flag is required.\n");
                return;
            }
            suanpan_debug("Internal flags are set.\n");
        }

        return_obj = std::make_unique<ConcreteK4>(tag, elastic_modulus, hardening, std::move(pool), density, enable_damage, enable_crack_closing, objective_damage);
    }

    void new_concretetable(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        std::string c_name, t_name;
        if(!get_input(command, t_name, c_name)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        std::error_code code;
        mat t_table, c_table;
        if(!fs::exists(t_name, code) || !t_table.load(t_name, raw_ascii) || t_table.n_cols != 2) {
            suanpan_error("Cannot load \"{}\".\n", t_name);
            return;
        }
        if(!fs::exists(c_name, code) || !c_table.load(c_name, raw_ascii) || c_table.n_cols != 2) {
            suanpan_error("Cannot load \"{}\".\n", c_name);
            return;
        }

        auto m_point = .2;
        if(!command.eof() && !get_input(command, m_point)) {
            suanpan_error("A valid transition switch is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<ConcreteTable>(tag, -abs(c_table), abs(t_table), m_point, density);
    }

    void new_concretetsai(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        const auto para = get_remaining<double>(command);

        if(para.size() == 8) return_obj = std::make_unique<ConcreteTsai>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7]);
        else if(para.size() == 9) return_obj = std::make_unique<ConcreteTsai>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8]);
        else
            suanpan_error("Eight or nine double inputs are required.\n");
    }

    void new_coulombfriction(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double max_friction, factor;
        if(!get_input(command, max_friction, factor)) {
            suanpan_error("Valid parameters are required.\n");
            return;
        }

        return_obj = std::make_unique<CoulombFriction>(tag, max_friction, factor);
    }

    void new_customelastic1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned expression_tag;
        if(!get_input(command, expression_tag)) {
            suanpan_error("A valid expression tag is required.\n");
            return;
        }

        auto density = 0.;
        if(!get_optional_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomElastic1D>(tag, expression_tag, density);
    }

    void new_custommises1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        unsigned k_tag, h_tag;
        if(!get_input(command, k_tag, h_tag)) {
            suanpan_error("A valid expression tag is required.\n");
            return;
        }

        auto density = 0.;
        if(!get_optional_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomMises1D>(tag, elastic_modulus, k_tag, h_tag, density);
    }

    void new_dhakal(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned mat_tag;
        if(!get_input(command, mat_tag)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        double y_strain, parameter;
        if(!get_input(command, y_strain)) {
            suanpan_error("A valid yield strain is required.\n");
            return;
        }
        if(!get_input(command, parameter)) {
            suanpan_error("A valid bar parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<Dhakal>(tag, mat_tag, y_strain, parameter);
    }

    void new_duncanselig(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto pool = vec{14.7, 400. * 14.7, .6, 300 * 14.7, .2, .7, .1, .7, .5};
        if(!get_input(command, pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<DuncanSelig>(tag, pool, density);
    }

    void new_sinh1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<Sinh1D>(tag, elastic_modulus, density);
    }

    void new_tanh1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<Tanh1D>(tag, elastic_modulus, density);
    }

    void new_timberpd(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec modulus(6);
        if(!get_input(command, modulus)) {
            suanpan_error("A valid modulus is required.\n");
            return;
        }

        vec poissons_ratio(3);
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        vec stress(9);
        if(!get_input(command, stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        vec para(7);
        if(!get_input(command, para)) {
            suanpan_error("A valid hardening parameter is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<TimberPD>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), std::move(para), density);
    }

    void new_asymmelastic1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double t_elastic_modulus, c_elastic_modulus;
        if(!get_input(command, t_elastic_modulus, c_elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<AsymmElastic1D>(tag, t_elastic_modulus, c_elastic_modulus, density);
    }

    void new_elastic1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<Elastic1D>(tag, elastic_modulus, density);
    }

    void new_elastic2d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        auto material_type = 0;
        if(command.eof())
            suanpan_debug("Plane stress assumed.\n");
        else if(!get_input(command, material_type)) {
            suanpan_error("A valid material type is required.\n");
            return;
        }

        return_obj = std::make_unique<Elastic2D>(tag, elastic_modulus, poissons_ratio, density, material_type == 0 ? PlaneType::S : PlaneType::E);
    }

    void new_expcc(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus, poissons_ratio, beta, m, pt, a0, e0, lambda, kappa;
        if(!get_input(command, elastic_modulus, poissons_ratio, beta, m, pt, a0, e0, lambda, kappa)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        if(std::fabs(lambda) < std::fabs(kappa)) {
            suanpan_error("The inelastic slope (lambda) must be greater than the elastic slope (kappa).\n");
            return;
        }

        return_obj = std::make_unique<ExpCC>(tag, std::fabs(elastic_modulus), std::fabs(poissons_ratio), std::fabs(beta), std::fabs(m), std::fabs(pt), std::fabs(a0), std::fabs(e0), std::fabs(lambda), std::fabs(kappa), density);
    }

    void new_customcc(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        double beta, m, pt;
        if(!get_input(command, beta)) {
            suanpan_error("A valid beta is required.\n");
            return;
        }
        if(!get_input(command, m)) {
            suanpan_error("A valid radius ratio is required.\n");
            return;
        }
        if(!get_input(command, pt)) {
            suanpan_error("A valid tensile yield strength is required.\n");
            return;
        }

        unsigned expression;
        if(!get_input(command, expression)) {
            suanpan_error("A valid expression tag is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomCC>(tag, elastic_modulus, poissons_ratio, beta, m, pt, expression, density);
    }

    void new_expdp(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        double eta_yield, eta_flow, xi, cohesion, cohesion_a, cohesion_b;
        if(!get_input(command, eta_yield)) {
            suanpan_error("A valid eta for yielding criterion is required.\n");
            return;
        }
        if(!get_input(command, eta_flow)) {
            suanpan_error("A valid eta for plasticity flow rule is required.\n");
            return;
        }
        if(!get_input(command, xi)) {
            suanpan_error("A valid xi is required.\n");
            return;
        }
        if(!get_input(command, cohesion, cohesion_a, cohesion_b)) {
            suanpan_error("A valid cohesion is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<ExpDP>(tag, elastic_modulus, poissons_ratio, eta_yield, eta_flow, xi, cohesion, cohesion_a, cohesion_b, density);
    }

    void new_expgurson(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{2E2, .3, .4, .2, 1., 1., 0., 1., 0., 0.};
        if(!get_optional_input(command, para_pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<ExpGurson>(tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7), para_pool(8), para_pool(9));
    }

    void new_customgurson(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag, expression_tag;
        if(!get_input(command, tag, expression_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{2E2, .3, 1., 1., 0., 1., 0., 0.};
        if(!get_optional_input(command, para_pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomGurson>(tag, expression_tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7));
    }

    void new_expgurson1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{2E2, .3, .4, .2, 1., 1., 0., 1., 0., 0.};
        if(!get_optional_input(command, para_pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<ExpGurson1D>(tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7), para_pool(8), para_pool(9));
    }

    void new_customgurson1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag, expression_tag;
        if(!get_input(command, tag, expression_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{2E2, .3, 1., 1., 0., 1., 0., 0.};
        if(!get_optional_input(command, para_pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomGurson1D>(tag, expression_tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7));
    }

    void new_exporthotropic(unique_ptr<Material>& return_obj, std::istringstream& command, const bool is_hoffman) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec modulus(6);
        if(!get_input(command, modulus)) {
            suanpan_error("A valid modulus is required.\n");
            return;
        }

        vec poissons_ratio(3);
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        vec stress(9);
        if(!get_input(command, stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        double a, b;
        if(!get_input(command, a)) {
            suanpan_error("A valid a is required.\n");
            return;
        }
        if(!get_input(command, b)) {
            suanpan_error("A valid b is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        if(is_hoffman) return_obj = std::make_unique<ExpHoffman>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), a, b, density);
        else return_obj = std::make_unique<ExpTsaiWu>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), a, b, density);
    }

    void new_customorthotropic(unique_ptr<Material>& return_obj, std::istringstream& command, const bool is_hoffman) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec modulus(6);
        if(!get_input(command, modulus)) {
            suanpan_error("A valid modulus is required.\n");
            return;
        }

        vec poissons_ratio(3);
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        vec stress(9);
        if(!get_input(command, stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        unsigned expression;
        if(!get_input(command, expression)) {
            suanpan_error("A valid expression tag is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        if(is_hoffman) return_obj = std::make_unique<CustomHoffman>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), expression, density);
        else return_obj = std::make_unique<CustomTsaiWu>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), expression, density);
    }

    void new_expj2(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        double yield_stress, a, b;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }
        if(!get_input(command, a, b)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<ExpJ2>(tag, elastic_modulus, poissons_ratio, yield_stress, a, b, density);
    }

    void new_customj2(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        unsigned k_tag, h_tag;
        if(!get_input(command, k_tag, h_tag)) {
            suanpan_error("A valid expression is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomJ2>(tag, elastic_modulus, poissons_ratio, k_tag, h_tag, density);
    }

    void new_expmises1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress, a, b, c;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }
        if(!get_input(command, a, b, c)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<ExpMises1D>(tag, elastic_modulus, yield_stress, a, b, c, density);
    }

    void new_dafaliasmanzari(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec p{125., .05, 1.25, .02, .9, .7, .01, 7., .1, .9, 1.1, -.7, 3.5, 4., 6E2, -130., .2, 0.};
        if(!get_optional_input(command, p)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<DafaliasManzari>(tag, p(0), p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), p(11), p(12), p(13), p(14), p(15), p(16), p(17));
    }

    void new_flag01(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        double residual;
        if(!get_input(command, residual)) {
            suanpan_error("A valid residual stress is required.\n");
            return;
        }

        auto hardening_ratio = 0.;
        if(command.eof())
            suanpan_debug("Zero hardening ratio assumed.\n");
        else if(!get_input(command, hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<Flag>(tag, elastic_modulus, yield_stress, residual, hardening_ratio, density);
    }

    void new_flag02(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double t_yield_stress;
        if(!get_input(command, t_yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        double t_residual;
        if(!get_input(command, t_residual)) {
            suanpan_error("A valid residual stress is required.\n");
            return;
        }

        double t_hardening_ratio;
        if(!get_input(command, t_hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        double c_yield_stress;
        if(!get_input(command, c_yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        double c_residual;
        if(!get_input(command, c_residual)) {
            suanpan_error("A valid residual stress is required.\n");
            return;
        }

        double c_hardening_ratio;
        if(!get_input(command, c_hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<Flag>(tag, elastic_modulus, t_yield_stress, t_residual, t_hardening_ratio, c_yield_stress, c_residual, c_hardening_ratio, density);
    }

    void new_fluid(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double bulk_modulus;
        if(!get_input(command, bulk_modulus)) {
            suanpan_error("A valid bulk modulus is required.\n");
            return;
        }

        double density;
        if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<Fluid>(tag, bulk_modulus, density);
    }

    void new_gap01(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto gap_strain = 0.;
        if(!command.eof() && !get_input(command, gap_strain)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<Gap01>(tag, elastic_modulus, yield_stress, gap_strain, density);
    }

    void new_isotropicelastic3d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<IsotropicElastic3D>(tag, elastic_modulus, poissons_ratio, density);
    }

    void new_elasticos(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<ElasticOS>(tag, elastic_modulus, poissons_ratio, density);
    }

    void new_kelvin(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned damper_tag, spring_tag;
        if(!get_input(command, damper_tag, spring_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Kelvin>(tag, damper_tag, spring_tag);
    }

    void new_lineardamage(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned mat_tag;
        if(!get_input(command, mat_tag)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        vec value(3, fill::zeros);
        if(!get_input(command, value(0))) {
            suanpan_error("A valid start strain is required.\n");
            return;
        }
        if(!get_input(command, value(1))) {
            suanpan_error("A valid end strain is required.\n");
            return;
        }
        if(!command.eof() && !get_input(command, value(2))) {
            suanpan_error("A valid end damage value is required.\n");
            return;
        }

        return_obj = std::make_unique<LinearDamage>(tag, mat_tag, value(0), value(1), value(2));
    }

    void new_laminated(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Laminated>(tag, get_remaining<uword>(command));
    }

    void new_maxwell(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag, damper_tag, spring_tag;
        if(!get_input(command, tag, damper_tag, spring_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        std::string matrix = "false";
        if(!command.eof() && !get_input(command, matrix)) {
            suanpan_error("A valid algorithm switch is required.\n");
            return;
        }

        unsigned proceed = 0;
        if(!command.eof() && !get_input(command, proceed)) {
            suanpan_error("A valid algorithm switch is required.\n");
            return;
        }

        auto beta = .5;
        if(!command.eof() && !get_input(command, beta)) {
            suanpan_error("A valid beta value is required.\n");
            return;
        }

        return_obj = std::make_unique<Maxwell>(tag, damper_tag, spring_tag, is_true(matrix), proceed, beta);
    }

    void new_mooneyrivlin(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double bulk_modulus;
        if(!get_input(command, bulk_modulus)) {
            suanpan_error("A valid bulk modulus is required.\n");
            return;
        }

        double a10;
        if(!get_input(command, a10)) {
            suanpan_error("A valid a10 is required.\n");
            return;
        }

        double a01;
        if(!get_input(command, a01)) {
            suanpan_error("A valid a01 is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<MooneyRivlin>(tag, bulk_modulus, a10, a01, density);
    }

    void new_mpf(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto hardening_ratio = .05;
        if(!command.eof() && !get_input(command, hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto R0 = 20.;
        if(!command.eof() && !get_input(command, R0)) {
            suanpan_error("A valid r0 is required.\n");
            return;
        }

        auto A1 = 18.5;
        if(!command.eof() && !get_input(command, A1)) {
            suanpan_error("A valid a1 is required.\n");
            return;
        }

        auto A2 = .15;
        if(!command.eof() && !get_input(command, A2)) {
            suanpan_error("A valid a2 is required.\n");
            return;
        }

        auto A3 = .01;
        if(!command.eof() && !get_input(command, A3)) {
            suanpan_error("A valid a3 is required.\n");
            return;
        }

        auto A4 = 7.;
        if(!command.eof() && !get_input(command, A4)) {
            suanpan_error("A valid a4 is required.\n");
            return;
        }

        std::string iso = "false";
        if(!command.eof() && !get_input(command, iso)) {
            suanpan_error("A valid isotropic hardening switch is required.\n");
            return;
        }

        std::string con = "false";
        if(!command.eof() && !get_input(command, con)) {
            suanpan_error("A valid constant radius switch is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<MPF>(tag, elastic_modulus, yield_stress, hardening_ratio, R0, A1, A2, A3, A4, is_true(iso), is_true(con), density);
    }

    void new_multilinearoo(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        mat t_backbone, c_backbone;

        std::string name;
        if(!get_input(command, name) || !t_backbone.load(name, raw_ascii) || t_backbone.empty()) {
            suanpan_error("A valid tension backbone file is required.\n");
            return;
        }
        if(!get_input(command, name) || !c_backbone.load(name, raw_ascii) || c_backbone.empty()) {
            suanpan_error("A valid compression backbone file is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        if(0. == t_backbone(0, 1)) t_backbone = t_backbone.tail_rows(t_backbone.n_rows - 1);
        if(0. == c_backbone(0, 1)) c_backbone = c_backbone.tail_rows(c_backbone.n_rows - 1);

        return_obj = std::make_unique<MultilinearOO>(tag, abs(t_backbone), -abs(c_backbone), density);
    }

    void new_multilinearpo(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        mat t_backbone, c_backbone;

        std::string name;
        if(!get_input(command, name) || !t_backbone.load(name, raw_ascii) || t_backbone.empty()) {
            suanpan_error("A valid tension backbone file is required.\n");
            return;
        }
        if(!get_input(command, name) || !c_backbone.load(name, raw_ascii) || c_backbone.empty()) {
            suanpan_error("A valid compression backbone file is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        if(0. == t_backbone(0, 1)) t_backbone = t_backbone.tail_rows(t_backbone.n_rows - 1);
        if(0. == c_backbone(0, 1)) c_backbone = c_backbone.tail_rows(c_backbone.n_rows - 1);

        return_obj = std::make_unique<MultilinearPO>(tag, abs(t_backbone), -abs(c_backbone), density);
    }

    void new_multilinearelastic1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        const auto all = get_remaining<double>(command);

        auto size = all.size();
        auto density = 0.;
        if(size % 2 == 1) {
            --size;
            density = all.back();
        }

        std::vector<double> e, s;
        for(size_t I = 0; I < size;) {
            e.emplace_back(all.at(I++));
            s.emplace_back(all.at(I++));
        }

        return_obj = std::make_unique<MultilinearElastic1D>(tag, join_rows(vec{e}, vec{s}), density);
    }

    void new_multilinearj2(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        double density;
        if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        const auto [p_strain, p_stress] = get_remaining<double, double>(command);

        return_obj = std::make_unique<MultilinearJ2>(tag, elastic_modulus, poissons_ratio, join_rows(vec{p_strain}, vec{p_stress}), density);
    }

    void new_multilinearmises1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double density;
        if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        const auto [p_strain, p_stress] = get_remaining<double, double>(command);

        return_obj = std::make_unique<MultilinearMises1D>(tag, elastic_modulus, join_rows(vec{p_strain}, vec{p_stress}), density);
    }

    void new_nle3d01(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec pool(4);
        if(!get_input(command, pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<NLE3D01>(tag, pool(0), pool(1), pool(2), pool(3), density);
    }

    void new_nonviscous01(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        const auto [m_r, m_i, s_r, s_i] = get_remaining<double, double, double, double>(command);

        auto m_imag = vec{m_i}, s_imag = vec{s_i};
        if(accu(m_imag) + accu(s_imag) > 1E-10) {
            suanpan_error("Parameters should be conjugate pairs.\n");
            return;
        }

        auto m = cx_vec{vec{m_r}, m_imag}, s = cx_vec{vec{s_r}, s_imag};

        if(const auto sum = accu(m % exp(-1E8 * s)); sum.real() * sum.real() + sum.imag() * sum.imag() > 1E-10) {
            suanpan_error("The provided kernel does not converge to zero.\n");
            return;
        }

        return_obj = std::make_unique<Nonviscous01>(tag, std::move(m), std::move(s));
    }

    void new_orthotropicelastic3d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec modulus(6);
        if(!get_input(command, modulus)) {
            suanpan_error("A valid modulus is required.\n");
            return;
        }

        vec poissons_ratio(3);
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        auto density = 0.;
        if(command.eof())
            suanpan_debug("Zero density assumed.\n");
        else if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<OrthotropicElastic3D>(tag, std::move(modulus), std::move(poissons_ratio), density);
    }

    void new_paraboliccc(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poisson's ratio is required.\n");
            return;
        }

        double beta, m, pt, a, a_slope;
        if(!get_input(command, beta)) {
            suanpan_error("A valid beta is required.\n");
            return;
        }
        if(!get_input(command, m)) {
            suanpan_error("A valid radius ratio is required.\n");
            return;
        }
        if(!get_input(command, pt)) {
            suanpan_error("A valid tensile yield strength is required.\n");
            return;
        }
        if(!get_input(command, a)) {
            suanpan_error("A valid initial size is required.\n");
            return;
        }
        if(!get_input(command, a_slope)) {
            suanpan_error("A valid hardening slope is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<ParabolicCC>(tag, elastic_modulus, poissons_ratio, beta, m, pt, a, a_slope, density);
    }

    void new_parallel(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Parallel>(tag, get_remaining<uword>(command));
    }

    void new_planestrain(unique_ptr<Material>& return_obj, std::istringstream& command, const unsigned type) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned full_tag;
        if(!get_input(command, full_tag)) {
            suanpan_error("A valid reference material tag is required.\n");
            return;
        }

        return_obj = std::make_unique<PlaneStrain>(tag, full_tag, type);
    }

    template<typename T> void new_wrapper(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned full_tag;
        if(!get_input(command, full_tag)) {
            suanpan_error("A valid reference material tag is required.\n");
            return;
        }

        auto max_iteration = 1;
        if(!command.eof() && !get_input(command, max_iteration)) {
            suanpan_error("A valid number of maximum iteration is required.\n");
            return;
        }

        return_obj = std::make_unique<T>(tag, full_tag, max_iteration);
    }

    void new_os146s(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned full_tag;
        if(!get_input(command, full_tag)) {
            suanpan_error("A valid reference material tag is required.\n");
            return;
        }

        double shear_modulus;
        if(!get_input(command, shear_modulus)) {
            suanpan_error("A valid shear modulus is required.\n");
            return;
        }

        return_obj = std::make_unique<OS146S>(tag, full_tag, shear_modulus);
    }

    void new_polyelastic1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<PolyElastic1D>(tag, get_remaining<double>(command), 0.);
    }

    void new_polyj2(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double poissons_ratio;
        if(!get_input(command, poissons_ratio)) {
            suanpan_error("A valid poissons ratio is required.\n");
            return;
        }

        double density;
        if(!get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        const auto pool = get_remaining<double>(command);

        if(pool.size() < 3) {
            suanpan_error("At least two valid parameters for hardening are required.\n");
            return;
        }

        return_obj = std::make_unique<PolyJ2>(tag, elastic_modulus, poissons_ratio, pool, density);
    }

    void new_prestrain(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag, base_tag, amplitude_tag;
        if(!get_input(command, tag, base_tag, amplitude_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid magnitude is required.\n");
            return;
        }

        return_obj = std::make_unique<Prestrain>(tag, base_tag, amplitude_tag, magnitude);
    }

    void new_rambergosgood(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto offset = 1.;
        if(!command.eof() && !get_input(command, offset)) {
            suanpan_error("A valid offset is required.\n");
            return;
        }

        auto n = 4.;
        if(!command.eof() && !get_input(command, n)) {
            suanpan_error("A valid n is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<RambergOsgood>(tag, elastic_modulus, yield_stress, offset, n, density);
    }

    void new_rebar2d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned major_tag, minor_tag;
        if(!get_input(command, major_tag, minor_tag)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        double major_ratio, minor_ratio;
        if(!get_input(command, major_ratio, minor_ratio)) {
            suanpan_error("A valid reinforcement ratio is required.\n");
            return;
        }

        return_obj = std::make_unique<Rebar2D>(tag, major_tag, minor_tag, major_ratio, minor_ratio);
    }

    void new_rebar3d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned tag_x, tag_y, tag_z;
        if(!get_input(command, tag_x, tag_y, tag_z)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        double ratio_x, ratio_y, ratio_z;
        if(!get_input(command, ratio_x, ratio_y, ratio_z)) {
            suanpan_error("A valid reinforcement ratio is required.\n");
            return;
        }

        return_obj = std::make_unique<Rebar3D>(tag, tag_x, tag_y, tag_z, ratio_x, ratio_y, ratio_z);
    }

    void new_sequential(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        const auto pool = get_remaining<uword>(command);

        if(1 == pool.size()) {
            suanpan_error("At least two material models are required.\n");
            return;
        }

        return_obj = std::make_unique<Sequential>(tag, pool);
    }

    void new_sliplock(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double elastic_modulus;
        if(!get_input(command, elastic_modulus)) {
            suanpan_error("A valid elastic modulus is required.\n");
            return;
        }

        double yield_strain;
        if(!get_input(command, yield_strain)) {
            suanpan_error("A valid yield strain is required.\n");
            return;
        }

        auto hardening_ratio = 100.;
        if(!command.eof() && !get_input(command, hardening_ratio)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        auto R0 = 20.;
        if(!command.eof() && !get_input(command, R0)) {
            suanpan_error("A valid r0 is required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        return_obj = std::make_unique<SlipLock>(tag, elastic_modulus, yield_strain, hardening_ratio, R0, density);
    }

    void new_stacked(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Stacked>(tag, get_remaining<uword>(command));
    }

    void new_subloading1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec p{2E5, 2E2, 0., 2E2, 1E1, 2E2, 0., 2E2, 1E1, 1E1, 1E1, 1E1, .7};
        if(!get_optional_input(command, p)) {
            suanpan_error("Valid inputs are required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        DataSubloading1D para{p(0), p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), 1., 0., 0., {{p(10), 1.}}, {{p(11), std::min(1. - datum::eps, p(12))}}};
        if(para.m_iso < 0. || para.m_kin < 0.) {
            suanpan_error("The evolution rate must be positive.\n");
            return;
        }
        if(para.cv < 1.) {
            suanpan_error("The viscous limit c_v must be greater than unity.\n");
            return;
        }

        return_obj = std::make_unique<Subloading1D>(tag, std::move(para), density);
    }

    void new_subloadingviscous1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec p{2E5, 2E2, 0., 2E2, 1E1, 2E2, 0., 2E2, 1E1, 1E1, 1., 0., 0., 1E1, 1E1, .7};
        if(!get_optional_input(command, p)) {
            suanpan_error("Valid inputs are required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        DataSubloading1D para{p(0), p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), p(11), p(12), {{p(13), 1.}}, {{p(14), std::min(1. - datum::eps, p(15))}}};
        if(para.m_iso < 0. || para.m_kin < 0.) {
            suanpan_error("The evolution rate must be positive.\n");
            return;
        }
        if(para.cv < 1.) {
            suanpan_error("The viscous limit c_v must be greater than unity.\n");
            return;
        }

        return_obj = std::make_unique<Subloading1D>(tag, std::move(para), density);
    }

    void new_multisubloading1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec p{2E5, 2E2, 0., 2E2, 1E1, 2E2, 0., 2E2, 1E1, 1E1, 0.};
        if(!get_input(command, p)) {
            suanpan_error("Valid inputs are required.\n");
            return;
        }

        std::vector<DataSubloading1D::Saturation> back, core;

        std::string token;
        while(!command.eof() && get_input(command, token)) {
            double a, b;
            if(is_equal("-back", token)) {
                if(!get_input(command, a, b)) {
                    suanpan_error("Valid saturation parameters are required.\n");
                    return;
                }
                back.emplace_back(a, b);
            }
            else if(is_equal("-core", token)) {
                if(!get_input(command, a, b)) {
                    suanpan_error("Valid saturation parameters are required.\n");
                    return;
                }
                core.emplace_back(a, b);
            }
            else {
                suanpan_error("Valid saturation type is required.\n");
                return;
            }
        }

        DataSubloading1D para{p(0), p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), 1., 0., 0., std::move(back), std::move(core)};
        if(para.m_iso < 0. || para.m_kin < 0.) {
            suanpan_error("The evolution rate must be positive.\n");
            return;
        }
        if(para.cv < 1.) {
            suanpan_error("The viscous limit c_v must be greater than unity.\n");
            return;
        }

        return_obj = std::make_unique<Subloading1D>(tag, std::move(para), p(10));
    }

    void new_subloading(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec p{2E5, .3, 2E2, 0., 2E2, 1E1, 2E2, 0., 2E2, 1E1, 1E1, 1E1, 1E1, .7};
        if(!get_optional_input(command, p)) {
            suanpan_error("Valid inputs are required.\n");
            return;
        }

        auto density = 0.;
        if(!command.eof() && !get_input(command, density)) {
            suanpan_error("A valid density is required.\n");
            return;
        }

        DataSubloading para{p(0), p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), {p(11), 1.}, {p(12), std::min(1. - datum::eps, p(13))}};
        if(para.m_iso < 0. || para.m_kin < 0.) {
            suanpan_error("The evolution rate must be positive.\n");
            return;
        }

        return_obj = std::make_unique<Subloading>(tag, std::move(para), density);
    }

    void new_substepping(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned material_tag;
        if(!get_input(command, material_tag)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto max_iteration = 20u;
        if(!get_optional_input(command, max_iteration)) {
            suanpan_error("A valid maximum iteration is required.\n");
            return;
        }

        return_obj = std::make_unique<Substepping>(tag, material_tag, max_iteration);
    }

    void new_rotation2d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned full_tag;
        if(!get_input(command, full_tag)) {
            suanpan_error("A valid reference material tag is required.\n");
            return;
        }

        double a;
        if(!get_input(command, a)) {
            suanpan_error("A valid angle is required.\n");
            return;
        }

        return_obj = std::make_unique<Rotation2D>(tag, full_tag, a);
    }

    void new_simplesand(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec pool{10E4, .2, .01, -.7, 5., 1.25, 1.1, 3.5, 1.915, -130., .02, 2., 0.};
        if(!get_optional_input(command, pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        return_obj = std::make_unique<SimpleSand>(tag, pool(0), pool(1), pool(2), pool(3), pool(4), pool(5), pool(6), pool(7), pool(8), pool(9), pool(10), pool(11), pool(12));
    }

    void new_steelbrb(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto pool = get_remaining<double>(command);

        if(6 == pool.size()) {
            pool.insert(pool.end(), pool.begin() + 3, pool.begin() + 6);
            pool.emplace_back(0.);
        }
        else if(7 == pool.size()) pool.insert(pool.end() - 1, pool.begin() + 3, pool.begin() + 6);
        else if(9 == pool.size()) pool.emplace_back(0.);

        if(10 != pool.size()) {
            suanpan_error("6 or 9 parameters are required.\n");
            return;
        }

        return_obj = std::make_unique<SteelBRB>(tag, pool);
    }

    void new_rotation3d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned full_tag;
        if(!get_input(command, full_tag)) {
            suanpan_error("A valid reference material tag is required.\n");
            return;
        }

        double a, b, c;
        if(!get_input(command, a, b, c)) {
            suanpan_error("A valid angle is required.\n");
            return;
        }

        return_obj = std::make_unique<Rotation3D>(tag, full_tag, a, b, c);
    }

    void new_tablecdp(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{3E4, .2, .2, 1.16, .5, 2400E-12};

        auto idx = 0;
        double para;
        while(!command.eof() && idx < 2)
            if(get_input(command, para)) para_pool(idx++) = para;

        mat c_table, t_table, dc_table, dt_table;

        auto check_file = [&](mat& table) {
            std::string table_name;
            if(!get_input(command, table_name)) {
                suanpan_error("A valid parameter is required.\n");
                return false;
            }
            if(std::error_code code; !fs::exists(table_name, code) || !table.load(table_name, raw_ascii) || table.n_cols < 2) {
                suanpan_error("Cannot load \"{}\".\n", table_name);
                return false;
            }
            if(0. != table(0)) {
                suanpan_error("Nonzero first plastic strain detected.\n");
                return false;
            }
            return true;
        };

        if(!check_file(c_table)) return;
        if(!check_file(t_table)) return;
        if(!check_file(dc_table)) return;
        if(!check_file(dt_table)) return;

        while(!command.eof() && idx < 6)
            if(get_input(command, para)) para_pool(idx++) = para;

        return_obj = std::make_unique<TableCDP>(tag, para_pool(0), para_pool(1), std::move(t_table), std::move(c_table), std::move(dt_table), std::move(dc_table), para_pool(2), para_pool(3), para_pool(4), para_pool(5));
    }

    void new_tablegurson(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec para_pool{2E2, .3, 1., 1., 0., 1., 0., 0.};

        auto idx = 0;
        double para;
        while(!command.eof() && idx < 2)
            if(get_input(command, para)) para_pool(idx++) = para;

        std::string table_name;
        if(!get_input(command, table_name)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        mat hardening_table;
        if(std::error_code code; !fs::exists(table_name, code) || !hardening_table.load(table_name, raw_ascii) || hardening_table.n_cols < 2) {
            suanpan_error("Cannot load \"{}\".\n", table_name);
            return;
        }

        while(!command.eof() && idx < 8)
            if(get_input(command, para)) para_pool(idx++) = para;

        return_obj = std::make_unique<TableGurson>(tag, para_pool(0), para_pool(1), std::move(hardening_table), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7));
    }

    void new_trilinearstraindegradation(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned mat_tag;
        if(!get_input(command, mat_tag)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        double s_strain, e_strain, e_damage;
        if(!get_input(command, s_strain)) {
            suanpan_error("A valid start strain is required.\n");
            return;
        }
        if(!get_input(command, e_strain)) {
            suanpan_error("A valid end strain is required.\n");
            return;
        }
        if(!get_input(command, e_damage)) {
            suanpan_error("A valid end damage is required.\n");
            return;
        }

        return_obj = std::make_unique<TrilinearStrainDegradation>(tag, mat_tag, s_strain, e_strain, e_damage);
    }

    void new_customdegradation(unique_ptr<Material>& return_obj, std::istringstream& command, const bool if_strain) {
        unsigned tag, mat_tag, p_expression_tag, n_expression_tag;
        if(!get_input(command, tag, mat_tag, p_expression_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        if(command.eof()) n_expression_tag = p_expression_tag;
        else if(!get_input(command, n_expression_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        if(if_strain) return_obj = std::make_unique<CustomStrainDegradation>(tag, mat_tag, p_expression_tag, n_expression_tag);
        else return_obj = std::make_unique<CustomStressDegradation>(tag, mat_tag, p_expression_tag, n_expression_tag);
    }

    void new_trivial(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Trivial>(tag);
    }

    void new_vafcrp(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec pool{2E5, .2, 4E2, 0., 5E2, 1E1, 0., 0.};
        if(!get_optional_input(command, pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        const auto all = get_remaining<double>(command);

        auto size = all.size();
        auto density = 0.;
        if(size % 2 == 1) {
            --size;
            density = all.back();
        }

        std::vector<double> ai, bi;
        for(size_t I = 0; I < size;) {
            ai.emplace_back(all.at(I++));
            bi.emplace_back(all.at(I++));
        }

        return_obj = std::make_unique<VAFCRP>(tag, DataVAFCRP{pool(0), pool(1), pool(2), pool(3), pool(4), pool(5), pool(6), pool(7), ai, bi}, density);
    }

    void new_vafcrp1d(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        vec pool{2E5, 4E2, 0., 1E2, 1E1, 0., 0.};
        if(!get_optional_input(command, pool)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        const auto all = get_remaining<double>(command);

        auto size = all.size();
        auto density = 0.;
        if(size % 2 == 1) {
            --size;
            density = all.back();
        }

        std::vector<double> ai, bi;
        for(size_t I = 0; I < size;) {
            ai.emplace_back(all.at(I++));
            bi.emplace_back(all.at(I++));
        }

        return_obj = std::make_unique<VAFCRP1D>(tag, DataVAFCRP1D{pool(0), pool(1), pool(2), pool(3), pool(4), pool(5), pool(6), ai, bi}, density);
    }

    void new_viscosity01(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double alpha;
        if(!get_input(command, alpha)) {
            suanpan_error("A valid alpha is required.\n");
            return;
        }

        double damping;
        if(!get_input(command, damping)) {
            suanpan_error("A valid damping coefficient is required.\n");
            return;
        }

        auto limit = 1.;
        if(!command.eof() && !get_input(command, limit)) {
            suanpan_error("A valid limit is required.\n");
            return;
        }

        return_obj = std::make_unique<Viscosity01>(tag, alpha, damping, limit);
    }

    void new_viscosity02(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double alpha;
        if(!get_input(command, alpha)) {
            suanpan_error("A valid alpha is required.\n");
            return;
        }

        double damping_a;
        if(!get_input(command, damping_a)) {
            suanpan_error("A valid damping coefficient for the first quadrant is required.\n");
            return;
        }

        auto damping_b = damping_a;
        if(!command.eof() && !get_input(command, damping_b)) {
            suanpan_error("A valid damping coefficient for the second quadrant is required.\n");
            return;
        }

        auto damping_c = damping_a;
        if(!command.eof() && !get_input(command, damping_c)) {
            suanpan_error("A valid damping coefficient for the third quadrant is required.\n");
            return;
        }

        auto damping_d = damping_a;
        if(!command.eof() && !get_input(command, damping_d)) {
            suanpan_error("A valid damping coefficient for the fourth quadrant is required.\n");
            return;
        }

        auto gap_a = 1E3;
        if(!command.eof() && !get_input(command, gap_a)) {
            suanpan_error("A valid gap size for strain axis is required.\n");
            return;
        }

        auto gap_b = 1E3;
        if(!command.eof() && !get_input(command, gap_b)) {
            suanpan_error("A valid gap size for strain rate axis is required.\n");
            return;
        }

        auto limit = 1.;
        if(!command.eof() && !get_input(command, limit)) {
            suanpan_error("A valid limit is required.\n");
            return;
        }

        return_obj = std::make_unique<Viscosity02>(tag, alpha, damping_a, damping_b, damping_c, damping_d, gap_a, gap_b, limit);
    }

    void new_bilinearviscosity(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        double damping;
        if(!get_input(command, damping)) {
            suanpan_error("A valid damping coefficient is required.\n");
            return;
        }

        double yield_stress;
        if(!get_input(command, yield_stress)) {
            suanpan_error("A valid yield stress is required.\n");
            return;
        }

        auto hardening = 0.;
        if(!command.eof() && !get_input(command, hardening)) {
            suanpan_error("A valid hardening ratio is required.\n");
            return;
        }

        return_obj = std::make_unique<BilinearViscosity>(tag, damping, yield_stress, hardening);
    }

    void new_customviscosity(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag, expression_tag;
        if(!get_input(command, tag, expression_tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<CustomViscosity>(tag, expression_tag);
    }

    void new_yeoh(unique_ptr<Material>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        const auto pool = get_remaining<double>(command);

        const auto t_size = static_cast<long long>(pool.size());
        const auto h_size = t_size / 2;

        return_obj = std::make_unique<Yeoh>(tag, std::vector(pool.begin(), pool.begin() + h_size), std::vector(pool.begin() + h_size, pool.begin() + 2 * h_size), t_size % 2 == 0 ? 0. : pool.back());
    }
} // namespace

int create_new_material(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string material_id;
    if(!get_input(command, material_id)) {
        suanpan_error("A valid material type is required.\n");
        return 0;
    }

    unique_ptr<Material> new_material = nullptr;

    if(is_equal(material_id, "AFC") || is_equal(material_id, "AFC01")) new_afc01(new_material, command);
    else if(is_equal(material_id, "AFC02") || is_equal(material_id, "AFCS")) new_afc02(new_material, command);
    else if(is_equal(material_id, "AFC03") || is_equal(material_id, "AFCN")) new_afc03(new_material, command);
    else if(is_equal(material_id, "ArmstrongFrederick")) new_armstrongfrederick(new_material, command);
    else if(is_equal(material_id, "ArmstrongFrederick1D")) new_armstrongfrederick1d(new_material, command);
    else if(is_equal(material_id, "AFCO1D")) new_armstrongfrederick1d(new_material, command, true);
    else if(is_equal(material_id, "AsymmElastic1D")) new_asymmelastic1d(new_material, command);
    else if(is_equal(material_id, "Axisymmetric")) new_axisymmetric(new_material, command);
    else if(is_equal(material_id, "AxisymmetricElastic")) new_axisymmetricelastic(new_material, command);
    else if(is_equal(material_id, "Bilinear1D")) new_bilinear1d(new_material, command);
    else if(is_equal(material_id, "BilinearCC")) new_bilinearcc(new_material, command);
    else if(is_equal(material_id, "BilinearDP")) new_bilineardp(new_material, command);
    else if(is_equal(material_id, "BilinearElastic1D")) new_bilinearelastic1d(new_material, command);
    else if(is_equal(material_id, "BilinearHoffman")) new_bilinearorthotropic(new_material, command, true);
    else if(is_equal(material_id, "BilinearTsaiWu")) new_bilinearorthotropic(new_material, command, false);
    else if(is_equal(material_id, "BilinearJ2")) new_bilinearj2(new_material, command);
    else if(is_equal(material_id, "BilinearMises1D")) new_bilinearmises1d(new_material, command);
    else if(is_equal(material_id, "BilinearOO")) new_bilinearoo(new_material, command);
    else if(is_equal(material_id, "BilinearPeric")) new_bilinearperic(new_material, command);
    else if(is_equal(material_id, "BilinearPO")) new_bilinearpo(new_material, command);
    else if(is_equal(material_id, "BilinearViscosity")) new_bilinearviscosity(new_material, command);
    else if(is_equal(material_id, "BlatzKo")) new_blatzko(new_material, command);
    else if(is_equal(material_id, "BoucWen")) new_boucwen(new_material, command);
    else if(is_equal(material_id, "BWBN")) new_bwbn(new_material, command);
    else if(is_equal(material_id, "CDP")) new_cdp(new_material, command);
    else if(is_equal(material_id, "CDPM2") || is_equal(material_id, "CDPM2ISO")) new_cdpm2(new_material, command, CDPM2::DamageType::ISOTROPIC);
    else if(is_equal(material_id, "CDPM2ANISO")) new_cdpm2(new_material, command, CDPM2::DamageType::ANISOTROPIC);
    else if(is_equal(material_id, "CDPM2NO")) new_cdpm2(new_material, command, CDPM2::DamageType::NODAMAGE);
    else if(is_equal(material_id, "Concrete21")) new_concrete21(new_material, command);
    else if(is_equal(material_id, "Concrete22")) new_concrete22(new_material, command);
    else if(is_equal(material_id, "ConcreteCM")) new_concretecm(new_material, command);
    else if(is_equal(material_id, "ConcreteExp")) new_concreteexp(new_material, command);
    else if(is_equal(material_id, "ConcreteK4")) new_concretek4(new_material, command);
    else if(is_equal(material_id, "ConcreteTable")) new_concretetable(new_material, command);
    else if(is_equal(material_id, "ConcreteTsai")) new_concretetsai(new_material, command);
    else if(is_equal(material_id, "CoulombFriction")) new_coulombfriction(new_material, command);
    else if(is_equal(material_id, "CustomCC")) new_customcc(new_material, command);
    else if(is_equal(material_id, "CustomCDP")) new_customcdp(new_material, command);
    else if(is_equal(material_id, "CustomDP")) new_customdp(new_material, command);
    else if(is_equal(material_id, "CustomElastic1D")) new_customelastic1d(new_material, command);
    else if(is_equal(material_id, "CustomGurson")) new_customgurson(new_material, command);
    else if(is_equal(material_id, "CustomGurson1D")) new_customgurson1d(new_material, command);
    else if(is_equal(material_id, "CustomHoffman")) new_customorthotropic(new_material, command, true);
    else if(is_equal(material_id, "CustomTsaiWu")) new_customorthotropic(new_material, command, false);
    else if(is_equal(material_id, "CustomJ2")) new_customj2(new_material, command);
    else if(is_equal(material_id, "CustomMises1D")) new_custommises1d(new_material, command);
    else if(is_equal(material_id, "CustomStrainDegradation")) new_customdegradation(new_material, command, true);
    else if(is_equal(material_id, "CustomStressDegradation")) new_customdegradation(new_material, command, false);
    else if(is_equal(material_id, "CustomViscosity")) new_customviscosity(new_material, command);
    else if(is_equal(material_id, "DafaliasManzari")) new_dafaliasmanzari(new_material, command);
    else if(is_equal(material_id, "Dhakal")) new_dhakal(new_material, command);
    else if(is_equal(material_id, "DuncanSelig")) new_duncanselig(new_material, command);
    else if(is_equal(material_id, "Elastic1D")) new_elastic1d(new_material, command);
    else if(is_equal(material_id, "Elastic2D")) new_elastic2d(new_material, command);
    else if(is_equal(material_id, "Elastic3D") || is_equal(material_id, "IsotropicElastic3D")) new_isotropicelastic3d(new_material, command);
    else if(is_equal(material_id, "ElasticOS")) new_elasticos(new_material, command);
    else if(is_equal(material_id, "ExpCC")) new_expcc(new_material, command);
    else if(is_equal(material_id, "ExpDP")) new_expdp(new_material, command);
    else if(is_equal(material_id, "ExpGurson")) new_expgurson(new_material, command);
    else if(is_equal(material_id, "ExpGurson1D")) new_expgurson1d(new_material, command);
    else if(is_equal(material_id, "ExpHoffman")) new_exporthotropic(new_material, command, true);
    else if(is_equal(material_id, "ExpTsaiWu")) new_exporthotropic(new_material, command, false);
    else if(is_equal(material_id, "ExpJ2")) new_expj2(new_material, command);
    else if(is_equal(material_id, "ExpMises1D")) new_expmises1d(new_material, command);
    else if(is_equal(material_id, "Flag01")) new_flag01(new_material, command);
    else if(is_equal(material_id, "Flag02")) new_flag02(new_material, command);
    else if(is_equal(material_id, "Fluid")) new_fluid(new_material, command);
    else if(is_equal(material_id, "Gap01")) new_gap01(new_material, command);
    else if(is_equal(material_id, "Kelvin")) new_kelvin(new_material, command);
    else if(is_equal(material_id, "Laminated")) new_laminated(new_material, command);
    else if(is_equal(material_id, "LinearDamage")) new_lineardamage(new_material, command);
    else if(is_equal(material_id, "Maxwell")) new_maxwell(new_material, command);
    else if(is_equal(material_id, "MooneyRivlin")) new_mooneyrivlin(new_material, command);
    else if(is_equal(material_id, "MPF")) new_mpf(new_material, command);
    else if(is_equal(material_id, "MultilinearElastic1D")) new_multilinearelastic1d(new_material, command);
    else if(is_equal(material_id, "MultilinearJ2")) new_multilinearj2(new_material, command);
    else if(is_equal(material_id, "MultilinearMises1D")) new_multilinearmises1d(new_material, command);
    else if(is_equal(material_id, "MultilinearOO")) new_multilinearoo(new_material, command);
    else if(is_equal(material_id, "MultilinearPO")) new_multilinearpo(new_material, command);
    else if(is_equal(material_id, "NLE1D01")) new_nle1d01(new_material, command);
    else if(is_equal(material_id, "NLE3D01")) new_nle3d01(new_material, command);
    else if(is_equal(material_id, "Nonviscous01")) new_nonviscous01(new_material, command);
    else if(is_equal(material_id, "OrthotropicElastic3D")) new_orthotropicelastic3d(new_material, command);
    else if(is_equal(material_id, "OS146")) new_wrapper<OS146>(new_material, command);
    else if(is_equal(material_id, "OS146S")) new_os146s(new_material, command);
    else if(is_equal(material_id, "ParabolicCC")) new_paraboliccc(new_material, command);
    else if(is_equal(material_id, "Parallel")) new_parallel(new_material, command);
    else if(is_equal(material_id, "PlaneStrain")) new_planestrain(new_material, command, 0);
    else if(is_equal(material_id, "PlaneStress")) new_wrapper<PlaneStress>(new_material, command);
    else if(is_equal(material_id, "PlaneSymmetric13")) new_planestrain(new_material, command, 1);
    else if(is_equal(material_id, "PlaneSymmetric23")) new_planestrain(new_material, command, 2);
    else if(is_equal(material_id, "PolyElastic1D")) new_polyelastic1d(new_material, command);
    else if(is_equal(material_id, "PolyJ2")) new_polyj2(new_material, command);
    else if(is_equal(material_id, "Prestrain")) new_prestrain(new_material, command);
    else if(is_equal(material_id, "RambergOsgood")) new_rambergosgood(new_material, command);
    else if(is_equal(material_id, "Rebar2D")) new_rebar2d(new_material, command);
    else if(is_equal(material_id, "Rebar3D")) new_rebar3d(new_material, command);
    else if(is_equal(material_id, "Rotation2D")) new_rotation2d(new_material, command);
    else if(is_equal(material_id, "Rotation3D")) new_rotation3d(new_material, command);
    else if(is_equal(material_id, "Sequential")) new_sequential(new_material, command);
    else if(is_equal(material_id, "SimpleSand")) new_simplesand(new_material, command);
    else if(is_equal(material_id, "Sinh1D")) new_sinh1d(new_material, command);
    else if(is_equal(material_id, "SlipLock")) new_sliplock(new_material, command);
    else if(is_equal(material_id, "Stacked")) new_stacked(new_material, command);
    else if(is_equal(material_id, "SteelBRB")) new_steelbrb(new_material, command);
    else if(is_equal(material_id, "Subloading1D")) new_subloading1d(new_material, command);
    else if(is_equal(material_id, "SubloadingViscous1D")) new_subloadingviscous1d(new_material, command);
    else if(is_equal(material_id, "MultiSubloading1D")) new_multisubloading1d(new_material, command);
    else if(is_equal(material_id, "Subloading")) new_subloading(new_material, command);
    else if(is_equal(material_id, "Substepping")) new_substepping(new_material, command);
    else if(is_equal(material_id, "TableCDP")) new_tablecdp(new_material, command);
    else if(is_equal(material_id, "TableGurson")) new_tablegurson(new_material, command);
    else if(is_equal(material_id, "Tanh1D")) new_tanh1d(new_material, command);
    else if(is_equal(material_id, "TimberPD")) new_timberpd(new_material, command);
    else if(is_equal(material_id, "TrilinearStrainDegradation")) new_trilinearstraindegradation(new_material, command);
    else if(is_equal(material_id, "Trivial")) new_trivial(new_material, command);
    else if(is_equal(material_id, "Uniaxial")) new_wrapper<Uniaxial>(new_material, command);
    else if(is_equal(material_id, "VAFCRP")) new_vafcrp(new_material, command);
    else if(is_equal(material_id, "VAFCRP1D")) new_vafcrp1d(new_material, command);
    else if(is_equal(material_id, "Viscosity01")) new_viscosity01(new_material, command);
    else if(is_equal(material_id, "Viscosity02")) new_viscosity02(new_material, command);
    else if(is_equal(material_id, "Yeoh")) new_yeoh(new_material, command);
    else external_module::object(new_material, domain, material_id, command);

    if(nullptr == new_material || !domain->insert(std::move(new_material)))
        suanpan_error("Fail to create new material via \"{}\".\n", command.str());

    return 0;
}
