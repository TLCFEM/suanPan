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
/**
 * @file TestNonlinearDruckerPrager.cpp
 * @brief Tests for the nonlinear Drucker--Prager material family.
 *
 * @author tlc
 * @date 15/11/2026
 * @version 0.1.0
 * @addtogroup UnitTest
 * @{
 */

#include "CatchHeader.h"

#include <Material/Material3D/DruckerPrager/BilinearDP.h>
#include <Material/Material3D/DruckerPrager/ExpDP.h>
#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

namespace {
    constexpr double elastic_modulus = 1E4;
    constexpr double poissons_ratio = .3;
    constexpr double eta_yield = .355;
    constexpr double eta_flow = .355;
    constexpr double xi = .974;
    constexpr double cohesion = .5;
    constexpr double cohesion_slope = 100.;

    unique_ptr<BilinearDP> make_bilinear(const double flow = eta_flow) {
        auto material = std::make_unique<BilinearDP>(1u, elastic_modulus, poissons_ratio, eta_yield, flow, xi, cohesion, cohesion_slope);

        material->initialize_base(nullptr);
        material->initialize(nullptr);
        material->set_initialized(true);

        return material;
    }

    unique_ptr<ExpDP> make_exp() {
        auto material = std::make_unique<ExpDP>(1u, elastic_modulus, poissons_ratio, eta_yield, eta_flow, xi, cohesion, 2., 200.);

        material->initialize_base(nullptr);
        material->initialize(nullptr);
        material->set_initialized(true);

        return material;
    }

    double plastic_strain(const Material& material) { return material.record(OutputType::HIST).front()(0); }
} // namespace

TEST_CASE("NonlinearDruckerPrager Material Properties", "[Material.NonlinearDruckerPrager]") {
    const auto material = make_bilinear();

    REQUIRE(material->get(BilinearDP::Parameter::ELASTIC) == Approx(elastic_modulus));
    REQUIRE(material->get(BilinearDP::Parameter::POISSON) == Approx(poissons_ratio));
    REQUIRE(material->get(BilinearDP::Parameter::SHEAR) == Approx(elastic_modulus / 2.6));
    REQUIRE(material->get(BilinearDP::Parameter::BULK) == Approx(elastic_modulus / 1.2));
    REQUIRE(material->get(BilinearDP::Parameter::PEAKSTRAIN) == Approx(0.));
    REQUIRE(material->get(BilinearDP::Parameter::CRACKSTRAIN) == Approx(0.));

    REQUIRE(material->get_initial_history().n_elem == 1);

    const auto copy = material->unique_copy();
    REQUIRE(copy != nullptr);
    REQUIRE(copy->get(BilinearDP::Parameter::ELASTIC) == Approx(elastic_modulus));
    REQUIRE(copy->get_initial_history().n_elem == 1);

    REQUIRE_NOTHROW(material->print());
    REQUIRE_NOTHROW(copy->print());
    // the derived classes override print, call the base implementation directly
    REQUIRE_NOTHROW(material->NonlinearDruckerPrager::print());
}

TEST_CASE("NonlinearDruckerPrager Zero Increment", "[Material.NonlinearDruckerPrager]") {
    const auto material = make_bilinear();

    const auto stress = material->get_current_stress();
    const auto stiffness = material->get_current_stiffness();

    REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(material->get_current_strain()));
    REQUIRE(norm(material->get_trial_stress() - stress) < 1E-12);
    REQUIRE(norm(material->get_trial_stiffness() - stiffness) < 1E-12);
    REQUIRE(0. == Approx(plastic_strain(*material)));
}

TEST_CASE("NonlinearDruckerPrager Elastic Response", "[Material.NonlinearDruckerPrager]") {
    const auto material = make_bilinear();

    const vec strain{-1E-5, -1E-5, -1E-5, 0., 0., 0.};

    REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(strain));
    REQUIRE(material->get_trial_stress().is_finite());
    REQUIRE(material->get_trial_stiffness().is_finite());
    REQUIRE(norm(material->get_trial_stress() - material->get_initial_stiffness() * strain) < 1E-10);
    REQUIRE(0. == Approx(plastic_strain(*material)));
}

TEST_CASE("NonlinearDruckerPrager Plastic Response", "[Material.NonlinearDruckerPrager]") {
    SECTION("Associated Flow") {
        const auto material = make_bilinear();

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(vec{-1E-3, 0., 0., 0., 0., 0.}));

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        material->commit_status();
        REQUIRE(plastic_strain(*material) > 0.);

        const auto previous = plastic_strain(*material);

        for(auto I = 0u; I < 10u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-5E-4, 0., 0., 0., 0., 0.}));
            material->commit_status();
        }

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        REQUIRE(plastic_strain(*material) > previous);
    }

    SECTION("Nonassociated Flow") {
        const auto material = make_bilinear(.2);

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(vec{-1E-3, 0., 0., 0., 0., 0.}));

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        material->commit_status();
        REQUIRE(plastic_strain(*material) > 0.);
    }

    SECTION("Exponential Hardening") {
        const auto material = make_exp();

        for(auto I = 0u; I < 10u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-2E-4, 0., 0., 0., 1E-4, 0.}));
            material->commit_status();
        }

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        REQUIRE(plastic_strain(*material) > 0.);
    }
}

TEST_CASE("NonlinearDruckerPrager Apex Return", "[Material.NonlinearDruckerPrager]") {
    SECTION("Associated Flow") {
        const auto material = make_bilinear();

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(vec{1E-3, 1E-3, 1E-3, 0., 0., 0.}));

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        REQUIRE(tensor::stress::norm(tensor::dev(material->get_trial_stress())) < 1E-10);
        material->commit_status();
        REQUIRE(plastic_strain(*material) > 0.);
    }

    SECTION("Nonassociated Flow") {
        const auto material = make_bilinear(.2);

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(vec{1E-3, 1E-3, 1E-3, 0., 0., 0.}));

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        REQUIRE(tensor::stress::norm(tensor::dev(material->get_trial_stress())) < 1E-10);
        material->commit_status();
        REQUIRE(plastic_strain(*material) > 0.);
    }
}
