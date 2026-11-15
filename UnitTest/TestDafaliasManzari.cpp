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
 * @file TestDafaliasManzari.cpp
 * @brief Tests for the Dafalias--Manzari sand plasticity model.
 *
 * @author tlc
 * @date 15/11/2026
 * @version 0.1.0
 * @addtogroup UnitTest
 * @{
 */

#include "CatchHeader.h"

#include <Material/Material3D/Sand/DafaliasManzari.h>
#include <Toolbox/tensor.h>

namespace {
    unique_ptr<DafaliasManzari> make_dm(const double a = -.7) {
        auto material = std::make_unique<DafaliasManzari>(1u, 125., .05, 1.25, .02, .9, .7, .01, 7., .1, .9, 1.1, a, 3.5, 4., 6E2, -130., .2);

        material->initialize_base(nullptr);
        material->initialize(nullptr);
        material->set_initialized(true);

        return material;
    }

    unique_ptr<DafaliasManzari> make_confined() {
        auto material = make_dm();
        for(auto I = 0u; I < 200u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-1.5E-5, -1.5E-5, 0., 0., 0., 0.}));
            material->commit_status();
        }
        return material;
    }

    void require_consistent_tangent(const vec& increment, const double relative_tolerance = 1E-4) {
        auto material = make_confined();
        const vec target = material->get_current_strain() + increment;

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(target));

        const mat tangent = material->get_trial_stiffness();
        const vec stress = material->get_trial_stress();

        mat numerical(6, 6, fill::none);

        for(auto I = 0u; I < 6u; ++I) {
            const auto delta = 1E-8 * std::max(1., std::fabs(target(I)));
            auto perturbed = target;
            perturbed(I) += delta;

            auto other = make_confined();
            REQUIRE(SUANPAN_SUCCESS == other->update_trial_status(perturbed));

            numerical.col(I) = (other->get_trial_stress() - stress) / delta;
        }

        REQUIRE(norm(tangent - numerical) < relative_tolerance * norm(tangent));
    }
} // namespace

TEST_CASE("DafaliasManzari Consistent Tangent", "[Material.DafaliasManzari]") {
    SECTION("Compression") { require_consistent_tangent(vec{-1.5E-3, -1.5E-3, 0., 0., 0., 0.}); }
    SECTION("Compression And Shear") { require_consistent_tangent(vec{-1E-3, -3E-4, -3E-4, 0., 2E-3, 0.}); }
    SECTION("Shear") { require_consistent_tangent(vec{-1E-4, 0., 0., 0., 2E-3, 0.}); }
}

TEST_CASE("DafaliasManzari Material Properties", "[Material.DafaliasManzari]") {
    const auto material = make_dm();

    constexpr auto gi = .2 * 125. * 130.;
    constexpr auto elastic_modulus = gi * (2. + 2. * .05);

    REQUIRE(material->get(DafaliasManzari::Parameter::ELASTIC) == Approx(elastic_modulus));
    REQUIRE(material->get(DafaliasManzari::Parameter::POISSON) == Approx(.05));
    REQUIRE(material->get(DafaliasManzari::Parameter::SHEAR) == Approx(gi));
    REQUIRE(material->get(DafaliasManzari::Parameter::BULK) == Approx(elastic_modulus / 2.7));
    REQUIRE(material->get(DafaliasManzari::Parameter::PEAKSTRAIN) == Approx(0.));
    REQUIRE(material->get(DafaliasManzari::Parameter::CRACKSTRAIN) == Approx(0.));

    REQUIRE(material->get_initial_history().n_elem == 18);

    const auto copy = material->unique_copy();
    REQUIRE(copy != nullptr);
    REQUIRE(copy->get(DafaliasManzari::Parameter::ELASTIC) == Approx(elastic_modulus));
    REQUIRE(copy->get_initial_history().n_elem == 18);

    REQUIRE_NOTHROW(material->print());
}

TEST_CASE("DafaliasManzari Zero Increment", "[Material.DafaliasManzari]") {
    const auto material = make_confined();

    const auto stress = material->get_current_stress();
    const auto stiffness = material->get_current_stiffness();

    REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(material->get_current_strain()));
    REQUIRE(norm(material->get_trial_stress() - stress) < 1E-12);
    REQUIRE(norm(material->get_trial_stiffness() - stiffness) < 1E-12);
}

TEST_CASE("DafaliasManzari Small Strain Response", "[Material.DafaliasManzari]") {
    SECTION("Elastic Compression") {
        const auto material = make_dm();

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(vec{-1E-6, -1E-6, -1E-6, 0., 0., 0.}));

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        REQUIRE(tensor::mean<3>(material->get_trial_stress()) < 0.);
    }

    SECTION("Elastic Shear") {
        const auto material = make_dm();

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(vec{0., 0., 0., 0., 1E-6, 0.}));

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        REQUIRE(tensor::stress::norm(tensor::dev(material->get_trial_stress())) > 0.);
    }
}

TEST_CASE("DafaliasManzari Drained Shearing", "[Material.DafaliasManzari]") {
    SECTION("Contractive") {
        const auto material = make_dm();

        // build up a dense state before shearing
        for(auto I = 0u; I < 1000u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-2E-5, -2E-5, 0., 0., 0., 0.}));
            material->commit_status();
        }

        // constant volume shearing mobilises the dilatancy memory
        for(auto I = 0u; I < 1000u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{2E-5, -1E-5, -1E-5, 0., 0., 0.}));
            material->commit_status();
        }

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
    }

    SECTION("Dilative") {
        // a positive dilatancy coefficient reverses the sign of the plastic volumetric
        // response, which exercises the dilative branch of the model
        const auto material = make_dm(.7);

        for(auto I = 0u; I < 200u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-1.5E-5, -1.5E-5, 0., 0., 0., 0.}));
            material->commit_status();
        }

        for(auto I = 0u; I < 400u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{1E-5, -.5E-5, -.5E-5, 0., 0., 0.}));
            material->commit_status();
        }

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
    }
}

TEST_CASE("DafaliasManzari Tension Cap", "[Material.DafaliasManzari]") {
    SECTION("Elastic Tension") {
        const auto material = make_dm();

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(vec{2E-3, 2E-3, 2E-3, 0., 0., 0.}));

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        REQUIRE(std::fabs(tensor::mean<3>(material->get_trial_stress())) < 1E-10);
    }

    SECTION("Drained Extension") {
        const auto material = make_confined();

        // the pressure is driven to zero by dilatancy, after which the cap governs the response
        for(auto I = 0u; I < 70u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{6E-5, 0., 0., 0., 0., 0.}));
            material->commit_status();
        }

        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        REQUIRE(std::fabs(tensor::mean<3>(material->get_trial_stress())) < 1E-10);
    }
}
