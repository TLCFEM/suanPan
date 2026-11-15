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
 * @file TestSimpleSand.cpp
 * @brief Tests for the simple sand plasticity model.
 *
 * @author tlc
 * @date 15/11/2026
 * @version 0.1.0
 * @addtogroup UnitTest
 * @{
 */

#include "CatchHeader.h"

#include <Material/Material3D/Sand/SimpleSand.h>
#include <Toolbox/tensor.h>

namespace {
    unique_ptr<SimpleSand> make_sand() {
        auto material = std::make_unique<SimpleSand>(1u, 10E4, .2, .01, -.7, 5., 1.25, 1.1, 3.5, 1.915, -130., .02, 2.);

        material->initialize_base(nullptr);
        material->initialize(nullptr);
        material->set_initialized(true);

        return material;
    }

    unique_ptr<SimpleSand> make_confined() {
        auto material = make_sand();
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

TEST_CASE("SimpleSand Consistent Tangent", "[Material.SimpleSand]") {
    SECTION("Compression") { require_consistent_tangent(vec{-1.5E-3, -1.5E-3, 0., 0., 0., 0.}); }
    SECTION("Compression And Shear") { require_consistent_tangent(vec{-1E-3, -3E-4, -3E-4, 0., 2E-3, 0.}); }
    SECTION("Shear") { require_consistent_tangent(vec{-1E-4, 0., 0., 0., 2E-3, 0.}); }
}

TEST_CASE("SimpleSand Tension Cap", "[Material.SimpleSand]") {
    SECTION("Elastic Tension") {
        const auto material = make_sand();

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
