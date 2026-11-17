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
 * @file TestDuncanSelig.cpp
 * @brief Tests for the Duncan--Selig soil model.
 *
 * @author tlc
 * @date 17/11/2026
 * @version 0.1.0
 * @addtogroup UnitTest
 * @{
 */

#include "CatchHeader.h"

#include <Material/Material2D/DuncanSelig.h>

namespace {
    unique_ptr<DuncanSelig> make_model() {
        auto material = std::make_unique<DuncanSelig>(1u, vec{14.7, 6000., .6, 4000., .2, .7, .1, .7, .5});

        material->initialize_base(nullptr);
        material->initialize(nullptr);
        material->set_initialized(true);

        return material;
    }

    constexpr auto ref_strain = 14.7 / 6000.;

    unique_ptr<DuncanSelig> make_confined() {
        auto material = make_model();
        for(auto I = 0u; I < 100u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-7.35E-5, -7.35E-5, 0.}));
            material->commit_status();
        }
        return material;
    }

    unique_ptr<DuncanSelig> make_uniaxial() {
        auto material = make_model();
        for(auto I = 0u; I < 200u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-3. * ref_strain / 200., 0., 0.}));
            material->commit_status();
        }
        return material;
    }

    unique_ptr<DuncanSelig> make_loaded() {
        auto material = make_confined();
        REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{0., 0., 1E-3}));
        material->commit_status();
        return material;
    }

    unique_ptr<DuncanSelig> make_unloaded() {
        auto material = make_loaded();
        REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{0., 0., -6E-4}));
        material->commit_status();
        return material;
    }

    void require_consistent_tangent(const unique_ptr<DuncanSelig>& material, const vec& increment, const double relative_tolerance = 1E-4) {
        const vec target = material->get_current_strain() + increment;

        auto reference = material->unique_copy();
        REQUIRE(SUANPAN_SUCCESS == reference->update_trial_status(target));

        const mat tangent = reference->get_trial_stiffness();
        const vec stress = reference->get_trial_stress();

        mat numerical(3, 3, fill::none);

        for(auto I = 0u; I < 3u; ++I) {
            const auto delta = 1E-8 * std::max(1., std::fabs(target(I)));
            auto perturbed = target;
            perturbed(I) += delta;

            auto other = material->unique_copy();
            REQUIRE(SUANPAN_SUCCESS == other->update_trial_status(perturbed));

            numerical.col(I) = (other->get_trial_stress() - stress) / delta;
        }

        REQUIRE(norm(tangent - numerical) < relative_tolerance * norm(tangent));
    }

    void require_finite_response(const unique_ptr<DuncanSelig>& material, const vec& increment) {
        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(material->get_current_strain() + increment));
        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
    }
} // namespace

TEST_CASE("DuncanSelig Consistent Tangent", "[Material.DuncanSelig]") {
    SECTION("Elastic Compression") { require_consistent_tangent(make_unloaded(), vec{-1E-4, -1E-4, 0.}); }
    SECTION("Elastic Shear") { require_consistent_tangent(make_unloaded(), vec{0., 0., 1E-4}); }
    SECTION("Plastic Compression") { require_consistent_tangent(make_loaded(), vec{-1E-2, -1E-2, 1E-4}); }
    SECTION("Plastic Compression And Shear") { require_consistent_tangent(make_loaded(), vec{-1E-3, -3E-4, 1E-3}); }
    SECTION("Uniaxial Compression") { require_consistent_tangent(make_confined(), vec{-3E-3, 0., 1E-4}); }
    SECTION("Tension") { require_consistent_tangent(make_loaded(), vec{1E-3, 1E-3, 1E-4}); }
}

TEST_CASE("DuncanSelig Robustness", "[Material.DuncanSelig]") {
    SECTION("Large Compression") { require_finite_response(make_confined(), vec{-1E-2, -1E-2, 0.}); }
    SECTION("Large Shear") { require_finite_response(make_confined(), vec{0., 0., 1E-2}); }
    SECTION("Large Tension") { require_finite_response(make_confined(), vec{1E-2, 1E-2, 0.}); }
    SECTION("Mixed") { require_finite_response(make_confined(), vec{-1E-2, 1E-2, 1E-2}); }
    SECTION("Neutral Loading") {
        const auto material = make_uniaxial();
        for(auto I = 0u; I < 100u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-0.02 * ref_strain, -0.02 * ref_strain, 0.}));
            material->commit_status();
        }
        REQUIRE(material->get_current_stress().is_finite());
        REQUIRE(material->get_current_stiffness().is_finite());
    }
    SECTION("Uniaxial Map Point") { require_finite_response(make_uniaxial(), vec{-1.067 * ref_strain, -2. * ref_strain, 0.}); }
    SECTION("Cyclic") {
        const auto material = make_confined();
        for(auto I = 0u; I < 20u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{1E-3, -1E-3, 1E-3}));
            material->commit_status();
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-1E-3, 1E-3, -1E-3}));
            material->commit_status();
        }
        REQUIRE(material->get_current_stress().is_finite());
        REQUIRE(material->get_current_stiffness().is_finite());
    }
}

//! @}
