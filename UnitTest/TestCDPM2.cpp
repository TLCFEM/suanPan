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
 * @file TestCDPM2.cpp
 * @brief Tests for the CDPM2 concrete damage plasticity model.
 *
 * @author tlc
 * @date 04/11/2026
 * @version 0.1.0
 * @addtogroup UnitTest
 * @{
 */

#include "CatchHeader.h"

#include <Material/Material3D/Concrete/CDPM2.h>

namespace {
    unique_ptr<CDPM2> make_cdpm2(const CDPM2::DamageType damage_type = CDPM2::DamageType::NODAMAGE) {
        auto material = std::make_unique<CDPM2>(1u, 3E4, .3, 3., 30., .3, .01, .85, .08, .003, 2., 1E-6, 5., 5E-4, 5E-4, damage_type, 0.);

        material->initialize_base(nullptr);
        material->initialize(nullptr);
        material->set_initialized(true);

        return material;
    }
} // namespace

TEST_CASE("CDPM2 Consistent Tangent", "[Material.CDPM2]") {
    const auto check = [](const vec& total_incre, const CDPM2::DamageType damage_type) {
        auto material = make_cdpm2(damage_type);
        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(total_incre));

        const mat tangent = material->get_trial_stiffness();
        const vec stress = material->get_trial_stress();

        mat numerical(6, 6, fill::zeros);

        for(auto I = 0u; I < 6u; ++I) {
            const auto delta = 1E-8 * std::max(1., std::fabs(total_incre(I)));
            auto perturbed = total_incre;
            perturbed(I) += delta;

            auto other = make_cdpm2(damage_type);
            REQUIRE(SUANPAN_SUCCESS == other->update_trial_status(perturbed));
            numerical.col(I) = (other->get_trial_stress() - stress) / delta;
        }

        REQUIRE(norm(tangent - numerical) <= 1E-3 * norm(numerical));
    };

    SECTION("Plastic") { check(vec{0., 0., -1E-2, 0., 0., 0.}, CDPM2::DamageType::NODAMAGE); }
    SECTION("Plastic With Damage") { check(vec{0., 0., -1E-2, 0., 0., 0.}, CDPM2::DamageType::ANISOTROPIC); }
    SECTION("Shear") { check(vec{0., 0., 0., 0., 0., -1E-2}, CDPM2::DamageType::NODAMAGE); }
}

TEST_CASE("CDPM2 Fully Damaged State", "[Material.CDPM2]") {
    auto material = make_cdpm2(CDPM2::DamageType::ANISOTROPIC);

    const vec incre{0., 0., 1E-3, 0., 0., 0.};

    // once fully damaged the stress vanishes, subsequent increments must still be integrated
    for(auto I = 0u; I < 100u; ++I) {
        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(material->get_current_strain() + incre));
        REQUIRE(material->get_trial_stress().is_finite());
        material->commit_status();
    }
}
