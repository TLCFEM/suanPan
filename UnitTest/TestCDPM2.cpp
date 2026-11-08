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
 * @version 0.2.0
 * @addtogroup UnitTest
 * @{
 */

#include "CatchHeader.h"

#include <Material/Material3D/Concrete/CDPM2.h>
#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

namespace {
    constexpr auto kp_index = 0u;
    constexpr auto kdt_index = 10u;
    constexpr auto kdc_index = 11u;
    constexpr auto kdt1_index = 12u;
    constexpr auto kdc1_index = 13u;
    constexpr auto kdt2_index = 14u;
    constexpr auto kdc2_index = 15u;
    constexpr auto omegat_index = 16u;
    constexpr auto omegac_index = 17u;

    constexpr double elastic_modulus = 3E4;
    constexpr double poissons_ratio = .3;
    constexpr double ft = 3.;
    constexpr double fc = 30.;
    constexpr double qh0 = .3;
    constexpr double hp = .01;
    constexpr double df = .85;
    constexpr double ah = .08;
    constexpr double bh = .003;
    constexpr double ch = 2.;
    constexpr double dh = 1E-6;
    constexpr double as = 5.;
    constexpr double eft = 5E-4;
    constexpr double efc = 5E-4;

    unique_ptr<CDPM2> make_cdpm2(const CDPM2::DamageType damage_type = CDPM2::DamageType::NODAMAGE) {
        auto material = std::make_unique<CDPM2>(1u, elastic_modulus, poissons_ratio, ft, fc, qh0, hp, df, ah, bh, ch, dh, as, eft, efc, damage_type, 0.);

        material->initialize_base(nullptr);
        material->initialize(nullptr);
        material->set_initialized(true);

        return material;
    }

    vec get_history(const CDPM2& material) { return material.record(OutputType::HIST).front(); }

    double qh1_value(const double kp) {
        if(kp >= 1.) return 1.;
        return qh0 + (1. - qh0) * kp * (kp * (kp - 3.) + 3.) - hp * kp * (kp - 1.) * (kp - 2.);
    }

    double qh2_value(const double kp) { return kp >= 1. ? 1. + hp * (kp - 1.) : 1.; }

    // Independent evaluation of the yield function (18) of Grassl et al. (2013).
    double yield_value(const vec& stress, const double kp) {
        const DataCDPM2 data;
        const auto deviator = tensor::dev(stress);
        const auto s = tensor::stress::norm(deviator);
        const auto p = tensor::mean<3>(stress);

        auto lode = 2.;
        if(s > datum::eps) lode = 2. * std::cos(std::acos(tensor::stress::lode(deviator)) / 3.);

        const auto square_lode = lode * lode;
        const auto square_term = data.ra * square_lode;
        const auto sqrt_term = std::sqrt(data.rb * square_term + data.rc);
        const auto r = (square_term + data.rb) / (data.ra * lode + sqrt_term);

        const auto hardening = qh1_value(kp);
        const auto hardening_two = qh2_value(kp);
        const auto g3 = (s / std::sqrt(6.) + p) / data.fc;
        const auto g1 = (1. - hardening) * g3 * g3 + std::sqrt(1.5) * s / data.fc;
        const auto g4 = (r * s / std::sqrt(6.) + p) / data.fc;

        return g1 * g1 + data.m0 * hardening * hardening * hardening_two * g4 - hardening * hardening * hardening_two * hardening_two;
    }

    std::pair<vec, vec> tensile_compressive_split(const vec& stress) {
        vec principal_stress;
        mat principal_direction;
        if(!eig_sym(principal_stress, principal_direction, tensor::stress::to_tensor(stress), "std")) throw std::runtime_error("eig_sym() failed");

        mat tensile(3, 3, fill::zeros);
        mat compressive(3, 3, fill::zeros);
        for(auto I = 0u; I < 3u; ++I) {
            const auto& direction = principal_direction.col(I);
            const mat dyad = direction * direction.t();
            if(principal_stress(I) > 0.) tensile += principal_stress(I) * dyad;
            else compressive += principal_stress(I) * dyad;
        }

        return {tensor::stress::to_voigt(tensile), tensor::stress::to_voigt(compressive)};
    }

    void require_damage_law(const CDPM2& material) {
        const auto history = get_history(material);
        const DataCDPM2 data;

        const auto check = [](const double kd, const double kd1, const double kd2, const double omega, const double ef, const double tensile_strength, const double modulus) {
            if(kd < tensile_strength / modulus) {
                REQUIRE(omega == Approx(0.).margin(1E-14));
                return;
            }

            const auto residual = tensile_strength * std::exp(-(kd1 + omega * kd2) / ef) - (1. - omega) * modulus * kd;
            REQUIRE(std::fabs(residual) <= 1E-8 * std::max({1., std::fabs(modulus * kd), tensile_strength}));
            REQUIRE(omega >= -1E-14);
            REQUIRE(omega <= 1. + 1E-14);
        };

        check(history(kdt_index), history(kdt1_index), history(kdt2_index), history(omegat_index), data.eft, data.ft, data.elastic_modulus);
        check(history(kdc_index), history(kdc1_index), history(kdc2_index), history(omegac_index), data.efc, data.ft, data.elastic_modulus);
    }

    void require_effective_stress(CDPM2& material, const CDPM2::DamageType damage_type) {
        const auto history = get_history(material);
        const vec plastic_strain = history.subvec(1, 6);
        const vec effective_stress = material.get_initial_stiffness() * (material.get_current_strain() - plastic_strain);
        const auto& nominal_stress = material.get_current_stress();

        if(CDPM2::DamageType::NODAMAGE == damage_type) {
            REQUIRE(norm(nominal_stress - effective_stress) <= 1E-9 * std::max(1., norm(effective_stress)));
            return;
        }

        const auto omegat = history(omegat_index);
        const auto omegac = history(omegac_index);

        if(CDPM2::DamageType::ISOTROPIC == damage_type) {
            const vec expected = (1. - omegat) * (1. - omegac) * effective_stress;
            REQUIRE(norm(nominal_stress - expected) <= 1E-9 * std::max(1., norm(expected)));
            return;
        }

        const auto [tensile_stress, compressive_stress] = tensile_compressive_split(effective_stress);
        const vec expected = (1. - omegat) * tensile_stress + (1. - omegac) * compressive_stress;
        REQUIRE(norm(nominal_stress - expected) <= 1E-9 * std::max(1., norm(expected)));
    }

    void require_consistent_tangent(const std::vector<vec>& preload, const vec& increment, const CDPM2::DamageType damage_type = CDPM2::DamageType::NODAMAGE, const double relative_tolerance = 1E-3) {
        const auto make_loaded = [&] {
            auto material = make_cdpm2(damage_type);
            for(const auto& preload_increment : preload) {
                REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(preload_increment));
                material->commit_status();
            }
            return material;
        };

        auto material = make_loaded();
        const vec target = material->get_current_strain() + increment;

        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(target));

        const mat tangent = material->get_trial_stiffness();
        const vec stress = material->get_trial_stress();

        mat numerical(6, 6, fill::none);

        for(auto I = 0u; I < 6u; ++I) {
            const auto delta = 1E-8 * std::max(1., std::fabs(target(I)));
            auto perturbed = target;
            perturbed(I) += delta;

            auto other = make_loaded();
            REQUIRE(SUANPAN_SUCCESS == other->update_trial_status(perturbed));

            numerical.col(I) = (other->get_trial_stress() - stress) / delta;
        }

        REQUIRE(norm(tangent - numerical) <= relative_tolerance * norm(numerical));
    }
} // namespace

TEST_CASE("CDPM2 Elastic Response", "[Material.CDPM2]") {
    const auto material = make_cdpm2();

    const vec strain{1E-6, -2E-6, 3E-7, 1E-7, -2E-7, 3E-7};

    REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(strain));
    REQUIRE(norm(material->get_trial_stress() - material->get_initial_stiffness() * strain) <= 1E-10 * norm(material->get_initial_stiffness() * strain));
    REQUIRE(norm(material->get_trial_stiffness() - material->get_initial_stiffness()) <= 1E-10 * norm(material->get_initial_stiffness()));

    material->commit_status();
    require_effective_stress(*material, CDPM2::DamageType::NODAMAGE);
}

TEST_CASE("CDPM2 Hydrostatic Elastic Response", "[Material.CDPM2]") {
    // a purely volumetric state has no unique deviatoric direction and must not generate NaNs
    const vec strain{-2E-6, -2E-6, -2E-6, 0., 0., 0.};

    auto material = make_cdpm2();

    REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(strain));
    REQUIRE(material->get_trial_stress().is_finite());
    REQUIRE(material->get_trial_stiffness().is_finite());
    REQUIRE(norm(material->get_trial_stress() - material->get_initial_stiffness() * strain) <= 1E-10 * norm(material->get_initial_stiffness() * strain));
    REQUIRE(norm(material->get_trial_stiffness() - material->get_initial_stiffness()) <= 1E-10 * norm(material->get_initial_stiffness()));
}

TEST_CASE("CDPM2 Hydrostatic Plastic Response", "[Material.CDPM2]") {
    const auto material = make_cdpm2(CDPM2::DamageType::ANISOTROPIC);

    for(auto I = 0u; I < 50u; ++I) {
        REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{-1E-4, -1E-4, -1E-4, 0., 0., 0.}));
        REQUIRE(material->get_trial_stress().is_finite());
        REQUIRE(material->get_trial_stiffness().is_finite());
        material->commit_status();
    }
}

TEST_CASE("CDPM2 Consistent Tangent", "[Material.CDPM2]") {
    SECTION("Uniaxial Compression") { require_consistent_tangent({}, vec{0., 0., -1E-2, 0., 0., 0.}); }
    SECTION("Uniaxial Tension") { require_consistent_tangent({}, vec{5E-4, 0., 0., 0., 0., 0.}); }
    SECTION("Shear") { require_consistent_tangent({}, vec{0., 0., 0., 0., 0., -1E-2}); }
    SECTION("Compression And Shear") { require_consistent_tangent({}, vec{-2E-3, -1E-3, -2E-3, 0., 1E-3, 0.}); }
    SECTION("Plastic With Isotropic Damage") { require_consistent_tangent({}, vec{0., 0., -1E-2, 0., 0., 0.}, CDPM2::DamageType::ISOTROPIC); }
    SECTION("Plastic With Anisotropic Damage") { require_consistent_tangent({}, vec{0., 0., -1E-2, 0., 0., 0.}, CDPM2::DamageType::ANISOTROPIC); }
}

TEST_CASE("CDPM2 Consistent Tangent After Damage", "[Material.CDPM2]") {
    const std::vector<vec> preload(40u, vec{1E-5, 0., 0., 0., 0., 0.});

    SECTION("Unloading With Isotropic Damage") { require_consistent_tangent(preload, vec{-5E-5, 0., 0., 0., 0., 0.}, CDPM2::DamageType::ISOTROPIC); }
    SECTION("Unloading With Anisotropic Damage") { require_consistent_tangent(preload, vec{-5E-5, 0., 0., 0., 0., 0.}, CDPM2::DamageType::ANISOTROPIC); }
    SECTION("Reloading With Anisotropic Damage") { require_consistent_tangent(preload, vec{6E-5, 0., 0., 0., 0., 0.}, CDPM2::DamageType::ANISOTROPIC); }
}

TEST_CASE("CDPM2 Damage Evolution", "[Material.CDPM2]") {
    const auto check = [](const vec& increment, const bool tensile, const unsigned steps) {
        auto material = make_cdpm2(CDPM2::DamageType::ANISOTROPIC);

        auto previous_omegat = 0.;
        auto previous_omegac = 0.;

        for(auto I = 0u; I < steps; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(increment));
            material->commit_status();

            const auto history = get_history(*material);
            const auto omegat = history(omegat_index);
            const auto omegac = history(omegac_index);

            REQUIRE(omegat >= previous_omegat - 1E-14);
            REQUIRE(omegac >= previous_omegac - 1E-14);
            REQUIRE(omegat >= -1E-14);
            REQUIRE(omegat <= 1. + 1E-14);
            REQUIRE(omegac >= -1E-14);
            REQUIRE(omegac <= 1. + 1E-14);

            REQUIRE(material->get_current_stress().is_finite());
            REQUIRE(material->get_current_stiffness().is_finite());

            require_damage_law(*material);
            require_effective_stress(*material, CDPM2::DamageType::ANISOTROPIC);

            previous_omegat = omegat;
            previous_omegac = omegac;
        }

        // compressive damage in CDPM2 is much weaker than tensile damage for these parameters
        const auto min_damage = tensile ? .1 : 1E-3;
        REQUIRE((tensile ? get_history(*material)(omegat_index) : get_history(*material)(omegac_index)) > min_damage);
    };

    SECTION("Monotonic Tension") { check(vec{1E-5, 0., 0., 0., 0., 0.}, true, 200u); }
    SECTION("Monotonic Compression") { check(vec{9E-5, 9E-5, -3E-4, 0., 0., 0.}, false, 300u); }
}

TEST_CASE("CDPM2 Return Mapping Satisfies Yield Condition", "[Material.CDPM2]") {
    const auto check = [](const vec& increment) {
        const auto material = make_cdpm2();

        auto previous_kp = 0.;
        auto yielded = false;

        for(auto I = 0u; I < 200u; ++I) {
            REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(increment));
            material->commit_status();

            const auto history = get_history(*material);
            const auto kp = history(kp_index);

            REQUIRE(kp >= previous_kp - 1E-14);

            if(kp > previous_kp + 1E-14) {
                const auto residual = yield_value(material->get_current_stress(), kp);
                REQUIRE(std::fabs(residual) <= 1E-8 * std::max(1., norm(material->get_current_stress())));
                yielded = true;
            }

            previous_kp = kp;
        }

        REQUIRE(yielded);
    };

    SECTION("Uniaxial Compression") { check(vec{0., 0., -2E-5, 0., 0., 0.}); }
    SECTION("Uniaxial Tension") { check(vec{1E-5, 0., 0., 0., 0., 0.}); }
}

TEST_CASE("CDPM2 Zero Increment", "[Material.CDPM2]") {
    const auto material = make_cdpm2(CDPM2::DamageType::ANISOTROPIC);

    for(auto I = 0u; I < 20u; ++I) {
        REQUIRE(SUANPAN_SUCCESS == material->update_incre_status(vec{1E-5, 0., 0., 0., 0., 0.}));
        material->commit_status();
    }

    REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(material->get_current_strain()));
    REQUIRE(norm(material->get_trial_stress() - material->get_current_stress()) <= 1E-10);
    REQUIRE(norm(material->get_trial_stiffness() - material->get_current_stiffness()) <= 1E-10);
}

TEST_CASE("CDPM2 Fully Damaged State", "[Material.CDPM2]") {
    const auto material = make_cdpm2(CDPM2::DamageType::ANISOTROPIC);

    const vec incre{0., 0., 1E-3, 0., 0., 0.};

    // once fully damaged the stress vanishes, subsequent increments must still be integrated
    for(auto I = 0u; I < 100u; ++I) {
        REQUIRE(SUANPAN_SUCCESS == material->update_trial_status(material->get_current_strain() + incre));
        REQUIRE(material->get_trial_stress().is_finite());
        material->commit_status();
    }
}

//! @}
