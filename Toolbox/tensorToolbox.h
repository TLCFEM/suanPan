/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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

#ifndef TENSORTOOLBOX_H
#define TENSORTOOLBOX_H

#include <suanPan.h>

template<typename T> class Quaternion;

namespace tensor {
	mat isotropic_stiffness(double, double);
	mat orthotropic_stiffness(const vec&, const vec&);

	mat unit_deviatoric_tensor4();
	mat unit_deviatoric_tensor4v2();
	mat unit_symmetric_tensor4();

	static const vec unit_tensor2{1., 1., 1., 0., 0., 0.};

	namespace stress {
		// applies to 3D tensor only, either principal or not
		double invariant1(const vec&);
		// applies to 3D tensor only, either principal or not
		double invariant2(const vec&);
		// applies to 3D tensor only, either principal or not
		double invariant3(const vec&);

		// compute load angle based on input deviatoric stress
		double lode(vec);

		static const vec norm_weight{1., 1., 1., 2., 2., 2.};
	} // namespace stress
	namespace strain {
		// applies to 3D tensor only, either principal or not
		double invariant1(const vec&);
		// applies to 3D tensor only, either principal or not
		double invariant2(const vec&);
		// applies to 3D tensor only, either principal or not
		double invariant3(const vec&);

		// compute load angle based on input deviatoric strain
		double lode(vec);

		static const vec norm_weight{1., 1., 1., .5, .5, .5};
	} // namespace strain
	double trace(const vec&);
	double mean(const vec&);
	vec dev(const vec&);
	vec dev(vec&&);

	mat dev(const mat&);
	mat dev(mat&&);

	namespace strain {
		mat to_green(mat&&);
		mat to_green(const mat&);
		mat to_tensor(const vec&);
		vec to_voigt(const mat&);
		double norm(const vec&);
		double norm(vec&&);
		double double_contraction(const vec&, const vec&);
		double double_contraction(vec&&, vec&&);
	} // namespace strain
	namespace stress {
		mat to_tensor(const vec&);
		vec to_voigt(const mat&);
		double norm(const vec&);
		double norm(vec&&);
		double double_contraction(const vec&, const vec&);
		double double_contraction(vec&&, vec&&);
	} // namespace stress

} // namespace tensor

namespace transform {
	double atan2(const vec&);
	mat compute_jacobian_nominal_to_principal(const mat&);
	mat compute_jacobian_principal_to_nominal(const mat&);

	template<typename T> Mat<T> skew_symm(const Mat<T>& R) {
		suanpan_debug([&]() { if(R.n_elem != 3) throw invalid_argument("need 3 element vector"); });

		Mat<T> S(3, 3, fill::zeros);

		S(0, 1) = -(S(1, 0) = R(2));
		S(2, 0) = -(S(0, 2) = R(1));
		S(1, 2) = -(S(2, 1) = R(0));

		return S;
	}

	template<typename T> Mat<T> rodrigues(const Col<T>& R) { return arma::expmat(transform::skew_symm(R)); }

	template<typename T> Quaternion<T> to_quaternion(const Col<T>& R) {
		const auto angle = arma::norm(R);

		return {std::cos(.5 * angle), std::sin(.5 * angle) / angle * R};
	}

	template<typename T> Col<T> to_pseudo(const Mat<T>& R) {
		const Mat<T> S = arma::real(arma::logmat(R));

		return {S(2, 1), S(0, 2), S(1, 0)};
	}

	namespace strain {
		double angle(const vec&);
		mat trans(double);
		vec principal(const vec&);
		vec rotate(const vec&, double);
	} // namespace strain
	namespace stress {
		double angle(const vec&);
		mat trans(double);
		vec principal(const vec&);
		vec rotate(const vec&, double);
	} // namespace stress
	namespace beam {
		mat global_to_local(double, double, double);
		mat global_to_local(const vec&, double);
	} // namespace beam
	namespace triangle {
		vec to_area_coordinate(const vec&, const mat&);
	}

} // namespace transform

namespace suanpan {
	template<typename T> T ramp(const T in) { return in > T(0) ? in : T(0); }
} // namespace suanpan

#endif
