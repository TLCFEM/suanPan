/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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
 * @class B3DC
 * @brief A B3DC class.
 * 
 * B3DC is a corotational transformation for 3D beam elements.
 * 
 * The implementation is mainly based on de Souza's thesis.
 * 
 * Force-based Finite Element for Large Displacement Inelastic Analysis of Frames
 * 
 * @author tlc
 * @date 16/12/2021
 * @version 0.1.0
 * @file B3DC.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef B3DC_H
#define B3DC_H

#include "B3DL.h"

class B3DC : public B3DL {
    mat basic{3, 3, fill::none};

    rowvec t6i, t6j;

protected:
    const span sa{0, 2}, sb{3, 5}, sc{6, 8}, sd{9, 11};

    double elongation = 0.;

    mat trial_rotation{3, 2, fill::zeros}, current_rotation{3, 2, fill::zeros};
    mat trial_n, current_n, reference;
    field<mat> sn{6}, se{3};

    vec theta;
    mat transformation;

    const double initial_length = 0.;

    [[nodiscard]] mat compute_a() const;
    [[nodiscard]] mat compute_l(const mat&, const subview_col<double>&) const;
    [[nodiscard]] mat compute_m(const mat&, const subview_col<double>&) const;
    [[nodiscard]] mat compute_g(const mat&, const subview_col<double>&, const subview_col<double>&) const;

    // some syntax sugar
    [[nodiscard]] subview_col<double> e(uword) const;
    [[nodiscard]] subview_col<double> r(uword) const;
    [[nodiscard]] subview_col<double> ni(uword) const;
    [[nodiscard]] subview_col<double> nj(uword) const;
    [[nodiscard]] const mat& sni(uword) const;
    [[nodiscard]] const mat& snj(uword) const;

    void update_direct_cosine(const vec&);
    void update_e(const vec&);
    void update_theta();

    void update_transformation() override;

    [[nodiscard]] virtual unsigned nodal_size() const;

public:
    using B3DL::B3DL;

    [[nodiscard]] bool is_nlgeom() const override;

    unique_ptr<Orientation> get_copy() override;

    void commit_status() override;
    void reset_status() override;
    void clear_status() override;

    [[nodiscard]] vec to_local_vec(const vec&) const override;
    [[nodiscard]] vec to_global_vec(const vec&) const override;
    [[nodiscard]] mat to_global_geometry_mat(const mat&) const override;
    [[nodiscard]] mat to_global_stiffness_mat(const mat&) const override;
};

#endif

//! @}
