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
 * @class Orientation
 * @brief A Orientation class.
 *
 * The Orientation class handles transformations between coordinate systems.
 * Be sure to call set_element_ptr() before updating anything.
 *
 * @author tlc
 * @date 27/06/2018
 * @version 0.1.0
 * @file Orientation.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef ORIENTATION_H
#define ORIENTATION_H

#include <Domain/Tag.h>

class Element;

class Orientation : public Tag {
protected:
    const Element* element_ptr = nullptr;

    vec z_axis;

    double length = 0., inclination = 0.;

    mat direction_cosine;

    void check_element_ptr() const;

    virtual void update_transformation() = 0;

public:
    explicit Orientation(unsigned = 0, vec&& = {});
    Orientation(const Orientation&) = default;           // copy allowed
    Orientation(Orientation&&) = delete;                 // move forbidden
    Orientation& operator=(const Orientation&) = delete; // copy assign forbidden
    Orientation& operator=(Orientation&&) = delete;      // move assign forbidden
    ~Orientation() override = default;

    void update_axis(const vec&);

    virtual void set_element_ptr(const Element*);

    [[nodiscard]] virtual bool is_nlgeom() const;

    /**
     * \return the size of nodal displacement vector
     */
    [[nodiscard]] virtual unsigned global_size() const = 0;
    /**
     * \return the size of displacement vector in the local system
     */
    [[nodiscard]] virtual unsigned local_size() const = 0;

    [[nodiscard]] double get_length() const;
    [[nodiscard]] double get_inclination() const;
    [[nodiscard]] const mat& get_transformation() const;

    virtual unique_ptr<Orientation> get_copy() = 0;

    virtual void update_status();
    virtual void commit_status();
    virtual void reset_status();
    virtual void clear_status();

    [[nodiscard]] virtual vec to_local_vec(double) const;
    [[nodiscard]] virtual vec to_global_vec(double) const;
    [[nodiscard]] virtual mat to_global_mass_mat(double) const;
    [[nodiscard]] virtual mat to_global_geometry_mat(double) const;
    [[nodiscard]] virtual mat to_global_stiffness_mat(double) const;

    [[nodiscard]] virtual vec to_local_vec(vec&&) const;
    [[nodiscard]] virtual vec to_global_vec(vec&&) const;
    [[nodiscard]] virtual mat to_global_mass_mat(mat&&) const;
    [[nodiscard]] virtual mat to_global_geometry_mat(mat&&) const;
    [[nodiscard]] virtual mat to_global_stiffness_mat(mat&&) const;

    /**
     * \brief transform anything from global to local system
     *        e.g., disp -> disp, vel -> vel, acc -> acc,
     *        not applicable to conversion such as disp -> strain
     * \return variable in local system
     */
    [[nodiscard]] virtual vec to_local_vec(const vec&) const = 0;
    /**
     * \brief transform anything from local to global system
     *        e.g., disp -> disp, vel -> vel, acc -> acc,
     *        not applicable to conversion such as disp -> strain
     * \return variable in global system
     */
    [[nodiscard]] virtual vec to_global_vec(const vec&) const = 0;
    /**
     * \brief transform anything from local to global system
     *        e.g., stiffness -> stiffness.
     * \return variable in global system
     */
    [[nodiscard]] virtual mat to_global_mass_mat(const mat&) const;
    [[nodiscard]] virtual mat to_global_geometry_mat(const mat&) const;
    [[nodiscard]] virtual mat to_global_stiffness_mat(const mat&) const = 0;
};

#endif

//! @}
