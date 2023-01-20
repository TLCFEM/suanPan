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
 * @class Patch
 * @brief A Patch class.
 * @author tlc
 * @date 14/11/2020
 * @version 0.1.0
 * @file Patch.h
 * @addtogroup Patch
 * @ingroup Element
 * @{
 */

#ifndef PATCH_H
#define PATCH_H

#include <Element/MaterialElement.h>
#include <Element/SectionElement.h>

class Patch {
protected:
    const field<vec> knot_pool;
    vector<uvec> element_span;

public:
    explicit Patch(field<vec>&&);

    [[nodiscard]] uvec get_number_of_control_points() const;
};

class MaterialPatch : public Patch, public MaterialElement {
public:
    MaterialPatch(unsigned,     // tag
                  unsigned,     // number of dofs
                  uvec&&,       // node encoding
                  uvec&&,       // material tags
                  field<vec>&&, // knot pool
                  bool,         // nonlinear geometry switch
                  MaterialType  // material type
    );
};

class MaterialPatch2D : public MaterialPatch {
public:
    MaterialPatch2D(unsigned,     // tag
                    unsigned,     // number of dofs
                    uvec&&,       // node encoding
                    uvec&&,       // material tags
                    field<vec>&&, // knot pool
                    bool          // nonlinear geometry switch
    );
};

class MaterialPatch3D : public MaterialPatch {
public:
    MaterialPatch3D(unsigned,     // tag
                    unsigned,     // number of dofs
                    uvec&&,       // node encoding
                    uvec&&,       // material tags
                    field<vec>&&, // knot pool
                    bool          // nonlinear geometry switch
    );
};

class SectionPatch : public Patch, public SectionElement {
public:
    SectionPatch(unsigned,     // tag
                 unsigned,     // number of dofs
                 uvec&&,       // node encoding
                 uvec&&,       // section tags
                 field<vec>&&, // knot pool
                 bool,         // nonlinear geometry switch
                 SectionType   // section type
    );
};

class SectionPatch2D : public SectionPatch {
public:
    SectionPatch2D(unsigned,     // tag
                   unsigned,     // number of dofs
                   uvec&&,       // node encoding
                   uvec&&,       // section tags
                   field<vec>&&, // knot pool
                   bool          // nonlinear geometry switch
    );
};

class SectionPatch3D : public SectionPatch {
public:
    SectionPatch3D(unsigned,     // tag
                   unsigned,     // number of dofs
                   uvec&&,       // node encoding
                   uvec&&,       // section tags
                   field<vec>&&, // knot pool
                   bool          // nonlinear geometry switch
    );
};

#endif

//! @}
