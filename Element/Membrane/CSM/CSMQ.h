/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
 * @class CSMQ
 * @brief The CSMQ class.
 *
 * @author tlc
 * @date 02/06/2021
 * @version 0.1.0
 * @file CSMQ.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CSMQ_H
#define CSMQ_H

#include <Element/MaterialElement.h>

class CSMQ : public MaterialElement2D {
    struct IntegrationPoint final {
        vec coor;
        double weight;
        unique_ptr<Material> m_material;
        mat b1, b2, b3;
        IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
    };

    static constexpr unsigned m_dof = 3;

    const unsigned m_node;

    const unsigned m_size = m_dof * m_node;

    const double thickness;

    vector<IntegrationPoint> int_pt;

    virtual const uvec& get_translation_dof() = 0;
    virtual const uvec& get_rotation_dof() = 0;

public:
    CSMQ(unsigned,    // tag
         uvec&&,      // node tag
         unsigned,    // material tag
         unsigned,    // number of nodes
         double = 1., // thickness
         double = -1. // length
    );

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;

    [[nodiscard]] mat compute_shape_function(const mat&, unsigned) const override;

    vector<vec> record(OutputType) override;

    void print() override;

#ifdef SUANPAN_VTK
    void Setup() override;
    void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
    mat GetData(OutputType) override;
    void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

class CSMQ5 final : public CSMQ {
    inline static const uvec t_dof{0, 1, 3, 4, 6, 7, 9, 10, 12, 13};
    inline static const uvec r_dof{2, 5, 8, 11, 14};

    const uvec& get_translation_dof() override;
    const uvec& get_rotation_dof() override;

public:
    CSMQ5(unsigned,    // tag
          uvec&&,      // node tag
          unsigned,    // material tag
          double = 1., // thickness
          double = -1. // length
    );
};

class CSMQ6 final : public CSMQ {
    inline static const uvec t_dof{0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16};
    inline static const uvec r_dof{2, 5, 8, 11, 14, 17};

    const uvec& get_translation_dof() override;
    const uvec& get_rotation_dof() override;

public:
    CSMQ6(unsigned,    // tag
          uvec&&,      // node tag
          unsigned,    // material tag
          double = 1., // thickness
          double = -1. // length
    );
};

class CSMQ7 final : public CSMQ {
    inline static const uvec t_dof{0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19};
    inline static const uvec r_dof{2, 5, 8, 11, 14, 17, 20};

    const uvec& get_translation_dof() override;
    const uvec& get_rotation_dof() override;

public:
    CSMQ7(unsigned,    // tag
          uvec&&,      // node tag
          unsigned,    // material tag
          double = 1., // thickness
          double = -1. // length
    );
};

#endif

//! @}
