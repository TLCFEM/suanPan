/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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

#include "ElementTemplate.h"

#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>

/**
 * @brief
 * Here we target our ElementTemplate class to fulfill the functionality of a
 * constant stress triangular element, viz., CPS3, in Abaqus notation.
 *
 * The example element is derived from the `MaterialElement2D` class, hence it
 * is a 2D element using material models (instead of sections).
 *
 * The constructor of `MaterialElement2D` class asks for six parameters:
 * - Unique Element Object Tag (T)
 * - Number of Nodes (NN)
 * - Number of DoFs (ND)
 * - Node Encoding Tags (NT)
 * - Material Tag (MT)
 * - Nonlinear Switch (F)
 *
 * In our example, CT and F will be constants, NN is 3 and ND is 2. So we have
 * three parameters plus thickness for our derived element. Except for
 * initializing private member variables, we do not have to do anything. All
 * other initializations will be handled by the Element and Domain class. As
 * simple as this.
 */
ElementTemplate::ElementTemplate(const unsigned T, uvec&& NT, const unsigned MT, const double TH)
    : MaterialElement2D(T, m_node, m_dof, std::move(NT), uvec{MT}, false)
    , thickness(TH) {}

/**
 * @brief
 * As explained before, this method get all necessary information, which
 * includes getting copies of Material objects and other operations, from the
 * Domain object.
 *
 * Please note that **we do not have to check the existence of any objects**
 * which are used in the element. The validity of the connected node objects and
 * the material models is checked in the base initialisation before calling
 * this method. The execution of this `initialize()` method automatically
 * implies that this is a valid Element object with valid material model.
 *
 * The displacement mode is
 * \f{gather}{\phi=\begin{bmatrix}1&x&y\end{bmatrix}.\f}
 *
 * The strain matrix is calculated as
 * \f{gather}{B=\partial{}N=\partial{}\left(\phi{}C^{-1}\right),\f}
 * where
 * \f{gather}{C=\begin{bmatrix}1&x_i&y_i\\1&x_j&y_j\\1&x_k&y_k\end{bmatrix}.\f}
 *
 * One can also initialize stiffness matrix and/or other build-in matrices from
 * Element class (check the definition for details) in the `initialize()` method.
 * However, this it not necessary, as the Solver will always call
 * update_status() method with a zero trial displacement to update current
 * stiffness and resistance before running iterations.
 */
int ElementTemplate::initialize(const shared_ptr<DomainBase>& D) {
    //! As CPS3 is a constant stress/strain element, one integration point at the
    //! center of the element is enough. Hence we only have one material model
    //! defined. First we get a reference of the Material object from the Domain
    //! and then call the `get_copy()` method to get a local copy. Direct
    //! assignment is allowed, the move semantics will automatically be invoked.
    //! There is no need to check if the material model is a 2D one. The validation
    //! is done in base Element class initialisation.
    m_material = D->get<Material>(material_tag(0))->get_copy();

    //! The node pointers are handled in the base Element class, we do not have to
    //! set it manually. Now we could fill in the `ele_coor` matrix. The
    //! area/natural coordinate is another version of implementation. Please refer
    //! to FEM textbooks for theories. This will be used for the computation of
    //! the shape function.
    mat ele_coor(m_node, m_node, fill::ones);
    for(unsigned i = 0; i < m_node; ++i) {
        auto& tmp_coor = node_ptr[i].lock()->get_coordinate();
        for(unsigned j = 0; j < m_dof; ++j) ele_coor(i, j + 1llu) = tmp_coor(j);
    }

    //! The area is half of the determinant of the above matrix.
    //! The area of 2D polygons can also be computed by the `shoelace` function.
    area = .5 * det(ele_coor);

    const mat inv_coor = inv(ele_coor);

    //! A standard way to construct the strain mat is to derive from the partial
    //! derivative of the shape function N. For CP3, as it is a constant
    //! stress/strain element, the derivatives are constants which can be directly
    //! obtained from above matrix.
    strain_mat.zeros(3, m_node * m_dof);
    for(unsigned i = 0, j = 0, k = 1; i < 3; ++i, j += m_dof, k += m_dof) {
        strain_mat(2, k) = strain_mat(0, j) = inv_coor(1, i);
        strain_mat(2, j) = strain_mat(1, k) = inv_coor(2, i);
    }

    trial_stiffness = current_stiffness = initial_stiffness = strain_mat.t() * m_material->get_initial_stiffness() * strain_mat * area * thickness;

    if(const auto t_density = area * thickness * m_material->get_density(); t_density > 0.) {
        initial_mass.zeros(m_size, m_size);
        const rowvec n = mean(ele_coor) * inv_coor;
        const mat t_mass = n.t() * n * t_density * area * thickness;
        initial_mass(uvec{1, 3, 5}, uvec{1, 3, 5}) = t_mass;
        initial_mass(uvec{0, 2, 4}, uvec{0, 2, 4}) = t_mass;
    }

    //! We use function `ConstantMass()` to indicate the mass matrix will not change
    //! so that `trial_mass`, `current_mass` and `initial_mass` matrices can point to
    //! the same memory location. This avoids unnecessary allocation of memory.
    //! It is not compulsory to call `ConstantMass()`, `ConstantStiffness()` and
    //! `ConstantDamping()` but highly recommended to do so when one or some matrices
    //! indeed remain unchanged for the whole analysis.
    ConstantMass(this);

    return SUANPAN_SUCCESS;
}

/**
 * @brief Now we handle the status update method. We get trial displacement via
 * build-in method and pass trial strain to the material model. Then get updated
 * stiffness and stress back to form element stiffness and resistance.
 *
 * For a static analysis, **stiffness** and **resistance** have to be
 * formulated. Apart from this, there is nothing you have to do. They will be
 * send to global assembler by methods in base Element class, which can also be
 * overridden to be customized.
 */
int ElementTemplate::update_status() {
    m_material->update_trial_status(strain_mat * get_trial_displacement());

    trial_stiffness = area * thickness * strain_mat.t() * m_material->get_trial_stiffness() * strain_mat;
    trial_resistance = area * thickness * strain_mat.t() * m_material->get_trial_stress();

    return 0;
}

/**
 * \brief Simply call corresponding methods in material objects. If the element
 * itself has history variables, they should also be updated/modified in
 * following methods.
 */
int ElementTemplate::commit_status() { return m_material->commit_status(); }

int ElementTemplate::clear_status() { return m_material->clear_status(); }

int ElementTemplate::reset_status() { return m_material->reset_status(); }
