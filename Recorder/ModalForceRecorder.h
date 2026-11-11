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
 * @class ModalForceRecorder
 * @brief A ModalForceRecorder class.
 *
 * The ModalForceRecorder extracts the modal elastic force from a generalised
 * eigenvalue analysis. For the generalised eigenproblem
 *
 *     K * phi = lambda * M * phi
 *
 * the modal elastic force of mode n is
 *
 *     f_n = K * phi_n.
 *
 * The mode shapes are normalised with respect to the mass matrix, i.e.
 * phi' * M * phi = 1, so that the force corresponds to a unit modal
 * displacement.
 *
 * The force is assembled element by element, f_n = sum_e K_e * phi_{n,e},
 * rather than from the global stiffness matrix. The global stiffness matrix is
 * modified when boundary conditions and constraints are processed, whilst the
 * element stiffness matrices are not. Consequently, the assembled force vector
 * contains the reaction components at restrained DoFs, which are generally
 * nonzero. The free-DoF components are the usual elastic restoring forces.
 *
 * Both the assembled nodal force vector and the element force vectors are
 * stored. The nodal force of node i is the sum of the contributions of all
 * elements attached to node i. The element force of element e contains the
 * element-local contribution K_e * phi_{e} in the element DoF order.
 *
 * @author tlc
 * @date 11/11/2026
 * @version 0.1.0
 * @file ModalForceRecorder.h
 * @addtogroup Recorder
 * @{
 */

#ifndef MODALFORCERECORDER_H
#define MODALFORCERECORDER_H

#include <Recorder/Recorder.h>

class ModalForceRecorder final : public Recorder {
    vec frequency;      // natural frequencies
    mat mode_shape;     // mass normalized mode shapes
    mat force_response; // assembled modal elastic force

    std::vector<std::pair<unsigned, uvec>> node_dof;     // node tag and corresponding DoF indices
    std::vector<std::pair<unsigned, mat>> element_force; // element tag and element force matrix

    void record_impl(const shared_ptr<DomainBase>&) override;

public:
    explicit ModalForceRecorder(unsigned = 0, bool = true);

    void save() override;

    void print() override;
};

#endif

//! @}
