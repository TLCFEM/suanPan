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
 * @class ModalRecorder
 * @brief A ModalRecorder class.
 *
 * The ModalRecorder extracts the modal properties from a generalised eigenvalue
 * analysis. For the generalised eigenproblem
 *
 *     K * phi = lambda * M * phi
 *
 * the following quantities are computed for each mode:
 *
 * 1. natural frequency f = sqrt(lambda) / (2 * pi),
 * 2. modal participation factor gamma = phi' * M * r / (phi' * M * phi),
 * 3. effective modal mass m_eff = gamma^2 * (phi' * M * phi),
 * 4. mode shape phi,
 * 5. force response f = lambda * M * phi, equivalent to K * phi.
 *
 * The mode shapes are normalised with respect to the mass matrix, i.e.
 * phi' * M * phi = 1. The influence vector r is a unit vector in each of the
 * requested degree of freedom directions. At least one direction must be
 * specified. Participation factors and effective masses are stored as matrices
 * with one row per mode and one column per direction.
 *
 * @author tlc
 * @date 10/11/2026
 * @version 0.1.0
 * @file ModalRecorder.h
 * @addtogroup Recorder
 * @{
 */

#ifndef MODALRECORDER_H
#define MODALRECORDER_H

#include <Domain/Node.h>
#include <Recorder/Recorder.h>

class ModalRecorder final : public Recorder {
    std::vector<Node::DOF> directions;

    vec frequency;            // natural frequencies
    mat participation_factor; // modal participation factors
    mat effective_mass;       // effective modal masses
    mat mode_shape;           // mass normalized mode shapes
    mat force_response;       // modal force response

    std::vector<std::pair<unsigned, uvec>> node_dof; // node tag and corresponding DoF indices

    void record_impl(const shared_ptr<DomainBase>&) override;

public:
    explicit ModalRecorder(unsigned = 0, std::vector<Node::DOF>&& = {}, bool = true);

    void save() override;

    void print() override;

    [[nodiscard]] const vec& get_frequency() const;
    [[nodiscard]] const mat& get_participation_factor() const;
    [[nodiscard]] const mat& get_effective_mass() const;
    [[nodiscard]] const mat& get_mode_shape() const;
    [[nodiscard]] const mat& get_force_response() const;
};

#endif

//! @}
