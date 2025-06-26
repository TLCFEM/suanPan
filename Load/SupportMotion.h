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
/**
 * @class SupportMotion
 * @brief A SupportAcceleration class.
 *
 * The SupportAcceleration class is in charge of handling displacement load.
 *
 * @author tlc
 * @date 01/04/2021
 * @version 0.1.0
 * @file SupportMotion.h
 * @addtogroup Load
 * @{
 */

#ifndef SUPPORTMOTION_H
#define SUPPORTMOTION_H

#include <Load/Load.h>

class SupportMotion : public Load {
protected:
    uvec encoding;

public:
    SupportMotion(
        unsigned, // tag
        double,   // magnitude
        uvec&&,   // node tags
        uvec&&,   // dof tags
        unsigned  // amplitude tag
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

class SupportDisplacement final : public SupportMotion {
public:
    using SupportMotion::SupportMotion;

    int process(const shared_ptr<DomainBase>&) override;
};

class SupportVelocity final : public SupportMotion {
public:
    using SupportMotion::SupportMotion;

    int process(const shared_ptr<DomainBase>&) override;
};

class SupportAcceleration final : public SupportMotion {
public:
    using SupportMotion::SupportMotion;

    int process(const shared_ptr<DomainBase>&) override;
};

#endif // SUPPORTMOTION_H

//! @}
