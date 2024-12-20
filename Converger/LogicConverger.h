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
 * @class LogicConverger
 * @brief The LogicConverger class that handles converger test to indicate if the
 * iteration converges.
 * @author tlc
 * @date 26/03/2022
 * @version 0.1.0
 * @file LogicConverger.h
 * @addtogroup Converger
 * @{
 */

#ifndef LOGICCONVERGER_H
#define LOGICCONVERGER_H

#include "Converger.h"

class LogicConverger : public Converger {
    const unsigned tag_a, tag_b;

protected:
    shared_ptr<Converger> converger_a, converger_b;

public:
    LogicConverger(unsigned, unsigned, unsigned);

    int initialize() override;
};

class LogicAND final : public LogicConverger {
public:
    using LogicConverger::LogicConverger;

    unique_ptr<Converger> get_copy() override;

    bool is_converged(unsigned) override;
};

class LogicOR final : public LogicConverger {
public:
    using LogicConverger::LogicConverger;

    unique_ptr<Converger> get_copy() override;

    bool is_converged(unsigned) override;
};

class LogicXOR final : public LogicConverger {
public:
    using LogicConverger::LogicConverger;

    unique_ptr<Converger> get_copy() override;

    bool is_converged(unsigned) override;
};

#endif

//! @}
