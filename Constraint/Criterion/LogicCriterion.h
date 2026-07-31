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
 * @class LogicCriterion
 * @brief A LogicCriterion class.
 *
 * The LogicCriterion class.
 *
 * @author tlc
 * @date 08/04/2022
 * @version 0.1.0
 * @file LogicCriterion.h
 * @addtogroup Criterion
 * @{
 */

#ifndef LOGICCRITERION_H
#define LOGICCRITERION_H

#include "Criterion.h"

class LogicCriterion : public Criterion {
    const std::vector<unsigned> tags;

    [[nodiscard]] virtual int check(const std::vector<int>&) const = 0;

protected:
    std::vector<shared_ptr<Criterion>> criteria;

public:
    LogicCriterion(unsigned, std::vector<unsigned>&&);

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) final;
};

class LogicCriterionAND final : public LogicCriterion {
    [[nodiscard]] int check(const std::vector<int>&) const override;

public:
    using LogicCriterion::LogicCriterion;

    unique_ptr<Criterion> unique_copy() override;
};

class LogicCriterionOR final : public LogicCriterion {
    [[nodiscard]] int check(const std::vector<int>&) const override;

public:
    using LogicCriterion::LogicCriterion;

    unique_ptr<Criterion> unique_copy() override;
};

#endif

//! @}
