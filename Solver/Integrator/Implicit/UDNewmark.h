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
 * @class UDNewmark
 *
 * @author tlc
 * @date 10/05/2026
 * @version 0.1.0
 * @file UDNewmark.h
 * @addtogroup Integrator
 * @{
 */

#ifndef UDNEWMARK_H
#define UDNEWMARK_H

#include "Newmark.h"

class UDNewmark : public Newmark {
    const cx_vec m, s;

protected:
    cx_vec s_para, m_para;

    cx_mat current_nonviscous;

    double aux_para;
    double accu_para{0.};

    void update_parameter(double) override;

    [[nodiscard]] virtual vec target_field() const = 0;

public:
    UDNewmark(unsigned, double, double, cx_vec&&, cx_vec&&);

    int initialize() override;

    void commit_status() override;
    void clear_status() override;

    void print() override;
};

class UDDNewmark final : public UDNewmark {
protected:
    [[nodiscard]] vec target_field() const override;

public:
    using UDNewmark::UDNewmark;

    void assemble_resistance() override;
    void assemble_effective_matrix() override;
};

class UDANewmark final : public UDNewmark {
protected:
    [[nodiscard]] vec target_field() const override;

public:
    using UDNewmark::UDNewmark;

    void assemble_resistance() override;
    void assemble_effective_matrix() override;
};

#endif

//! @}
