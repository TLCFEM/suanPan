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
 * @class Converger
 * @brief The Converger class handles converger test to indicate if the
 * iteration converges according to various rules.
 *
 * The class stores a pointer `factory` pointed to the Workshop and get
 * information from this Workshop. The `tolerance` and `error` are stored
 * independently so that the Workshop will not be modified.
 *
 * The class further provides a `print_flag` to indicate if the test information
 * should be printed out.
 *
 * @author tlc
 * @date 15/09/2017
 * @version 0.2.1
 * @file Converger.h
 * @addtogroup Converger
 * @{
 */

#ifndef CONVERGER_H
#define CONVERGER_H

#include <Domain/Tag.h>

class DomainBase;

class Converger : public CopiableTag {
    std::weak_ptr<DomainBase> database; /**< pointer to DomainBase */

    double tolerance; /**< tolerance */

    unsigned max_iteration; /**< max iteration number */

    const bool print_flag; /**< print flag */

    double error = 0.; /**< current error */

    bool conv_flag = false; /**< converger flag */
protected:
    [[nodiscard]] vec get_residual() const;

    [[nodiscard]] bool is_print() const;

public:
    explicit Converger(unsigned = 0, double = 1E-8, unsigned = 10, bool = false);

    virtual int initialize();

    virtual unique_ptr<Converger> get_copy() = 0;

    void set_tolerance(double);
    [[nodiscard]] double get_tolerance() const;

    void set_max_iteration(unsigned);
    [[nodiscard]] unsigned get_max_iteration() const;

    void set_domain(const std::weak_ptr<DomainBase>&);
    [[nodiscard]] const std::weak_ptr<DomainBase>& get_domain() const;

    virtual void set_error(double);
    [[nodiscard]] double get_error() const;

    virtual void set_conv_flag(bool);
    [[nodiscard]] bool get_conv_flag() const;

    virtual bool is_converged(unsigned) = 0;
};

#endif

//! @}
