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
 * @class Amplitude
 * @brief An Amplitude class that can generate Amplitude pattern.
 *
 * Similar to the strategy used in ABAQUS, the load pattern and the magnitude of
 * that load can be separated. The Amplitude class is designed to generate
 * different Amplitude patterns. The time and amplitude are stored in two
 * different vectors: `TIME` and `AMP`. By calling the method
 *
 *     double get_amplitude(const double& NEW_TIME_POINT);
 *
 * The Amplitude is interpolated. If the given `NEW_TIME_POINT` is beyond the
 * right bound, the Amplitude will be set to the maximum one. Such a case may
 * happen to the auto step size control scheme, in which the (pseudo) time may
 * exceed the given bound due to the machine error.
 *
 * Currently, there are five predefined amplitude patterns.
 * - Tabular
 * - Linear/Ramp
 *   \f{gather}{a=A_0+A\left(t-t_0\right)/t_s\f}
 * - Periodic
 *   \f{gather}{a=A_0+\sum_{n=1}^N\left(A_n\cos{}n\omega\left(t-t_0\right)+B_n\sin{}n\omega\left(t-t_0\right)\right)\f}
 * - Modulated
 *   \f{gather}{a=A_0+A\sin\omega_1\left(t-t_0\right)\sin\omega_2\left(t-t_0\right)\f}
 * - Decay
 *   \f{gather}{a=A_0+A\exp\left(-\left(t-t_0\right)/t_d\right)\f}
 *
 * @author tlc
 * @date 17/09/2017
 * @version 0.1.0
 * @file Amplitude.h
 * @addtogroup Amplitude
 * @{
 */

#ifndef AMPLITUDE_H
#define AMPLITUDE_H

enum class AmplitudeType {
    RAMP,
    TABULAR,
    PERIODIC,
    MODULATED,
    DECAY
};

#include <Domain/Tag.h>

class DomainBase;

class Amplitude : public Tag {
protected:
    unsigned start_step;
    double start_time = 0.; // T0
public:
    explicit Amplitude(unsigned = 0, unsigned = 0);
    ~Amplitude() override = default;

    virtual void initialize(const shared_ptr<DomainBase>&);

    virtual double get_amplitude(double);

    void set_start_step(unsigned);
    [[nodiscard]] unsigned get_start_step() const;

    void set_start_time(double);
    [[nodiscard]] double set_start_time() const;
};

#endif

//! @}
