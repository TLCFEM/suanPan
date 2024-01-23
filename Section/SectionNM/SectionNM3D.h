/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class SectionNM3D
 * @brief A SectionNM3D class.
 * @author tlc
 * @date 26/11/2021
 * @version 0.1.0
 * @file SectionNM3D.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef SECTIONNM3D_H
#define SECTIONNM3D_H

#include <Section/SectionNM/SectionNM.h>

struct DataSectionNM3D {
    const double EA;
    const double EIS;
    const double EIW;
};

class SectionNM3D : protected DataSectionNM3D, public SectionNM {
public:
    SectionNM3D(
        unsigned, // tag
        double, double, double, double
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
