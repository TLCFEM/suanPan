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
 * @class SectionNM2D
 * @brief A SectionNM2D class.
 * @author tlc
 * @date 26/11/2021
 * @version 0.1.0
 * @file SectionNM2D.h
 * @addtogroup Section-NM
 * @ingroup Section
 * @{
 */

#ifndef SECTIONNM2D_H
#define SECTIONNM2D_H

#include <Section/SectionNM/SectionNM.h>

struct DataSectionNM2D {
    const double EA;
    const double EIS;
};

class SectionNM2D : protected DataSectionNM2D, public SectionNM {
public:
    SectionNM2D(
        unsigned, // tag
        double,
        double,
        double
    );

    int initialize(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
