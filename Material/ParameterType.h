/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#ifndef PARAMETERTYPE_H
#define PARAMETERTYPE_H

enum class ParameterType {
    ELASTICMODULUS,
    POISSONSRATIO,
    SHEARMODULUS,
    BULKMODULUS,
    PEAKSTRAIN,
    CRACKSTRAIN
};

struct material_property {
    const double elastic_modulus, poissons_ratio;

    material_property(const double E, const double P)
        : elastic_modulus(E)
        , poissons_ratio(P) {}

    double operator()(const ParameterType P) const {
        switch(P) {
        case ParameterType::ELASTICMODULUS:
            return elastic_modulus;
        case ParameterType::POISSONSRATIO:
            return poissons_ratio;
        case ParameterType::SHEARMODULUS:
            return elastic_modulus / (2. + 2. * poissons_ratio);
        case ParameterType::BULKMODULUS:
            return elastic_modulus / (3. - 6. * poissons_ratio);
        case ParameterType::PEAKSTRAIN:
        case ParameterType::CRACKSTRAIN:
        default:
            return 0.;
        }
    }
};

#endif
