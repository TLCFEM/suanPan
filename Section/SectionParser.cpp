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

#include "SectionParser.h"

#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Section/Section>
#include <Toolbox/TokenRegistry.h>

namespace {
    void new_cell2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double area;
        if(!get_input(command, area)) {
            suanpan_error("A valid area is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Cell2D>(tag, area, material_id, eccentricity);
    }

    void new_cell3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double area;
        if(!get_input(command, area)) {
            suanpan_error("A valid area is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto eccentricity_a = 0., eccentricity_b = 0.;
        if(!command.eof() && !get_input(command, eccentricity_a)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }
        if(!command.eof() && !get_input(command, eccentricity_b)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Cell3D>(tag, area, material_id, eccentricity_a, eccentricity_b);
    }

    void new_box2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double width;
        if(!get_input(command, width)) {
            suanpan_error("A valid width is required.\n");
            return;
        }

        double height;
        if(!get_input(command, height)) {
            suanpan_error("A valid height is required.\n");
            return;
        }

        double thickness;
        if(!get_input(command, thickness)) {
            suanpan_error("A valid thickness is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(command.eof())
            suanpan_debug("Six integration points assumed.\n");
        else if(!get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Box2D>(tag, width, height, thickness, material_id, int_pt, eccentricity);
    }

    void new_box3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double width;
        if(!get_input(command, width)) {
            suanpan_error("A valid width is required.\n");
            return;
        }

        double height;
        if(!get_input(command, height)) {
            suanpan_error("A valid height is required.\n");
            return;
        }

        double thickness;
        if(!get_input(command, thickness)) {
            suanpan_error("A valid thickness is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(command.eof())
            suanpan_debug("Six integration points assumed.\n");
        else if(!get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_a = 0., eccentricity_b = 0.;
        if(!command.eof() && !get_input(command, eccentricity_a)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }
        if(!command.eof() && !get_input(command, eccentricity_b)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Box3D>(tag, width, height, thickness, material_id, int_pt, eccentricity_a, eccentricity_b);
    }

    void new_cell3dos(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double area, omega, py, pz;
        if(!get_input(command, area, omega, py, pz)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto eccentricity_a = 0., eccentricity_b = 0.;
        if(!command.eof() && !get_input(command, eccentricity_a)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }
        if(!command.eof() && !get_input(command, eccentricity_b)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Cell3DOS>(tag, area, omega, py, pz, material_id, eccentricity_a, eccentricity_b);
    }

    void new_circle1d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Circle1D>(tag, radius, material_id);
    }

    void new_circle2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Circle2D>(tag, radius, material_id, int_pt, eccentricity);
    }

    void new_circle3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_y = 0.;
        if(!command.eof() && !get_input(command, eccentricity_y)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }
        auto eccentricity_z = 0.;
        if(!command.eof() && !get_input(command, eccentricity_z)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Circle3D>(tag, radius, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
    }

    void new_circularhollow2D(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        double thickness;
        if(!get_input(command, thickness)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<CircularHollow2D>(tag, radius, thickness, material_id, int_pt, eccentricity);
    }

    void new_circularhollow3D(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double radius;
        if(!get_input(command, radius)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        double thickness;
        if(!get_input(command, thickness)) {
            suanpan_error("A valid radius is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_y = 0.;
        if(!command.eof() && !get_input(command, eccentricity_y)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }
        auto eccentricity_z = 0.;
        if(!command.eof() && !get_input(command, eccentricity_z)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<CircularHollow3D>(tag, radius, thickness, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
    }

    void new_fibre1d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        return_obj = std::make_unique<Fibre1D>(tag, get_remaining<uword>(command));
    }

    void new_fibre2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        return_obj = std::make_unique<Fibre2D>(tag, get_remaining<uword>(command));
    }

    template<bool if_os> void new_fibre3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        if constexpr(if_os) return_obj = std::make_unique<Fibre3DOS>(tag, get_remaining<uword>(command));
        else return_obj = std::make_unique<Fibre3D>(tag, get_remaining<uword>(command));
    }

    void new_hsection2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        vec dim(6);
        if(!get_input(command, dim)) {
            suanpan_error("A valid dimension is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<HSection2D>(tag, dim(0), dim(1), dim(2), dim(3), dim(4), dim(5), material_id, int_pt, eccentricity);
    }

    double barycenter(const vec& dim) {
        if(4u == dim.n_elem) {
            // dim(0): flange width
            // dim(1): flange thickness
            // dim(2): web height
            // dim(3): web thickness
            const auto flange_area = dim(0) * dim(1);
            return -.5 * flange_area * (dim(1) + dim(2)) / (flange_area + dim(2) * dim(3));
        }

        // dim(0): top flange width
        // dim(1): top flange thickness
        // dim(2): bottom flange width
        // dim(3): bottom flange thickness
        // dim(4): web height
        // dim(5): web thickness
        const auto top_flange_area = dim(0) * dim(1);
        const auto bottom_flange_area = dim(2) * dim(3);
        return -.5 * (top_flange_area * (dim(1) + dim(4)) + bottom_flange_area * (dim(3) + dim(4))) / (top_flange_area + bottom_flange_area + dim(4) * dim(5));
    }

    template<bool recenter> void new_isection2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        vec dim(6);
        if(!get_input(command, dim)) {
            suanpan_error("A valid dimension is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if constexpr(recenter) eccentricity = barycenter(dim);
        else if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<ISection2D>(tag, dim(0), dim(1), dim(2), dim(3), dim(4), dim(5), material_id, int_pt, eccentricity);
    }

    template<bool recenter> void new_isection3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        vec dim(6);
        if(!get_input(command, dim)) {
            suanpan_error("A valid dimension is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_y = 0., eccentricity_z = 0.;
        if constexpr(recenter) eccentricity_y = barycenter(dim);
        else {
            if(!command.eof() && !get_input(command, eccentricity_y)) {
                suanpan_error("A valid eccentricity is required.\n");
                return;
            }
            if(!command.eof() && !get_input(command, eccentricity_z)) {
                suanpan_error("A valid eccentricity is required.\n");
                return;
            }
        }

        return_obj = std::make_unique<ISection3D>(tag, std::move(dim), material_id, int_pt, vec{eccentricity_y, eccentricity_z});
    }

    void new_rectangle1d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double width;
        if(!get_input(command, width)) {
            suanpan_error("A valid width is required.\n");
            return;
        }

        double height;
        if(!get_input(command, height)) {
            suanpan_error("A valid height is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Rectangle1D>(tag, width, height, material_id);
    }

    void new_rectangle2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double width;
        if(!get_input(command, width)) {
            suanpan_error("A valid width is required.\n");
            return;
        }

        double height;
        if(!get_input(command, height)) {
            suanpan_error("A valid height is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(command.eof())
            suanpan_debug("Six integration points assumed.\n");
        else if(!get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Rectangle2D>(tag, width, height, material_id, int_pt, eccentricity);
    }

    void new_rectangle3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double width;
        if(!get_input(command, width)) {
            suanpan_error("A valid width is required.\n");
            return;
        }

        double height;
        if(!get_input(command, height)) {
            suanpan_error("A valid height is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(command.eof())
            suanpan_debug("Six integration points assumed.\n");
        else if(!get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_a = 0., eccentricity_b = 0.;
        if(!command.eof() && !get_input(command, eccentricity_a)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }
        if(!command.eof() && !get_input(command, eccentricity_b)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<Rectangle3D>(tag, width, height, material_id, int_pt, eccentricity_a, eccentricity_b);
    }

    void new_trusssection(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        double area;
        if(!get_input(command, area)) {
            suanpan_error("A valid area is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        return_obj = std::make_unique<TrussSection>(tag, area, material_id);
    }

    template<bool recenter> void new_tsection2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        vec dim(4);
        if(!get_input(command, dim)) {
            suanpan_error("A valid dimension is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if constexpr(recenter) eccentricity = barycenter(dim);
        else if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        return_obj = std::make_unique<TSection2D>(tag, dim(0), dim(1), dim(2), dim(3), material_id, int_pt, eccentricity);
    }

    template<bool recenter> void new_tsection3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        vec dim(4);
        if(!get_input(command, dim)) {
            suanpan_error("A valid dimension is required.\n");
            return;
        }

        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_y = 0., eccentricity_z = 0.;
        if constexpr(recenter) eccentricity_y = barycenter(dim);
        else {
            if(!command.eof() && !get_input(command, eccentricity_y)) {
                suanpan_error("A valid eccentricity is required.\n");
                return;
            }
            if(!command.eof() && !get_input(command, eccentricity_z)) {
                suanpan_error("A valid eccentricity is required.\n");
                return;
            }
        }

        return_obj = std::make_unique<TSection3D>(tag, std::move(dim), material_id, int_pt, vec{eccentricity_y, eccentricity_z});
    }

    template<unsigned size> void new_nm2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        vec P(size);
        if(!get_input(command, P)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        if constexpr(3u == size) {
            return_obj = std::make_unique<NM2D1>(tag, P(0), P(1), P(2));
            return;
        }

        mat poly_set;
        if(const auto para_set = get_remaining<double>(command); para_set.size() > 0) {
            if(para_set.size() % 3 != 0) {
                suanpan_error("A valid parameter set is required.\n");
                return;
            }

            inplace_trans((poly_set = para_set).reshape(3, poly_set.n_elem / 3));
        }

        if constexpr(8u == size) return_obj = std::make_unique<NM2D2>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), std::move(poly_set));
        else if constexpr(11u == size) return_obj = std::make_unique<NM2D3>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), vec{P(8), P(8)}, vec{P(9), P(9)}, P(10), std::move(poly_set));
    }

    template<unsigned size> void new_nm3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        vec P(size);
        if(!get_input(command, P)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        if constexpr(4u == size) {
            return_obj = std::make_unique<NM3D1>(tag, P(0), P(1), P(2), P(3));
            return;
        }

        mat poly_set;
        if(const auto para_set = get_remaining<double>(command); para_set.size() > 0) {
            if(para_set.size() % 4 != 0) {
                suanpan_error("A valid parameter set is required.\n");
                return;
            }

            inplace_trans((poly_set = para_set).reshape(4, poly_set.n_elem / 4));
        }

        if constexpr(10u == size) return_obj = std::make_unique<NM3D2>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), std::move(poly_set));
        else if constexpr(13u == size) return_obj = std::make_unique<NM3D3>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), vec{P(10), P(10), P(10)}, vec{P(11), P(11), P(11)}, P(12), std::move(poly_set));
    }

    template<unsigned size> void new_nmk(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string&) {
        vec P(size);
        if(!get_input(command, P)) {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

        mat poly_set;
        if(const auto para_set = get_remaining<double>(command); para_set.size() > 0) {
            const auto p_size = 13u == size ? 3 : 4;

            if(para_set.size() % p_size != 0) {
                suanpan_error("A valid parameter set is required.\n");
                return;
            }

            inplace_trans((poly_set = para_set).reshape(p_size, poly_set.n_elem / p_size));
        }

        if constexpr(13u == size) return_obj = std::make_unique<NM2D3>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), vec{P(8), P(9)}, vec{P(10), P(11)}, P(12), std::move(poly_set));
        else if constexpr(17u == size) return_obj = std::make_unique<NM3D3>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), vec{P(10), P(11), P(12)}, vec{P(13), P(14), P(15)}, P(16), std::move(poly_set));
    }

    vec euisection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"HEA100", {100.0, 8.0, 100.0, 8.0, 80.0, 5.0}},
            {"HEA120", {120.0, 8.0, 120.0, 8.0, 98.0, 5.0}},
            {"HEA140", {140.0, 8.5, 140.0, 8.5, 116.0, 5.5}},
            {"HEA160", {160.0, 9.0, 160.0, 9.0, 134.0, 6.0}},
            {"HEA180", {180.0, 9.5, 180.0, 9.5, 152.0, 6.0}},
            {"HEA200", {200.0, 10.0, 200.0, 10.0, 170.0, 6.5}},
            {"HEA220", {220.0, 11.0, 220.0, 11.0, 188.0, 7.0}},
            {"HEA240", {240.0, 12.0, 240.0, 12.0, 206.0, 7.5}},
            {"HEA260", {260.0, 12.5, 260.0, 12.5, 225.0, 7.5}},
            {"HEA280", {280.0, 13.0, 280.0, 13.0, 244.0, 8.0}},
            {"HEA300", {300.0, 14.0, 300.0, 14.0, 262.0, 8.5}},
            {"HEA320", {300.0, 15.5, 300.0, 15.5, 279.0, 9.0}},
            {"HEA340", {300.0, 16.5, 300.0, 16.5, 297.0, 9.5}},
            {"HEA360", {300.0, 17.5, 300.0, 17.5, 315.0, 10.0}},
            {"HEA400", {300.0, 19.0, 300.0, 19.0, 352.0, 11.0}},
            {"HEA450", {300.0, 21.0, 300.0, 21.0, 398.0, 11.5}},
            {"HEA500", {300.0, 23.0, 300.0, 23.0, 444.0, 12.0}},
            {"HEA550", {300.0, 24.0, 300.0, 24.0, 492.0, 12.5}},
            {"HEA600", {300.0, 25.0, 300.0, 25.0, 540.0, 13.0}},
            {"HEA650", {300.0, 26.0, 300.0, 26.0, 588.0, 13.5}},
            {"HEA700", {300.0, 27.0, 300.0, 27.0, 636.0, 14.5}},
            {"HEA800", {300.0, 28.0, 300.0, 28.0, 734.0, 15.0}},
            {"HEA900", {300.0, 30.0, 300.0, 30.0, 830.0, 16.0}},
            {"HEA1000", {300.0, 31.0, 300.0, 31.0, 928.0, 16.5}},
            {"HEB100", {100.0, 10.0, 100.0, 10.0, 80.0, 6.0}},
            {"HEB120", {120.0, 11.0, 120.0, 11.0, 98.0, 6.5}},
            {"HEB140", {140.0, 12.0, 140.0, 12.0, 116.0, 7.0}},
            {"HEB160", {160.0, 13.0, 160.0, 13.0, 134.0, 8.0}},
            {"HEB180", {180.0, 14.0, 180.0, 14.0, 152.0, 8.5}},
            {"HEB200", {200.0, 15.0, 200.0, 15.0, 170.0, 9.0}},
            {"HEB220", {220.0, 16.0, 220.0, 16.0, 188.0, 9.5}},
            {"HEB240", {240.0, 17.0, 240.0, 17.0, 206.0, 10.0}},
            {"HEB260", {260.0, 17.5, 260.0, 17.5, 225.0, 10.0}},
            {"HEB280", {280.0, 18.0, 280.0, 18.0, 244.0, 10.5}},
            {"HEB300", {300.0, 19.0, 300.0, 19.0, 262.0, 11.0}},
            {"HEB320", {300.0, 20.5, 300.0, 20.5, 279.0, 11.5}},
            {"HEB340", {300.0, 21.5, 300.0, 21.5, 297.0, 12.0}},
            {"HEB360", {300.0, 22.5, 300.0, 22.5, 315.0, 12.5}},
            {"HEB400", {300.0, 24.0, 300.0, 24.0, 352.0, 13.5}},
            {"HEB450", {300.0, 26.0, 300.0, 26.0, 398.0, 14.0}},
            {"HEB500", {300.0, 28.0, 300.0, 28.0, 444.0, 14.5}},
            {"HEB550", {300.0, 29.0, 300.0, 29.0, 492.0, 15.0}},
            {"HEB600", {300.0, 30.0, 300.0, 30.0, 540.0, 15.5}},
            {"HEB650", {300.0, 31.0, 300.0, 31.0, 588.0, 16.0}},
            {"HEB700", {300.0, 32.0, 300.0, 32.0, 636.0, 17.0}},
            {"HEB800", {300.0, 33.0, 300.0, 33.0, 734.0, 17.5}},
            {"HEB900", {300.0, 35.0, 300.0, 35.0, 830.0, 18.5}},
            {"HEB1000", {300.0, 36.0, 300.0, 36.0, 928.0, 19.0}},
            {"IPE80", {46.0, 5.2, 46.0, 5.2, 69.6, 3.8}},
            {"IPE100", {55.0, 5.7, 55.0, 5.7, 88.6, 4.1}},
            {"IPE120", {64.0, 6.3, 64.0, 6.3, 107.4, 4.4}},
            {"IPE140", {73.0, 6.9, 73.0, 6.9, 126.2, 4.7}},
            {"IPE160", {82.0, 7.4, 82.0, 7.4, 145.2, 5.0}},
            {"IPE180", {91.0, 8.0, 91.0, 8.0, 164.0, 5.3}},
            {"IPE200", {100.0, 8.5, 100.0, 8.5, 183.0, 5.6}},
            {"IPE220", {110.0, 9.2, 110.0, 9.2, 201.6, 5.9}},
            {"IPE240", {120.0, 9.8, 120.0, 9.8, 220.4, 6.2}},
            {"IPE270", {135.0, 10.2, 135.0, 10.2, 249.6, 6.6}},
            {"IPE300", {150.0, 10.7, 150.0, 10.7, 278.6, 7.1}},
            {"IPE330", {160.0, 11.5, 160.0, 11.5, 307.0, 7.5}},
            {"IPE360", {170.0, 12.7, 170.0, 12.7, 334.6, 8.0}},
            {"IPE400", {180.0, 13.5, 180.0, 13.5, 373.0, 8.6}},
            {"IPE450", {190.0, 14.6, 190.0, 14.6, 420.8, 9.4}},
            {"IPE500", {200.0, 16.0, 200.0, 16.0, 468.0, 10.2}},
            {"IPE550", {210.0, 17.2, 210.0, 17.2, 515.6, 11.1}},
            {"IPE600", {220.0, 19.0, 220.0, 19.0, 562.0, 12.0}},
            {"HEM100", {106.0, 20.0, 106.0, 20.0, 80.0, 12.0}},
            {"HEM120", {126.0, 21.0, 126.0, 21.0, 98.0, 12.5}},
            {"HEM140", {146.0, 22.0, 146.0, 22.0, 116.0, 13.0}},
            {"HEM160", {166.0, 23.0, 166.0, 23.0, 134.0, 14.0}},
            {"HEM180", {186.0, 24.0, 186.0, 24.0, 152.0, 14.5}},
            {"HEM200", {206.0, 25.0, 206.0, 25.0, 170.0, 15.0}},
            {"HEM220", {226.0, 26.0, 226.0, 26.0, 188.0, 15.5}},
            {"HEM240", {248.0, 32.0, 248.0, 32.0, 206.0, 18.0}},
            {"HEM260", {268.0, 32.5, 268.0, 32.5, 225.0, 18.0}},
            {"HEM280", {288.0, 33.0, 288.0, 33.0, 244.0, 18.5}},
            {"HEM300", {310.0, 39.0, 310.0, 39.0, 262.0, 21.0}},
            {"HEM320", {309.0, 40.0, 309.0, 40.0, 279.0, 21.0}},
            {"HEM340", {309.0, 40.0, 309.0, 40.0, 297.0, 21.0}},
            {"HEM360", {308.0, 40.0, 308.0, 40.0, 315.0, 21.0}},
            {"HEM400", {307.0, 40.0, 307.0, 40.0, 352.0, 21.0}},
            {"HEM450", {307.0, 40.0, 307.0, 40.0, 398.0, 21.0}},
            {"HEM500", {306.0, 40.0, 306.0, 40.0, 444.0, 21.0}},
            {"HEM550", {306.0, 40.0, 306.0, 40.0, 492.0, 21.0}},
            {"HEM600", {305.0, 40.0, 305.0, 40.0, 540.0, 21.0}},
            {"HEM650", {305.0, 40.0, 305.0, 40.0, 588.0, 21.0}},
            {"HEM700", {304.0, 40.0, 304.0, 40.0, 636.0, 21.0}},
            {"HEM800", {303.0, 40.0, 303.0, 40.0, 734.0, 21.0}},
            {"HEM900", {302.0, 40.0, 302.0, 40.0, 830.0, 21.0}},
            {"HEM1000", {302.0, 40.0, 302.0, 40.0, 928.0, 21.0}},
            {"UB127X76X13", {76, 7.6, 76, 7.6, 111.8, 4}},
            {"UB152X89X16", {88.7, 7.7, 88.7, 7.7, 137, 4.5}},
            {"UB178X102X19", {101.2, 7.9, 101.2, 7.9, 162, 4.8}},
            {"UB203X102X23", {101.8, 9.3, 101.8, 9.3, 184.6, 5.4}},
            {"UB203X133X25", {133.2, 7.8, 133.2, 7.8, 187.6, 5.7}},
            {"UB203X133X30", {133.9, 9.6, 133.9, 9.6, 187.6, 6.4}},
            {"UB254X102X22", {101.6, 6.8, 101.6, 6.8, 240.4, 5.7}},
            {"UB254X102X25", {101.9, 8.4, 101.9, 8.4, 240.4, 6}},
            {"UB254X102X28", {102.2, 10, 102.2, 10, 240.4, 6.3}},
            {"UB254X146X31", {146.1, 8.6, 146.1, 8.6, 234.2, 6}},
            {"UB254X146X37", {146.4, 10.9, 146.4, 10.9, 234.2, 6.3}},
            {"UB254X146X43", {147.3, 12.7, 147.3, 12.7, 234.2, 7.2}},
            {"UB305X102X25", {101.6, 7, 101.6, 7, 291.1, 5.8}},
            {"UB305X102X28", {101.8, 8.8, 101.8, 8.8, 291.1, 6}},
            {"UB305X102X33", {102.4, 10.8, 102.4, 10.8, 291.1, 6.6}},
            {"UB305X127X37", {123.4, 10.7, 123.4, 10.7, 283, 7.1}},
            {"UB305X127X42", {124.3, 12.1, 124.3, 12.1, 283, 8}},
            {"UB305X127X48", {125.3, 14, 125.3, 14, 283, 9}},
            {"UB305X165X40", {165, 10.2, 165, 10.2, 283, 6}},
            {"UB305X165X46", {165.7, 11.8, 165.7, 11.8, 283, 6.7}},
            {"UB305X165X54", {166.9, 13.7, 166.9, 13.7, 283, 7.9}},
            {"UB356X127X33", {125.4, 8.5, 125.4, 8.5, 332, 6}},
            {"UB356X127X39", {126, 10.7, 126, 10.7, 332, 6.6}},
            {"UB356X171X45", {171.1, 9.7, 171.1, 9.7, 332, 7}},
            {"UB356X171X51", {171.5, 11.5, 171.5, 11.5, 332, 7.4}},
            {"UB356X171X57", {172.2, 13, 172.2, 13, 332, 8.1}},
            {"UB356X171X67", {173.2, 15.7, 173.2, 15.7, 332, 9.1}},
            {"UB406X140X39", {141.8, 8.6, 141.8, 8.6, 380.8, 6.4}},
            {"UB406X140X46", {142.2, 11.2, 142.2, 11.2, 380.8, 6.8}},
            {"UB406X140X53", {143.3, 12.9, 143.3, 12.9, 380.8, 7.9}},
            {"UB406X178X54", {177.7, 10.9, 177.7, 10.9, 380.8, 7.7}},
            {"UB406X178X60", {177.9, 12.8, 177.9, 12.8, 380.8, 7.9}},
            {"UB406X178X67", {178.8, 14.3, 178.8, 14.3, 380.8, 8.8}},
            {"UB406X178X74", {179.5, 16, 179.5, 16, 380.8, 9.5}},
            {"UB406X178X85", {181, 18.2, 181, 18.2, 380.6, 10.9}},
            {"UB457X152X52", {152.4, 10.9, 152.4, 10.9, 428, 7.6}},
            {"UB457X152X60", {152.9, 13.3, 152.9, 13.3, 428, 8.1}},
            {"UB457X152X67", {153.8, 15, 153.8, 15, 428, 9}},
            {"UB457X152X74", {154.4, 17, 154.4, 17, 428, 9.6}},
            {"UB457X152X82", {155.3, 18.9, 155.3, 18.9, 428, 10.5}},
            {"UB457X191X67", {189.9, 12.7, 189.9, 12.7, 428, 8.5}},
            {"UB457X191X74", {190.4, 14.5, 190.4, 14.5, 428, 9}},
            {"UB457X191X82", {191.3, 16, 191.3, 16, 428, 9.9}},
            {"UB457X191X89", {191.9, 17.7, 191.9, 17.7, 428, 10.5}},
            {"UB457X191X98", {192.8, 19.6, 192.8, 19.6, 428, 11.4}},
            {"UB457X191X106", {194, 20.6, 194, 20.6, 427.8, 12.6}},
            {"UB533X165X66", {165, 11.4, 165, 11.4, 502.2, 8.9}},
            {"UB533X165X74", {166, 13.6, 166, 13.6, 501.8, 9.7}},
            {"UB533X165X85", {166, 16.5, 166, 16.5, 502, 10.3}},
            {"UB533X210X82", {208.8, 13.2, 208.8, 13.2, 501.9, 9.6}},
            {"UB533X210X92", {209.3, 15.6, 209.3, 15.6, 501.9, 10.1}},
            {"UB533X210X101", {210, 17.4, 210, 17.4, 501.9, 10.8}},
            {"UB533X210X109", {210.8, 18.8, 210.8, 18.8, 501.9, 11.6}},
            {"UB533X210X122", {211.9, 21.3, 211.9, 21.3, 501.9, 12.7}},
            {"UB533X210X138", {214, 23.6, 214, 23.6, 501.8, 14.7}},
            {"UB610X178X82", {178, 12.8, 178, 12.8, 573.4, 10}},
            {"UB610X178X92", {179, 15, 179, 15, 573, 10.9}},
            {"UB610X229X101", {227.6, 14.8, 227.6, 14.8, 573, 10.5}},
            {"UB610X229X113", {228.2, 17.3, 228.2, 17.3, 573, 11.1}},
            {"UB610X229X125", {229, 19.6, 229, 19.6, 573, 11.9}},
            {"UB610X229X140", {230.2, 22.1, 230.2, 22.1, 573, 13.1}},
            {"UB610X229X153", {229, 24.9, 229, 24.9, 573.2, 14}},
            {"UB610X305X149", {304.8, 19.7, 304.8, 19.7, 573, 11.8}},
            {"UB610X305X179", {307.1, 23.6, 307.1, 23.6, 573, 14.1}},
            {"UB610X305X238", {311.4, 31.4, 311.4, 31.4, 573, 18.4}},
            {"UB686X254X125", {253, 16.2, 253, 16.2, 645.5, 11.7}},
            {"UB686X254X140", {253.7, 19, 253.7, 19, 645.5, 12.4}},
            {"UB686X254X152", {254.5, 21, 254.5, 21, 645.5, 13.2}},
            {"UB686X254X170", {255.8, 23.7, 255.8, 23.7, 645.5, 14.5}},
            {"UB686X254X192", {254, 27.9, 254, 27.9, 646.2, 15.5}},
            {"UB762X267X134", {264.4, 15.5, 264.4, 15.5, 719, 12}},
            {"UB762X267X147", {265.2, 17.5, 265.2, 17.5, 719, 12.8}},
            {"UB762X267X173", {266.7, 21.6, 266.7, 21.6, 719, 14.3}},
            {"UB762X267X197", {268, 25.4, 268, 25.4, 719, 15.6}},
            {"UB762X267X220", {266, 30, 266, 30, 719, 16.5}},
            {"UB838X292X176", {291.7, 18.8, 291.7, 18.8, 797.3, 14}},
            {"UB838X292X194", {292.4, 21.7, 292.4, 21.7, 797.3, 14.7}},
            {"UB838X292X226", {293.8, 26.8, 293.8, 26.8, 797.3, 16.1}},
            {"UB838X292X251", {292, 31, 292, 31, 797, 17}},
            {"UB914X305X201", {303.3, 20.2, 303.3, 20.2, 862.6, 15.1}},
            {"UB914X305X224", {304.1, 23.9, 304.1, 23.9, 862.6, 15.9}},
            {"UB914X305X238", {305, 25.9, 305, 25.9, 863.2, 16.5}},
            {"UB914X305X253", {305.5, 27.9, 305.5, 27.9, 862.6, 17.3}},
            {"UB914X305X271", {307, 30, 307, 30, 863, 18.4}},
            {"UB914X305X289", {307.7, 32, 307.7, 32, 862.6, 19.5}},
            {"UB914X305X313", {309, 34.5, 309, 34.5, 863, 21.1}},
            {"UB914X305X345", {308, 39.9, 308, 39.9, 863.2, 22.1}},
            {"UB914X305X381", {310, 43.9, 310, 43.9, 863.2, 24.4}},
            {"UB914X305X425", {313, 49, 313, 49, 863, 26.9}},
            {"UB914X305X474", {316, 54.1, 316, 54.1, 862.8, 30}},
            {"UB914X305X521", {319, 58.9, 319, 58.9, 863.2, 33}},
            {"UB914X305X576", {322, 65, 322, 65, 863, 36.1}},
            {"UB920X420X344", {418, 32, 418, 32, 863, 19.3}},
            {"UB920X420X368", {419, 34.3, 419, 34.3, 862.4, 20.3}},
            {"UB920X420X390", {420, 36.6, 420, 36.6, 862.8, 21.3}},
            {"UB920X420X420", {422, 39.9, 422, 39.9, 863.2, 22.5}},
            {"UB920X420X449", {423, 42.7, 423, 42.7, 862.6, 24}},
            {"UB920X420X491", {422, 47, 422, 47, 863, 25.9}},
            {"UB920X420X537", {425, 51.1, 425, 51.1, 862.8, 28.4}},
            {"UB920X420X588", {427, 55.9, 427, 55.9, 863.2, 31}},
            {"UB920X420X656", {431, 62, 431, 62, 863, 34.5}},
            {"UB920X420X725", {434, 68.1, 434, 68.1, 862.8, 38.1}},
            {"UB920X420X787", {437, 73.9, 437, 73.9, 863.2, 40.9}},
            {"UB920X420X970", {446, 89.9, 446, 89.9, 863.2, 50}},
            {"UB920X420X1077", {451, 99.1, 451, 99.1, 862.8, 55}},
            {"UB920X420X1194", {457, 109, 457, 109, 863, 60.5}},
            {"UB920X420X1269", {461, 115.1, 461, 115.1, 862.8, 64}},
            {"UB920X420X1377", {473, 115.1, 473, 115.1, 862.8, 76.7}},
            {"UB1000X400X296", {400, 27.1, 400, 27.1, 927.8, 16.5}},
            {"UB1000X400X321", {400, 31, 400, 31, 928, 16.5}},
            {"UB1000X400X371", {400, 36.1, 400, 36.1, 927.8, 19}},
            {"UB1000X400X412", {402, 40, 402, 40, 928, 21.1}},
            {"UB1000X400X443", {402, 41.9, 402, 41.9, 928.2, 23.6}},
            {"UB1000X400X483", {404, 46, 404, 46, 928, 25.4}},
            {"UB1000X400X539", {407, 51.1, 407, 51.1, 927.8, 28.4}},
            {"UB1000X400X554", {408, 52, 408, 52, 928, 29.5}},
            {"UB1000X400X591", {409, 55.9, 409, 55.9, 928.2, 31}},
            {"UB1000X400X642", {412, 60, 412, 60, 928, 34}},
            {"UB1000X400X748", {417, 70, 417, 70, 928, 39}},
            {"UB1000X400X883", {424, 82, 424, 82, 928, 45.5}},
            {"UB1000X400X976", {428, 89.9, 428, 89.9, 928.2, 50}},
            {"UB1016X305X222", {300, 21.1, 300, 21.1, 927.8, 16}},
            {"UB1016X305X249", {300, 26, 300, 26, 928, 16.5}},
            {"UB1016X305X272", {300, 31, 300, 31, 928, 16.5}},
            {"UB1016X305X314", {300, 35.9, 300, 35.9, 928.2, 19.1}},
            {"UB1016X305X350", {302, 40, 302, 40, 928, 21.1}},
            {"UB1016X305X393", {303, 43.9, 303, 43.9, 928.2, 24.4}},
            {"UB1016X305X415", {304, 46, 304, 46, 928, 26}},
            {"UB1016X305X438", {305, 49, 305, 49, 928, 26.9}},
            {"UB1016X305X494", {309, 54, 309, 54, 928, 31}},
            {"UB1016X305X584", {314, 64, 314, 64, 928, 36}},
            {"UB1100X400X343", {400, 31, 400, 31, 1028, 18}},
            {"UB1100X400X390", {400, 36, 400, 36, 1028, 20}},
            {"UB1100X400X433", {402, 40, 402, 40, 1028, 22}},
            {"UB1100X400X499", {405, 45, 405, 45, 1028, 26}},
            {"UB1100X400X548", {407, 50, 407, 50, 1028, 28}},
            {"UB1100X400X607", {410, 55, 410, 55, 1028, 31}},
            {"UC152X152X23", {152.2, 6.8, 152.2, 6.8, 138.8, 5.8}},
            {"UC152X152X30", {152.9, 9.4, 152.9, 9.4, 138.8, 6.5}},
            {"UC152X152X37", {154.4, 11.5, 154.4, 11.5, 138.8, 8}},
            {"UC152X152X44", {155.9, 13.6, 155.9, 13.6, 138.8, 9.5}},
            {"UC152X152X51", {157.4, 15.7, 157.4, 15.7, 138.8, 11}},
            {"UC203X203X46", {203.6, 11, 203.6, 11, 181.2, 7.2}},
            {"UC203X203X52", {204.3, 12.5, 204.3, 12.5, 181.2, 7.9}},
            {"UC203X203X60", {205.8, 14.2, 205.8, 14.2, 181.2, 9.4}},
            {"UC203X203X71", {206.4, 17.3, 206.4, 17.3, 181.2, 10}},
            {"UC203X203X86", {209.1, 20.5, 209.1, 20.5, 181.2, 12.7}},
            {"UC203X203X100", {210, 23.7, 210, 23.7, 181.6, 14.5}},
            {"UC254X254X73", {254.6, 14.2, 254.6, 14.2, 225.7, 8.6}},
            {"UC254X254X89", {256.3, 17.3, 256.3, 17.3, 225.7, 10.3}},
            {"UC254X254X107", {258.8, 20.5, 258.8, 20.5, 225.7, 12.8}},
            {"UC254X254X132", {261.3, 25.3, 261.3, 25.3, 225.7, 15.3}},
            {"UC254X254X167", {265.2, 31.7, 265.2, 31.7, 225.7, 19.2}},
            {"UC305X305X97", {305.3, 15.4, 305.3, 15.4, 277.1, 9.9}},
            {"UC305X305X118", {307.4, 18.7, 307.4, 18.7, 277.1, 12}},
            {"UC305X305X137", {309.2, 21.7, 309.2, 21.7, 277.1, 13.8}},
            {"UC305X305X158", {311.2, 25, 311.2, 25, 277.1, 15.8}},
            {"UC305X305X198", {314.5, 31.4, 314.5, 31.4, 277.1, 19.1}},
            {"UC305X305X240", {318.4, 37.7, 318.4, 37.7, 277.1, 23}},
            {"UC305X305X283", {322.2, 44.1, 322.2, 44.1, 277.1, 26.8}},
            {"UC305X305X313", {325, 48.3, 325, 48.3, 277.4, 30}},
            {"UC305X305X342", {328, 52.6, 328, 52.6, 276.8, 32.6}},
            {"UC356X368X129", {368.6, 17.5, 368.6, 17.5, 320.6, 10.4}},
            {"UC356X368X153", {370.5, 20.7, 370.5, 20.7, 320.6, 12.3}},
            {"UC356X368X177", {372.6, 23.8, 372.6, 23.8, 320.6, 14.4}},
            {"UC356X368X202", {374.7, 27, 374.7, 27, 320.6, 16.5}},
            {"UC356X406X235", {394.8, 30.2, 394.8, 30.2, 320.6, 18.4}},
            {"UC356X406X287", {399, 36.5, 399, 36.5, 320.6, 22.6}},
            {"UC356X406X340", {403, 42.9, 403, 42.9, 320.6, 26.6}},
            {"UC356X406X393", {407, 49.2, 407, 49.2, 320.6, 30.6}},
            {"UC356X406X467", {412.2, 58, 412.2, 58, 320.6, 35.8}},
            {"UC356X406X509", {416, 62.7, 416, 62.7, 320.6, 39.1}},
            {"UC356X406X551", {418.5, 67.5, 418.5, 67.5, 320.6, 42.1}},
            {"UC356X406X592", {421, 72.3, 421, 72.3, 320.4, 45}},
            {"UC356X406X634", {424, 77, 424, 77, 320.6, 47.6}},
            {"UC356X406X677", {428, 81.5, 428, 81.5, 320, 51.2}},
            {"UC356X406X744", {432, 88.9, 432, 88.9, 320.2, 55.6}},
            {"UC356X406X818", {437, 97, 437, 97, 320, 60.5}},
            {"UC356X406X900", {442, 106, 442, 106, 319, 65.9}},
            {"UC356X406X990", {448, 115, 448, 115, 320, 71.9}},
            {"UC356X406X1086", {454, 125, 454, 125, 319, 78}},
            {"UC356X406X1202", {471, 130, 471, 130, 320, 95}},
            {"UC356X406X1299", {476, 140, 476, 140, 320, 100}},
            {"UBP203X203X45", {205.9, 9.5, 205.9, 9.5, 181.2, 9.5}},
            {"UBP203X203X54", {207.7, 11.4, 207.7, 11.4, 181.2, 11.3}},
            {"UBP254X254X63", {256.6, 10.7, 256.6, 10.7, 225.7, 10.6}},
            {"UBP254X254X71", {258, 12, 258, 12, 225.7, 12}},
            {"UBP254X254X85", {260.4, 14.3, 260.4, 14.3, 225.7, 14.4}},
            {"UBP305X305X79", {306.4, 11.1, 306.4, 11.1, 277.1, 11}},
            {"UBP305X305X88", {307.8, 12.3, 307.8, 12.3, 277.1, 12.4}},
            {"UBP305X305X95", {308.7, 13.3, 308.7, 13.3, 277.1, 13.3}},
            {"UBP305X305X110", {310.7, 15.4, 310.7, 15.4, 277.1, 15.3}},
            {"UBP305X305X126", {312.9, 17.6, 312.9, 17.6, 277.1, 17.5}},
            {"UBP305X305X149", {316, 20.7, 316, 20.7, 277.1, 20.6}},
            {"UBP305X305X186", {320.9, 25.6, 320.9, 25.6, 277.1, 25.5}},
            {"UBP305X305X223", {325.7, 30.4, 325.7, 30.4, 277.1, 30.3}},
            {"UBP356X368X109", {371, 12.9, 371, 12.9, 320.6, 12.8}},
            {"UBP356X368X133", {373.8, 15.7, 373.8, 15.7, 320.6, 15.6}},
            {"UBP356X368X152", {376, 17.9, 376, 17.9, 320.6, 17.8}},
            {"UBP356X368X174", {378.5, 20.4, 378.5, 20.4, 320.6, 20.3}},
        };

        return registry.fetch(type);
    }

    // Adjusted registry: dimensions optimized against the EN 10365:2017
    // tabulated section properties with eu_registry_generator.py.
    vec euisection_adjusted(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"HEA100", {97.691, 8.568, 97.691, 8.568, 79.6896, 5.5}},
            {"HEA120", {117.9505, 8.4099, 117.9505, 8.4099, 98.1606, 5.5}},
            {"HEA140", {138.417, 8.7843, 138.417, 8.7843, 116.4839, 6.05}},
            {"HEA160", {157.1275, 9.4804, 157.1275, 9.4804, 134.385, 6.6}},
            {"HEA180", {177.4459, 9.8953, 177.4459, 9.8953, 152.5743, 6.6}},
            {"HEA200", {196.3389, 10.5449, 196.3389, 10.5449, 170.7469, 7.15}},
            {"HEA220", {217.0883, 11.4293, 217.0883, 11.4293, 189.0169, 7.7}},
            {"HEA240", {235.8553, 12.612, 235.8553, 12.612, 206.9081, 8.25}},
            {"HEA260", {254.1675, 13.3329, 254.1675, 13.3329, 225.7152, 8.25}},
            {"HEA280", {274.7888, 13.7119, 274.7888, 13.7119, 245.0676, 8.8}},
            {"HEA300", {293.5974, 14.8848, 293.5974, 14.8848, 263.0602, 9.35}},
            {"HEA320", {294.5305, 16.3334, 294.5305, 16.3334, 280.2272, 9.9}},
            {"HEA340", {295.327, 17.2591, 295.327, 17.2591, 298.5394, 10.45}},
            {"HEA360", {296.03, 18.1865, 296.03, 18.1865, 316.8691, 11}},
            {"HEA400", {297.272, 19.5197, 297.272, 19.5197, 354.6969, 12.1}},
            {"HEA450", {298.3663, 21.354, 298.3663, 21.354, 401.5175, 12.65}},
            {"HEA500", {298.4433, 23.3776, 298.4433, 23.3776, 447.4821, 13.0017}},
            {"HEA550", {298.5961, 24.3522, 298.5961, 24.3522, 495.7862, 13.432}},
            {"HEA600", {298.5439, 25.3765, 298.5439, 25.3765, 543.9308, 13.8233}},
            {"HEA650", {298.553, 26.3841, 298.553, 26.3841, 592.2004, 14.2303}},
            {"HEA700", {298.7392, 27.3603, 298.7392, 27.3603, 640.3927, 15.1855}},
            {"HEA800", {298.5747, 28.4297, 298.5747, 28.4297, 740.0272, 15.7171}},
            {"HEA900", {298.769, 30.4044, 298.769, 30.4044, 836.5271, 16.6276}},
            {"HEA1000", {298.4732, 31.4976, 298.4732, 31.4976, 934.3554, 17.0193}},
            {"HEB100", {98.3149, 10.5214, 98.3149, 10.5214, 79.6808, 6.6}},
            {"HEB120", {118.723, 11.353, 118.723, 11.353, 98.0659, 7.15}},
            {"HEB140", {139.0868, 12.2349, 139.0868, 12.2349, 116.3071, 7.7}},
            {"HEB160", {158.368, 13.3974, 158.368, 13.3974, 134.2476, 8.8}},
            {"HEB180", {178.8755, 14.2625, 178.8755, 14.2625, 152.5877, 9.35}},
            {"HEB200", {198.0148, 15.441, 198.0148, 15.441, 170.4936, 9.9}},
            {"HEB220", {218.5923, 16.3049, 218.5923, 16.3049, 188.8673, 10.45}},
            {"HEB240", {237.6463, 17.4952, 237.6463, 17.4952, 206.7062, 11}},
            {"HEB260", {256.5056, 18.1978, 256.5056, 18.1978, 225.6217, 11}},
            {"HEB280", {276.93, 18.583, 276.93, 18.583, 244.8957, 11.55}},
            {"HEB300", {295.8913, 19.7685, 295.8913, 19.7685, 262.8302, 12.1}},
            {"HEB320", {296.6303, 21.183, 296.6303, 21.183, 280.1067, 12.65}},
            {"HEB340", {297.1553, 22.1076, 297.1553, 22.1076, 298.4484, 13.2}},
            {"HEB360", {297.7069, 23.0123, 297.7069, 23.0123, 316.8539, 13.75}},
            {"HEB400", {298.5746, 24.3565, 298.5746, 24.3565, 354.53, 14.8271}},
            {"HEB450", {298.5356, 26.3925, 298.5356, 26.3925, 400.5413, 15.1335}},
            {"HEB500", {298.5444, 28.4145, 298.5444, 28.4145, 446.669, 15.4803}},
            {"HEB550", {298.9226, 29.3338, 298.9226, 29.3338, 495.1858, 15.9425}},
            {"HEB600", {298.862, 30.3568, 298.862, 30.3568, 543.2262, 16.3409}},
            {"HEB650", {298.8843, 31.3538, 298.8843, 31.3538, 591.5389, 16.7531}},
            {"HEB700", {298.9326, 32.3589, 298.9326, 32.3589, 639.6524, 17.6892}},
            {"HEB800", {298.6181, 33.4794, 298.6181, 33.4794, 738.823, 18.2072}},
            {"HEB900", {298.9461, 35.4115, 298.9461, 35.4115, 835.4178, 19.1301}},
            {"HEB1000", {298.8989, 36.4395, 298.8989, 36.4395, 933.7528, 19.5424}},
            {"IPE80", {45.7955, 5.2815, 45.7955, 5.2815, 70.1448, 3.9639}},
            {"IPE100", {53.9583, 6.0481, 53.9583, 6.0481, 88.6381, 4.2675}},
            {"IPE120", {63.5475, 6.4532, 63.5475, 6.4532, 108.1512, 4.6318}},
            {"IPE140", {72.3342, 7.097, 72.3342, 7.097, 126.636, 4.8523}},
            {"IPE160", {81.3965, 7.5744, 81.3965, 7.5744, 146.4131, 5.3257}},
            {"IPE180", {90.6581, 8.112, 90.6581, 8.112, 165.4307, 5.5649}},
            {"IPE200", {98.8244, 8.801, 98.8244, 8.801, 184.7715, 6.0364}},
            {"IPE220", {109.3756, 9.3775, 109.3756, 9.3775, 203.6655, 6.3416}},
            {"IPE240", {118.9326, 10.1011, 118.9326, 10.1011, 223.2791, 6.7843}},
            {"IPE270", {133.8222, 10.4908, 133.8222, 10.4908, 252.375, 7.0887}},
            {"IPE300", {149.394, 10.8456, 149.394, 10.8456, 281.9668, 7.6118}},
            {"IPE330", {158.1958, 11.9168, 158.1958, 11.9168, 310.2552, 8.0501}},
            {"IPE360", {168.7851, 12.9867, 168.7851, 12.9867, 338.3886, 8.5571}},
            {"IPE400", {178.4041, 13.8935, 178.4041, 13.8935, 377.6917, 9.2802}},
            {"IPE450", {188.8299, 14.8964, 188.8299, 14.8964, 426.0294, 10.0199}},
            {"IPE500", {198.6215, 16.3589, 198.6215, 16.3589, 472.5746, 10.7171}},
            {"IPE550", {208.0608, 17.7206, 208.0608, 17.7206, 520.8236, 11.6783}},
            {"IPE600", {218.3108, 19.4718, 218.3108, 19.4718, 567.0712, 12.5455}},
            {"HEM100", {105.5019, 20.3133, 105.5019, 20.3133, 79.8257, 12.9802}},
            {"HEM120", {125.6042, 21.2127, 125.6042, 21.2127, 98.0436, 13.4019}},
            {"HEM140", {145.6138, 22.1756, 145.6138, 22.1756, 116.095, 13.8232}},
            {"HEM160", {165.4564, 23.2411, 165.4564, 23.2411, 134.1825, 15.0849}},
            {"HEM180", {185.5054, 24.1978, 185.5054, 24.1978, 152.2794, 15.4807}},
            {"HEM200", {205.2962, 25.2657, 205.2962, 25.2657, 170.3865, 16.2204}},
            {"HEM220", {225.3201, 26.242, 225.3201, 26.242, 188.4152, 16.5708}},
            {"HEM240", {247.2405, 32.3065, 247.2405, 32.3065, 206.4009, 19.3584}},
            {"HEM260", {267.1516, 32.8245, 267.1516, 32.8245, 225.7206, 19.6576}},
            {"HEM280", {286.9407, 33.3683, 286.9407, 33.3683, 244.5611, 19.9738}},
            {"HEM300", {308.9561, 39.4015, 308.9561, 39.4015, 262.6004, 22.7783}},
            {"HEM320", {307.9445, 40.4264, 307.9445, 40.4264, 279.5918, 22.6043}},
            {"HEM340", {307.9217, 40.4327, 307.9217, 40.4327, 297.682, 22.5018}},
            {"HEM360", {307.0607, 40.3775, 307.0607, 40.3775, 315.9547, 22.4766}},
            {"HEM400", {306.119, 40.3686, 306.119, 40.3686, 353.1491, 22.3274}},
            {"HEM450", {306.0943, 40.373, 306.0943, 40.373, 399.4615, 22.1385}},
            {"HEM500", {304.9735, 40.4116, 304.9735, 40.4116, 445.5218, 21.9985}},
            {"HEM550", {305.1661, 40.3487, 305.1661, 40.3487, 494.046, 21.9296}},
            {"HEM600", {304.2427, 40.3259, 304.2427, 40.3259, 542.3115, 21.8577}},
            {"HEM650", {304.0394, 40.4015, 304.0394, 40.4015, 590.3673, 21.7198}},
            {"HEM700", {303.0963, 40.3841, 303.0963, 40.3841, 638.6716, 21.6711}},
            {"HEM800", {301.9789, 40.4432, 301.9789, 40.4432, 737.9423, 21.7248}},
            {"HEM900", {300.8195, 40.5028, 300.8195, 40.5028, 834.1924, 21.5979}},
            {"HEM1000", {300.8277, 40.5054, 300.8277, 40.5054, 932.7207, 21.525}},
            {"UB127X76X13", {75.9571, 7.6524, 75.9571, 7.6524, 112.8212, 4.3387}},
            {"UB152X89X16", {88.5804, 7.755, 88.5804, 7.755, 138.0714, 4.7679}},
            {"UB178X102X19", {101.1588, 7.9281, 101.1588, 7.9281, 163.1115, 5.0749}},
            {"UB203X102X23", {101.6269, 9.357, 101.6269, 9.357, 185.5984, 5.6039}},
            {"UB203X133X25", {133.1503, 7.8186, 133.1503, 7.8186, 188.6285, 5.9355}},
            {"UB203X133X30", {133.797, 9.6309, 133.797, 9.6309, 188.4496, 6.606}},
            {"UB254X102X22", {100.9484, 6.9164, 100.9484, 6.9164, 241.7426, 5.8136}},
            {"UB254X102X25", {101.8137, 8.4425, 101.8137, 8.4425, 241.9236, 6.1304}},
            {"UB254X102X28", {102.2714, 10.0078, 102.2714, 10.0078, 241.7678, 6.4724}},
            {"UB254X146X31", {146.0431, 8.6191, 146.0431, 8.6191, 235.2706, 6.1821}},
            {"UB254X146X37", {146.3021, 10.9288, 146.3021, 10.9288, 235.0015, 6.4839}},
            {"UB254X146X43", {147.044, 12.76, 147.044, 12.76, 234.7559, 7.3643}},
            {"UB305X102X25", {101.315, 7.0663, 101.315, 7.0663, 292.9078, 5.9065}},
            {"UB305X102X28", {101.2985, 8.9139, 101.2985, 8.9139, 292.3051, 6.1082}},
            {"UB305X102X33", {102.0982, 10.8942, 102.0982, 10.8942, 292.233, 6.698}},
            {"UB305X127X37", {123.0561, 10.7879, 123.0561, 10.7879, 284.321, 7.2704}},
            {"UB305X127X42", {124.1138, 12.1654, 124.1138, 12.1654, 284.2673, 8.171}},
            {"UB305X127X48", {124.98, 14.111, 124.98, 14.111, 283.9119, 9.1415}},
            {"UB305X165X40", {164.6817, 10.2543, 164.6817, 10.2543, 284.2759, 6.1751}},
            {"UB305X165X46", {165.5371, 11.839, 165.5371, 11.839, 284.1885, 6.8736}},
            {"UB305X165X54", {166.7541, 13.7369, 166.7541, 13.7369, 283.983, 8.1018}},
            {"UB356X127X33", {124.7752, 8.6261, 124.7752, 8.6261, 334.4901, 6.1615}},
            {"UB356X127X39", {125.6851, 10.7904, 125.6851, 10.7904, 334.0697, 6.7972}},
            {"UB356X171X45", {170.6649, 9.7748, 170.6649, 9.7748, 333.8763, 7.179}},
            {"UB356X171X51", {171.1301, 11.5726, 171.1301, 11.5726, 333.59, 7.5913}},
            {"UB356X171X57", {171.8676, 13.0747, 171.8676, 13.0747, 333.2899, 8.3061}},
            {"UB356X171X67", {172.9119, 15.7788, 172.9119, 15.7788, 333.0377, 9.2988}},
            {"UB406X140X39", {141.4732, 8.6673, 141.4732, 8.6673, 383.502, 6.5726}},
            {"UB406X140X46", {141.8737, 11.2793, 141.8737, 11.2793, 383.0765, 6.9519}},
            {"UB406X140X53", {143.1754, 12.9449, 143.1754, 12.9449, 382.722, 8.0908}},
            {"UB406X178X54", {177.3479, 10.9639, 177.3479, 10.9639, 382.5319, 7.8814}},
            {"UB406X178X60", {177.6661, 12.8513, 177.6661, 12.8513, 382.4841, 8.0692}},
            {"UB406X178X67", {178.4844, 14.3778, 178.4844, 14.3778, 382.131, 8.953}},
            {"UB406X178X74", {179.1807, 16.082, 179.1807, 16.082, 381.9492, 9.6621}},
            {"UB406X178X85", {180.7295, 18.2797, 180.7295, 18.2797, 381.5856, 11.0495}},
            {"UB457X152X52", {152.0464, 10.98, 152.0464, 10.98, 430.282, 7.7271}},
            {"UB457X152X60", {152.6167, 13.3834, 152.6167, 13.3834, 429.8311, 8.2316}},
            {"UB457X152X67", {153.6325, 15.0598, 153.6325, 15.0598, 429.6667, 9.1592}},
            {"UB457X152X74", {154.266, 17.0547, 154.266, 17.0547, 429.5065, 9.7602}},
            {"UB457X152X82", {155.1303, 18.974, 155.1303, 18.974, 429.3043, 10.6359}},
            {"UB457X191X67", {189.6189, 12.7563, 189.6189, 12.7563, 429.6725, 8.6482}},
            {"UB457X191X74", {190.122, 14.5633, 190.122, 14.5633, 429.4684, 9.1427}},
            {"UB457X191X82", {191.1134, 16.0496, 191.1134, 16.0496, 429.3108, 10.0588}},
            {"UB457X191X89", {191.6522, 17.7666, 191.6522, 17.7666, 429.1275, 10.6573}},
            {"UB457X191X98", {192.6407, 19.6501, 192.6407, 19.6501, 429.0563, 11.564}},
            {"UB457X191X106", {193.9316, 20.6271, 193.9316, 20.6271, 428.7478, 12.7611}},
            {"UB533X165X66", {164.3772, 11.5331, 164.3772, 11.5331, 505.5143, 9.0861}},
            {"UB533X165X74", {165.6, 13.7113, 165.6, 13.7113, 504.7683, 9.8933}},
            {"UB533X165X85", {165.6431, 16.6205, 165.6431, 16.6205, 504.5127, 10.4845}},
            {"UB533X210X82", {208.3251, 13.2909, 208.3251, 13.2909, 504.372, 9.7905}},
            {"UB533X210X92", {208.9309, 15.6845, 208.9309, 15.6845, 504.0456, 10.2987}},
            {"UB533X210X101", {209.7159, 17.4732, 209.7159, 17.4732, 503.8796, 11.0071}},
            {"UB533X210X109", {210.5263, 18.878, 210.5263, 18.878, 503.6587, 11.8053}},
            {"UB533X210X122", {211.6512, 21.3806, 211.6512, 21.3806, 503.4647, 12.9011}},
            {"UB533X210X138", {213.6507, 23.721, 213.6507, 23.721, 503.2034, 14.8864}},
            {"UB610X178X82", {177.5014, 12.9154, 177.5014, 12.9154, 576.5473, 10.1646}},
            {"UB610X178X92", {178.578, 15.1118, 178.578, 15.1118, 575.8173, 11.0604}},
            {"UB610X229X101", {227.2959, 14.8626, 227.2959, 14.8626, 575.493, 10.6677}},
            {"UB610X229X113", {227.8443, 17.3827, 227.8443, 17.3827, 575.0646, 11.2593}},
            {"UB610X229X125", {228.7045, 19.6771, 228.7045, 19.6771, 574.8515, 12.0646}},
            {"UB610X229X140", {229.9591, 22.1705, 229.9591, 22.1705, 574.6775, 13.2738}},
            {"UB610X229X153", {228.7199, 24.9922, 228.7199, 24.9922, 574.6922, 14.1703}},
            {"UB610X305X149", {304.324, 19.7937, 304.324, 19.7937, 575.504, 12.0973}},
            {"UB610X305X179", {306.7058, 23.6954, 306.7058, 23.6954, 574.897, 14.4091}},
            {"UB610X305X238", {311.1327, 31.4892, 311.1327, 31.4892, 574.4791, 18.7015}},
            {"UB686X254X125", {252.5416, 16.2916, 252.5416, 16.2916, 648.6499, 11.9151}},
            {"UB686X254X140", {253.2568, 19.1025, 253.2568, 19.1025, 648.2564, 12.6055}},
            {"UB686X254X152", {254.1343, 21.0944, 254.1343, 21.0944, 648.013, 13.4174}},
            {"UB686X254X170", {255.4304, 23.8049, 255.4304, 23.8049, 647.5816, 14.7124}},
            {"UB686X254X192", {253.6859, 28.0083, 253.6859, 28.0083, 647.8938, 15.7082}},
            {"UB762X267X134", {263.8659, 15.5992, 263.8659, 15.5992, 723.0507, 12.2213}},
            {"UB762X267X147", {264.6616, 17.6099, 264.6616, 17.6099, 722.5388, 13.021}},
            {"UB762X267X173", {266.2722, 21.7092, 266.2722, 21.7092, 721.9037, 14.5277}},
            {"UB762X267X197", {267.6014, 25.5188, 267.6014, 25.5188, 721.5134, 15.8158}},
            {"UB762X267X220", {265.6591, 30.116, 265.6591, 30.116, 721.2825, 16.7462}},
            {"UB838X292X176", {291.1272, 18.9154, 291.1272, 18.9154, 801.1563, 14.2267}},
            {"UB838X292X194", {291.8826, 21.8186, 291.8826, 21.8186, 800.7314, 14.9282}},
            {"UB838X292X226", {293.2845, 26.9467, 293.2845, 26.9467, 799.8448, 16.3327}},
            {"UB838X292X251", {291.5324, 31.1532, 291.5324, 31.1532, 799.2662, 17.2312}},
            {"UB914X305X201", {302.684, 20.3292, 302.684, 20.3292, 866.9094, 15.3361}},
            {"UB914X305X224", {303.8061, 23.9824, 303.8061, 23.9824, 866.6151, 16.157}},
            {"UB914X305X238", {304.718, 25.9868, 304.718, 25.9868, 866.8248, 16.7577}},
            {"UB914X305X253", {304.878, 28.0747, 304.878, 28.0747, 865.4506, 17.5301}},
            {"UB914X305X271", {306.6145, 30.1233, 306.6145, 30.1233, 865.8843, 18.6478}},
            {"UB914X305X289", {307.406, 32.1027, 307.406, 32.1027, 865.4473, 19.7642}},
            {"UB914X305X313", {308.5368, 34.6632, 308.5368, 34.6632, 865.241, 21.338}},
            {"UB914X305X345", {307.6642, 40.0463, 307.6642, 40.0463, 865.2398, 22.3528}},
            {"UB914X305X381", {309.6698, 44.0476, 309.6698, 44.0476, 865.0275, 24.6451}},
            {"UB914X305X425", {312.6254, 49.1793, 312.6254, 49.1793, 864.4827, 27.1396}},
            {"UB914X305X474", {315.7882, 54.2219, 315.7882, 54.2219, 864.2906, 30.2614}},
            {"UB914X305X521", {318.746, 59.0516, 318.746, 59.0516, 864.4704, 33.2442}},
            {"UB914X305X576", {321.6714, 65.2075, 321.6714, 65.2075, 864.011, 36.3172}},
            {"UB920X420X344", {417.6769, 32.0742, 417.6769, 32.0742, 865.2707, 19.5783}},
            {"UB920X420X368", {418.7278, 34.3688, 418.7278, 34.3688, 864.4813, 20.5761}},
            {"UB920X420X390", {419.6215, 36.6996, 419.6215, 36.6996, 864.5803, 21.5606}},
            {"UB920X420X420", {421.7192, 39.9837, 421.7192, 39.9837, 864.8939, 22.7753}},
            {"UB920X420X449", {422.7444, 42.7761, 422.7444, 42.7761, 864.2307, 24.2831}},
            {"UB920X420X491", {421.7309, 47.089, 421.7309, 47.089, 864.3928, 26.1706}},
            {"UB920X420X537", {424.6831, 51.2127, 424.6831, 51.2127, 864.1659, 28.6541}},
            {"UB920X420X588", {426.7336, 56.006, 426.7336, 56.006, 864.1731, 31.2666}},
            {"UB920X420X656", {430.7785, 62.0961, 430.7785, 62.0961, 863.858, 34.7761}},
            {"UB920X420X725", {433.7455, 68.2227, 433.7455, 68.2227, 863.4699, 38.3606}},
            {"UB920X420X787", {436.8241, 73.9947, 436.8241, 73.9947, 864.0282, 41.158}},
            {"UB920X420X970", {446.0741, 89.8871, 446.0741, 89.8871, 864.0012, 50.3661}},
            {"UB920X420X1077", {450.7827, 99.2186, 450.7827, 99.2186, 863.6192, 55.2995}},
            {"UB920X420X1194", {457.0228, 108.9702, 457.0228, 108.9702, 863.9832, 60.8855}},
            {"UB920X420X1269", {460.7681, 115.3044, 460.7681, 115.3044, 863.2935, 64.224}},
            {"UB920X420X1377", {472.6957, 115.3332, 472.6957, 115.3332, 863.1919, 76.8896}},
            {"UB1000X400X296", {398.6039, 27.3923, 398.6039, 27.3923, 934.2465, 17.078}},
            {"UB1000X400X321", {398.7522, 31.2977, 398.7522, 31.2977, 933.6837, 17.084}},
            {"UB1000X400X371", {398.9975, 36.3791, 398.9975, 36.3791, 932.705, 19.5968}},
            {"UB1000X400X412", {400.9456, 40.3192, 400.9456, 40.3192, 932.1719, 21.6753}},
            {"UB1000X400X443", {401.1324, 42.1857, 401.1324, 42.1857, 932.2506, 24.1949}},
            {"UB1000X400X483", {403.0866, 46.324, 403.0866, 46.324, 931.2582, 25.9796}},
            {"UB1000X400X539", {406.2029, 51.4155, 406.2029, 51.4155, 930.9929, 28.9698}},
            {"UB1000X400X554", {407.1611, 52.3372, 407.1611, 52.3372, 930.7674, 30.0736}},
            {"UB1000X400X591", {408.2034, 56.2404, 408.2034, 56.2404, 930.8293, 31.5646}},
            {"UB1000X400X642", {411.2536, 60.3413, 411.2536, 60.3413, 930.5202, 34.5546}},
            {"UB1000X400X748", {416.3896, 70.3273, 416.3896, 70.3273, 929.9945, 39.5668}},
            {"UB1000X400X883", {423.4938, 82.34, 423.4938, 82.34, 929.4604, 46.0315}},
            {"UB1000X400X976", {427.711, 90.084, 427.711, 90.084, 929.7939, 50.8215}},
            {"UB1016X305X222", {298.2029, 21.5148, 298.2029, 21.5148, 937.2657, 16.512}},
            {"UB1016X305X249", {298.1759, 26.5, 298.1759, 26.5, 935.4911, 17.0029}},
            {"UB1016X305X272", {298.4732, 31.4976, 298.4732, 31.4976, 934.3554, 17.0193}},
            {"UB1016X305X314", {298.6882, 36.4041, 298.6882, 36.4041, 933.5699, 19.6225}},
            {"UB1016X305X350", {300.8097, 40.5108, 300.8097, 40.5108, 932.7164, 21.6185}},
            {"UB1016X305X393", {302.0048, 44.3843, 302.0048, 44.3843, 932.4592, 24.9208}},
            {"UB1016X305X415", {302.9716, 46.519, 302.9716, 46.519, 931.7769, 26.5148}},
            {"UB1016X305X438", {304.0305, 49.5234, 304.0305, 49.5234, 931.5197, 27.4208}},
            {"UB1016X305X494", {308.0624, 54.5367, 308.0624, 54.5367, 930.9929, 31.5042}},
            {"UB1016X305X584", {313.1208, 64.5907, 313.1208, 64.5907, 930.1174, 36.4937}},
            {"UB1100X400X343", {398.8177, 31.2781, 398.8177, 31.2781, 1034.3391, 18.5161}},
            {"UB1100X400X390", {399.0617, 36.2618, 399.0617, 36.2618, 1033.5009, 20.5322}},
            {"UB1100X400X433", {401.0933, 40.2847, 401.0933, 40.2847, 1032.7817, 22.5169}},
            {"UB1100X400X499", {404.2517, 45.2672, 404.2517, 45.2672, 1032.1865, 26.5281}},
            {"UB1100X400X548", {406.313, 50.2709, 406.313, 50.2709, 1031.7298, 28.5241}},
            {"UB1100X400X607", {409.3642, 55.2699, 409.3642, 55.2699, 1031.2904, 31.528}},
            {"UC152X152X23", {151.831, 6.8487, 151.831, 6.8487, 139.5233, 6.0387}},
            {"UC152X152X30", {152.6927, 9.4384, 152.6927, 9.4384, 139.2727, 6.8161}},
            {"UC152X152X37", {154.1423, 11.5563, 154.1423, 11.5563, 139.1227, 8.2616}},
            {"UC152X152X44", {155.6319, 13.6703, 155.6319, 13.6703, 139.0847, 9.7557}},
            {"UC152X152X51", {157.1331, 15.777, 157.1331, 15.777, 139.0038, 11.2496}},
            {"UC203X203X46", {203.1727, 11.067, 203.1727, 11.067, 181.9317, 7.5647}},
            {"UC203X203X52", {204.0837, 12.5413, 204.0837, 12.5413, 181.9066, 8.3244}},
            {"UC203X203X60", {205.6036, 14.2426, 205.6036, 14.2426, 181.7972, 9.8266}},
            {"UC203X203X71", {206.0987, 17.3729, 206.0987, 17.3729, 181.5939, 10.3649}},
            {"UC203X203X86", {208.8532, 20.5692, 208.8532, 20.5692, 181.5189, 13.0633}},
            {"UC203X203X100", {209.7951, 23.766, 209.7951, 23.766, 181.8789, 14.8488}},
            {"UC254X254X73", {254.1924, 14.2677, 254.1924, 14.2677, 226.6429, 9.0953}},
            {"UC254X254X89", {255.9089, 17.3762, 255.9089, 17.3762, 226.417, 10.7837}},
            {"UC254X254X107", {258.5582, 20.5585, 258.5582, 20.5585, 226.2917, 13.314}},
            {"UC254X254X132", {261.0219, 25.3798, 261.0219, 25.3798, 226.1077, 15.7649}},
            {"UC254X254X167", {265.0566, 31.7505, 265.0566, 31.7505, 226.0702, 19.7439}},
            {"UC305X305X97", {304.7749, 15.4787, 304.7749, 15.4787, 278.3949, 10.4619}},
            {"UC305X305X118", {306.8896, 18.7912, 306.8896, 18.7912, 278.0536, 12.5675}},
            {"UC305X305X137", {308.6205, 21.8218, 308.6205, 21.8218, 277.7603, 14.3229}},
            {"UC305X305X158", {310.913, 25.0697, 310.913, 25.0697, 277.8268, 16.4021}},
            {"UC305X305X198", {314.1354, 31.5105, 314.1354, 31.5105, 277.5083, 19.6374}},
            {"UC305X305X240", {317.9706, 37.845, 317.9706, 37.845, 277.3012, 23.5066}},
            {"UC305X305X283", {321.9044, 44.2133, 321.9044, 44.2133, 277.3246, 27.3359}},
            {"UC305X305X313", {324.6954, 48.4318, 324.6954, 48.4318, 277.5185, 30.4994}},
            {"UC305X305X342", {327.6905, 52.7454, 327.6905, 52.7454, 276.8442, 33.0781}},
            {"UC356X368X129", {368.1835, 17.5557, 368.1835, 17.5557, 321.7823, 10.9074}},
            {"UC356X368X153", {369.9659, 20.7839, 369.9659, 20.7839, 321.4312, 12.7826}},
            {"UC356X368X177", {372.2277, 23.8706, 372.2277, 23.8706, 321.3415, 14.8941}},
            {"UC356X368X202", {374.2938, 27.0879, 374.2938, 27.0879, 321.1644, 16.9683}},
            {"UC356X406X235", {394.3497, 30.2969, 394.3497, 30.2969, 321.1233, 18.7329}},
            {"UC356X406X287", {398.8362, 36.5459, 398.8362, 36.5459, 321.0743, 23.1246}},
            {"UC356X406X340", {402.6698, 43.0007, 402.6698, 43.0007, 320.7707, 27.0524}},
            {"UC356X406X393", {406.7922, 49.2768, 406.7922, 49.2768, 320.7993, 31.095}},
            {"UC356X406X467", {412.0674, 58.0534, 412.0674, 58.0534, 320.8129, 36.3135}},
            {"UC356X406X509", {415.9006, 62.7445, 415.9006, 62.7445, 320.7708, 39.6317}},
            {"UC356X406X551", {418.2305, 67.6273, 418.2305, 67.6273, 320.5908, 42.5219}},
            {"UC356X406X592", {420.8652, 72.3672, 420.8652, 72.3672, 320.5476, 45.4912}},
            {"UC356X406X634", {423.7824, 77.1189, 423.7824, 77.1189, 320.5833, 48.0283}},
            {"UC356X406X677", {427.9341, 81.5571, 427.9341, 81.5571, 320.129, 51.6838}},
            {"UC356X406X744", {431.6653, 89.0847, 431.6653, 89.0847, 320.0127, 55.9552}},
            {"UC356X406X818", {436.6966, 97.1833, 436.6966, 97.1833, 319.8254, 60.8362}},
            {"UC356X406X900", {441.5606, 106.3007, 441.5606, 106.3007, 318.5405, 66.0695}},
            {"UC356X406X990", {448.1227, 114.9289, 448.1227, 114.9289, 320.2834, 72.5575}},
            {"UC356X406X1086", {453.517, 125.3919, 453.517, 125.3919, 318.325, 78.0497}},
            {"UC356X406X1202", {470.7864, 130.1704, 470.7864, 130.1704, 319.798, 95.334}},
            {"UC356X406X1299", {475.8747, 140.1326, 475.8747, 140.1326, 319.859, 100.3596}},
            {"UBP203X203X45", {205.4937, 9.557, 205.4937, 9.557, 182.0479, 9.8648}},
            {"UBP203X203X54", {207.3064, 11.4638, 207.3064, 11.4638, 181.8645, 11.6589}},
            {"UBP254X254X63", {256.0631, 10.7659, 256.0631, 10.7659, 226.902, 11.0708}},
            {"UBP254X254X71", {257.5438, 12.0631, 257.5438, 12.0631, 226.7363, 12.4893}},
            {"UBP254X254X85", {260.0113, 14.3631, 260.0113, 14.3631, 226.4902, 14.9029}},
            {"UBP305X305X79", {305.6881, 11.1763, 305.6881, 11.1763, 278.7828, 11.5666}},
            {"UBP305X305X88", {307.1941, 12.3719, 307.1941, 12.3719, 278.593, 12.9779}},
            {"UBP305X305X95", {308.088, 13.38, 308.088, 13.38, 278.4725, 13.8363}},
            {"UBP305X305X110", {310.1867, 15.4761, 310.1867, 15.4761, 278.2621, 15.8711}},
            {"UBP305X305X126", {312.4071, 17.6838, 312.4071, 17.6838, 278.1041, 18.0461}},
            {"UBP305X305X149", {315.4455, 20.8099, 315.4455, 20.8099, 277.7977, 21.1202}},
            {"UBP305X305X186", {320.5349, 25.6831, 320.5349, 25.6831, 277.7278, 26.0414}},
            {"UBP305X305X223", {325.3769, 30.4985, 325.3769, 30.4985, 277.5135, 30.8406}},
            {"UBP356X368X109", {370.5287, 12.952, 370.5287, 12.952, 322.126, 13.2863}},
            {"UBP356X368X133", {373.1768, 15.7782, 373.1768, 15.7782, 321.6612, 16.0779}},
            {"UBP356X368X152", {375.7099, 17.9439, 375.7099, 17.9439, 321.6695, 18.3248}},
            {"UBP356X368X174", {378.0484, 20.4696, 378.0484, 20.4696, 321.4349, 20.7836}},
        };

        return registry.fetch(type);
    }

    vec nzchsection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"610.0X12.7CHS", {305.00, 12.7}},
            {"610.0X9.5CHS", {305.00, 9.5}},
            {"610.0X6.4CHS", {305.00, 6.4}},
            {"508.0X12.7CHS", {254.00, 12.7}},
            {"508.0X9.5CHS", {254.00, 9.5}},
            {"508.0X6.4CHS", {254.00, 6.4}},
            {"165.1X5.4CHS", {82.55, 5.4}},
            {"165.1X5.0CHS", {82.55, 5.0}},
            {"139.7X5.4CHS", {69.85, 5.4}},
            {"139.7X5.0CHS", {69.85, 5.0}},
            {"114.3X5.4CHS", {57.15, 5.4}},
            {"114.3X4.5CHS", {57.15, 4.5}},
            {"101.6X5.0CHS", {50.80, 5.0}},
            {"101.6X4.0CHS", {50.80, 4.0}},
            {"88.9X5.9CHS", {44.45, 5.9}},
            {"88.9X5.0CHS", {44.45, 5.0}},
            {"88.9X4.0CHS", {44.45, 4.0}},
            {"76.1X5.9CHS", {38.05, 5.9}},
            {"76.1X4.5CHS", {38.05, 4.5}},
            {"76.1X3.6CHS", {38.05, 3.6}},
            {"60.3X5.4CHS", {30.15, 5.4}},
            {"60.3X4.5CHS", {30.15, 4.5}},
            {"60.3X3.6CHS", {30.15, 3.6}},
            {"48.3X5.4CHS", {24.15, 5.4}},
            {"48.3X4.0CHS", {24.15, 4.0}},
            {"48.3X3.2CHS", {24.15, 3.2}},
            {"42.4X4.9CHS", {21.20, 4.9}},
            {"42.4X4.0CHS", {21.20, 4.0}},
            {"42.4X3.2CHS", {21.20, 3.2}},
            {"457.0X12.7CHS", {228.50, 12.7}},
            {"457.0X9.5CHS", {228.50, 9.5}},
            {"457.0X6.4CHS", {228.50, 6.4}},
            {"406.4X12.7CHS", {203.20, 12.7}},
            {"406.4X9.5CHS", {203.20, 9.5}},
            {"406.4X6.4CHS", {203.20, 6.4}},
            {"355.6X12.7CHS", {177.80, 12.7}},
            {"355.6X9.5CHS", {177.80, 9.5}},
            {"355.6X6.4CHS", {177.80, 6.4}},
            {"323.9X12.7CHS", {161.95, 12.7}},
            {"323.9X9.5CHS", {161.95, 9.5}},
            {"323.9X6.4CHS", {161.95, 6.4}},
            {"273.1X9.3CHS", {136.55, 9.3}},
            {"273.1X6.4CHS", {136.55, 6.4}},
            {"273.1X4.8CHS", {136.55, 4.8}},
            {"219.1X8.2CHS", {109.55, 8.2}},
            {"219.1X6.4CHS", {109.55, 6.4}},
            {"219.1X4.8CHS", {109.55, 4.8}},
            {"168.3X7.1CHS", {84.15, 7.1}},
            {"168.3X6.4CHS", {84.15, 6.4}},
            {"168.3X4.8CHS", {84.15, 4.8}},
            {"165.1X3.5CHS", {82.55, 3.5}},
            {"165.1X3.0CHS", {82.55, 3.0}},
            {"139.7X3.5CHS", {69.85, 3.5}},
            {"139.7X3.0CHS", {69.85, 3.0}},
            {"114.3X6.0CHS", {57.15, 6.0}},
            {"114.3X4.8CHS", {57.15, 4.8}},
            {"114.3X3.6CHS", {57.15, 3.6}},
            {"114.3X3.2CHS", {57.15, 3.2}},
            {"101.6X3.2CHS", {50.80, 3.2}},
            {"101.6X2.6CHS", {50.80, 2.6}},
            {"88.9X5.5CHS", {44.45, 5.5}},
            {"88.9X4.8CHS", {44.45, 4.8}},
            {"88.9X3.2CHS", {44.45, 3.2}},
            {"88.9X2.6CHS", {44.45, 2.6}},
            {"76.1X3.2CHS", {38.05, 3.2}},
            {"76.1X2.3CHS", {38.05, 2.3}},
        };

        return registry.fetch(type);
    }

    vec nzisection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"1200WB455", {500.0, 40.0, 500.0, 40.0, 1120.0, 16.0}},
            {"1200WB423", {500.0, 36.0, 500.0, 36.0, 1120.0, 16.0}},
            {"1200WB392", {500.0, 32.0, 500.0, 32.0, 1120.0, 16.0}},
            {"1200WB342", {400.0, 32.0, 400.0, 32.0, 1120.0, 16.0}},
            {"1200WB317", {400.0, 28.0, 400.0, 28.0, 1120.0, 16.0}},
            {"1200WB278", {350.0, 25.0, 350.0, 25.0, 1120.0, 16.0}},
            {"1200WB249", {275.0, 25.0, 275.0, 25.0, 1120.0, 16.0}},
            {"1000WB322", {400.0, 32.0, 400.0, 32.0, 960.0, 16.0}},
            {"1000WB296", {400.0, 28.0, 400.0, 28.0, 960.0, 16.0}},
            {"1000WB258", {350.0, 25.0, 350.0, 25.0, 960.0, 16.0}},
            {"1000WB215", {300.0, 20.0, 300.0, 20.0, 960.0, 16.0}},
            {"900WB282", {400.0, 32.0, 400.0, 32.0, 860.0, 12.0}},
            {"900WB257", {400.0, 28.0, 400.0, 28.0, 860.0, 12.0}},
            {"900WB218", {350.0, 25.0, 350.0, 25.0, 860.0, 12.0}},
            {"900WB175", {300.0, 20.0, 300.0, 20.0, 860.0, 12.0}},
            {"800WB192", {300.0, 28.0, 300.0, 28.0, 760.0, 10.0}},
            {"800WB168", {275.0, 25.0, 275.0, 25.0, 760.0, 10.0}},
            {"800WB146", {275.0, 20.0, 275.0, 20.0, 760.0, 10.0}},
            {"800WB122", {250.0, 16.0, 250.0, 16.0, 760.0, 10.0}},
            {"700WB173", {275.0, 28.0, 275.0, 28.0, 660.0, 10.0}},
            {"700WB150", {250.0, 25.0, 250.0, 25.0, 660.0, 10.0}},
            {"700WB130", {250.0, 20.0, 250.0, 20.0, 660.0, 10.0}},
            {"700WB115", {250.0, 16.0, 250.0, 16.0, 660.0, 10.0}},
            {"500WC440", {500.0, 40.0, 500.0, 40.0, 400.0, 40.0}},
            {"500WC414", {500.0, 40.0, 500.0, 40.0, 400.0, 32.0}},
            {"500WC383", {500.0, 36.0, 500.0, 36.0, 400.0, 32.0}},
            {"500WC340", {500.0, 32.0, 500.0, 32.0, 450.0, 25.0}},
            {"500WC290", {500.0, 28.0, 500.0, 28.0, 450.0, 20.0}},
            {"500WC267", {500.0, 25.0, 500.0, 25.0, 450.0, 20.0}},
            {"500WC228", {500.0, 20.0, 500.0, 20.0, 450.0, 20.0}},
            {"400WC361", {400.0, 40.0, 400.0, 40.0, 350.0, 40.0}},
            {"400WC328", {400.0, 40.0, 400.0, 40.0, 350.0, 28.0}},
            {"400WC303", {400.0, 36.0, 400.0, 36.0, 350.0, 28.0}},
            {"400WC270", {400.0, 32.0, 400.0, 32.0, 350.0, 25.0}},
            {"400WC212", {400.0, 25.0, 400.0, 25.0, 350.0, 20.0}},
            {"400WC181", {400.0, 20.0, 400.0, 20.0, 350.0, 20.0}},
            {"400WC144", {400.0, 16.0, 400.0, 16.0, 350.0, 16.0}},
            {"350WC280", {350.0, 40.0, 350.0, 40.0, 275.0, 28.0}},
            {"350WC258", {350.0, 36.0, 350.0, 36.0, 275.0, 28.0}},
            {"350WC230", {350.0, 32.0, 350.0, 32.0, 275.0, 25.0}},
            {"350WC197", {350.0, 28.0, 350.0, 28.0, 275.0, 20.0}},
            {"610UB125", {229.0, 19.6, 229.0, 19.6, 572.0, 11.9}},
            {"610UB113", {228.0, 17.3, 228.0, 17.3, 572.0, 11.2}},
            {"610UB101", {228.0, 14.8, 228.0, 14.8, 572.0, 10.6}},
            {"530UB92.4", {209.0, 15.6, 209.0, 15.6, 502.0, 10.2}},
            {"530UB82.0", {209.0, 13.2, 209.0, 13.2, 502.0, 9.6}},
            {"460UB82.1", {191.0, 16.0, 191.0, 16.0, 428.0, 9.9}},
            {"460UB74.6", {190.0, 14.5, 190.0, 14.5, 428.0, 9.1}},
            {"460UB67.1", {190.0, 12.7, 190.0, 12.7, 428.0, 8.5}},
            {"410UB59.7", {178.0, 12.8, 178.0, 12.8, 381.0, 7.8}},
            {"410UB53.7", {178.0, 10.9, 178.0, 10.9, 381.0, 7.6}},
            {"360UB56.7", {172.0, 13.0, 172.0, 13.0, 333.0, 8.0}},
            {"360UB50.7", {171.0, 11.5, 171.0, 11.5, 333.0, 7.3}},
            {"360UB44.7", {171.0, 9.7, 171.0, 9.7, 333.0, 6.9}},
            {"310UB46.2", {166.0, 11.8, 166.0, 11.8, 284.0, 6.7}},
            {"310UB40.4", {165.0, 10.2, 165.0, 10.2, 284.0, 6.1}},
            {"310UB32.0", {149.0, 8.0, 149.0, 8.0, 282.0, 5.5}},
            {"250UB37.3", {146.0, 10.9, 146.0, 10.9, 234.0, 6.4}},
            {"250UB31.4", {146.0, 8.6, 146.0, 8.6, 234.0, 6.1}},
            {"250UB25.7", {124.0, 8.0, 124.0, 8.0, 232.0, 5.0}},
            {"200UB29.8", {134.0, 9.6, 134.0, 9.6, 188.0, 6.3}},
            {"200UB25.4", {133.0, 7.8, 133.0, 7.8, 188.0, 5.8}},
            {"200UB22.3", {133.0, 7.0, 133.0, 7.0, 188.0, 5.0}},
            {"200UB18.2", {99.0, 7.0, 99.0, 7.0, 184.0, 4.5}},
            {"180UB22.2", {90.0, 10.0, 90.0, 10.0, 159.0, 6.0}},
            {"180UB18.1", {90.0, 8.0, 90.0, 8.0, 159.0, 5.0}},
            {"180UB16.1", {90.0, 7.0, 90.0, 7.0, 159.0, 4.5}},
            {"150UB18.0", {75.0, 9.5, 75.0, 9.5, 136.0, 6.0}},
            {"150UB14.0", {75.0, 7.0, 75.0, 7.0, 136.0, 5.0}},
            {"320UC158", {311.0, 25.0, 311.0, 25.0, 277.0, 15.7}},
            {"310UC137", {309.0, 21.7, 309.0, 21.7, 277.0, 13.8}},
            {"310UC118", {307.0, 18.7, 307.0, 18.7, 277.0, 11.9}},
            {"310UC96.8", {305.0, 15.4, 305.0, 15.4, 277.0, 9.9}},
            {"250UC89.5", {256.0, 17.3, 256.0, 17.3, 225.0, 10.5}},
            {"250UC72.9", {254.0, 14.2, 254.0, 14.2, 225.0, 8.6}},
            {"200UC59.5", {205.0, 14.2, 205.0, 14.2, 181.0, 9.3}},
            {"200UC52.2", {204.0, 12.5, 204.0, 12.5, 181.0, 8.0}},
            {"200UC46.2", {203.0, 11.0, 203.0, 11.0, 181.0, 7.3}},
            {"150UC37.2", {154.0, 11.5, 154.0, 11.5, 139.0, 8.1}},
            {"150UC30.0", {153.0, 9.4, 153.0, 9.4, 139.0, 6.6}},
            {"150UC23.4", {152.0, 6.8, 152.0, 6.8, 139.0, 6.1}},
            {"100UC14.8", {99.0, 7.0, 99.0, 7.0, 83.0, 5.0}},
        };

        return registry.fetch(type);
    }

    vec nzrhsection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"250X150X9.0RHS", {150.0, 250.0, 9.0}},
            {"250X150X6.0RHS", {150.0, 250.0, 6.0}},
            {"250X150X5.0RHS", {150.0, 250.0, 5.0}},
            {"200X100X9.0RHS", {100.0, 200.0, 9.0}},
            {"200X100X6.0RHS", {100.0, 200.0, 6.0}},
            {"200X100X5.0RHS", {100.0, 200.0, 5.0}},
            {"200X100X4.0RHS", {100.0, 200.0, 4.0}},
            {"150X100X6.0RHS", {100.0, 150.0, 6.0}},
            {"150X100X5.0RHS", {100.0, 150.0, 5.0}},
            {"150X100X4.0RHS", {100.0, 150.0, 4.0}},
            {"150X50X5.0RHS", {50.0, 150.0, 5.0}},
            {"150X50X4.0RHS", {50.0, 150.0, 4.0}},
            {"150X50X3.0RHS", {50.0, 150.0, 3.0}},
            {"125X75X5.0RHS", {75.0, 125.0, 5.0}},
            {"125X75X4.0RHS", {75.0, 125.0, 4.0}},
            {"125X75X3.0RHS", {75.0, 125.0, 3.0}},
            {"100X50X6.0RHS", {50.0, 100.0, 6.0}},
            {"100X50X5.0RHS", {50.0, 100.0, 5.0}},
            {"100X50X4.0RHS", {50.0, 100.0, 4.0}},
            {"100X50X3.5RHS", {50.0, 100.0, 3.5}},
            {"100X50X3.0RHS", {50.0, 100.0, 3.0}},
            {"100X50X2.5RHS", {50.0, 100.0, 2.5}},
            {"100X50X2.0RHS", {50.0, 100.0, 2.0}},
            {"75X50X4.0RHS", {50.0, 75.0, 4.0}},
            {"75X50X3.0RHS", {50.0, 75.0, 3.0}},
            {"75X50X2.5RHS", {50.0, 75.0, 2.5}},
            {"75X50X2.0RHS", {50.0, 75.0, 2.0}},
            {"75X25X2.5RHS", {25.0, 75.0, 2.5}},
            {"75X25X2.0RHS", {25.0, 75.0, 2.0}},
            {"75X25X1.6RHS", {25.0, 75.0, 1.6}},
            {"65X35X3.0RHS", {35.0, 65.0, 3.0}},
            {"65X35X2.5RHS", {35.0, 65.0, 2.5}},
            {"65X35X2.0RHS", {35.0, 65.0, 2.0}},
            {"50X25X3.0RHS", {25.0, 50.0, 3.0}},
            {"50X25X2.5RHS", {25.0, 50.0, 2.5}},
            {"50X25X2.0RHS", {25.0, 50.0, 2.0}},
            {"50X25X1.6RHS", {25.0, 50.0, 1.6}},
            {"50X20X3.0RHS", {20.0, 50.0, 3.0}},
            {"50X20X2.5RHS", {20.0, 50.0, 2.5}},
            {"50X20X2.0RHS", {20.0, 50.0, 2.0}},
            {"50X20X1.6RHS", {20.0, 50.0, 1.6}},
        };

        return registry.fetch(type);
    }

    vec nzshsection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"250X9.0SHS", {250.0, 250.0, 9.0}},
            {"250X6.0SHS", {250.0, 250.0, 6.0}},
            {"200X9.0SHS", {200.0, 200.0, 9.0}},
            {"200X6.0SHS", {200.0, 200.0, 6.0}},
            {"200X5.0SHS", {200.0, 200.0, 5.0}},
            {"150X9.0SHS", {150.0, 150.0, 9.0}},
            {"150X6.0SHS", {150.0, 150.0, 6.0}},
            {"150X5.0SHS", {150.0, 150.0, 5.0}},
            {"125X9.0SHS", {125.0, 125.0, 9.0}},
            {"125X6.0SHS", {125.0, 125.0, 6.0}},
            {"125X5.0SHS", {125.0, 125.0, 5.0}},
            {"125X4.0SHS", {125.0, 125.0, 4.0}},
            {"100X9.0SHS", {100.0, 100.0, 9.0}},
            {"100X6.0SHS", {100.0, 100.0, 6.0}},
            {"100X5.0SHS", {100.0, 100.0, 5.0}},
            {"100X4.0SHS", {100.0, 100.0, 4.0}},
            {"100X3.0SHS", {100.0, 100.0, 3.0}},
            {"89X6.0SHS", {89.0, 89.0, 6.0}},
            {"89X5.0SHS", {89.0, 89.0, 5.0}},
            {"89X3.5SHS", {89.0, 89.0, 3.5}},
            {"75X6.0SHS", {75.0, 75.0, 6.0}},
            {"75X5.0SHS", {75.0, 75.0, 5.0}},
            {"75X4.0SHS", {75.0, 75.0, 4.0}},
            {"75X3.5SHS", {75.0, 75.0, 3.5}},
            {"75X3.0SHS", {75.0, 75.0, 3.0}},
            {"75X2.5SHS", {75.0, 75.0, 2.5}},
            {"65X3.0SHS", {65.0, 65.0, 3.0}},
            {"65X2.5SHS", {65.0, 65.0, 2.5}},
            {"65X2.0SHS", {65.0, 65.0, 2.0}},
            {"50X4.0SHS", {50.0, 50.0, 4.0}},
            {"50X3.0SHS", {50.0, 50.0, 3.0}},
            {"50X2.5SHS", {50.0, 50.0, 2.5}},
            {"50X2.0SHS", {50.0, 50.0, 2.0}},
            {"50X1.6SHS", {50.0, 50.0, 1.6}},
            {"40X4.0SHS", {40.0, 40.0, 4.0}},
            {"40X2.5SHS", {40.0, 40.0, 2.5}},
            {"40X2.0SHS", {40.0, 40.0, 2.0}},
            {"40X1.6SHS", {40.0, 40.0, 1.6}},
            {"35X3.0SHS", {35.0, 35.0, 3.0}},
            {"35X2.5SHS", {35.0, 35.0, 2.5}},
            {"35X2.0SHS", {35.0, 35.0, 2.0}},
            {"35X1.6SHS", {35.0, 35.0, 1.6}},
            {"30X2.0SHS", {30.0, 30.0, 2.0}},
            {"30X1.6SHS", {30.0, 30.0, 1.6}},
            {"25X3.0SHS", {25.0, 25.0, 3.0}},
            {"25X2.5SHS", {25.0, 25.0, 2.5}},
            {"25X2.0SHS", {25.0, 25.0, 2.0}},
            {"25X1.6SHS", {25.0, 25.0, 1.6}},
            {"20X1.6SHS", {20.0, 20.0, 1.6}},
        };

        return registry.fetch(type);
    }

    vec usisection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"W44X408", {16.100, 2.170, 16.100, 2.170, 40.460, 1.220}},
            {"W44X368", {16.000, 1.970, 16.000, 1.970, 40.460, 1.100}},
            {"W44X335", {15.900, 1.770, 15.900, 1.770, 40.460, 1.030}},
            {"W44X290", {15.800, 1.580, 15.800, 1.580, 40.440, 0.865}},
            {"W44X262", {15.800, 1.420, 15.800, 1.420, 40.460, 0.785}},
            {"W44X230", {15.800, 1.220, 15.800, 1.220, 40.460, 0.710}},
            {"W40X655", {16.900, 3.540, 16.900, 3.540, 36.520, 1.970}},
            {"W40X593", {16.700, 3.230, 16.700, 3.230, 36.540, 1.790}},
            {"W40X503", {16.400, 2.760, 16.400, 2.760, 36.580, 1.540}},
            {"W40X431", {16.200, 2.360, 16.200, 2.360, 36.580, 1.340}},
            {"W40X397", {16.100, 2.200, 16.100, 2.200, 36.600, 1.220}},
            {"W40X372", {16.100, 2.050, 16.100, 2.050, 36.500, 1.160}},
            {"W40X362", {16.000, 2.010, 16.000, 2.010, 36.580, 1.120}},
            {"W40X324", {15.900, 1.810, 15.900, 1.810, 36.580, 1.000}},
            {"W40X297", {15.800, 1.650, 15.800, 1.650, 36.500, 0.930}},
            {"W40X277", {15.800, 1.580, 15.800, 1.580, 36.540, 0.830}},
            {"W40X249", {15.800, 1.420, 15.800, 1.420, 36.560, 0.750}},
            {"W40X215", {15.800, 1.220, 15.800, 1.220, 36.560, 0.650}},
            {"W40X199", {15.800, 1.070, 15.800, 1.070, 36.560, 0.650}},
            {"W40X392", {12.400, 2.520, 12.400, 2.520, 36.560, 1.420}},
            {"W40X331", {12.200, 2.130, 12.200, 2.130, 36.540, 1.220}},
            {"W40X327", {12.100, 2.130, 12.100, 2.130, 36.540, 1.180}},
            {"W40X294", {12.000, 1.930, 12.000, 1.930, 36.540, 1.060}},
            {"W40X278", {12.000, 1.810, 12.000, 1.810, 36.580, 1.030}},
            {"W40X264", {11.900, 1.730, 11.900, 1.730, 36.540, 0.960}},
            {"W40X235", {11.900, 1.580, 11.900, 1.580, 36.540, 0.830}},
            {"W40X211", {11.800, 1.420, 11.800, 1.420, 36.560, 0.750}},
            {"W40X183", {11.800, 1.200, 11.800, 1.200, 36.600, 0.650}},
            {"W40X167", {11.800, 1.030, 11.800, 1.030, 36.540, 0.650}},
            {"W40X149", {11.800, 0.830, 11.800, 0.830, 36.540, 0.630}},
            {"W36X925", {18.600, 4.530, 18.600, 4.530, 34.040, 3.020}},
            {"W36X853", {18.200, 4.530, 18.200, 4.530, 34.040, 2.520}},
            {"W36X802", {18.000, 4.290, 18.000, 4.290, 34.020, 2.380}},
            {"W36X723", {17.800, 3.900, 17.800, 3.900, 34.000, 2.170}},
            {"W36X652", {17.600, 3.540, 17.600, 3.540, 34.020, 1.970}},
            {"W36X529", {17.200, 2.910, 17.200, 2.910, 33.980, 1.610}},
            {"W36X487", {17.100, 2.680, 17.100, 2.680, 33.940, 1.500}},
            {"W36X441", {17.000, 2.440, 17.000, 2.440, 34.020, 1.360}},
            {"W36X395", {16.800, 2.200, 16.800, 2.200, 34.000, 1.220}},
            {"W36X361", {16.700, 2.010, 16.700, 2.010, 33.980, 1.120}},
            {"W36X330", {16.600, 1.850, 16.600, 1.850, 34.000, 1.020}},
            {"W36X302", {16.700, 1.680, 16.700, 1.680, 33.940, 0.945}},
            {"W36X282", {16.600, 1.570, 16.600, 1.570, 33.960, 0.885}},
            {"W36X262", {16.600, 1.440, 16.600, 1.440, 34.020, 0.840}},
            {"W36X247", {16.500, 1.350, 16.500, 1.350, 34.000, 0.800}},
            {"W36X231", {16.500, 1.260, 16.500, 1.260, 33.980, 0.760}},
            {"W36X387", {12.700, 2.560, 12.700, 2.560, 33.980, 1.420}},
            {"W36X350", {12.600, 2.320, 12.600, 2.320, 33.960, 1.300}},
            {"W36X318", {12.400, 2.130, 12.400, 2.130, 33.940, 1.180}},
            {"W36X286", {12.300, 1.930, 12.300, 1.930, 33.940, 1.060}},
            {"W36X256", {12.200, 1.730, 12.200, 1.730, 33.940, 0.960}},
            {"W36X232", {12.100, 1.570, 12.100, 1.570, 33.960, 0.870}},
            {"W36X210", {12.200, 1.360, 12.200, 1.360, 33.980, 0.830}},
            {"W36X194", {12.100, 1.260, 12.100, 1.260, 33.980, 0.765}},
            {"W36X182", {12.100, 1.180, 12.100, 1.180, 33.940, 0.725}},
            {"W36X170", {12.000, 1.100, 12.000, 1.100, 34.000, 0.680}},
            {"W36X160", {12.000, 1.020, 12.000, 1.020, 33.960, 0.650}},
            {"W36X150", {12.000, 0.940, 12.000, 0.940, 34.020, 0.625}},
            {"W36X135", {12.000, 0.790, 12.000, 0.790, 34.020, 0.600}},
            {"W33X387", {16.200, 2.280, 16.200, 2.280, 31.440, 1.260}},
            {"W33X354", {16.100, 2.090, 16.100, 2.090, 31.420, 1.160}},
            {"W33X318", {16.000, 1.890, 16.000, 1.890, 31.420, 1.040}},
            {"W33X291", {15.900, 1.730, 15.900, 1.730, 31.340, 0.960}},
            {"W33X263", {15.800, 1.570, 15.800, 1.570, 31.360, 0.870}},
            {"W33X241", {15.900, 1.400, 15.900, 1.400, 31.400, 0.830}},
            {"W33X221", {15.800, 1.280, 15.800, 1.280, 31.340, 0.775}},
            {"W33X201", {15.700, 1.150, 15.700, 1.150, 31.400, 0.715}},
            {"W33X169", {11.500, 1.220, 11.500, 1.220, 31.360, 0.670}},
            {"W33X152", {11.600, 1.060, 11.600, 1.060, 31.380, 0.635}},
            {"W33X141", {11.500, 0.960, 11.500, 0.960, 31.380, 0.605}},
            {"W33X130", {11.500, 0.855, 11.500, 0.855, 31.390, 0.580}},
            {"W33X118", {11.500, 0.740, 11.500, 0.740, 31.420, 0.550}},
            {"W30X391", {15.600, 2.440, 15.600, 2.440, 28.320, 1.360}},
            {"W30X357", {15.500, 2.240, 15.500, 2.240, 28.320, 1.240}},
            {"W30X326", {15.400, 2.050, 15.400, 2.050, 28.300, 1.140}},
            {"W30X292", {15.300, 1.850, 15.300, 1.850, 28.300, 1.020}},
            {"W30X261", {15.200, 1.650, 15.200, 1.650, 28.300, 0.930}},
            {"W30X235", {15.100, 1.500, 15.100, 1.500, 28.300, 0.830}},
            {"W30X211", {15.100, 1.320, 15.100, 1.320, 28.260, 0.775}},
            {"W30X191", {15.000, 1.190, 15.000, 1.190, 28.320, 0.710}},
            {"W30X173", {15.000, 1.070, 15.000, 1.070, 28.260, 0.655}},
            {"W30X148", {10.500, 1.180, 10.500, 1.180, 28.340, 0.650}},
            {"W30X132", {10.500, 1.000, 10.500, 1.000, 28.300, 0.615}},
            {"W30X124", {10.500, 0.930, 10.500, 0.930, 28.340, 0.585}},
            {"W30X116", {10.500, 0.850, 10.500, 0.850, 28.300, 0.565}},
            {"W30X108", {10.500, 0.760, 10.500, 0.760, 28.280, 0.545}},
            {"W30X99", {10.500, 0.670, 10.500, 0.670, 28.360, 0.520}},
            {"W30X90", {10.400, 0.610, 10.400, 0.610, 28.280, 0.470}},
            {"W27X539", {15.300, 3.540, 15.300, 3.540, 25.420, 1.970}},
            {"W27X368", {14.700, 2.480, 14.700, 2.480, 25.440, 1.380}},
            {"W27X336", {14.600, 2.280, 14.600, 2.280, 25.440, 1.260}},
            {"W27X307", {14.400, 2.090, 14.400, 2.090, 25.420, 1.160}},
            {"W27X281", {14.400, 1.930, 14.400, 1.930, 25.440, 1.060}},
            {"W27X258", {14.300, 1.770, 14.300, 1.770, 25.460, 0.980}},
            {"W27X235", {14.200, 1.610, 14.200, 1.610, 25.480, 0.910}},
            {"W27X217", {14.100, 1.500, 14.100, 1.500, 25.400, 0.830}},
            {"W27X194", {14.000, 1.340, 14.000, 1.340, 25.420, 0.750}},
            {"W27X178", {14.100, 1.190, 14.100, 1.190, 25.420, 0.725}},
            {"W27X161", {14.000, 1.080, 14.000, 1.080, 25.440, 0.660}},
            {"W27X146", {14.000, 0.975, 14.000, 0.975, 25.450, 0.605}},
            {"W27X129", {10.000, 1.100, 10.000, 1.100, 25.400, 0.610}},
            {"W27X114", {10.100, 0.930, 10.100, 0.930, 25.440, 0.570}},
            {"W27X102", {10.000, 0.830, 10.000, 0.830, 25.440, 0.515}},
            {"W27X94", {10.000, 0.745, 10.000, 0.745, 25.410, 0.490}},
            {"W27X84", {10.000, 0.640, 10.000, 0.640, 25.420, 0.460}},
            {"W24X370", {13.700, 2.720, 13.700, 2.720, 22.560, 1.520}},
            {"W24X335", {13.500, 2.480, 13.500, 2.480, 22.540, 1.380}},
            {"W24X306", {13.400, 2.280, 13.400, 2.280, 22.540, 1.260}},
            {"W24X279", {13.300, 2.090, 13.300, 2.090, 22.520, 1.160}},
            {"W24X250", {13.200, 1.890, 13.200, 1.890, 22.520, 1.040}},
            {"W24X229", {13.100, 1.730, 13.100, 1.730, 22.540, 0.960}},
            {"W24X207", {13.000, 1.570, 13.000, 1.570, 22.560, 0.870}},
            {"W24X192", {13.000, 1.460, 13.000, 1.460, 22.580, 0.810}},
            {"W24X176", {12.900, 1.340, 12.900, 1.340, 22.520, 0.750}},
            {"W24X162", {13.000, 1.220, 13.000, 1.220, 22.560, 0.705}},
            {"W24X146", {12.900, 1.090, 12.900, 1.090, 22.520, 0.650}},
            {"W24X131", {12.900, 0.960, 12.900, 0.960, 22.580, 0.605}},
            {"W24X117", {12.800, 0.850, 12.800, 0.850, 22.600, 0.550}},
            {"W24X104", {12.800, 0.750, 12.800, 0.750, 22.600, 0.500}},
            {"W24X103", {9.000, 0.980, 9.000, 0.980, 22.540, 0.550}},
            {"W24X94", {9.070, 0.875, 9.070, 0.875, 22.550, 0.515}},
            {"W24X84", {9.020, 0.770, 9.020, 0.770, 22.560, 0.470}},
            {"W24X76", {8.990, 0.680, 8.990, 0.680, 22.540, 0.440}},
            {"W24X68", {8.970, 0.585, 8.970, 0.585, 22.530, 0.415}},
            {"W24X62", {7.040, 0.590, 7.040, 0.590, 22.520, 0.430}},
            {"W24X55", {7.010, 0.505, 7.010, 0.505, 22.590, 0.395}},
            {"W21X275", {12.900, 2.190, 12.900, 2.190, 19.720, 1.220}},
            {"W21X248", {12.800, 1.990, 12.800, 1.990, 19.720, 1.100}},
            {"W21X223", {12.700, 1.790, 12.700, 1.790, 19.820, 1.000}},
            {"W21X201", {12.600, 1.630, 12.600, 1.630, 19.740, 0.910}},
            {"W21X182", {12.500, 1.480, 12.500, 1.480, 19.740, 0.830}},
            {"W21X166", {12.400, 1.360, 12.400, 1.360, 19.780, 0.750}},
            {"W21X147", {12.500, 1.150, 12.500, 1.150, 19.800, 0.720}},
            {"W21X132", {12.400, 1.040, 12.400, 1.040, 19.720, 0.650}},
            {"W21X122", {12.400, 0.960, 12.400, 0.960, 19.780, 0.600}},
            {"W21X111", {12.300, 0.875, 12.300, 0.875, 19.750, 0.550}},
            {"W21X101", {12.300, 0.800, 12.300, 0.800, 19.800, 0.500}},
            {"W21X93", {8.420, 0.930, 8.420, 0.930, 19.740, 0.580}},
            {"W21X83", {8.360, 0.835, 8.360, 0.835, 19.730, 0.515}},
            {"W21X73", {8.300, 0.740, 8.300, 0.740, 19.720, 0.455}},
            {"W21X68", {8.270, 0.685, 8.270, 0.685, 19.730, 0.430}},
            {"W21X62", {8.240, 0.615, 8.240, 0.615, 19.770, 0.400}},
            {"W21X55", {8.220, 0.522, 8.220, 0.522, 19.756, 0.375}},
            {"W21X48", {8.140, 0.430, 8.140, 0.430, 19.740, 0.350}},
            {"W21X57", {6.560, 0.650, 6.560, 0.650, 19.800, 0.405}},
            {"W21X50", {6.530, 0.535, 6.530, 0.535, 19.730, 0.380}},
            {"W21X44", {6.500, 0.450, 6.500, 0.450, 19.800, 0.350}},
            {"W18X311", {12.000, 2.740, 12.000, 2.740, 16.820, 1.520}},
            {"W18X283", {11.900, 2.500, 11.900, 2.500, 16.900, 1.400}},
            {"W18X258", {11.800, 2.300, 11.800, 2.300, 16.900, 1.280}},
            {"W18X234", {11.700, 2.110, 11.700, 2.110, 16.880, 1.160}},
            {"W18X211", {11.600, 1.910, 11.600, 1.910, 16.880, 1.060}},
            {"W18X192", {11.500, 1.750, 11.500, 1.750, 16.900, 0.960}},
            {"W18X175", {11.400, 1.590, 11.400, 1.590, 16.820, 0.890}},
            {"W18X158", {11.300, 1.440, 11.300, 1.440, 16.820, 0.810}},
            {"W18X143", {11.200, 1.320, 11.200, 1.320, 16.860, 0.730}},
            {"W18X130", {11.200, 1.200, 11.200, 1.200, 16.900, 0.670}},
            {"W18X119", {11.300, 1.060, 11.300, 1.060, 16.880, 0.655}},
            {"W18X106", {11.200, 0.940, 11.200, 0.940, 16.820, 0.590}},
            {"W18X97", {11.100, 0.870, 11.100, 0.870, 16.860, 0.535}},
            {"W18X86", {11.100, 0.770, 11.100, 0.770, 16.860, 0.480}},
            {"W18X76", {11.000, 0.680, 11.000, 0.680, 16.840, 0.425}},
            {"W18X71", {7.640, 0.810, 7.640, 0.810, 16.880, 0.495}},
            {"W18X65", {7.590, 0.750, 7.590, 0.750, 16.900, 0.450}},
            {"W18X60", {7.560, 0.695, 7.560, 0.695, 16.810, 0.415}},
            {"W18X55", {7.530, 0.630, 7.530, 0.630, 16.840, 0.390}},
            {"W18X50", {7.500, 0.570, 7.500, 0.570, 16.860, 0.355}},
            {"W18X46", {6.060, 0.605, 6.060, 0.605, 16.890, 0.360}},
            {"W18X40", {6.020, 0.525, 6.020, 0.525, 16.850, 0.315}},
            {"W18X35", {6.000, 0.425, 6.000, 0.425, 16.850, 0.300}},
            {"W16X100", {10.400, 0.985, 10.400, 0.985, 15.030, 0.585}},
            {"W16X89", {10.400, 0.875, 10.400, 0.875, 15.050, 0.525}},
            {"W16X77", {10.300, 0.760, 10.300, 0.760, 14.980, 0.455}},
            {"W16X67", {10.200, 0.665, 10.200, 0.665, 14.970, 0.395}},
            {"W16X57", {7.120, 0.715, 7.120, 0.715, 14.970, 0.430}},
            {"W16X50", {7.070, 0.630, 7.070, 0.630, 15.040, 0.380}},
            {"W16X45", {7.040, 0.565, 7.040, 0.565, 14.970, 0.345}},
            {"W16X40", {7.000, 0.505, 7.000, 0.505, 14.990, 0.305}},
            {"W16X36", {6.990, 0.430, 6.990, 0.430, 15.040, 0.295}},
            {"W16X31", {5.530, 0.440, 5.530, 0.440, 15.020, 0.275}},
            {"W16X26", {5.500, 0.345, 5.500, 0.345, 15.010, 0.250}},
            {"W14X873", {18.800, 5.510, 18.800, 5.510, 12.580, 3.940}},
            {"W14X808", {18.600, 5.120, 18.600, 5.120, 12.560, 3.740}},
            {"W14X730", {17.900, 4.910, 17.900, 4.910, 12.580, 3.070}},
            {"W14X665", {17.700, 4.520, 17.700, 4.520, 12.560, 2.830}},
            {"W14X605", {17.400, 4.160, 17.400, 4.160, 12.580, 2.600}},
            {"W14X550", {17.200, 3.820, 17.200, 3.820, 12.560, 2.380}},
            {"W14X500", {17.000, 3.500, 17.000, 3.500, 12.600, 2.190}},
            {"W14X455", {16.800, 3.210, 16.800, 3.210, 12.580, 2.020}},
            {"W14X426", {16.700, 3.040, 16.700, 3.040, 12.620, 1.880}},
            {"W14X398", {16.600, 2.850, 16.600, 2.850, 12.600, 1.770}},
            {"W14X370", {16.500, 2.660, 16.500, 2.660, 12.580, 1.660}},
            {"W14X342", {16.400, 2.470, 16.400, 2.470, 12.560, 1.540}},
            {"W14X311", {16.200, 2.260, 16.200, 2.260, 12.580, 1.410}},
            {"W14X283", {16.100, 2.070, 16.100, 2.070, 12.560, 1.290}},
            {"W14X257", {16.000, 1.890, 16.000, 1.890, 12.620, 1.180}},
            {"W14X233", {15.900, 1.720, 15.900, 1.720, 12.560, 1.070}},
            {"W14X211", {15.800, 1.560, 15.800, 1.560, 12.580, 0.980}},
            {"W14X193", {15.700, 1.440, 15.700, 1.440, 12.620, 0.890}},
            {"W14X176", {15.700, 1.310, 15.700, 1.310, 12.580, 0.830}},
            {"W14X159", {15.600, 1.190, 15.600, 1.190, 12.620, 0.745}},
            {"W14X145", {15.500, 1.090, 15.500, 1.090, 12.620, 0.680}},
            {"W14X132", {14.700, 1.030, 14.700, 1.030, 12.640, 0.645}},
            {"W14X120", {14.700, 0.940, 14.700, 0.940, 12.620, 0.590}},
            {"W14X109", {14.600, 0.860, 14.600, 0.860, 12.580, 0.525}},
            {"W14X99", {14.600, 0.780, 14.600, 0.780, 12.640, 0.485}},
            {"W14X90", {14.500, 0.710, 14.500, 0.710, 12.580, 0.440}},
            {"W14X82", {10.100, 0.855, 10.100, 0.855, 12.590, 0.510}},
            {"W14X74", {10.100, 0.785, 10.100, 0.785, 12.630, 0.450}},
            {"W14X68", {10.000, 0.720, 10.000, 0.720, 12.560, 0.415}},
            {"W14X61", {10.000, 0.645, 10.000, 0.645, 12.610, 0.375}},
            {"W14X53", {8.060, 0.660, 8.060, 0.660, 12.580, 0.370}},
            {"W14X48", {8.030, 0.595, 8.030, 0.595, 12.610, 0.340}},
            {"W14X43", {8.000, 0.530, 8.000, 0.530, 12.640, 0.305}},
            {"W14X38", {6.770, 0.515, 6.770, 0.515, 13.070, 0.310}},
            {"W14X34", {6.750, 0.455, 6.750, 0.455, 13.090, 0.285}},
            {"W14X30", {6.730, 0.385, 6.730, 0.385, 13.030, 0.270}},
            {"W14X26", {5.030, 0.420, 5.030, 0.420, 13.060, 0.255}},
            {"W14X22", {5.000, 0.335, 5.000, 0.335, 13.030, 0.230}},
            {"W12X336", {13.400, 2.960, 13.400, 2.960, 10.880, 1.780}},
            {"W12X305", {13.200, 2.710, 13.200, 2.710, 10.880, 1.630}},
            {"W12X279", {13.100, 2.470, 13.100, 2.470, 10.960, 1.530}},
            {"W12X252", {13.000, 2.250, 13.000, 2.250, 10.900, 1.400}},
            {"W12X230", {12.900, 2.070, 12.900, 2.070, 10.960, 1.290}},
            {"W12X210", {12.800, 1.900, 12.800, 1.900, 10.900, 1.180}},
            {"W12X190", {12.700, 1.740, 12.700, 1.740, 10.920, 1.060}},
            {"W12X170", {12.600, 1.560, 12.600, 1.560, 10.880, 0.960}},
            {"W12X152", {12.500, 1.400, 12.500, 1.400, 10.900, 0.870}},
            {"W12X136", {12.400, 1.250, 12.400, 1.250, 10.900, 0.790}},
            {"W12X120", {12.300, 1.110, 12.300, 1.110, 10.880, 0.710}},
            {"W12X106", {12.200, 0.990, 12.200, 0.990, 10.920, 0.610}},
            {"W12X96", {12.200, 0.900, 12.200, 0.900, 10.900, 0.550}},
            {"W12X87", {12.100, 0.810, 12.100, 0.810, 10.880, 0.515}},
            {"W12X79", {12.100, 0.735, 12.100, 0.735, 10.930, 0.470}},
            {"W12X72", {12.000, 0.670, 12.000, 0.670, 10.960, 0.430}},
            {"W12X65", {12.000, 0.605, 12.000, 0.605, 10.890, 0.390}},
            {"W12X58", {10.000, 0.640, 10.000, 0.640, 10.920, 0.360}},
            {"W12X53", {10.000, 0.575, 10.000, 0.575, 10.950, 0.345}},
            {"W12X50", {8.080, 0.640, 8.080, 0.640, 10.920, 0.370}},
            {"W12X45", {8.050, 0.575, 8.050, 0.575, 10.950, 0.335}},
            {"W12X40", {8.010, 0.515, 8.010, 0.515, 10.870, 0.295}},
            {"W12X35", {6.560, 0.520, 6.560, 0.520, 11.460, 0.300}},
            {"W12X30", {6.520, 0.440, 6.520, 0.440, 11.420, 0.260}},
            {"W12X26", {6.490, 0.380, 6.490, 0.380, 11.440, 0.230}},
            {"W12X22", {4.030, 0.425, 4.030, 0.425, 11.450, 0.260}},
            {"W12X19", {4.010, 0.350, 4.010, 0.350, 11.500, 0.235}},
            {"W12X16", {3.990, 0.265, 3.990, 0.265, 11.470, 0.220}},
            {"W12X14", {3.970, 0.225, 3.970, 0.225, 11.450, 0.200}},
            {"W10X112", {10.400, 1.250, 10.400, 1.250, 8.900, 0.755}},
            {"W10X100", {10.300, 1.120, 10.300, 1.120, 8.860, 0.680}},
            {"W10X88", {10.300, 0.990, 10.300, 0.990, 8.820, 0.605}},
            {"W10X77", {10.200, 0.870, 10.200, 0.870, 8.860, 0.530}},
            {"W10X68", {10.100, 0.770, 10.100, 0.770, 8.860, 0.470}},
            {"W10X60", {10.100, 0.680, 10.100, 0.680, 8.840, 0.420}},
            {"W10X54", {10.000, 0.615, 10.000, 0.615, 8.870, 0.370}},
            {"W10X49", {10.000, 0.560, 10.000, 0.560, 8.880, 0.340}},
            {"W10X45", {8.020, 0.620, 8.020, 0.620, 8.860, 0.350}},
            {"W10X39", {7.990, 0.530, 7.990, 0.530, 8.860, 0.315}},
            {"W10X33", {7.960, 0.435, 7.960, 0.435, 8.860, 0.290}},
            {"W10X30", {5.810, 0.510, 5.810, 0.510, 9.480, 0.300}},
            {"W10X26", {5.770, 0.440, 5.770, 0.440, 9.420, 0.260}},
            {"W10X22", {5.750, 0.360, 5.750, 0.360, 9.480, 0.240}},
            {"W10X19", {4.020, 0.395, 4.020, 0.395, 9.410, 0.250}},
            {"W10X17", {4.010, 0.330, 4.010, 0.330, 9.440, 0.240}},
            {"W10X15", {4.000, 0.270, 4.000, 0.270, 9.450, 0.230}},
            {"W10X12", {3.960, 0.210, 3.960, 0.210, 9.450, 0.190}},
            {"W8X67", {8.280, 0.935, 8.280, 0.935, 7.130, 0.570}},
            {"W8X58", {8.220, 0.810, 8.220, 0.810, 7.130, 0.510}},
            {"W8X48", {8.110, 0.685, 8.110, 0.685, 7.130, 0.400}},
            {"W8X40", {8.070, 0.560, 8.070, 0.560, 7.130, 0.360}},
            {"W8X35", {8.020, 0.495, 8.020, 0.495, 7.130, 0.310}},
            {"W8X31", {8.000, 0.435, 8.000, 0.435, 7.130, 0.285}},
            {"W8X28", {6.540, 0.465, 6.540, 0.465, 7.130, 0.285}},
            {"W8X24", {6.500, 0.400, 6.500, 0.400, 7.130, 0.245}},
            {"W8X21", {5.270, 0.400, 5.270, 0.400, 7.480, 0.250}},
            {"W8X18", {5.250, 0.330, 5.250, 0.330, 7.480, 0.230}},
            {"W8X15", {4.020, 0.315, 4.020, 0.315, 7.480, 0.245}},
            {"W8X13", {4.000, 0.255, 4.000, 0.255, 7.480, 0.230}},
            {"W8X10", {3.940, 0.205, 3.940, 0.205, 7.480, 0.170}},
            {"W6X25", {6.080, 0.455, 6.080, 0.455, 5.470, 0.320}},
            {"W6X20", {6.020, 0.365, 6.020, 0.365, 5.470, 0.260}},
            {"W6X15", {5.990, 0.260, 5.990, 0.260, 5.470, 0.230}},
            {"W6X16", {4.030, 0.405, 4.030, 0.405, 5.470, 0.260}},
            {"W6X12", {4.000, 0.280, 4.000, 0.280, 5.470, 0.230}},
            {"W6X9", {3.940, 0.215, 3.940, 0.215, 5.470, 0.170}},
            {"W6X8.5", {3.940, 0.195, 3.940, 0.195, 5.440, 0.170}},
            {"W5X19", {5.030, 0.430, 5.030, 0.430, 4.290, 0.270}},
            {"W5X16", {5.000, 0.360, 5.000, 0.360, 4.290, 0.240}},
            {"W4X13", {4.060, 0.345, 4.060, 0.345, 3.470, 0.280}},
            {"M12.5X12.4", {3.750, 0.228, 3.750, 0.228, 12.044, 0.155}},
            {"M12.5X11.6", {3.500, 0.211, 3.500, 0.211, 12.078, 0.155}},
            {"M12X11.8", {3.070, 0.225, 3.070, 0.225, 11.550, 0.177}},
            {"M12X10.8", {3.070, 0.210, 3.070, 0.210, 11.580, 0.160}},
            {"M12X10", {3.250, 0.180, 3.250, 0.180, 11.640, 0.149}},
            {"M10X9", {2.690, 0.206, 2.690, 0.206, 9.588, 0.157}},
            {"M10X8", {2.690, 0.182, 2.690, 0.182, 9.586, 0.141}},
            {"M10X7.5", {2.690, 0.173, 2.690, 0.173, 9.644, 0.130}},
            {"M8X6.5", {2.280, 0.189, 2.280, 0.189, 7.622, 0.135}},
            {"M8X6.2", {2.280, 0.177, 2.280, 0.177, 7.646, 0.129}},
            {"M6X4.4", {1.840, 0.171, 1.840, 0.171, 5.658, 0.114}},
            {"M6X3.7", {2.000, 0.129, 2.000, 0.129, 5.662, 0.098}},
            {"M5X18.9", {5.000, 0.416, 5.000, 0.416, 4.168, 0.316}},
            {"M4X6", {3.800, 0.160, 3.800, 0.160, 3.480, 0.130}},
            {"M4X4.08", {2.250, 0.170, 2.250, 0.170, 3.660, 0.115}},
            {"M3X2.9", {2.250, 0.130, 2.250, 0.130, 2.740, 0.090}},
            {"S24X121", {8.050, 1.090, 8.050, 1.090, 22.320, 0.800}},
            {"S24X106", {7.870, 1.090, 7.870, 1.090, 22.320, 0.620}},
            {"S24X100", {7.250, 0.870, 7.250, 0.870, 22.260, 0.745}},
            {"S24X90", {7.130, 0.870, 7.130, 0.870, 22.260, 0.625}},
            {"S24X80", {7.000, 0.870, 7.000, 0.870, 22.260, 0.500}},
            {"S20X96", {7.200, 0.920, 7.200, 0.920, 18.460, 0.800}},
            {"S20X86", {7.060, 0.920, 7.060, 0.920, 18.460, 0.660}},
            {"S20X75", {6.390, 0.795, 6.390, 0.795, 18.410, 0.635}},
            {"S20X66", {6.260, 0.795, 6.260, 0.795, 18.410, 0.505}},
            {"S18X70", {6.250, 0.691, 6.250, 0.691, 16.618, 0.711}},
            {"S18X54.7", {6.000, 0.691, 6.000, 0.691, 16.618, 0.461}},
            {"S15X50", {5.640, 0.622, 5.640, 0.622, 13.756, 0.550}},
            {"S15X42.9", {5.500, 0.622, 5.500, 0.622, 13.756, 0.411}},
            {"S12X50", {5.480, 0.659, 5.480, 0.659, 10.682, 0.687}},
            {"S12X40.8", {5.250, 0.659, 5.250, 0.659, 10.682, 0.462}},
            {"S12X35", {5.080, 0.544, 5.080, 0.544, 10.912, 0.428}},
            {"S12X31.8", {5.000, 0.544, 5.000, 0.544, 10.912, 0.350}},
            {"S10X35", {4.940, 0.491, 4.940, 0.491, 9.018, 0.594}},
            {"S10X25.4", {4.660, 0.491, 4.660, 0.491, 9.018, 0.311}},
            {"S8X23", {4.170, 0.425, 4.170, 0.425, 7.150, 0.441}},
            {"S8X18.4", {4.000, 0.425, 4.000, 0.425, 7.150, 0.271}},
            {"S6X17.25", {3.570, 0.359, 3.570, 0.359, 5.282, 0.465}},
            {"S6X12.5", {3.330, 0.359, 3.330, 0.359, 5.282, 0.232}},
            {"S5X10", {3.000, 0.326, 3.000, 0.326, 4.348, 0.214}},
            {"S4X9.5", {2.800, 0.293, 2.800, 0.293, 3.414, 0.326}},
            {"S4X7.7", {2.660, 0.293, 2.660, 0.293, 3.414, 0.193}},
            {"S3X7.5", {2.510, 0.260, 2.510, 0.260, 2.480, 0.349}},
            {"S3X5.7", {2.330, 0.260, 2.330, 0.260, 2.480, 0.170}},
            {"HP18X204", {18.100, 1.130, 18.100, 1.130, 16.040, 1.130}},
            {"HP18X181", {18.000, 1.000, 18.000, 1.000, 16.000, 1.000}},
            {"HP18X157", {17.900, 0.870, 17.900, 0.870, 15.960, 0.870}},
            {"HP18X135", {17.800, 0.750, 17.800, 0.750, 16.000, 0.750}},
            {"HP16X183", {16.300, 1.130, 16.300, 1.130, 14.240, 1.130}},
            {"HP16X162", {16.100, 1.000, 16.100, 1.000, 14.300, 1.000}},
            {"HP16X141", {16.000, 0.875, 16.000, 0.875, 14.250, 0.875}},
            {"HP16X121", {15.900, 0.750, 15.900, 0.750, 14.300, 0.750}},
            {"HP16X101", {15.800, 0.625, 15.800, 0.625, 14.250, 0.625}},
            {"HP16X88", {15.700, 0.540, 15.700, 0.540, 14.220, 0.540}},
            {"HP14X117", {14.900, 0.805, 14.900, 0.805, 12.590, 0.805}},
            {"HP14X102", {14.800, 0.705, 14.800, 0.705, 12.590, 0.705}},
            {"HP14X89", {14.700, 0.615, 14.700, 0.615, 12.570, 0.615}},
            {"HP14X73", {14.600, 0.505, 14.600, 0.505, 12.590, 0.505}},
            {"HP12X89", {12.300, 0.720, 12.300, 0.720, 10.960, 0.720}},
            {"HP12X84", {12.300, 0.685, 12.300, 0.685, 10.930, 0.685}},
            {"HP12X74", {12.200, 0.610, 12.200, 0.610, 10.880, 0.605}},
            {"HP12X63", {12.100, 0.515, 12.100, 0.515, 10.870, 0.515}},
            {"HP12X53", {12.000, 0.435, 12.000, 0.435, 10.930, 0.435}},
            {"HP10X57", {10.200, 0.565, 10.200, 0.565, 8.860, 0.565}},
            {"HP10X42", {10.100, 0.420, 10.100, 0.420, 8.860, 0.415}},
            {"HP8X36", {8.160, 0.445, 8.160, 0.445, 7.130, 0.445}},
        };

        return registry.fetch(type);
    }

    // Adjusted registry: dimensions optimized against the AISC tabulated
    // section properties with section_optimizer.py.  The M and S entries use
    // a 7 % dimension margin; strong-axis properties stay within 5 %.
    vec usisection_adjusted(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"W44X408", {16.1301, 2.1638, 16.1301, 2.1638, 40.5299, 1.239}},
            {"W44X368", {15.9694, 1.982, 15.9694, 1.982, 40.5136, 1.1042}},
            {"W44X335", {15.9467, 1.7701, 15.9467, 1.7701, 40.5744, 1.0355}},
            {"W44X290", {15.7845, 1.5833, 15.7845, 1.5833, 40.5109, 0.8739}},
            {"W44X262", {15.7412, 1.4167, 15.7412, 1.4167, 40.5553, 0.8053}},
            {"W44X230", {15.7156, 1.2286, 15.7156, 1.2286, 40.4905, 0.7203}},
            {"W40X655", {16.8778, 3.5502, 16.8778, 3.5502, 36.5831, 2.0023}},
            {"W40X593", {16.6138, 3.2745, 16.6138, 3.2745, 36.5135, 1.7857}},
            {"W40X503", {16.3479, 2.7844, 16.3479, 2.7844, 36.5662, 1.5613}},
            {"W40X431", {16.233, 2.3588, 16.233, 2.3588, 36.6562, 1.3774}},
            {"W40X397", {16.1221, 2.1964, 16.1221, 2.1964, 36.6978, 1.2579}},
            {"W40X372", {16.0974, 2.0344, 16.0974, 2.0344, 36.7185, 1.2132}},
            {"W40X362", {15.9464, 2.0345, 15.9464, 2.0345, 36.6758, 1.1244}},
            {"W40X324", {15.9002, 1.8155, 15.9002, 1.8155, 36.684, 1.0261}},
            {"W40X297", {15.7657, 1.6644, 15.7657, 1.6644, 36.6924, 0.951}},
            {"W40X277", {15.7788, 1.5851, 15.7788, 1.5851, 36.7342, 0.8582}},
            {"W40X249", {15.7214, 1.4276, 15.7214, 1.4276, 36.7642, 0.7781}},
            {"W40X215", {15.8282, 1.2127, 15.8282, 1.2127, 36.8184, 0.6825}},
            {"W40X199", {15.6822, 1.0797, 15.6822, 1.0797, 36.7183, 0.6791}},
            {"W40X392", {12.3849, 2.5056, 12.3849, 2.5056, 36.6755, 1.4705}},
            {"W40X331", {12.1273, 2.1462, 12.1273, 2.1462, 36.6461, 1.2466}},
            {"W40X327", {12.1556, 2.12, 12.1556, 2.12, 36.767, 1.2066}},
            {"W40X294", {12.0139, 1.9304, 12.0139, 1.9304, 36.8203, 1.0838}},
            {"W40X278", {11.9783, 1.8065, 11.9783, 1.8065, 36.7477, 1.0604}},
            {"W40X264", {11.9343, 1.7296, 11.9343, 1.7296, 36.8037, 0.9825}},
            {"W40X235", {11.9335, 1.5606, 11.9335, 1.5606, 36.9031, 0.863}},
            {"W40X211", {11.7004, 1.4557, 11.7004, 1.4557, 36.6398, 0.7653}},
            {"W40X183", {11.7613, 1.2172, 11.7613, 1.2172, 36.9012, 0.6689}},
            {"W40X167", {11.7844, 1.034, 11.7844, 1.034, 36.9078, 0.6758}},
            {"W40X149", {11.7293, 0.8482, 11.7293, 0.8482, 36.9112, 0.6486}},
            {"W36X925", {18.6235, 4.5143, 18.6235, 4.5143, 34.0914, 3.0461}},
            {"W36X853", {18.1883, 4.541, 18.1883, 4.541, 34.0362, 2.5227}},
            {"W36X802", {17.9816, 4.3036, 17.9816, 4.3036, 34.0009, 2.3911}},
            {"W36X723", {17.8242, 3.8881, 17.8242, 3.8881, 34.0811, 2.1849}},
            {"W36X652", {17.5833, 3.5411, 17.5833, 3.5411, 34.0273, 1.9814}},
            {"W36X529", {17.2283, 2.906, 17.2283, 2.906, 34.0372, 1.6438}},
            {"W36X487", {17.091, 2.6917, 17.091, 2.6917, 34.0491, 1.5005}},
            {"W36X441", {16.9065, 2.4619, 16.9065, 2.4619, 33.9317, 1.3772}},
            {"W36X395", {16.7352, 2.2329, 16.7352, 2.2329, 33.9966, 1.2167}},
            {"W36X361", {16.6286, 2.043, 16.6286, 2.043, 33.9727, 1.1219}},
            {"W36X330", {16.6018, 1.8575, 16.6018, 1.8575, 34.0591, 1.0356}},
            {"W36X302", {16.7174, 1.6663, 16.7174, 1.6663, 34.132, 0.9745}},
            {"W36X282", {16.6481, 1.5579, 16.6481, 1.5579, 34.1767, 0.9067}},
            {"W36X262", {16.5083, 1.4511, 16.5083, 1.4511, 34.0218, 0.8616}},
            {"W36X247", {16.3989, 1.3719, 16.3989, 1.3719, 33.9734, 0.8103}},
            {"W36X231", {16.4877, 1.2564, 16.4877, 1.2564, 34.1629, 0.784}},
            {"W36X387", {12.6862, 2.5664, 12.6862, 2.5664, 34.006, 1.4394}},
            {"W36X350", {12.5804, 2.3317, 12.5804, 2.3317, 33.9536, 1.3053}},
            {"W36X318", {12.3928, 2.135, 12.3928, 2.135, 34.004, 1.1904}},
            {"W36X286", {12.2836, 1.9384, 12.2836, 1.9384, 33.994, 1.065}},
            {"W36X256", {12.2183, 1.7282, 12.2183, 1.7282, 34.1366, 0.969}},
            {"W36X232", {12.1032, 1.577, 12.1032, 1.577, 34.0574, 0.8769}},
            {"W36X210", {12.1613, 1.3653, 12.1613, 1.3653, 34.0653, 0.8425}},
            {"W36X194", {12.1021, 1.2647, 12.1021, 1.2647, 34.0635, 0.7757}},
            {"W36X182", {12.0457, 1.1872, 12.0457, 1.1872, 34.062, 0.7346}},
            {"W36X170", {11.9953, 1.1091, 11.9953, 1.1091, 34.0777, 0.6868}},
            {"W36X160", {11.9854, 1.0251, 11.9854, 1.0251, 34.1168, 0.658}},
            {"W36X150", {11.9796, 0.9397, 11.9796, 0.9397, 34.1342, 0.6384}},
            {"W36X135", {11.9161, 0.7956, 11.9161, 0.7956, 34.0936, 0.614}},
            {"W33X387", {16.1871, 2.284, 16.1871, 2.284, 31.4021, 1.2757}},
            {"W33X354", {16.0879, 2.0976, 16.0879, 2.0976, 31.4374, 1.1621}},
            {"W33X318", {15.9995, 1.8847, 15.9995, 1.8847, 31.4404, 1.0637}},
            {"W33X291", {15.8628, 1.7393, 15.8628, 1.7393, 31.4564, 0.9703}},
            {"W33X263", {15.9092, 1.5468, 15.9092, 1.5468, 31.6308, 0.8915}},
            {"W33X241", {15.8413, 1.4058, 15.8413, 1.4058, 31.431, 0.8451}},
            {"W33X221", {15.8248, 1.2698, 15.8248, 1.2698, 31.5356, 0.7962}},
            {"W33X201", {15.7053, 1.1587, 15.7053, 1.1587, 31.5167, 0.7201}},
            {"W33X169", {11.4791, 1.2263, 11.4791, 1.2263, 31.4759, 0.6787}},
            {"W33X152", {11.574, 1.0536, 11.574, 1.0536, 31.5041, 0.6514}},
            {"W33X141", {11.5045, 0.9668, 11.5045, 0.9668, 31.4983, 0.6123}},
            {"W33X130", {11.4852, 0.8611, 11.4852, 0.8611, 31.5066, 0.5886}},
            {"W33X118", {11.4566, 0.7442, 11.4566, 0.7442, 31.5284, 0.5603}},
            {"W30X391", {15.607, 2.4363, 15.607, 2.4363, 28.3817, 1.3737}},
            {"W30X357", {15.5074, 2.2277, 15.5074, 2.2277, 28.4611, 1.2647}},
            {"W30X326", {15.2933, 2.0745, 15.2933, 2.0745, 28.27, 1.1462}},
            {"W30X292", {15.3074, 1.836, 15.3074, 1.836, 28.4253, 1.0463}},
            {"W30X261", {15.1388, 1.6549, 15.1388, 1.6549, 28.3625, 0.9484}},
            {"W30X235", {15.0693, 1.4963, 15.0693, 1.4963, 28.3919, 0.8535}},
            {"W30X211", {15.0571, 1.3286, 15.0571, 1.3286, 28.3352, 0.7864}},
            {"W30X191", {14.997, 1.1954, 14.997, 1.1954, 28.3813, 0.7143}},
            {"W30X173", {14.9421, 1.0742, 14.9421, 1.0742, 28.3901, 0.6628}},
            {"W30X148", {10.4625, 1.1856, 10.4625, 1.1856, 28.3687, 0.6625}},
            {"W30X132", {10.5322, 1.0036, 10.5322, 1.0036, 28.4332, 0.6214}},
            {"W30X124", {10.5163, 0.9311, 10.5163, 0.9311, 28.4286, 0.5953}},
            {"W30X116", {10.4683, 0.8554, 10.4683, 0.8554, 28.4086, 0.5737}},
            {"W30X108", {10.4646, 0.7622, 10.4646, 0.7622, 28.4709, 0.5539}},
            {"W30X99", {10.4486, 0.6713, 10.4486, 0.6713, 28.4808, 0.5266}},
            {"W30X90", {10.3972, 0.6124, 10.3972, 0.6124, 28.4875, 0.4773}},
            {"W27X539", {15.2707, 3.5263, 15.2707, 3.5263, 25.5109, 2.0104}},
            {"W27X368", {14.7461, 2.4406, 14.7461, 2.4406, 25.5749, 1.4434}},
            {"W27X336", {14.6172, 2.258, 14.6172, 2.258, 25.5548, 1.2987}},
            {"W27X307", {14.334, 2.131, 14.334, 2.131, 25.359, 1.153}},
            {"W27X281", {14.3617, 1.9244, 14.3617, 1.9244, 25.4854, 1.0929}},
            {"W27X258", {14.2297, 1.7846, 14.2297, 1.7846, 25.4498, 0.9939}},
            {"W27X235", {14.1542, 1.6235, 14.1542, 1.6235, 25.4572, 0.9213}},
            {"W27X217", {14.0973, 1.5047, 14.0973, 1.5047, 25.5111, 0.8427}},
            {"W27X194", {14.0062, 1.3496, 14.0062, 1.3496, 25.496, 0.7574}},
            {"W27X178", {14.0178, 1.207, 14.0178, 1.207, 25.4466, 0.7336}},
            {"W27X161", {14.0303, 1.078, 14.0303, 1.078, 25.5617, 0.68}},
            {"W27X146", {13.9335, 0.9813, 13.9335, 0.9813, 25.5153, 0.6219}},
            {"W27X129", {9.9677, 1.1115, 9.9677, 1.1115, 25.4728, 0.615}},
            {"W27X114", {10.0915, 0.9256, 10.0915, 0.9256, 25.5349, 0.5849}},
            {"W27X102", {9.9853, 0.8358, 9.9853, 0.8358, 25.5154, 0.5218}},
            {"W27X94", {9.9645, 0.7502, 9.9645, 0.7502, 25.5416, 0.496}},
            {"W27X84", {9.9764, 0.6392, 9.9764, 0.6392, 25.6106, 0.4666}},
            {"W24X370", {13.6658, 2.7125, 13.6658, 2.7125, 22.5779, 1.54}},
            {"W24X335", {13.5146, 2.4911, 13.5146, 2.4911, 22.557, 1.375}},
            {"W24X306", {13.3913, 2.2867, 13.3913, 2.2867, 22.5762, 1.2605}},
            {"W24X279", {13.2616, 2.1098, 13.2616, 2.1098, 22.5153, 1.1521}},
            {"W24X250", {13.1285, 1.9144, 13.1285, 1.9144, 22.5033, 1.032}},
            {"W24X229", {13.0961, 1.7346, 13.0961, 1.7346, 22.5773, 0.964}},
            {"W24X207", {13.0481, 1.5574, 13.0481, 1.5574, 22.6528, 0.8859}},
            {"W24X192", {13.0051, 1.4428, 13.0051, 1.4428, 22.6414, 0.8376}},
            {"W24X176", {12.8237, 1.3606, 12.8237, 1.3606, 22.4993, 0.7471}},
            {"W24X162", {12.9988, 1.208, 12.9988, 1.208, 22.6461, 0.7247}},
            {"W24X146", {12.8993, 1.0913, 12.8993, 1.0913, 22.5929, 0.6581}},
            {"W24X131", {12.8294, 0.9647, 12.8294, 0.9647, 22.5815, 0.6137}},
            {"W24X117", {12.7799, 0.8528, 12.7799, 0.8528, 22.6158, 0.5572}},
            {"W24X104", {12.7385, 0.751, 12.7385, 0.751, 22.5983, 0.5123}},
            {"W24X103", {8.9553, 0.9915, 8.9553, 0.9915, 22.5756, 0.5555}},
            {"W24X94", {9.069, 0.8745, 9.069, 0.8745, 22.6487, 0.5233}},
            {"W24X84", {9.0256, 0.7687, 9.0256, 0.7687, 22.6707, 0.4775}},
            {"W24X76", {8.9988, 0.6779, 8.9988, 0.6779, 22.6664, 0.4497}},
            {"W24X68", {8.9714, 0.5836, 8.9714, 0.5836, 22.69, 0.4251}},
            {"W24X62", {7.0454, 0.5892, 7.0454, 0.5892, 22.6963, 0.4364}},
            {"W24X55", {7.0221, 0.5021, 7.0221, 0.5021, 22.7656, 0.4017}},
            {"W21X275", {12.8208, 2.2312, 12.8208, 2.2312, 19.6925, 1.2497}},
            {"W21X248", {12.7655, 2.0086, 12.7655, 2.0086, 19.7742, 1.1404}},
            {"W21X223", {12.6884, 1.797, 12.6884, 1.797, 19.9306, 1.05}},
            {"W21X201", {12.5975, 1.6225, 12.5975, 1.6225, 19.8211, 0.93}},
            {"W21X182", {12.5159, 1.4749, 12.5159, 1.4749, 19.8059, 0.8431}},
            {"W21X166", {12.3842, 1.372, 12.3842, 1.372, 19.7555, 0.7498}},
            {"W21X147", {12.5137, 1.1492, 12.5137, 1.1492, 19.8006, 0.7296}},
            {"W21X132", {12.4428, 1.0356, 12.4428, 1.0356, 19.7865, 0.6589}},
            {"W21X122", {12.3894, 0.9611, 12.3894, 0.9611, 19.7915, 0.6109}},
            {"W21X111", {12.299, 0.8825, 12.299, 0.8825, 19.7771, 0.5523}},
            {"W21X101", {12.3043, 0.798, 12.3043, 0.798, 19.8195, 0.5131}},
            {"W21X93", {8.4083, 0.934, 8.4083, 0.934, 19.8179, 0.586}},
            {"W21X83", {8.3592, 0.8336, 8.3592, 0.8336, 19.8254, 0.5281}},
            {"W21X73", {8.2732, 0.7463, 8.2732, 0.7463, 19.7819, 0.4626}},
            {"W21X68", {8.2488, 0.69, 8.2488, 0.69, 19.8312, 0.4354}},
            {"W21X62", {8.2574, 0.6116, 8.2574, 0.6116, 19.8714, 0.4123}},
            {"W21X55", {8.1882, 0.5277, 8.1882, 0.5277, 19.815, 0.3827}},
            {"W21X48", {8.1132, 0.4338, 8.1132, 0.4338, 19.8899, 0.3557}},
            {"W21X57", {6.5429, 0.6527, 6.5429, 0.6527, 19.9, 0.4112}},
            {"W21X50", {6.4889, 0.5447, 6.4889, 0.5447, 19.8345, 0.385}},
            {"W21X44", {6.4647, 0.458, 6.4647, 0.458, 19.8529, 0.3568}},
            {"W18X311", {12.0253, 2.7247, 12.0253, 2.7247, 16.903, 1.5433}},
            {"W18X283", {11.8982, 2.4934, 11.8982, 2.4934, 16.884, 1.4194}},
            {"W18X258", {11.8051, 2.2784, 11.8051, 2.2784, 16.9238, 1.3131}},
            {"W18X234", {11.6151, 2.1284, 11.6151, 2.1284, 16.8067, 1.1406}},
            {"W18X211", {11.612, 1.8821, 11.612, 1.8821, 16.9417, 1.097}},
            {"W18X192", {11.4173, 1.769, 11.4173, 1.769, 16.815, 0.9403}},
            {"W18X175", {11.4146, 1.5729, 11.4146, 1.5729, 16.9296, 0.9155}},
            {"W18X158", {11.2891, 1.444, 11.2891, 1.444, 16.8568, 0.8124}},
            {"W18X143", {11.198, 1.3265, 11.198, 1.3265, 16.8548, 0.7296}},
            {"W18X130", {11.1456, 1.203, 11.1456, 1.203, 16.8465, 0.6809}},
            {"W18X119", {11.2644, 1.0601, 11.2644, 1.0601, 16.8831, 0.6652}},
            {"W18X106", {11.1663, 0.9468, 11.1663, 0.9468, 16.8433, 0.591}},
            {"W18X97", {11.136, 0.8723, 11.136, 0.8723, 16.895, 0.5372}},
            {"W18X86", {11.0577, 0.7758, 11.0577, 0.7758, 16.8758, 0.483}},
            {"W18X76", {10.9919, 0.686, 10.9919, 0.686, 16.8502, 0.4292}},
            {"W18X71", {7.6251, 0.8133, 7.6251, 0.8133, 16.8659, 0.505}},
            {"W18X65", {7.5852, 0.7516, 7.5852, 0.7516, 16.8837, 0.4559}},
            {"W18X60", {7.5449, 0.6982, 7.5449, 0.6982, 16.906, 0.4188}},
            {"W18X55", {7.5325, 0.6289, 7.5325, 0.6289, 16.9264, 0.3983}},
            {"W18X50", {7.4808, 0.5736, 7.4808, 0.5736, 16.8945, 0.3626}},
            {"W18X46", {6.0344, 0.6125, 6.0344, 0.6125, 16.9004, 0.3619}},
            {"W18X40", {5.9743, 0.5361, 5.9743, 0.5361, 16.8512, 0.32}},
            {"W18X35", {5.9703, 0.4302, 5.9703, 0.4302, 16.9252, 0.3053}},
            {"W16X100", {10.3953, 0.9922, 10.3953, 0.9922, 15.0242, 0.5834}},
            {"W16X89", {10.3799, 0.8732, 10.3799, 0.8732, 15.0403, 0.5379}},
            {"W16X77", {10.2599, 0.766, 10.2599, 0.766, 15.0255, 0.4582}},
            {"W16X67", {10.2078, 0.6706, 10.2078, 0.6706, 15.0245, 0.3943}},
            {"W16X57", {7.097, 0.7218, 7.097, 0.7218, 15.0114, 0.436}},
            {"W16X50", {7.0704, 0.6301, 7.0704, 0.6301, 15.0757, 0.3847}},
            {"W16X45", {6.9987, 0.5731, 6.9987, 0.5731, 15.006, 0.3517}},
            {"W16X40", {7.022, 0.5, 7.022, 0.5, 15.1032, 0.3169}},
            {"W16X36", {7.0238, 0.4235, 7.0238, 0.4235, 15.1251, 0.308}},
            {"W16X31", {5.5125, 0.4431, 5.5125, 0.4431, 15.0862, 0.2816}},
            {"W16X26", {5.4825, 0.3483, 5.4825, 0.3483, 15.0989, 0.2561}},
            {"W14X873", {18.903, 5.4136, 18.903, 5.4136, 12.7625, 4.1018}},
            {"W14X808", {18.6488, 5.0793, 18.6488, 5.0793, 12.6258, 3.8455}},
            {"W14X730", {17.9145, 4.891, 17.9145, 4.891, 12.6153, 3.1526}},
            {"W14X665", {17.6595, 4.5137, 17.6595, 4.5137, 12.5861, 2.9088}},
            {"W14X605", {17.4057, 4.164, 17.4057, 4.164, 12.5826, 2.6298}},
            {"W14X550", {17.1547, 3.8434, 17.1547, 3.8434, 12.5708, 2.4027}},
            {"W14X500", {16.9512, 3.5331, 16.9512, 3.5331, 12.5396, 2.1747}},
            {"W14X455", {16.8322, 3.2093, 16.8322, 3.2093, 12.6036, 2.0599}},
            {"W14X426", {16.6792, 3.0423, 16.6792, 3.0423, 12.6063, 1.8671}},
            {"W14X398", {16.5574, 2.8607, 16.5574, 2.8607, 12.5737, 1.7708}},
            {"W14X370", {16.503, 2.6496, 16.503, 2.6496, 12.6273, 1.7064}},
            {"W14X342", {16.4055, 2.4525, 16.4055, 2.4525, 12.6472, 1.617}},
            {"W14X311", {16.194, 2.2703, 16.194, 2.2703, 12.5959, 1.4192}},
            {"W14X283", {16.0574, 2.0837, 16.0574, 2.0837, 12.5835, 1.3012}},
            {"W14X257", {15.9977, 1.8874, 15.9977, 1.8874, 12.6297, 1.2065}},
            {"W14X233", {15.8587, 1.728, 15.8587, 1.728, 12.606, 1.0864}},
            {"W14X211", {15.8455, 1.5516, 15.8455, 1.5516, 12.6512, 1.0149}},
            {"W14X193", {15.7328, 1.4329, 15.7328, 1.4329, 12.6542, 0.9267}},
            {"W14X176", {15.6263, 1.3167, 15.6263, 1.3167, 12.6027, 0.8452}},
            {"W14X159", {15.5433, 1.1941, 15.5433, 1.1941, 12.615, 0.761}},
            {"W14X145", {15.4365, 1.1038, 15.4365, 1.1038, 12.5697, 0.6857}},
            {"W14X132", {14.7231, 1.0296, 14.7231, 1.0296, 12.6329, 0.6719}},
            {"W14X120", {14.7086, 0.9325, 14.7086, 0.9325, 12.6925, 0.6195}},
            {"W14X109", {14.5943, 0.8622, 14.5943, 0.8622, 12.6469, 0.5417}},
            {"W14X99", {14.552, 0.7823, 14.552, 0.7823, 12.6371, 0.5018}},
            {"W14X90", {14.4703, 0.7164, 14.4703, 0.7164, 12.6251, 0.4588}},
            {"W14X82", {10.076, 0.8667, 10.076, 0.8667, 12.6214, 0.5199}},
            {"W14X74", {10.0663, 0.7869, 10.0663, 0.7869, 12.6598, 0.4725}},
            {"W14X68", {9.9724, 0.7311, 9.9724, 0.7311, 12.6259, 0.4311}},
            {"W14X61", {9.9205, 0.6571, 9.9205, 0.6571, 12.6255, 0.3853}},
            {"W14X53", {8.0239, 0.6692, 8.0239, 0.6692, 12.6486, 0.3851}},
            {"W14X48", {8.013, 0.5986, 8.013, 0.5986, 12.6952, 0.3564}},
            {"W14X43", {7.9646, 0.536, 7.9646, 0.536, 12.6977, 0.3202}},
            {"W14X38", {6.7992, 0.5088, 6.7992, 0.5088, 13.1636, 0.3255}},
            {"W14X34", {6.7656, 0.4508, 6.7656, 0.4508, 13.1877, 0.2964}},
            {"W14X30", {6.7189, 0.3872, 6.7189, 0.3872, 13.1535, 0.2776}},
            {"W14X26", {5.0171, 0.4223, 5.0171, 0.4223, 13.1549, 0.2628}},
            {"W14X22", {4.9779, 0.3397, 4.9779, 0.3397, 13.168, 0.2366}},
            {"W12X336", {13.4065, 2.9494, 13.4065, 2.9494, 10.9258, 1.8145}},
            {"W12X305", {13.2665, 2.6873, 13.2665, 2.6873, 10.9705, 1.6593}},
            {"W12X279", {13.1184, 2.4815, 13.1184, 2.4815, 10.8908, 1.542}},
            {"W12X252", {13.0011, 2.2538, 13.0011, 2.2538, 10.9156, 1.4192}},
            {"W12X230", {12.8847, 2.076, 12.8847, 2.076, 10.9194, 1.2993}},
            {"W12X210", {12.8156, 1.8873, 12.8156, 1.8873, 10.9508, 1.2287}},
            {"W12X190", {12.6144, 1.7575, 12.6144, 1.7575, 10.8553, 1.073}},
            {"W12X170", {12.547, 1.5678, 12.547, 1.5678, 10.9287, 0.9753}},
            {"W12X152", {12.4797, 1.3987, 12.4797, 1.3987, 10.9392, 0.8985}},
            {"W12X136", {12.3716, 1.2589, 12.3716, 1.2589, 10.8983, 0.8064}},
            {"W12X120", {12.3103, 1.1084, 12.3103, 1.1084, 10.9301, 0.7244}},
            {"W12X106", {12.1812, 0.9981, 12.1812, 0.9981, 10.9182, 0.6319}},
            {"W12X96", {12.1546, 0.9016, 12.1546, 0.9016, 10.9438, 0.5741}},
            {"W12X87", {12.116, 0.8124, 12.116, 0.8124, 10.9438, 0.5408}},
            {"W12X79", {12.0605, 0.7381, 12.0605, 0.7381, 10.9476, 0.4935}},
            {"W12X72", {12.004, 0.6759, 12.004, 0.6759, 10.9501, 0.4468}},
            {"W12X65", {11.9474, 0.6118, 11.9474, 0.6118, 10.9439, 0.4095}},
            {"W12X58", {9.9766, 0.6459, 9.9766, 0.6459, 10.9602, 0.3769}},
            {"W12X53", {9.9408, 0.5838, 9.9408, 0.5838, 10.9441, 0.3622}},
            {"W12X50", {8.0702, 0.6419, 8.0702, 0.6419, 10.9587, 0.3879}},
            {"W12X45", {8.0231, 0.5804, 8.0231, 0.5804, 10.9588, 0.3461}},
            {"W12X40", {7.9779, 0.5204, 7.9779, 0.5204, 10.951, 0.3097}},
            {"W12X35", {6.5237, 0.5288, 6.5237, 0.5288, 11.4578, 0.2973}},
            {"W12X30", {6.4977, 0.4435, 6.4977, 0.4435, 11.4777, 0.264}},
            {"W12X26", {6.47, 0.3829, 6.47, 0.3829, 11.4836, 0.235}},
            {"W12X22", {4.0253, 0.427, 4.0253, 0.427, 11.5081, 0.2645}},
            {"W12X19", {3.9908, 0.3537, 3.9908, 0.3537, 11.5252, 0.2384}},
            {"W12X16", {3.9817, 0.267, 3.9817, 0.267, 11.5574, 0.224}},
            {"W12X14", {3.9662, 0.2262, 3.9662, 0.2262, 11.5558, 0.2047}},
            {"W10X112", {10.4172, 1.2509, 10.4172, 1.2509, 8.8717, 0.7703}},
            {"W10X100", {10.3334, 1.1239, 10.3334, 1.1239, 8.8808, 0.6859}},
            {"W10X88", {10.2652, 0.9915, 10.2652, 0.9915, 8.8757, 0.6352}},
            {"W10X77", {10.1948, 0.8709, 10.1948, 0.8709, 8.8757, 0.5565}},
            {"W10X68", {10.1377, 0.7709, 10.1377, 0.7709, 8.8992, 0.4813}},
            {"W10X60", {10.0459, 0.6858, 10.0459, 0.6858, 8.8754, 0.441}},
            {"W10X54", {9.9663, 0.624, 9.9663, 0.624, 8.8759, 0.3793}},
            {"W10X49", {9.981, 0.5631, 9.981, 0.5631, 8.8923, 0.357}},
            {"W10X45", {7.9885, 0.6277, 7.9885, 0.6277, 8.8742, 0.3675}},
            {"W10X39", {7.9404, 0.5386, 7.9404, 0.5386, 8.8768, 0.3308}},
            {"W10X33", {7.9252, 0.4404, 7.9252, 0.4404, 8.924, 0.3045}},
            {"W10X30", {5.8025, 0.5121, 5.8025, 0.5121, 9.4773, 0.306}},
            {"W10X26", {5.7573, 0.4427, 5.7573, 0.4427, 9.4711, 0.2659}},
            {"W10X22", {5.734, 0.3624, 5.734, 0.3624, 9.4744, 0.2466}},
            {"W10X19", {4.0146, 0.3965, 4.0146, 0.3965, 9.5029, 0.2568}},
            {"W10X17", {4.003, 0.3318, 4.003, 0.3318, 9.5134, 0.2459}},
            {"W10X15", {3.9808, 0.2739, 3.9808, 0.2739, 9.5064, 0.2349}},
            {"W10X12", {3.9571, 0.2105, 3.9571, 0.2105, 9.5395, 0.1965}},
            {"W8X67", {8.2742, 0.9372, 8.2742, 0.9372, 7.1409, 0.5866}},
            {"W8X58", {8.2114, 0.8129, 8.2114, 0.8129, 7.1438, 0.5247}},
            {"W8X48", {8.075, 0.6935, 8.075, 0.6935, 7.1286, 0.4072}},
            {"W8X40", {8.0514, 0.5638, 8.0514, 0.5638, 7.1422, 0.3685}},
            {"W8X35", {8.01, 0.4967, 8.01, 0.4967, 7.1686, 0.3255}},
            {"W8X31", {7.9779, 0.4381, 7.9779, 0.4381, 7.156, 0.2992}},
            {"W8X28", {6.5251, 0.4678, 6.5251, 0.4678, 7.1623, 0.2992}},
            {"W8X24", {6.4743, 0.404, 6.4743, 0.404, 7.1592, 0.2572}},
            {"W8X21", {5.2563, 0.4031, 5.2563, 0.4031, 7.5052, 0.2565}},
            {"W8X18", {5.2435, 0.3313, 5.2435, 0.3313, 7.5153, 0.2377}},
            {"W8X15", {3.9927, 0.3204, 3.9927, 0.3204, 7.508, 0.2513}},
            {"W8X13", {3.9925, 0.2565, 3.9925, 0.2565, 7.5409, 0.2381}},
            {"W8X10", {3.9032, 0.2105, 3.9032, 0.2105, 7.5291, 0.1754}},
            {"W6X25", {6.0945, 0.4528, 6.0945, 0.4528, 5.4902, 0.3317}},
            {"W6X20", {6.0205, 0.3654, 6.0205, 0.3654, 5.4867, 0.2681}},
            {"W6X15", {5.9751, 0.2619, 5.9751, 0.2619, 5.4851, 0.2376}},
            {"W6X16", {4.0236, 0.4071, 4.0236, 0.4071, 5.4844, 0.2679}},
            {"W6X12", {3.9899, 0.2819, 3.9899, 0.2819, 5.5018, 0.2366}},
            {"W6X9", {3.9288, 0.2174, 3.9288, 0.2174, 5.499, 0.177}},
            {"W6X8.5", {3.9302, 0.1964, 3.9302, 0.1964, 5.484, 0.1782}},
            {"W5X19", {5.0224, 0.4319, 5.0224, 0.4319, 4.3068, 0.2835}},
            {"W5X16", {4.9819, 0.364, 4.9819, 0.364, 4.2984, 0.252}},
            {"W4X13", {4.0539, 0.3467, 4.0539, 0.3467, 3.4747, 0.294}},
            {"M12.5X12.4", {3.7486, 0.2284, 3.7486, 0.2284, 12.1624, 0.158}},
            {"M12.5X11.6", {3.4902, 0.2125, 3.4902, 0.2125, 12.1531, 0.1577}},
            {"M12X11.8", {3.0811, 0.2224, 3.0811, 0.2224, 11.6656, 0.1802}},
            {"M12X10.8", {3.0327, 0.2164, 3.0327, 0.2164, 11.5984, 0.1612}},
            {"M12X10", {3.2311, 0.1826, 3.2311, 0.1826, 11.6887, 0.1516}},
            {"M10X9", {2.6804, 0.2083, 2.6804, 0.2083, 9.6481, 0.1591}},
            {"M10X8", {2.6774, 0.1846, 2.6774, 0.1846, 9.6475, 0.1433}},
            {"M10X7.5", {2.6741, 0.1757, 2.6741, 0.1757, 9.7099, 0.132}},
            {"M8X6.5", {2.2744, 0.1907, 2.2744, 0.1907, 7.6584, 0.1376}},
            {"M8X6.2", {2.2724, 0.1791, 2.2724, 0.1791, 7.6997, 0.1307}},
            {"M6X4.4", {1.8296, 0.1751, 1.8296, 0.1751, 5.6773, 0.1144}},
            {"M6X3.7", {1.9922, 0.1306, 1.9922, 0.1306, 5.6975, 0.1}},
            {"M5X18.9", {4.9899, 0.4194, 4.9899, 0.4194, 4.1827, 0.3297}},
            {"M4X6", {3.7425, 0.1675, 3.7425, 0.1675, 3.5063, 0.1391}},
            {"M4X4.08", {2.2093, 0.1819, 2.2093, 0.1819, 3.7016, 0.123}},
            {"M3X2.9", {2.2373, 0.1391, 2.2373, 0.1391, 2.7965, 0.0963}},
            {"S24X121", {7.5429, 1.1462, 7.5429, 1.1462, 22.4048, 0.8133}},
            {"S24X106", {7.3437, 1.1562, 7.3437, 1.1562, 22.3715, 0.6314}},
            {"S24X100", {6.7425, 0.9126, 6.7425, 0.9126, 22.3861, 0.7604}},
            {"S24X90", {6.6309, 0.9113, 6.6309, 0.9113, 22.421, 0.6438}},
            {"S24X80", {6.51, 0.9087, 6.51, 0.9087, 22.4288, 0.5213}},
            {"S20X96", {6.7562, 0.9533, 6.7562, 0.9533, 18.5912, 0.8258}},
            {"S20X86", {6.5912, 0.9659, 6.5912, 0.9659, 18.5366, 0.6791}},
            {"S20X75", {5.9427, 0.8323, 5.9427, 0.8323, 18.5672, 0.6522}},
            {"S20X66", {5.8218, 0.8313, 5.8218, 0.8313, 18.5603, 0.5244}},
            {"S18X70", {5.8235, 0.7127, 5.8235, 0.7127, 16.795, 0.7266}},
            {"S18X54.7", {5.58, 0.7121, 5.58, 0.7121, 16.8247, 0.4806}},
            {"S15X50", {5.2452, 0.6413, 5.2452, 0.6413, 13.9134, 0.5752}},
            {"S15X42.9", {5.115, 0.6401, 5.115, 0.6401, 13.9405, 0.4365}},
            {"S12X50", {5.1377, 0.6749, 5.1377, 0.6749, 10.7871, 0.7209}},
            {"S12X40.8", {4.8825, 0.6912, 4.8825, 0.6912, 10.7428, 0.4813}},
            {"S12X35", {4.7244, 0.5573, 4.7244, 0.5573, 11.0648, 0.4488}},
            {"S12X31.8", {4.65, 0.5571, 4.65, 0.5571, 11.0699, 0.3745}},
            {"S10X35", {4.5942, 0.5024, 4.5942, 0.5024, 9.1453, 0.623}},
            {"S10X25.4", {4.3338, 0.5009, 4.3338, 0.5009, 9.1954, 0.3328}},
            {"S8X23", {3.8781, 0.4338, 3.8781, 0.4338, 7.2594, 0.4697}},
            {"S8X18.4", {3.72, 0.4348, 3.72, 0.4348, 7.3007, 0.29}},
            {"S6X17.25", {3.3201, 0.3676, 3.3201, 0.3676, 5.3509, 0.4891}},
            {"S6X12.5", {3.0969, 0.3685, 3.0969, 0.3685, 5.3869, 0.2482}},
            {"S5X10", {2.79, 0.3347, 2.79, 0.3347, 4.4518, 0.229}},
            {"S4X9.5", {2.604, 0.3005, 2.604, 0.3005, 3.4733, 0.3488}},
            {"S4X7.7", {2.4738, 0.3035, 2.4738, 0.3035, 3.4778, 0.2065}},
            {"S3X7.5", {2.3343, 0.2694, 2.3343, 0.2694, 2.5077, 0.3734}},
            {"S3X5.7", {2.1669, 0.2727, 2.1669, 0.2727, 2.5009, 0.1819}},
            {"HP18X204", {18.0658, 1.1369, 18.0658, 1.1369, 16.1371, 1.1865}},
            {"HP18X181", {17.9227, 1.013, 17.9227, 1.013, 16.0881, 1.05}},
            {"HP18X157", {17.7846, 0.886, 17.7846, 0.886, 16.048, 0.9135}},
            {"HP18X135", {17.6241, 0.7719, 17.6241, 0.7719, 16.0742, 0.7875}},
            {"HP16X183", {16.1867, 1.1541, 16.1867, 1.1541, 14.2764, 1.1753}},
            {"HP16X162", {15.9948, 1.0189, 15.9948, 1.0189, 14.3747, 1.05}},
            {"HP16X141", {15.8485, 0.9009, 15.8485, 0.9009, 14.2849, 0.9188}},
            {"HP16X121", {15.7171, 0.7761, 15.7171, 0.7761, 14.3627, 0.7875}},
            {"HP16X101", {15.5437, 0.6553, 15.5437, 0.6553, 14.3224, 0.6562}},
            {"HP16X88", {15.4376, 0.567, 15.4376, 0.567, 14.3582, 0.567}},
            {"HP14X117", {14.8657, 0.808, 14.8657, 0.808, 12.6109, 0.8234}},
            {"HP14X102", {14.7443, 0.7102, 14.7443, 0.7102, 12.6103, 0.7286}},
            {"HP14X89", {14.6915, 0.6161, 14.6915, 0.6161, 12.6598, 0.6339}},
            {"HP14X73", {14.5722, 0.5059, 14.5722, 0.5059, 12.6651, 0.5249}},
            {"HP12X89", {12.2973, 0.7212, 12.2973, 0.7212, 11.0129, 0.7435}},
            {"HP12X84", {12.2899, 0.6871, 12.2899, 0.6871, 10.9623, 0.706}},
            {"HP12X74", {12.2529, 0.606, 12.2529, 0.606, 10.9816, 0.6322}},
            {"HP12X63", {12.0832, 0.5197, 12.0832, 0.5197, 10.9579, 0.5341}},
            {"HP12X53", {12.0226, 0.438, 12.0226, 0.438, 10.9801, 0.4537}},
            {"HP10X57", {10.2237, 0.5661, 10.2237, 0.5661, 8.914, 0.5766}},
            {"HP10X42", {10.0388, 0.4247, 10.0388, 0.4247, 8.8801, 0.4358}},
            {"HP8X36", {8.1453, 0.4467, 8.1453, 0.4467, 7.1434, 0.4658}},
        };

        return registry.fetch(type);
    }

    vec ustsection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"WT22X167.5", {15.900, 1.770, 20.230, 1.030}},
            {"WT22X145", {15.800, 1.580, 20.220, 0.865}},
            {"WT22X131", {15.800, 1.420, 20.280, 0.785}},
            {"WT22X115", {15.800, 1.220, 20.280, 0.710}},
            {"WT20X327.5", {16.900, 3.540, 18.260, 1.970}},
            {"WT20X296.5", {16.700, 3.230, 18.270, 1.790}},
            {"WT20X251.5", {16.400, 2.760, 18.240, 1.540}},
            {"WT20X215.5", {16.200, 2.360, 18.240, 1.340}},
            {"WT20X198.5", {16.100, 2.200, 18.300, 1.220}},
            {"WT20X186", {16.100, 2.050, 18.250, 1.160}},
            {"WT20X181", {16.000, 2.010, 18.290, 1.120}},
            {"WT20X162", {15.900, 1.810, 18.290, 1.000}},
            {"WT20X148.5", {15.800, 1.650, 18.250, 0.930}},
            {"WT20X138.5", {15.800, 1.580, 18.220, 0.830}},
            {"WT20X124.5", {15.800, 1.420, 18.280, 0.750}},
            {"WT20X107.5", {15.800, 1.220, 18.280, 0.650}},
            {"WT20X99.5", {15.800, 1.070, 18.230, 0.650}},
            {"WT20X196", {12.400, 2.520, 18.280, 1.420}},
            {"WT20X165.5", {12.200, 2.130, 18.270, 1.220}},
            {"WT20X163.5", {12.100, 2.130, 18.270, 1.180}},
            {"WT20X147", {12.000, 1.930, 18.270, 1.060}},
            {"WT20X139", {12.000, 1.810, 18.290, 1.030}},
            {"WT20X132", {11.900, 1.730, 18.270, 0.960}},
            {"WT20X117.5", {11.900, 1.580, 18.220, 0.830}},
            {"WT20X105.5", {11.800, 1.420, 18.280, 0.750}},
            {"WT20X91.5", {11.800, 1.200, 18.300, 0.650}},
            {"WT20X83.5", {11.800, 1.030, 18.270, 0.650}},
            {"WT20X74.5", {11.800, 0.830, 18.270, 0.630}},
            {"WT18X462.5", {18.600, 4.530, 17.070, 3.020}},
            {"WT18X426.5", {18.200, 4.530, 17.070, 2.520}},
            {"WT18X401", {18.000, 4.290, 17.010, 2.380}},
            {"WT18X361.5", {17.800, 3.900, 17.000, 2.170}},
            {"WT18X326", {17.600, 3.540, 16.960, 1.970}},
            {"WT18X264.5", {17.200, 2.910, 16.990, 1.610}},
            {"WT18X243.5", {17.100, 2.680, 17.020, 1.500}},
            {"WT18X220.5", {17.000, 2.440, 16.960, 1.360}},
            {"WT18X197.5", {16.800, 2.200, 17.000, 1.220}},
            {"WT18X180.5", {16.700, 2.010, 16.990, 1.120}},
            {"WT18X165", {16.600, 1.850, 16.950, 1.020}},
            {"WT18X151", {16.700, 1.680, 17.020, 0.945}},
            {"WT18X141", {16.600, 1.570, 17.030, 0.885}},
            {"WT18X131", {16.600, 1.440, 16.960, 0.840}},
            {"WT18X123.5", {16.500, 1.350, 16.950, 0.800}},
            {"WT18X115.5", {16.500, 1.260, 16.940, 0.760}},
            {"WT18X128", {12.200, 1.730, 16.970, 0.960}},
            {"WT18X116", {12.100, 1.570, 17.030, 0.870}},
            {"WT18X105", {12.200, 1.360, 16.940, 0.830}},
            {"WT18X97", {12.100, 1.260, 16.940, 0.765}},
            {"WT18X91", {12.100, 1.180, 17.020, 0.725}},
            {"WT18X85", {12.000, 1.100, 17.000, 0.680}},
            {"WT18X80", {12.000, 1.020, 16.980, 0.650}},
            {"WT18X75", {12.000, 0.940, 16.960, 0.625}},
            {"WT18X67.5", {12.000, 0.790, 17.010, 0.600}},
            {"WT16.5X193.5", {16.200, 2.280, 15.720, 1.260}},
            {"WT16.5X177", {16.100, 2.090, 15.710, 1.160}},
            {"WT16.5X159", {16.000, 1.890, 15.710, 1.040}},
            {"WT16.5X145.5", {15.900, 1.730, 15.670, 0.960}},
            {"WT16.5X131.5", {15.800, 1.570, 15.730, 0.870}},
            {"WT16.5X120.5", {15.900, 1.400, 15.700, 0.830}},
            {"WT16.5X110.5", {15.800, 1.280, 15.720, 0.775}},
            {"WT16.5X100.5", {15.700, 1.150, 15.650, 0.715}},
            {"WT16.5X84.5", {11.500, 1.220, 15.680, 0.670}},
            {"WT16.5X76", {11.600, 1.060, 15.640, 0.635}},
            {"WT16.5X70.5", {11.500, 0.960, 15.740, 0.605}},
            {"WT16.5X65", {11.500, 0.855, 15.645, 0.580}},
            {"WT16.5X59", {11.500, 0.740, 15.660, 0.550}},
            {"WT15X195.5", {15.600, 2.440, 14.160, 1.360}},
            {"WT15X178.5", {15.500, 2.240, 14.160, 1.240}},
            {"WT15X163", {15.400, 2.050, 14.150, 1.140}},
            {"WT15X146", {15.300, 1.850, 14.150, 1.020}},
            {"WT15X130.5", {15.200, 1.650, 14.150, 0.930}},
            {"WT15X117.5", {15.100, 1.500, 14.200, 0.830}},
            {"WT15X105.5", {15.100, 1.320, 14.180, 0.775}},
            {"WT15X95.5", {15.000, 1.190, 14.110, 0.710}},
            {"WT15X86.5", {15.000, 1.070, 14.130, 0.655}},
            {"WT15X74", {10.500, 1.180, 14.120, 0.650}},
            {"WT15X66", {10.500, 1.000, 14.200, 0.615}},
            {"WT15X62", {10.500, 0.930, 14.170, 0.585}},
            {"WT15X58", {10.500, 0.850, 14.150, 0.565}},
            {"WT15X54", {10.500, 0.760, 14.140, 0.545}},
            {"WT15X49.5", {10.500, 0.670, 14.130, 0.520}},
            {"WT15X45", {10.400, 0.610, 14.190, 0.470}},
            {"WT13.5X269.5", {15.300, 3.540, 12.760, 1.970}},
            {"WT13.5X184", {14.700, 2.480, 12.720, 1.380}},
            {"WT13.5X168", {14.600, 2.280, 12.720, 1.260}},
            {"WT13.5X153.5", {14.400, 2.090, 12.710, 1.160}},
            {"WT13.5X140.5", {14.400, 1.930, 12.670, 1.060}},
            {"WT13.5X129", {14.300, 1.770, 12.730, 0.980}},
            {"WT13.5X117.5", {14.200, 1.610, 12.690, 0.910}},
            {"WT13.5X108.5", {14.100, 1.500, 12.700, 0.830}},
            {"WT13.5X97", {14.000, 1.340, 12.760, 0.750}},
            {"WT13.5X89", {14.100, 1.190, 12.710, 0.725}},
            {"WT13.5X80.5", {14.000, 1.080, 12.720, 0.660}},
            {"WT13.5X73", {14.000, 0.975, 12.725, 0.605}},
            {"WT13.5X64.5", {10.000, 1.100, 12.700, 0.610}},
            {"WT13.5X57", {10.100, 0.930, 12.670, 0.570}},
            {"WT13.5X51", {10.000, 0.830, 12.670, 0.515}},
            {"WT13.5X47", {10.000, 0.745, 12.755, 0.490}},
            {"WT13.5X42", {10.000, 0.640, 12.760, 0.460}},
            {"WT12X185", {13.700, 2.720, 11.280, 1.520}},
            {"WT12X167.5", {13.500, 2.480, 11.320, 1.380}},
            {"WT12X153", {13.400, 2.280, 11.320, 1.260}},
            {"WT12X139.5", {13.300, 2.090, 11.310, 1.160}},
            {"WT12X125", {13.200, 1.890, 11.310, 1.040}},
            {"WT12X114.5", {13.100, 1.730, 11.270, 0.960}},
            {"WT12X103.5", {13.000, 1.570, 11.330, 0.870}},
            {"WT12X96", {13.000, 1.460, 11.240, 0.810}},
            {"WT12X88", {12.900, 1.340, 11.260, 0.750}},
            {"WT12X81", {13.000, 1.220, 11.280, 0.705}},
            {"WT12X73", {12.900, 1.090, 11.310, 0.650}},
            {"WT12X65.5", {12.900, 0.960, 11.240, 0.605}},
            {"WT12X58.5", {12.800, 0.850, 11.250, 0.550}},
            {"WT12X52", {12.800, 0.750, 11.250, 0.500}},
            {"WT12X51.5", {9.000, 0.980, 11.320, 0.550}},
            {"WT12X47", {9.070, 0.875, 11.325, 0.515}},
            {"WT12X42", {9.020, 0.770, 11.330, 0.470}},
            {"WT12X38", {8.990, 0.680, 11.320, 0.440}},
            {"WT12X34", {8.970, 0.585, 11.315, 0.415}},
            {"WT12X31", {7.040, 0.590, 11.310, 0.430}},
            {"WT12X27.5", {7.010, 0.505, 11.295, 0.395}},
            {"WT10.5X137.5", {12.900, 2.190, 9.910, 1.220}},
            {"WT10.5X124", {12.800, 1.990, 9.910, 1.100}},
            {"WT10.5X111.5", {12.700, 1.790, 9.910, 1.000}},
            {"WT10.5X100.5", {12.600, 1.630, 9.870, 0.910}},
            {"WT10.5X91", {12.500, 1.480, 9.920, 0.830}},
            {"WT10.5X83", {12.400, 1.360, 9.840, 0.750}},
            {"WT10.5X73.5", {12.500, 1.150, 9.850, 0.720}},
            {"WT10.5X66", {12.400, 1.040, 9.860, 0.650}},
            {"WT10.5X61", {12.400, 0.960, 9.840, 0.600}},
            {"WT10.5X55.5", {12.300, 0.875, 9.925, 0.550}},
            {"WT10.5X50.5", {12.300, 0.800, 9.900, 0.500}},
            {"WT10.5X46.5", {8.420, 0.930, 9.870, 0.580}},
            {"WT10.5X41.5", {8.360, 0.835, 9.865, 0.515}},
            {"WT10.5X36.5", {8.300, 0.740, 9.860, 0.455}},
            {"WT10.5X34", {8.270, 0.685, 9.915, 0.430}},
            {"WT10.5X31", {8.240, 0.615, 9.885, 0.400}},
            {"WT10.5X27.5", {8.220, 0.522, 9.878, 0.375}},
            {"WT10.5X24", {8.140, 0.430, 9.870, 0.350}},
            {"WT10.5X28.5", {6.560, 0.650, 9.850, 0.405}},
            {"WT10.5X25", {6.530, 0.535, 9.865, 0.380}},
            {"WT10.5X22", {6.500, 0.450, 9.850, 0.350}},
            {"WT9X155.5", {12.000, 2.740, 8.460, 1.520}},
            {"WT9X141.5", {11.900, 2.500, 8.400, 1.400}},
            {"WT9X129", {11.800, 2.300, 8.400, 1.280}},
            {"WT9X117", {11.700, 2.110, 8.390, 1.160}},
            {"WT9X105.5", {11.600, 1.910, 8.390, 1.060}},
            {"WT9X96", {11.500, 1.750, 8.450, 0.960}},
            {"WT9X87.5", {11.400, 1.590, 8.410, 0.890}},
            {"WT9X79", {11.300, 1.440, 8.420, 0.810}},
            {"WT9X71.5", {11.200, 1.320, 8.430, 0.730}},
            {"WT9X65", {11.200, 1.200, 8.430, 0.670}},
            {"WT9X59.5", {11.300, 1.060, 8.430, 0.655}},
            {"WT9X53", {11.200, 0.940, 8.430, 0.590}},
            {"WT9X48.5", {11.100, 0.870, 8.430, 0.535}},
            {"WT9X43", {11.100, 0.770, 8.430, 0.480}},
            {"WT9X38", {11.000, 0.680, 8.430, 0.425}},
            {"WT9X35.5", {7.640, 0.810, 8.430, 0.495}},
            {"WT9X32.5", {7.590, 0.750, 8.430, 0.450}},
            {"WT9X30", {7.560, 0.695, 8.425, 0.415}},
            {"WT9X27.5", {7.530, 0.630, 8.430, 0.390}},
            {"WT9X25", {7.500, 0.570, 8.430, 0.355}},
            {"WT9X23", {6.060, 0.605, 8.425, 0.360}},
            {"WT9X20", {6.020, 0.525, 8.425, 0.315}},
            {"WT9X17.5", {6.000, 0.425, 8.425, 0.300}},
            {"WT8X50", {10.400, 0.985, 7.505, 0.585}},
            {"WT8X44.5", {10.400, 0.875, 7.505, 0.525}},
            {"WT8X38.5", {10.300, 0.760, 7.500, 0.455}},
            {"WT8X33.5", {10.200, 0.665, 7.505, 0.395}},
            {"WT8X28.5", {7.120, 0.715, 7.505, 0.430}},
            {"WT8X25", {7.070, 0.630, 7.500, 0.380}},
            {"WT8X22.5", {7.040, 0.565, 7.505, 0.345}},
            {"WT8X20", {7.000, 0.505, 7.505, 0.305}},
            {"WT8X18", {6.990, 0.430, 7.500, 0.295}},
            {"WT8X15.5", {5.530, 0.440, 7.500, 0.275}},
            {"WT8X13", {5.500, 0.345, 7.505, 0.250}},
            {"WT7X436.5", {18.800, 5.510, 6.290, 3.940}},
            {"WT7X404", {18.600, 5.120, 6.280, 3.740}},
            {"WT7X365", {17.900, 4.910, 6.290, 3.070}},
            {"WT7X332.5", {17.700, 4.520, 6.280, 2.830}},
            {"WT7X302.5", {17.400, 4.160, 6.340, 2.600}},
            {"WT7X275", {17.200, 3.820, 6.280, 2.380}},
            {"WT7X250", {17.000, 3.500, 6.300, 2.190}},
            {"WT7X227.5", {16.800, 3.210, 6.300, 2.020}},
            {"WT7X213", {16.700, 3.040, 6.300, 1.880}},
            {"WT7X199", {16.600, 2.850, 6.300, 1.770}},
            {"WT7X185", {16.500, 2.660, 6.300, 1.660}},
            {"WT7X171", {16.400, 2.470, 6.300, 1.540}},
            {"WT7X155.5", {16.200, 2.260, 6.300, 1.410}},
            {"WT7X141.5", {16.100, 2.070, 6.300, 1.290}},
            {"WT7X128.5", {16.000, 1.890, 6.300, 1.180}},
            {"WT7X116.5", {15.900, 1.720, 6.300, 1.070}},
            {"WT7X105.5", {15.800, 1.560, 6.300, 0.980}},
            {"WT7X96.5", {15.700, 1.440, 6.300, 0.890}},
            {"WT7X88", {15.700, 1.310, 6.300, 0.830}},
            {"WT7X79.5", {15.600, 1.190, 6.300, 0.745}},
            {"WT7X72.5", {15.500, 1.090, 6.300, 0.680}},
            {"WT7X66", {14.700, 1.030, 6.300, 0.645}},
            {"WT7X60", {14.700, 0.940, 6.300, 0.590}},
            {"WT7X54.5", {14.600, 0.860, 6.300, 0.525}},
            {"WT7X49.5", {14.600, 0.780, 6.300, 0.485}},
            {"WT7X45", {14.500, 0.710, 6.300, 0.440}},
            {"WT7X41", {10.100, 0.855, 6.305, 0.510}},
            {"WT7X37", {10.100, 0.785, 6.305, 0.450}},
            {"WT7X34", {10.000, 0.720, 6.300, 0.415}},
            {"WT7X30.5", {10.000, 0.645, 6.305, 0.375}},
            {"WT7X26.5", {8.060, 0.660, 6.300, 0.370}},
            {"WT7X24", {8.030, 0.595, 6.305, 0.340}},
            {"WT7X21.5", {8.000, 0.530, 6.300, 0.305}},
            {"WT7X19", {6.770, 0.515, 6.535, 0.310}},
            {"WT7X17", {6.750, 0.455, 6.535, 0.285}},
            {"WT7X15", {6.730, 0.385, 6.535, 0.270}},
            {"WT7X13", {5.030, 0.420, 6.540, 0.255}},
            {"WT7X11", {5.000, 0.335, 6.535, 0.230}},
            {"WT6X168", {13.400, 2.960, 5.450, 1.780}},
            {"WT6X152.5", {13.200, 2.710, 5.450, 1.630}},
            {"WT6X139.5", {13.100, 2.470, 5.460, 1.530}},
            {"WT6X126", {13.000, 2.250, 5.460, 1.400}},
            {"WT6X115", {12.900, 2.070, 5.460, 1.290}},
            {"WT6X105", {12.800, 1.900, 5.460, 1.180}},
            {"WT6X95", {12.700, 1.740, 5.450, 1.060}},
            {"WT6X85", {12.600, 1.560, 5.460, 0.960}},
            {"WT6X76", {12.500, 1.400, 5.460, 0.870}},
            {"WT6X68", {12.400, 1.250, 5.460, 0.790}},
            {"WT6X60", {12.300, 1.110, 5.450, 0.710}},
            {"WT6X53", {12.200, 0.990, 5.460, 0.610}},
            {"WT6X48", {12.200, 0.900, 5.460, 0.550}},
            {"WT6X43.5", {12.100, 0.810, 5.460, 0.515}},
            {"WT6X39.5", {12.100, 0.735, 5.455, 0.470}},
            {"WT6X36", {12.000, 0.670, 5.460, 0.430}},
            {"WT6X32.5", {12.000, 0.605, 5.455, 0.390}},
            {"WT6X29", {10.000, 0.640, 5.460, 0.360}},
            {"WT6X26.5", {10.000, 0.575, 5.455, 0.345}},
            {"WT6X25", {8.080, 0.640, 5.460, 0.370}},
            {"WT6X22.5", {8.050, 0.575, 5.455, 0.335}},
            {"WT6X20", {8.010, 0.515, 5.455, 0.295}},
            {"WT6X17.5", {6.560, 0.520, 5.730, 0.300}},
            {"WT6X15", {6.520, 0.440, 5.730, 0.260}},
            {"WT6X13", {6.490, 0.380, 5.730, 0.230}},
            {"WT6X11", {4.030, 0.425, 5.735, 0.260}},
            {"WT6X9.5", {4.010, 0.350, 5.730, 0.235}},
            {"WT6X8", {3.990, 0.265, 5.735, 0.220}},
            {"WT6X7", {3.970, 0.225, 5.735, 0.200}},
            {"WT5X56", {10.400, 1.250, 4.430, 0.755}},
            {"WT5X50", {10.300, 1.120, 4.430, 0.680}},
            {"WT5X44", {10.300, 0.990, 4.430, 0.605}},
            {"WT5X38.5", {10.200, 0.870, 4.430, 0.530}},
            {"WT5X34", {10.100, 0.770, 4.430, 0.470}},
            {"WT5X30", {10.100, 0.680, 4.430, 0.420}},
            {"WT5X27", {10.000, 0.615, 4.435, 0.370}},
            {"WT5X24.5", {10.000, 0.560, 4.430, 0.340}},
            {"WT5X22.5", {8.020, 0.620, 4.430, 0.350}},
            {"WT5X19.5", {7.990, 0.530, 4.430, 0.315}},
            {"WT5X16.5", {7.960, 0.435, 4.435, 0.290}},
            {"WT5X15", {5.810, 0.510, 4.730, 0.300}},
            {"WT5X13", {5.770, 0.440, 4.730, 0.260}},
            {"WT5X11", {5.750, 0.360, 4.730, 0.240}},
            {"WT5X9.5", {4.020, 0.395, 4.725, 0.250}},
            {"WT5X8.5", {4.010, 0.330, 4.730, 0.240}},
            {"WT5X7.5", {4.000, 0.270, 4.730, 0.230}},
            {"WT5X6", {3.960, 0.210, 4.730, 0.190}},
            {"WT4X33.5", {8.280, 0.935, 3.565, 0.570}},
            {"WT4X29", {8.220, 0.810, 3.570, 0.510}},
            {"WT4X24", {8.110, 0.685, 3.565, 0.400}},
            {"WT4X20", {8.070, 0.560, 3.570, 0.360}},
            {"WT4X17.5", {8.020, 0.495, 3.565, 0.310}},
            {"WT4X15.5", {8.000, 0.435, 3.565, 0.285}},
            {"WT4X14", {6.540, 0.465, 3.565, 0.285}},
            {"WT4X12", {6.500, 0.400, 3.570, 0.245}},
            {"WT4X10.5", {5.270, 0.400, 3.740, 0.250}},
            {"WT4X9", {5.250, 0.330, 3.740, 0.230}},
            {"WT4X7.5", {4.020, 0.315, 3.745, 0.245}},
            {"WT4X6.5", {4.000, 0.255, 3.745, 0.230}},
            {"WT4X5", {3.940, 0.205, 3.745, 0.170}},
            {"WT3X12.5", {6.080, 0.455, 2.735, 0.320}},
            {"WT3X10", {6.020, 0.365, 2.735, 0.260}},
            {"WT3X7.5", {5.990, 0.260, 2.740, 0.230}},
            {"WT3X8", {4.030, 0.405, 2.735, 0.260}},
            {"WT3X6", {4.000, 0.280, 2.740, 0.230}},
            {"WT3X4.5", {3.940, 0.215, 2.735, 0.170}},
            {"WT3X4.25", {3.940, 0.195, 2.725, 0.170}},
            {"WT2.5X9.5", {5.030, 0.430, 2.150, 0.270}},
            {"WT2.5X8", {5.000, 0.360, 2.150, 0.240}},
            {"WT2X6.5", {4.060, 0.345, 1.735, 0.280}},
            {"MT6.25X6.2", {3.750, 0.228, 6.042, 0.155}},
            {"MT6.25X5.8", {3.500, 0.211, 6.039, 0.155}},
            {"MT6X5.9", {3.070, 0.225, 5.775, 0.177}},
            {"MT6X5.4", {3.070, 0.210, 5.780, 0.160}},
            {"MT6X5", {3.250, 0.180, 5.810, 0.149}},
            {"MT5X4.5", {2.690, 0.206, 4.794, 0.157}},
            {"MT5X4", {2.690, 0.182, 4.798, 0.141}},
            {"MT5X3.75", {2.690, 0.173, 4.827, 0.130}},
            {"MT4X3.25", {2.280, 0.189, 3.811, 0.135}},
            {"MT4X3.1", {2.280, 0.177, 3.823, 0.129}},
            {"MT3X2.2", {1.840, 0.171, 2.829, 0.114}},
            {"MT3X1.85", {2.000, 0.129, 2.831, 0.098}},
            {"MT2.5X9.45", {5.000, 0.416, 2.084, 0.316}},
            {"MT2X3", {3.800, 0.160, 1.740, 0.130}},
            {"ST12X60.5", {8.050, 1.090, 11.210, 0.800}},
            {"ST12X53", {7.870, 1.090, 11.210, 0.620}},
            {"ST12X50", {7.250, 0.870, 11.130, 0.745}},
            {"ST12X45", {7.130, 0.870, 11.130, 0.625}},
            {"ST12X40", {7.000, 0.870, 11.130, 0.500}},
            {"ST10X48", {7.200, 0.920, 9.280, 0.800}},
            {"ST10X43", {7.060, 0.920, 9.280, 0.660}},
            {"ST10X37.5", {6.390, 0.795, 9.205, 0.635}},
            {"ST10X33", {6.260, 0.795, 9.205, 0.505}},
            {"ST9X35", {6.250, 0.691, 8.309, 0.711}},
            {"ST9X27.35", {6.000, 0.691, 8.309, 0.461}},
            {"ST7.5X25", {5.640, 0.622, 6.878, 0.550}},
            {"ST7.5X21.45", {5.500, 0.622, 6.878, 0.411}},
            {"ST6X25", {5.480, 0.659, 5.341, 0.687}},
            {"ST6X20.4", {5.250, 0.659, 5.341, 0.462}},
            {"ST6X17.5", {5.080, 0.544, 5.456, 0.428}},
            {"ST6X15.9", {5.000, 0.544, 5.456, 0.350}},
            {"ST5X17.5", {4.940, 0.491, 4.509, 0.594}},
            {"ST5X12.7", {4.660, 0.491, 4.509, 0.311}},
            {"ST4X11.5", {4.170, 0.425, 3.575, 0.441}},
            {"ST4X9.2", {4.000, 0.425, 3.575, 0.271}},
            {"ST3X8.6", {3.570, 0.359, 2.641, 0.465}},
            {"ST3X6.25", {3.330, 0.359, 2.641, 0.232}},
            {"ST2.5X5", {3.000, 0.326, 2.174, 0.214}},
            {"ST2X4.75", {2.800, 0.293, 1.707, 0.326}},
            {"ST2X3.85", {2.660, 0.293, 1.707, 0.193}},
            {"ST1.5X3.75", {2.510, 0.260, 1.240, 0.349}},
            {"ST1.5X2.85", {2.330, 0.260, 1.240, 0.170}},
        };

        return registry.fetch(type);
    }

    // Adjusted registry: dimensions optimized against the AISC tabulated
    // section properties with section_optimizer.py.  The ST entries use a
    // 7 % dimension margin; the two smallest ST shapes relax the weak axis.
    vec ustsection_adjusted(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"WT22X204", {16.0701, 2.1804, 20.1925, 1.225}},
            {"WT22X184", {15.9781, 1.9768, 20.2141, 1.1033}},
            {"WT22X167.5", {15.8865, 1.7891, 20.2453, 1.0231}},
            {"WT22X145", {15.804, 1.5781, 20.2356, 0.8669}},
            {"WT22X131", {15.6767, 1.434, 20.2275, 0.7863}},
            {"WT22X115", {15.6778, 1.2358, 20.1625, 0.7164}},
            {"WT20X327.5", {16.8625, 3.571, 18.2079, 1.9761}},
            {"WT20X296.5", {16.6037, 3.2774, 18.2241, 1.7902}},
            {"WT20X251.5", {16.2729, 2.822, 18.2379, 1.5303}},
            {"WT20X215.5", {16.1176, 2.4035, 18.2191, 1.3425}},
            {"WT20X198.5", {16.0241, 2.2374, 18.2146, 1.2241}},
            {"WT20X186", {15.9578, 2.0815, 18.1879, 1.1679}},
            {"WT20X181", {15.906, 2.051, 18.2309, 1.1204}},
            {"WT20X162", {15.8108, 1.8397, 18.205, 1.0107}},
            {"WT20X148.5", {15.7562, 1.6673, 18.17, 0.9438}},
            {"WT20X138.5", {15.7015, 1.6116, 18.2019, 0.8365}},
            {"WT20X124.5", {15.5904, 1.46, 18.1677, 0.7576}},
            {"WT20X107.5", {15.5822, 1.257, 18.0871, 0.666}},
            {"WT20X99.5", {15.5419, 1.1044, 18.2258, 0.6521}},
            {"WT20X196", {12.2245, 2.6016, 18.1861, 1.4196}},
            {"WT20X165.5", {12.0629, 2.1777, 18.1926, 1.2276}},
            {"WT20X163.5", {12.0438, 2.1773, 18.2403, 1.1798}},
            {"WT20X147", {11.9096, 1.9793, 18.2278, 1.062}},
            {"WT20X139", {11.8708, 1.8562, 18.2248, 1.0308}},
            {"WT20X132", {11.8059, 1.7808, 18.234, 0.9617}},
            {"WT20X117.5", {11.7401, 1.6338, 18.2111, 0.8336}},
            {"WT20X105.5", {11.649, 1.4701, 18.1719, 0.7579}},
            {"WT20X91.5", {11.6018, 1.26, 18.22, 0.6531}},
            {"WT20X83.5", {11.6419, 1.0655, 18.2221, 0.6554}},
            {"WT20X74.5", {11.5961, 0.8715, 18.2152, 0.6351}},
            {"WT18X462.5", {18.5928, 4.5374, 17.061, 3.024}},
            {"WT18X426.5", {18.151, 4.5677, 17.0191, 2.5257}},
            {"WT18X401", {17.9241, 4.3364, 17.0031, 2.3668}},
            {"WT18X361.5", {17.7385, 3.9412, 16.9694, 2.1682}},
            {"WT18X326", {17.5042, 3.5755, 16.9299, 1.9777}},
            {"WT18X264.5", {17.1385, 2.9413, 16.924, 1.6174}},
            {"WT18X243.5", {17.0219, 2.7127, 16.9225, 1.5065}},
            {"WT18X220.5", {16.8862, 2.4745, 16.9756, 1.3576}},
            {"WT18X197.5", {16.7807, 2.2177, 16.9709, 1.2223}},
            {"WT18X180.5", {16.6425, 2.0384, 16.9978, 1.1146}},
            {"WT18X165", {16.5924, 1.8604, 16.9547, 1.0248}},
            {"WT18X151", {16.5852, 1.6973, 16.8948, 0.9571}},
            {"WT18X141", {16.4523, 1.6099, 16.976, 0.8802}},
            {"WT18X131", {16.4212, 1.4728, 16.9841, 0.8379}},
            {"WT18X123.5", {16.3802, 1.3793, 16.9383, 0.8011}},
            {"WT18X115.5", {16.3557, 1.2846, 16.9263, 0.766}},
            {"WT18X193.5", {12.6532, 2.5875, 16.922, 1.4194}},
            {"WT18X175", {12.5485, 2.3484, 16.9963, 1.2923}},
            {"WT18X159", {12.3514, 2.1548, 16.941, 1.1804}},
            {"WT18X143", {12.2427, 1.955, 16.9818, 1.0579}},
            {"WT18X128", {12.1542, 1.7545, 17.0318, 0.9514}},
            {"WT18X116", {12.0692, 1.5889, 17.0021, 0.8672}},
            {"WT18X105", {12.1415, 1.3736, 16.9597, 0.8335}},
            {"WT18X97", {12.0354, 1.2815, 16.9498, 0.7683}},
            {"WT18X91", {12.0252, 1.1951, 16.953, 0.7285}},
            {"WT18X85", {11.9425, 1.1221, 16.9788, 0.6782}},
            {"WT18X80", {11.9048, 1.0415, 16.9637, 0.6513}},
            {"WT18X75", {11.9049, 0.9559, 16.9601, 0.6276}},
            {"WT18X67.5", {11.8687, 0.8068, 16.9618, 0.6029}},
            {"WT16.5X193.5", {16.1527, 2.2977, 15.636, 1.2685}},
            {"WT16.5X177", {16.0319, 2.1158, 15.6848, 1.1552}},
            {"WT16.5X159", {15.9603, 1.8979, 15.6085, 1.0537}},
            {"WT16.5X145.5", {15.8603, 1.7425, 15.6998, 0.9616}},
            {"WT16.5X131.5", {15.7566, 1.5819, 15.657, 0.876}},
            {"WT16.5X120.5", {15.7305, 1.4312, 15.682, 0.826}},
            {"WT16.5X110.5", {15.6981, 1.2989, 15.6736, 0.7734}},
            {"WT16.5X100.5", {15.647, 1.1704, 15.638, 0.7204}},
            {"WT16.5X84.5", {11.4535, 1.2333, 15.682, 0.6708}},
            {"WT16.5X76", {11.4393, 1.0849, 15.6312, 0.6393}},
            {"WT16.5X70.5", {11.4808, 0.9721, 15.6861, 0.6057}},
            {"WT16.5X65", {11.4475, 0.8687, 15.6789, 0.5809}},
            {"WT16.5X59", {11.3615, 0.761, 15.6473, 0.5534}},
            {"WT15X195.5", {15.5402, 2.4634, 14.1255, 1.3613}},
            {"WT15X178.5", {15.4194, 2.2612, 14.0887, 1.253}},
            {"WT15X163", {15.3019, 2.076, 14.101, 1.1474}},
            {"WT15X146", {15.1981, 1.868, 14.1287, 1.0213}},
            {"WT15X130.5", {15.1085, 1.6649, 14.093, 0.9416}},
            {"WT15X117.5", {14.9563, 1.5265, 14.1027, 0.8343}},
            {"WT15X105.5", {15.0227, 1.3342, 14.126, 0.7777}},
            {"WT15X95.5", {14.9853, 1.1963, 14.1451, 0.711}},
            {"WT15X86.5", {14.9264, 1.0768, 14.1343, 0.6574}},
            {"WT15X74", {10.4475, 1.1933, 14.1125, 0.6545}},
            {"WT15X66", {10.4591, 1.0229, 14.1038, 0.6189}},
            {"WT15X62", {10.4522, 0.9465, 14.1571, 0.5843}},
            {"WT15X58", {10.4177, 0.8678, 14.1305, 0.567}},
            {"WT15X54", {10.393, 0.7761, 14.1195, 0.5492}},
            {"WT15X49.5", {10.3776, 0.6833, 14.1438, 0.521}},
            {"WT15X45", {10.3022, 0.6258, 14.1166, 0.4735}},
            {"WT13.5X269.5", {15.3061, 3.5157, 12.7113, 1.9914}},
            {"WT13.5X184", {14.6113, 2.5091, 12.6826, 1.3833}},
            {"WT13.5X168", {14.4878, 2.3069, 12.6863, 1.264}},
            {"WT13.5X153.5", {14.4066, 2.1045, 12.7051, 1.1586}},
            {"WT13.5X140.5", {14.2913, 1.9531, 12.6857, 1.0627}},
            {"WT13.5X129", {14.194, 1.7974, 12.6735, 0.9848}},
            {"WT13.5X117.5", {14.1171, 1.6327, 12.6878, 0.9131}},
            {"WT13.5X108.5", {14.0427, 1.521, 12.6711, 0.8347}},
            {"WT13.5X97", {13.962, 1.3584, 12.6728, 0.7545}},
            {"WT13.5X89", {14.004, 1.2106, 12.6588, 0.7315}},
            {"WT13.5X80.5", {13.9183, 1.1005, 12.674, 0.6638}},
            {"WT13.5X73", {13.8703, 0.995, 12.6759, 0.6088}},
            {"WT13.5X64.5", {9.9623, 1.1151, 12.6879, 0.6117}},
            {"WT13.5X57", {10.0107, 0.9442, 12.6678, 0.5752}},
            {"WT13.5X51", {9.9435, 0.8465, 12.7125, 0.5147}},
            {"WT13.5X47", {9.9244, 0.7587, 12.7025, 0.4913}},
            {"WT13.5X42", {9.8605, 0.6581, 12.6801, 0.4623}},
            {"WT12X185", {13.6535, 2.7201, 11.2498, 1.5312}},
            {"WT12X167.5", {13.475, 2.5051, 11.3165, 1.3581}},
            {"WT12X153", {13.3802, 2.294, 11.2893, 1.2546}},
            {"WT12X139.5", {13.2938, 2.0958, 11.2598, 1.1633}},
            {"WT12X125", {13.1539, 1.9018, 11.2645, 1.042}},
            {"WT12X114.5", {13.0975, 1.7361, 11.2713, 0.9617}},
            {"WT12X103.5", {13.0067, 1.5722, 11.2674, 0.8729}},
            {"WT12X96", {12.9157, 1.4717, 11.2544, 0.8124}},
            {"WT12X88", {12.8799, 1.345, 11.2863, 0.749}},
            {"WT12X81", {12.8931, 1.2344, 11.2611, 0.7057}},
            {"WT12X73", {12.8481, 1.1014, 11.2723, 0.6505}},
            {"WT12X65.5", {12.8031, 0.9698, 11.2472, 0.609}},
            {"WT12X58.5", {12.7639, 0.8578, 11.2616, 0.5517}},
            {"WT12X52", {12.7207, 0.7561, 11.2717, 0.5009}},
            {"WT12X51.5", {8.9773, 0.9867, 11.2752, 0.5515}},
            {"WT12X47", {9.0532, 0.8783, 11.2645, 0.5166}},
            {"WT12X42", {8.955, 0.7854, 11.2398, 0.4732}},
            {"WT12X38", {8.9473, 0.6893, 11.2325, 0.4449}},
            {"WT12X34", {8.9076, 0.5962, 11.2571, 0.4158}},
            {"WT12X31", {6.9904, 0.6008, 11.2471, 0.4343}},
            {"WT12X27.5", {6.9472, 0.516, 11.2749, 0.3976}},
            {"WT10.5X137.5", {12.824, 2.2294, 9.8217, 1.2398}},
            {"WT10.5X124", {12.6838, 2.0417, 9.8183, 1.116}},
            {"WT10.5X111.5", {12.605, 1.8304, 9.8283, 1.0166}},
            {"WT10.5X100.5", {12.5611, 1.636, 9.8658, 0.9137}},
            {"WT10.5X91", {12.4656, 1.4893, 9.8535, 0.8332}},
            {"WT10.5X83", {12.3634, 1.3753, 9.885, 0.7469}},
            {"WT10.5X73.5", {12.4845, 1.1567, 9.8509, 0.7242}},
            {"WT10.5X66", {12.3932, 1.0444, 9.8562, 0.6529}},
            {"WT10.5X61", {12.3281, 0.9722, 9.8952, 0.5971}},
            {"WT10.5X55.5", {12.2901, 0.8843, 9.8956, 0.5476}},
            {"WT10.5X50.5", {12.242, 0.8088, 9.8731, 0.5015}},
            {"WT10.5X46.5", {8.3683, 0.9448, 9.85, 0.5827}},
            {"WT10.5X41.5", {8.3074, 0.8474, 9.8714, 0.5172}},
            {"WT10.5X36.5", {8.2578, 0.7504, 9.8548, 0.4565}},
            {"WT10.5X34", {8.2287, 0.6955, 9.8735, 0.4309}},
            {"WT10.5X31", {8.1655, 0.6309, 9.8492, 0.4027}},
            {"WT10.5X27.5", {8.1478, 0.5347, 9.8609, 0.3763}},
            {"WT10.5X24", {8.0593, 0.4425, 9.8633, 0.351}},
            {"WT10.5X28.5", {6.5018, 0.6643, 9.8561, 0.4079}},
            {"WT10.5X25", {6.4666, 0.5511, 9.8688, 0.3802}},
            {"WT10.5X22", {6.4227, 0.4639, 9.8567, 0.3531}},
            {"WT9X155.5", {11.9748, 2.7652, 8.4121, 1.5115}},
            {"WT9X141.5", {11.8617, 2.5154, 8.397, 1.406}},
            {"WT9X129", {11.7426, 2.3146, 8.4215, 1.2782}},
            {"WT9X117", {11.6604, 2.104, 8.4295, 1.1608}},
            {"WT9X105.5", {11.5013, 1.9317, 8.3945, 1.0624}},
            {"WT9X96", {11.4688, 1.7457, 8.4228, 0.9615}},
            {"WT9X87.5", {11.3758, 1.5921, 8.3953, 0.8986}},
            {"WT9X79", {11.2901, 1.4466, 8.4056, 0.8128}},
            {"WT9X71.5", {11.2164, 1.3237, 8.4384, 0.7269}},
            {"WT9X65", {11.0958, 1.2174, 8.4021, 0.6716}},
            {"WT9X59.5", {11.1862, 1.0771, 8.4184, 0.6546}},
            {"WT9X53", {11.1495, 0.9501, 8.3942, 0.5931}},
            {"WT9X48.5", {11.1122, 0.8739, 8.4317, 0.5334}},
            {"WT9X43", {11.0312, 0.7814, 8.3991, 0.4826}},
            {"WT9X38", {11.0241, 0.682, 8.427, 0.4249}},
            {"WT9X35.5", {7.625, 0.8117, 8.4073, 0.499}},
            {"WT9X32.5", {7.5693, 0.755, 8.4103, 0.4522}},
            {"WT9X30", {7.5152, 0.7048, 8.4043, 0.4177}},
            {"WT9X27.5", {7.5056, 0.6365, 8.4018, 0.3931}},
            {"WT9X25", {7.4467, 0.5795, 8.3948, 0.3583}},
            {"WT9X23", {6.0394, 0.6126, 8.4077, 0.3624}},
            {"WT9X20", {5.971, 0.5361, 8.4046, 0.3164}},
            {"WT9X17.5", {5.9477, 0.4353, 8.4089, 0.3017}},
            {"WT8X50", {10.407, 0.9893, 7.4953, 0.5857}},
            {"WT8X44.5", {10.3247, 0.8841, 7.486, 0.526}},
            {"WT8X38.5", {10.2695, 0.765, 7.486, 0.4565}},
            {"WT8X33.5", {10.2115, 0.6694, 7.5025, 0.3941}},
            {"WT8X28.5", {7.0915, 0.724, 7.4801, 0.4321}},
            {"WT8X25", {7.0335, 0.6393, 7.4851, 0.3811}},
            {"WT8X22.5", {6.996, 0.5731, 7.4711, 0.3488}},
            {"WT8X20", {6.941, 0.5152, 7.4726, 0.3074}},
            {"WT8X18", {6.9118, 0.4419, 7.4834, 0.2963}},
            {"WT8X15.5", {5.4706, 0.4524, 7.492, 0.2757}},
            {"WT8X13", {5.4357, 0.356, 7.4829, 0.2515}},
            {"WT7X436.5", {18.7264, 5.5686, 6.2582, 3.8937}},
            {"WT7X404", {18.568, 5.1398, 6.2461, 3.7645}},
            {"WT7X365", {17.9118, 4.8974, 6.3159, 3.0709}},
            {"WT7X332.5", {17.5555, 4.5892, 6.3467, 2.7036}},
            {"WT7X302.5", {17.378, 4.1857, 6.315, 2.5555}},
            {"WT7X275", {17.2286, 3.807, 6.2326, 2.4601}},
            {"WT7X250", {16.9645, 3.5257, 6.2923, 2.1707}},
            {"WT7X227.5", {16.8493, 3.1997, 6.2262, 2.0914}},
            {"WT7X213", {16.6662, 3.0485, 6.2852, 1.8793}},
            {"WT7X199", {16.6533, 2.824, 6.2523, 1.8198}},
            {"WT7X185", {16.4634, 2.6661, 6.2622, 1.6749}},
            {"WT7X171", {16.3314, 2.4811, 6.2771, 1.5478}},
            {"WT7X155.5", {16.2215, 2.2638, 6.2715, 1.4261}},
            {"WT7X141.5", {16.0905, 2.0761, 6.2636, 1.3059}},
            {"WT7X128.5", {15.9705, 1.8965, 6.2573, 1.193}},
            {"WT7X116.5", {15.9015, 1.7151, 6.2582, 1.0935}},
            {"WT7X105.5", {15.7477, 1.5736, 6.2961, 0.9799}},
            {"WT7X96.5", {15.667, 1.4515, 6.2903, 0.8901}},
            {"WT7X88", {15.6173, 1.3178, 6.2607, 0.8423}},
            {"WT7X79.5", {15.4939, 1.2046, 6.2814, 0.7465}},
            {"WT7X72.5", {15.4462, 1.0993, 6.2813, 0.6834}},
            {"WT7X66", {14.6582, 1.0422, 6.2681, 0.6512}},
            {"WT7X60", {14.5578, 0.959, 6.2915, 0.5877}},
            {"WT7X54.5", {14.5283, 0.8713, 6.2851, 0.5265}},
            {"WT7X49.5", {14.4602, 0.7959, 6.292, 0.4837}},
            {"WT7X45", {14.4532, 0.7188, 6.2664, 0.4466}},
            {"WT7X41", {10.0819, 0.8662, 6.2553, 0.5198}},
            {"WT7X37", {10.0188, 0.7958, 6.2396, 0.462}},
            {"WT7X34", {9.9721, 0.7323, 6.2342, 0.4263}},
            {"WT7X30.5", {9.8995, 0.6627, 6.2614, 0.3791}},
            {"WT7X26.5", {7.9623, 0.6828, 6.262, 0.3737}},
            {"WT7X24", {7.9373, 0.6144, 6.2606, 0.3446}},
            {"WT7X21.5", {7.896, 0.5486, 6.2423, 0.3111}},
            {"WT7X19", {6.7072, 0.5275, 6.5229, 0.3108}},
            {"WT7X17", {6.6699, 0.4679, 6.5139, 0.2865}},
            {"WT7X15", {6.6665, 0.3952, 6.5256, 0.271}},
            {"WT7X13", {4.9715, 0.4322, 6.5022, 0.2581}},
            {"WT7X11", {4.9345, 0.3477, 6.4988, 0.2331}},
            {"WT6X168", {13.331, 2.9902, 5.4208, 1.7701}},
            {"WT6X152.5", {13.2735, 2.683, 5.4074, 1.6809}},
            {"WT6X139.5", {13.1111, 2.4873, 5.417, 1.5401}},
            {"WT6X126", {12.9654, 2.2709, 5.44, 1.3936}},
            {"WT6X115", {12.8807, 2.0775, 5.4548, 1.2884}},
            {"WT6X105", {12.75, 1.9168, 5.4384, 1.1807}},
            {"WT6X95", {12.633, 1.7504, 5.4182, 1.0723}},
            {"WT6X85", {12.5577, 1.5656, 5.4158, 0.9776}},
            {"WT6X76", {12.4188, 1.4189, 5.4407, 0.8692}},
            {"WT6X68", {12.3516, 1.2638, 5.4245, 0.7977}},
            {"WT6X60", {12.2541, 1.1202, 5.442, 0.7093}},
            {"WT6X53", {12.1844, 0.9995, 5.399, 0.6256}},
            {"WT6X48", {12.105, 0.9114, 5.4188, 0.5591}},
            {"WT6X43.5", {12.0198, 0.8279, 5.4133, 0.5218}},
            {"WT6X39.5", {12.0043, 0.7472, 5.4134, 0.4778}},
            {"WT6X36", {11.9261, 0.6879, 5.4234, 0.4343}},
            {"WT6X32.5", {11.9113, 0.6175, 5.4135, 0.3966}},
            {"WT6X29", {9.9249, 0.6545, 5.4178, 0.366}},
            {"WT6X26.5", {9.9124, 0.5878, 5.3992, 0.3534}},
            {"WT6X25", {8.032, 0.6502, 5.4253, 0.3748}},
            {"WT6X22.5", {7.9821, 0.588, 5.4253, 0.3389}},
            {"WT6X20", {7.9215, 0.5294, 5.4153, 0.2993}},
            {"WT6X17.5", {6.511, 0.5294, 5.7295, 0.2995}},
            {"WT6X15", {6.4924, 0.4461, 5.7256, 0.2603}},
            {"WT6X13", {6.4522, 0.386, 5.7264, 0.2304}},
            {"WT6X11", {3.9987, 0.435, 5.719, 0.2607}},
            {"WT6X9.5", {3.9638, 0.3604, 5.7109, 0.2368}},
            {"WT6X8", {3.9513, 0.2727, 5.712, 0.2225}},
            {"WT6X7", {3.941, 0.23, 5.7104, 0.2033}},
            {"WT5X56", {10.3911, 1.2594, 4.3842, 0.773}},
            {"WT5X50", {10.2692, 1.1394, 4.4225, 0.674}},
            {"WT5X44", {10.1973, 1.0078, 4.4101, 0.6068}},
            {"WT5X38.5", {10.1643, 0.8761, 4.3928, 0.5404}},
            {"WT5X34", {10.0705, 0.7819, 4.401, 0.4759}},
            {"WT5X30", {10.0058, 0.6943, 4.4156, 0.4218}},
            {"WT5X27", {9.9789, 0.6229, 4.3916, 0.3776}},
            {"WT5X24.5", {9.9424, 0.5685, 4.3919, 0.3473}},
            {"WT5X22.5", {7.9737, 0.63, 4.3737, 0.3606}},
            {"WT5X19.5", {7.9145, 0.5429, 4.3935, 0.3206}},
            {"WT5X16.5", {7.8785, 0.4476, 4.388, 0.2966}},
            {"WT5X15", {5.7854, 0.516, 4.7099, 0.3021}},
            {"WT5X13", {5.7292, 0.4487, 4.7126, 0.2608}},
            {"WT5X11", {5.7138, 0.3664, 4.7093, 0.2417}},
            {"WT5X9.5", {3.9953, 0.4022, 4.7063, 0.2525}},
            {"WT5X8.5", {3.9669, 0.3403, 4.7089, 0.2418}},
            {"WT5X7.5", {3.9488, 0.2808, 4.7111, 0.231}},
            {"WT5X6", {3.892, 0.2205, 4.7139, 0.1907}},
            {"WT4X33.5", {8.2615, 0.9405, 3.5522, 0.576}},
            {"WT4X29", {8.1937, 0.8162, 3.5422, 0.5173}},
            {"WT4X24", {8.0818, 0.6918, 3.5415, 0.4058}},
            {"WT4X20", {8.0184, 0.5692, 3.5416, 0.3648}},
            {"WT4X17.5", {7.9704, 0.5036, 3.55, 0.3123}},
            {"WT4X15.5", {7.9249, 0.445, 3.5487, 0.287}},
            {"WT4X14", {6.4779, 0.4755, 3.5391, 0.2897}},
            {"WT4X12", {6.4358, 0.4102, 3.5294, 0.2508}},
            {"WT4X10.5", {5.231, 0.408, 3.7275, 0.2513}},
            {"WT4X9", {5.1987, 0.3389, 3.7322, 0.2302}},
            {"WT4X7.5", {3.9651, 0.3256, 3.7244, 0.2469}},
            {"WT4X6.5", {3.9442, 0.2643, 3.7239, 0.2323}},
            {"WT4X5", {3.8778, 0.215, 3.7304, 0.1707}},
            {"WT3X12.5", {6.0478, 0.4618, 2.7382, 0.3174}},
            {"WT3X10", {5.9888, 0.3703, 2.7191, 0.2636}},
            {"WT3X7.5", {5.9556, 0.2641, 2.7291, 0.2312}},
            {"WT3X8", {3.9974, 0.4136, 2.7317, 0.259}},
            {"WT3X6", {3.9726, 0.2859, 2.7159, 0.2346}},
            {"WT3X4.5", {3.9007, 0.2215, 2.717, 0.1728}},
            {"WT3X4.25", {3.8865, 0.2026, 2.7005, 0.1728}},
            {"WT2.5X9.5", {4.9994, 0.4368, 2.1183, 0.2772}},
            {"WT2.5X8", {4.9683, 0.3659, 2.1208, 0.2466}},
            {"WT2X6.5", {4.0424, 0.3496, 1.7169, 0.2869}},
            {"MT6.25X6.2", {3.7024, 0.2355, 6.0232, 0.1562}},
            {"MT6.25X5.8", {3.4716, 0.2158, 6.0232, 0.157}},
            {"MT6X5.9", {3.0357, 0.231, 5.7558, 0.1795}},
            {"MT6X5.4", {3.0317, 0.2164, 5.7661, 0.1613}},
            {"MT6X5", {3.2051, 0.1871, 5.7818, 0.1513}},
            {"MT5X4.5", {2.6434, 0.2158, 4.7759, 0.1584}},
            {"MT5X4", {2.6416, 0.1911, 4.7682, 0.1428}},
            {"MT5X3.75", {2.6401, 0.1816, 4.8099, 0.1307}},
            {"MT4X3.25", {2.2397, 0.1985, 3.7967, 0.1353}},
            {"MT4X3.1", {2.2399, 0.1858, 3.809, 0.1297}},
            {"MT3X2.2", {1.8093, 0.1796, 2.8229, 0.114}},
            {"MT3X1.85", {1.9647, 0.1355, 2.8234, 0.0984}},
            {"MT2.5X9.45", {4.9783, 0.4219, 2.0567, 0.3266}},
            {"MT2X3", {3.7383, 0.168, 1.718, 0.1365}},
            {"ST12X60.5", {7.4953, 1.1663, 11.0837, 0.8078}},
            {"ST12X53", {7.3229, 1.1663, 11.0628, 0.6303}},
            {"ST12X50", {6.7425, 0.9137, 11.0559, 0.7649}},
            {"ST12X45", {6.6309, 0.914, 11.0623, 0.6378}},
            {"ST12X40", {6.51, 0.9144, 11.0529, 0.5097}},
            {"ST10X48", {6.696, 0.9839, 9.1733, 0.8096}},
            {"ST10X43", {6.5658, 0.9814, 9.1228, 0.6756}},
            {"ST10X37.5", {5.9427, 0.8374, 9.1178, 0.6519}},
            {"ST10X33", {5.8218, 0.8349, 9.1189, 0.5191}},
            {"ST9X35", {5.8125, 0.7174, 8.2444, 0.7362}},
            {"ST9X27.35", {5.58, 0.7203, 8.2481, 0.4716}},
            {"ST7.5X25", {5.2452, 0.6448, 6.8191, 0.5674}},
            {"ST7.5X21.45", {5.115, 0.6449, 6.7991, 0.4272}},
            {"ST6X25", {5.0964, 0.6944, 5.2797, 0.7071}},
            {"ST6X20.4", {4.8825, 0.6942, 5.2802, 0.4754}},
            {"ST6X17.5", {4.7244, 0.562, 5.4096, 0.4415}},
            {"ST6X15.9", {4.65, 0.5614, 5.3795, 0.366}},
            {"ST5X17.5", {4.5942, 0.5066, 4.4658, 0.6165}},
            {"ST5X12.7", {4.3338, 0.5032, 4.4209, 0.3309}},
            {"ST4X11.5", {3.8781, 0.4386, 3.5231, 0.4609}},
            {"ST4X9.2", {3.72, 0.4369, 3.4938, 0.29}},
            {"ST3X8.6", {3.3201, 0.3692, 2.6054, 0.4856}},
            {"ST3X6.25", {3.0969, 0.3706, 2.5768, 0.2482}},
            {"ST2.5X5", {2.79, 0.3388, 2.1162, 0.229}},
            {"ST2X4.75", {2.604, 0.3063, 1.6592, 0.3488}},
            {"ST2X3.85", {2.4738, 0.3105, 1.6539, 0.2065}},
            {"ST1.5X3.75", {2.3343, 0.2747, 1.2034, 0.3734}},
            {"ST1.5X2.85", {2.1669, 0.2772, 1.2056, 0.1819}},
        };

        return registry.fetch(type);
    }

    vec usrhssection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"HSS34X10X1", {8.605, 32.600, 0.930}},
            {"HSS34X10X7/8", {8.780, 32.800, 0.814}},
            {"HSS34X10X3/4", {8.955, 32.950, 0.698}},
            {"HSS34X10X5/8", {9.130, 33.150, 0.581}},
            {"HSS30X10X1", {8.605, 28.600, 0.930}},
            {"HSS30X10X7/8", {8.780, 28.800, 0.814}},
            {"HSS30X10X3/4", {8.955, 28.950, 0.698}},
            {"HSS30X10X5/8", {9.130, 29.150, 0.581}},
            {"HSS30X10X1/2", {9.305, 29.300, 0.465}},
            {"HSS24X20X3/4", {18.950, 22.950, 0.698}},
            {"HSS24X20X5/8", {19.150, 23.150, 0.581}},
            {"HSS24X20X1/2", {19.300, 23.300, 0.465}},
            {"HSS24X20X3/8", {19.500, 23.500, 0.349}},
            {"HSS24X20X5/16", {19.550, 23.550, 0.291}},
            {"HSS24X18X3/4", {16.950, 22.950, 0.698}},
            {"HSS24X18X5/8", {17.150, 23.150, 0.581}},
            {"HSS24X18X1/2", {17.300, 23.300, 0.465}},
            {"HSS24X18X3/8", {17.500, 23.500, 0.349}},
            {"HSS24X18X5/16", {17.550, 23.550, 0.291}},
            {"HSS24X16X3/4", {14.950, 22.950, 0.698}},
            {"HSS24X16X5/8", {15.150, 23.150, 0.581}},
            {"HSS24X16X1/2", {15.300, 23.300, 0.465}},
            {"HSS24X16X3/8", {15.500, 23.500, 0.349}},
            {"HSS24X16X5/16", {15.550, 23.550, 0.291}},
            {"HSS24X14X3/4", {12.950, 22.950, 0.698}},
            {"HSS24X14X5/8", {13.150, 23.150, 0.581}},
            {"HSS24X14X1/2", {13.300, 23.300, 0.465}},
            {"HSS24X14X3/8", {13.500, 23.500, 0.349}},
            {"HSS24X14X5/16", {13.550, 23.550, 0.291}},
            {"HSS24X14X1/4", {13.650, 23.650, 0.233}},
            {"HSS24X12X1", {10.605, 22.600, 0.930}},
            {"HSS24X12X7/8", {10.780, 22.800, 0.814}},
            {"HSS24X12X3/4", {10.955, 22.950, 0.698}},
            {"HSS24X12X5/8", {11.150, 23.150, 0.581}},
            {"HSS24X12X1/2", {11.300, 23.300, 0.465}},
            {"HSS24X12X3/8", {11.500, 23.500, 0.349}},
            {"HSS24X12X5/16", {11.550, 23.550, 0.291}},
            {"HSS24X12X1/4", {11.650, 23.650, 0.233}},
            {"HSS24X8X1/2", {7.305, 23.300, 0.465}},
            {"HSS24X8X3/8", {7.475, 23.500, 0.349}},
            {"HSS24X8X5/16", {7.565, 23.550, 0.291}},
            {"HSS24X8X1/4", {7.650, 23.650, 0.233}},
            {"HSS22X22X1", {20.600, 20.600, 0.930}},
            {"HSS22X22X7/8", {20.800, 20.800, 0.814}},
            {"HSS22X22X3/4", {20.950, 20.950, 0.698}},
            {"HSS22X22X5/8", {21.150, 21.150, 0.581}},
            {"HSS22X22X1/2", {21.300, 21.300, 0.465}},
            {"HSS22X20X3/4", {18.950, 20.950, 0.698}},
            {"HSS22X20X5/8", {19.150, 21.150, 0.581}},
            {"HSS22X20X1/2", {19.300, 21.300, 0.465}},
            {"HSS22X20X3/8", {19.500, 21.500, 0.349}},
            {"HSS22X20X5/16", {19.550, 21.550, 0.291}},
            {"HSS22X18X3/4", {16.950, 20.950, 0.698}},
            {"HSS22X18X5/8", {17.150, 21.150, 0.581}},
            {"HSS22X18X1/2", {17.300, 21.300, 0.465}},
            {"HSS22X18X3/8", {17.500, 21.500, 0.349}},
            {"HSS22X18X5/16", {17.550, 21.550, 0.291}},
            {"HSS22X16X3/4", {14.950, 20.950, 0.698}},
            {"HSS22X16X5/8", {15.150, 21.150, 0.581}},
            {"HSS22X16X1/2", {15.300, 21.300, 0.465}},
            {"HSS22X16X3/8", {15.500, 21.500, 0.349}},
            {"HSS22X16X5/16", {15.550, 21.550, 0.291}},
            {"HSS22X16X1/4", {15.650, 21.650, 0.233}},
            {"HSS22X14X3/4", {12.950, 20.950, 0.698}},
            {"HSS22X14X5/8", {13.150, 21.150, 0.581}},
            {"HSS22X14X1/2", {13.300, 21.300, 0.465}},
            {"HSS22X14X3/8", {13.500, 21.500, 0.349}},
            {"HSS22X14X5/16", {13.550, 21.550, 0.291}},
            {"HSS22X14X1/4", {13.650, 21.650, 0.233}},
            {"HSS22X10X5/8", {9.130, 21.150, 0.581}},
            {"HSS22X10X1/2", {9.305, 21.300, 0.465}},
            {"HSS22X10X3/8", {9.475, 21.500, 0.349}},
            {"HSS22X10X5/16", {9.565, 21.550, 0.291}},
            {"HSS22X10X1/4", {9.650, 21.650, 0.233}},
            {"HSS20X20X1", {18.600, 18.600, 0.930}},
            {"HSS20X20X7/8", {18.800, 18.800, 0.814}},
            {"HSS20X20X3/4", {18.950, 18.950, 0.698}},
            {"HSS20X20X5/8", {19.150, 19.150, 0.581}},
            {"HSS20X20X1/2", {19.300, 19.300, 0.465}},
            {"HSS20X20X3/8", {19.500, 19.500, 0.349}},
            {"HSS20X20X5/16", {19.550, 19.550, 0.291}},
            {"HSS20X16X3/4", {14.950, 18.950, 0.698}},
            {"HSS20X16X5/8", {15.150, 19.150, 0.581}},
            {"HSS20X16X1/2", {15.300, 19.300, 0.465}},
            {"HSS20X16X3/8", {15.500, 19.500, 0.349}},
            {"HSS20X16X5/16", {15.550, 19.550, 0.291}},
            {"HSS20X16X1/4", {15.650, 19.650, 0.233}},
            {"HSS20X12X1", {10.605, 18.600, 0.930}},
            {"HSS20X12X7/8", {10.780, 18.800, 0.814}},
            {"HSS20X12X3/4", {10.955, 18.950, 0.698}},
            {"HSS20X12X5/8", {11.150, 19.150, 0.581}},
            {"HSS20X12X1/2", {11.300, 19.300, 0.465}},
            {"HSS20X12X3/8", {11.500, 19.500, 0.349}},
            {"HSS20X12X5/16", {11.550, 19.550, 0.291}},
            {"HSS20X8X1", {6.605, 18.600, 0.930}},
            {"HSS20X8X7/8", {6.780, 18.800, 0.814}},
            {"HSS20X8X3/4", {6.955, 18.950, 0.698}},
            {"HSS20X8X5/8", {7.130, 19.150, 0.581}},
            {"HSS20X8X1/2", {7.300, 19.300, 0.465}},
            {"HSS20X8X3/8", {7.475, 19.500, 0.349}},
            {"HSS20X8X5/16", {7.565, 19.550, 0.291}},
            {"HSS20X6X5/8", {5.130, 19.150, 0.581}},
            {"HSS20X6X1/2", {5.305, 19.300, 0.465}},
            {"HSS20X6X3/8", {5.475, 19.500, 0.349}},
            {"HSS20X6X5/16", {5.565, 19.550, 0.291}},
            {"HSS20X6X1/4", {5.650, 19.650, 0.233}},
            {"HSS20X4X1/2", {3.300, 19.300, 0.465}},
            {"HSS20X4X3/8", {3.475, 19.500, 0.349}},
            {"HSS20X4X5/16", {3.565, 19.550, 0.291}},
            {"HSS20X4X1/4", {3.650, 19.650, 0.233}},
            {"HSS18X18X1", {16.600, 16.600, 0.930}},
            {"HSS18X18X7/8", {16.800, 16.800, 0.814}},
            {"HSS18X18X3/4", {16.950, 16.950, 0.698}},
            {"HSS18X18X5/8", {17.150, 17.150, 0.581}},
            {"HSS18X18X1/2", {17.300, 17.300, 0.465}},
            {"HSS18X18X3/8", {17.500, 17.500, 0.349}},
            {"HSS18X18X5/16", {17.550, 17.550, 0.291}},
            {"HSS18X18X1/4", {17.650, 17.650, 0.233}},
            {"HSS18X10X5/8", {9.130, 17.150, 0.581}},
            {"HSS18X10X1/2", {9.305, 17.300, 0.465}},
            {"HSS18X10X3/8", {9.475, 17.500, 0.349}},
            {"HSS18X10X5/16", {9.565, 17.550, 0.291}},
            {"HSS18X10X1/4", {9.650, 17.650, 0.233}},
            {"HSS18X8X5/8", {7.130, 17.150, 0.581}},
            {"HSS18X8X1/2", {7.305, 17.300, 0.465}},
            {"HSS18X8X3/8", {7.475, 17.500, 0.349}},
            {"HSS18X8X5/16", {7.565, 17.550, 0.291}},
            {"HSS18X8X1/4", {7.650, 17.650, 0.233}},
            {"HSS18X6X3/4", {4.955, 16.950, 0.698}},
            {"HSS18X6X5/8", {5.130, 17.150, 0.581}},
            {"HSS18X6X1/2", {5.305, 17.300, 0.465}},
            {"HSS18X6X3/8", {5.475, 17.500, 0.349}},
            {"HSS18X6X5/16", {5.565, 17.550, 0.291}},
            {"HSS18X6X1/4", {5.650, 17.650, 0.233}},
            {"HSS16X16X1", {14.600, 14.600, 0.930}},
            {"HSS16X16X7/8", {14.800, 14.800, 0.814}},
            {"HSS16X16X3/4", {14.950, 14.950, 0.698}},
            {"HSS16X16X5/8", {15.150, 15.150, 0.581}},
            {"HSS16X16X1/2", {15.300, 15.300, 0.465}},
            {"HSS16X16X3/8", {15.500, 15.500, 0.349}},
            {"HSS16X16X5/16", {15.550, 15.550, 0.291}},
            {"HSS16X16X1/4", {15.650, 15.650, 0.233}},
            {"HSS16X12X1", {10.605, 14.600, 0.930}},
            {"HSS16X12X7/8", {10.780, 14.800, 0.814}},
            {"HSS16X12X3/4", {10.955, 14.950, 0.698}},
            {"HSS16X12X5/8", {11.150, 15.150, 0.581}},
            {"HSS16X12X1/2", {11.300, 15.300, 0.465}},
            {"HSS16X12X3/8", {11.500, 15.500, 0.349}},
            {"HSS16X12X5/16", {11.550, 15.550, 0.291}},
            {"HSS16X10X5/8", {9.130, 15.150, 0.581}},
            {"HSS16X10X1/2", {9.305, 15.300, 0.465}},
            {"HSS16X10X3/8", {9.475, 15.500, 0.349}},
            {"HSS16X10X5/16", {9.565, 15.550, 0.291}},
            {"HSS16X10X1/4", {9.650, 15.650, 0.233}},
            {"HSS16X8X7/8", {6.780, 14.800, 0.814}},
            {"HSS16X8X3/4", {6.955, 14.950, 0.698}},
            {"HSS16X8X5/8", {7.130, 15.150, 0.581}},
            {"HSS16X8X1/2", {7.300, 15.300, 0.465}},
            {"HSS16X8X3/8", {7.475, 15.500, 0.349}},
            {"HSS16X8X5/16", {7.565, 15.550, 0.291}},
            {"HSS16X8X1/4", {7.650, 15.650, 0.233}},
            {"HSS16X6X5/8", {5.130, 15.150, 0.581}},
            {"HSS16X6X1/2", {5.305, 15.300, 0.465}},
            {"HSS16X6X3/8", {5.475, 15.500, 0.349}},
            {"HSS16X6X5/16", {5.565, 15.550, 0.291}},
            {"HSS16X6X1/4", {5.650, 15.650, 0.233}},
            {"HSS16X6X3/16", {5.740, 15.750, 0.174}},
            {"HSS16X4X5/8", {3.130, 15.150, 0.581}},
            {"HSS16X4X1/2", {3.300, 15.300, 0.465}},
            {"HSS16X4X3/8", {3.475, 15.500, 0.349}},
            {"HSS16X4X5/16", {3.565, 15.550, 0.291}},
            {"HSS16X4X1/4", {3.650, 15.650, 0.233}},
            {"HSS16X4X3/16", {3.740, 15.750, 0.174}},
            {"HSS14X14X1", {12.600, 12.600, 0.930}},
            {"HSS14X14X7/8", {12.800, 12.800, 0.814}},
            {"HSS14X14X3/4", {12.950, 12.950, 0.698}},
            {"HSS14X14X5/8", {13.150, 13.150, 0.581}},
            {"HSS14X14X1/2", {13.300, 13.300, 0.465}},
            {"HSS14X14X3/8", {13.500, 13.500, 0.349}},
            {"HSS14X14X5/16", {13.550, 13.550, 0.291}},
            {"HSS14X14X1/4", {13.650, 13.650, 0.233}},
            {"HSS14X12X5/8", {11.150, 13.150, 0.581}},
            {"HSS14X12X1/2", {11.300, 13.300, 0.465}},
            {"HSS14X12X3/8", {11.500, 13.500, 0.349}},
            {"HSS14X12X5/16", {11.550, 13.550, 0.291}},
            {"HSS14X12X1/4", {11.650, 13.650, 0.233}},
            {"HSS14X10X7/8", {8.780, 12.800, 0.814}},
            {"HSS14X10X3/4", {8.955, 12.950, 0.698}},
            {"HSS14X10X5/8", {9.130, 13.150, 0.581}},
            {"HSS14X10X1/2", {9.300, 13.300, 0.465}},
            {"HSS14X10X3/8", {9.475, 13.500, 0.349}},
            {"HSS14X10X5/16", {9.565, 13.550, 0.291}},
            {"HSS14X10X1/4", {9.650, 13.650, 0.233}},
            {"HSS14X8X5/8", {7.130, 13.150, 0.581}},
            {"HSS14X8X1/2", {7.305, 13.300, 0.465}},
            {"HSS14X8X3/8", {7.475, 13.500, 0.349}},
            {"HSS14X8X5/16", {7.565, 13.550, 0.291}},
            {"HSS14X8X1/4", {7.650, 13.650, 0.233}},
            {"HSS14X8X3/16", {7.740, 13.750, 0.174}},
            {"HSS14X6X5/8", {5.130, 13.150, 0.581}},
            {"HSS14X6X1/2", {5.305, 13.300, 0.465}},
            {"HSS14X6X3/8", {5.475, 13.500, 0.349}},
            {"HSS14X6X5/16", {5.565, 13.550, 0.291}},
            {"HSS14X6X1/4", {5.650, 13.650, 0.233}},
            {"HSS14X6X3/16", {5.740, 13.750, 0.174}},
            {"HSS14X4X5/8", {3.130, 13.150, 0.581}},
            {"HSS14X4X1/2", {3.300, 13.300, 0.465}},
            {"HSS14X4X3/8", {3.475, 13.500, 0.349}},
            {"HSS14X4X5/16", {3.565, 13.550, 0.291}},
            {"HSS14X4X1/4", {3.650, 13.650, 0.233}},
            {"HSS14X4X3/16", {3.740, 13.750, 0.174}},
            {"HSS12X12X1", {10.605, 10.605, 0.930}},
            {"HSS12X12X7/8", {10.780, 10.780, 0.814}},
            {"HSS12X12X3/4", {10.955, 10.955, 0.698}},
            {"HSS12X12X5/8", {11.150, 11.150, 0.581}},
            {"HSS12X12X1/2", {11.300, 11.300, 0.465}},
            {"HSS12X12X3/8", {11.500, 11.500, 0.349}},
            {"HSS12X12X5/16", {11.550, 11.550, 0.291}},
            {"HSS12X12X1/4", {11.650, 11.650, 0.233}},
            {"HSS12X12X3/16", {11.750, 11.750, 0.174}},
            {"HSS12X10X5/8", {9.130, 11.150, 0.581}},
            {"HSS12X10X1/2", {9.300, 11.300, 0.465}},
            {"HSS12X10X3/8", {9.475, 11.500, 0.349}},
            {"HSS12X10X5/16", {9.565, 11.550, 0.291}},
            {"HSS12X10X1/4", {9.650, 11.650, 0.233}},
            {"HSS12X10X3/16", {9.740, 11.750, 0.174}},
            {"HSS12X8X5/8", {7.130, 11.150, 0.581}},
            {"HSS12X8X1/2", {7.300, 11.300, 0.465}},
            {"HSS12X8X3/8", {7.475, 11.500, 0.349}},
            {"HSS12X8X5/16", {7.565, 11.550, 0.291}},
            {"HSS12X8X1/4", {7.650, 11.650, 0.233}},
            {"HSS12X8X3/16", {7.740, 11.750, 0.174}},
            {"HSS12X6X5/8", {5.130, 11.150, 0.581}},
            {"HSS12X6X1/2", {5.305, 11.300, 0.465}},
            {"HSS12X6X3/8", {5.475, 11.500, 0.349}},
            {"HSS12X6X5/16", {5.565, 11.550, 0.291}},
            {"HSS12X6X1/4", {5.650, 11.650, 0.233}},
            {"HSS12X6X3/16", {5.740, 11.750, 0.174}},
            {"HSS12X4X5/8", {3.130, 11.150, 0.581}},
            {"HSS12X4X1/2", {3.300, 11.300, 0.465}},
            {"HSS12X4X3/8", {3.475, 11.500, 0.349}},
            {"HSS12X4X5/16", {3.565, 11.550, 0.291}},
            {"HSS12X4X1/4", {3.650, 11.650, 0.233}},
            {"HSS12X4X3/16", {3.740, 11.750, 0.174}},
            {"HSS12X3X5/16", {2.565, 11.550, 0.291}},
            {"HSS12X3X1/4", {2.650, 11.650, 0.233}},
            {"HSS12X3X3/16", {2.740, 11.750, 0.174}},
            {"HSS12X2X5/16", {1.565, 11.550, 0.291}},
            {"HSS12X2X1/4", {1.650, 11.650, 0.233}},
            {"HSS12X2X3/16", {1.740, 11.750, 0.174}},
            {"HSS10X10X3/4", {8.955, 8.955, 0.698}},
            {"HSS10X10X5/8", {9.130, 9.130, 0.581}},
            {"HSS10X10X1/2", {9.300, 9.300, 0.465}},
            {"HSS10X10X3/8", {9.475, 9.475, 0.349}},
            {"HSS10X10X5/16", {9.565, 9.565, 0.291}},
            {"HSS10X10X1/4", {9.650, 9.650, 0.233}},
            {"HSS10X10X3/16", {9.740, 9.740, 0.174}},
            {"HSS10X8X5/8", {7.130, 9.130, 0.581}},
            {"HSS10X8X1/2", {7.300, 9.300, 0.465}},
            {"HSS10X8X3/8", {7.475, 9.475, 0.349}},
            {"HSS10X8X5/16", {7.565, 9.565, 0.291}},
            {"HSS10X8X1/4", {7.650, 9.650, 0.233}},
            {"HSS10X8X3/16", {7.740, 9.740, 0.174}},
            {"HSS10X6X5/8", {5.130, 9.130, 0.581}},
            {"HSS10X6X1/2", {5.305, 9.300, 0.465}},
            {"HSS10X6X3/8", {5.475, 9.475, 0.349}},
            {"HSS10X6X5/16", {5.565, 9.565, 0.291}},
            {"HSS10X6X1/4", {5.650, 9.650, 0.233}},
            {"HSS10X6X3/16", {5.740, 9.740, 0.174}},
            {"HSS10X5X3/8", {4.475, 9.475, 0.349}},
            {"HSS10X5X5/16", {4.565, 9.565, 0.291}},
            {"HSS10X5X1/4", {4.650, 9.650, 0.233}},
            {"HSS10X4X5/8", {3.130, 9.130, 0.581}},
            {"HSS10X4X1/2", {3.300, 9.300, 0.465}},
            {"HSS10X4X3/8", {3.475, 9.475, 0.349}},
            {"HSS10X4X5/16", {3.565, 9.565, 0.291}},
            {"HSS10X4X1/4", {3.650, 9.650, 0.233}},
            {"HSS10X4X3/16", {3.740, 9.740, 0.174}},
            {"HSS10X4X1/8", {3.825, 9.825, 0.116}},
            {"HSS10X3-1/2X3/8", {2.975, 9.475, 0.349}},
            {"HSS10X3-1/2X5/16", {3.065, 9.565, 0.291}},
            {"HSS10X3-1/2X1/4", {3.150, 9.650, 0.233}},
            {"HSS10X3-1/2X3/16", {3.240, 9.740, 0.174}},
            {"HSS10X3X3/8", {2.475, 9.475, 0.349}},
            {"HSS10X3X5/16", {2.565, 9.565, 0.291}},
            {"HSS10X3X1/4", {2.650, 9.650, 0.233}},
            {"HSS10X3X3/16", {2.740, 9.740, 0.174}},
            {"HSS10X3X1/8", {2.825, 9.825, 0.116}},
            {"HSS10X2X3/8", {1.476, 9.475, 0.349}},
            {"HSS10X2X5/16", {1.565, 9.565, 0.291}},
            {"HSS10X2X1/4", {1.650, 9.650, 0.233}},
            {"HSS10X2X3/16", {1.740, 9.740, 0.174}},
            {"HSS10X2X1/8", {1.825, 9.825, 0.116}},
            {"HSS9X9X5/8", {8.130, 8.130, 0.581}},
            {"HSS9X9X1/2", {8.300, 8.300, 0.465}},
            {"HSS9X9X3/8", {8.475, 8.475, 0.349}},
            {"HSS9X9X5/16", {8.565, 8.565, 0.291}},
            {"HSS9X9X1/4", {8.650, 8.650, 0.233}},
            {"HSS9X9X3/16", {8.740, 8.740, 0.174}},
            {"HSS9X9X1/8", {8.825, 8.825, 0.116}},
            {"HSS9X7X5/8", {6.130, 8.130, 0.581}},
            {"HSS9X7X1/2", {6.300, 8.300, 0.465}},
            {"HSS9X7X3/8", {6.475, 8.475, 0.349}},
            {"HSS9X7X5/16", {6.565, 8.565, 0.291}},
            {"HSS9X7X1/4", {6.650, 8.650, 0.233}},
            {"HSS9X7X3/16", {6.740, 8.740, 0.174}},
            {"HSS9X5X5/8", {4.130, 8.130, 0.581}},
            {"HSS9X5X1/2", {4.300, 8.300, 0.465}},
            {"HSS9X5X3/8", {4.475, 8.475, 0.349}},
            {"HSS9X5X5/16", {4.565, 8.565, 0.291}},
            {"HSS9X5X1/4", {4.650, 8.650, 0.233}},
            {"HSS9X5X3/16", {4.740, 8.740, 0.174}},
            {"HSS9X3X1/2", {2.300, 8.300, 0.465}},
            {"HSS9X3X3/8", {2.475, 8.475, 0.349}},
            {"HSS9X3X5/16", {2.565, 8.565, 0.291}},
            {"HSS9X3X1/4", {2.650, 8.650, 0.233}},
            {"HSS9X3X3/16", {2.740, 8.740, 0.174}},
            {"HSS8X8X5/8", {7.130, 7.130, 0.581}},
            {"HSS8X8X1/2", {7.300, 7.300, 0.465}},
            {"HSS8X8X3/8", {7.475, 7.475, 0.349}},
            {"HSS8X8X5/16", {7.565, 7.565, 0.291}},
            {"HSS8X8X1/4", {7.650, 7.650, 0.233}},
            {"HSS8X8X3/16", {7.740, 7.740, 0.174}},
            {"HSS8X8X1/8", {7.825, 7.825, 0.116}},
            {"HSS8X6X5/8", {5.130, 7.130, 0.581}},
            {"HSS8X6X1/2", {5.305, 7.300, 0.465}},
            {"HSS8X6X3/8", {5.475, 7.475, 0.349}},
            {"HSS8X6X5/16", {5.565, 7.565, 0.291}},
            {"HSS8X6X1/4", {5.650, 7.650, 0.233}},
            {"HSS8X6X3/16", {5.740, 7.740, 0.174}},
            {"HSS8X4X5/8", {3.130, 7.130, 0.581}},
            {"HSS8X4X1/2", {3.300, 7.300, 0.465}},
            {"HSS8X4X3/8", {3.475, 7.475, 0.349}},
            {"HSS8X4X5/16", {3.565, 7.565, 0.291}},
            {"HSS8X4X1/4", {3.650, 7.650, 0.233}},
            {"HSS8X4X3/16", {3.740, 7.740, 0.174}},
            {"HSS8X4X1/8", {3.825, 7.825, 0.116}},
            {"HSS8X3X1/2", {2.300, 7.300, 0.465}},
            {"HSS8X3X3/8", {2.475, 7.475, 0.349}},
            {"HSS8X3X5/16", {2.565, 7.565, 0.291}},
            {"HSS8X3X1/4", {2.650, 7.650, 0.233}},
            {"HSS8X3X3/16", {2.740, 7.740, 0.174}},
            {"HSS8X3X1/8", {2.825, 7.825, 0.116}},
            {"HSS8X2X1/2", {1.302, 7.305, 0.465}},
            {"HSS8X2X3/8", {1.476, 7.475, 0.349}},
            {"HSS8X2X5/16", {1.565, 7.565, 0.291}},
            {"HSS8X2X1/4", {1.650, 7.650, 0.233}},
            {"HSS8X2X3/16", {1.740, 7.740, 0.174}},
            {"HSS8X2X1/8", {1.825, 7.825, 0.116}},
            {"HSS7X7X5/8", {6.130, 6.130, 0.581}},
            {"HSS7X7X1/2", {6.300, 6.300, 0.465}},
            {"HSS7X7X3/8", {6.475, 6.475, 0.349}},
            {"HSS7X7X5/16", {6.565, 6.565, 0.291}},
            {"HSS7X7X1/4", {6.650, 6.650, 0.233}},
            {"HSS7X7X3/16", {6.740, 6.740, 0.174}},
            {"HSS7X7X1/8", {6.825, 6.825, 0.116}},
            {"HSS7X5X1/2", {4.300, 6.300, 0.465}},
            {"HSS7X5X3/8", {4.475, 6.475, 0.349}},
            {"HSS7X5X5/16", {4.565, 6.565, 0.291}},
            {"HSS7X5X1/4", {4.650, 6.650, 0.233}},
            {"HSS7X5X3/16", {4.740, 6.740, 0.174}},
            {"HSS7X5X1/8", {4.825, 6.825, 0.116}},
            {"HSS7X4X1/2", {3.300, 6.300, 0.465}},
            {"HSS7X4X3/8", {3.475, 6.475, 0.349}},
            {"HSS7X4X5/16", {3.565, 6.565, 0.291}},
            {"HSS7X4X1/4", {3.650, 6.650, 0.233}},
            {"HSS7X4X3/16", {3.740, 6.740, 0.174}},
            {"HSS7X4X1/8", {3.825, 6.825, 0.116}},
            {"HSS7X3X1/2", {2.300, 6.300, 0.465}},
            {"HSS7X3X3/8", {2.475, 6.475, 0.349}},
            {"HSS7X3X5/16", {2.565, 6.565, 0.291}},
            {"HSS7X3X1/4", {2.650, 6.650, 0.233}},
            {"HSS7X3X3/16", {2.740, 6.740, 0.174}},
            {"HSS7X3X1/8", {2.825, 6.825, 0.116}},
            {"HSS7X2X1/4", {1.650, 6.650, 0.233}},
            {"HSS7X2X3/16", {1.740, 6.740, 0.174}},
            {"HSS7X2X1/8", {1.825, 6.825, 0.116}},
            {"HSS6X6X5/8", {5.130, 5.130, 0.581}},
            {"HSS6X6X1/2", {5.305, 5.305, 0.465}},
            {"HSS6X6X3/8", {5.475, 5.475, 0.349}},
            {"HSS6X6X5/16", {5.565, 5.565, 0.291}},
            {"HSS6X6X1/4", {5.650, 5.650, 0.233}},
            {"HSS6X6X3/16", {5.740, 5.740, 0.174}},
            {"HSS6X6X1/8", {5.825, 5.825, 0.116}},
            {"HSS6X5X1/2", {4.300, 5.305, 0.465}},
            {"HSS6X5X3/8", {4.475, 5.475, 0.349}},
            {"HSS6X5X5/16", {4.565, 5.565, 0.291}},
            {"HSS6X5X1/4", {4.650, 5.650, 0.233}},
            {"HSS6X5X3/16", {4.740, 5.740, 0.174}},
            {"HSS6X5X1/8", {4.825, 5.825, 0.116}},
            {"HSS6X4X1/2", {3.300, 5.305, 0.465}},
            {"HSS6X4X3/8", {3.475, 5.475, 0.349}},
            {"HSS6X4X5/16", {3.565, 5.565, 0.291}},
            {"HSS6X4X1/4", {3.650, 5.650, 0.233}},
            {"HSS6X4X3/16", {3.740, 5.740, 0.174}},
            {"HSS6X4X1/8", {3.825, 5.825, 0.116}},
            {"HSS6X3X1/2", {2.300, 5.305, 0.465}},
            {"HSS6X3X3/8", {2.475, 5.475, 0.349}},
            {"HSS6X3X5/16", {2.565, 5.565, 0.291}},
            {"HSS6X3X1/4", {2.650, 5.650, 0.233}},
            {"HSS6X3X3/16", {2.740, 5.740, 0.174}},
            {"HSS6X3X1/8", {2.825, 5.825, 0.116}},
            {"HSS6X2X3/8", {1.476, 5.475, 0.349}},
            {"HSS6X2X5/16", {1.565, 5.565, 0.291}},
            {"HSS6X2X1/4", {1.650, 5.650, 0.233}},
            {"HSS6X2X3/16", {1.740, 5.740, 0.174}},
            {"HSS6X2X1/8", {1.825, 5.825, 0.116}},
            {"HSS5-1/2X5-1/2X3/8", {4.975, 4.975, 0.349}},
            {"HSS5-1/2X5-1/2X5/16", {5.065, 5.065, 0.291}},
            {"HSS5-1/2X5-1/2X1/4", {5.150, 5.150, 0.233}},
            {"HSS5-1/2X5-1/2X3/16", {5.240, 5.240, 0.174}},
            {"HSS5-1/2X5-1/2X1/8", {5.325, 5.325, 0.116}},
            {"HSS5X5X1/2", {4.300, 4.300, 0.465}},
            {"HSS5X5X3/8", {4.475, 4.475, 0.349}},
            {"HSS5X5X5/16", {4.565, 4.565, 0.291}},
            {"HSS5X5X1/4", {4.650, 4.650, 0.233}},
            {"HSS5X5X3/16", {4.740, 4.740, 0.174}},
            {"HSS5X5X1/8", {4.825, 4.825, 0.116}},
            {"HSS5X4X1/2", {3.300, 4.300, 0.465}},
            {"HSS5X4X3/8", {3.475, 4.475, 0.349}},
            {"HSS5X4X5/16", {3.565, 4.565, 0.291}},
            {"HSS5X4X1/4", {3.650, 4.650, 0.233}},
            {"HSS5X4X3/16", {3.740, 4.740, 0.174}},
            {"HSS5X4X1/8", {3.825, 4.825, 0.116}},
            {"HSS5X3X1/2", {2.300, 4.300, 0.465}},
            {"HSS5X3X3/8", {2.475, 4.475, 0.349}},
            {"HSS5X3X5/16", {2.565, 4.565, 0.291}},
            {"HSS5X3X1/4", {2.650, 4.650, 0.233}},
            {"HSS5X3X3/16", {2.740, 4.740, 0.174}},
            {"HSS5X3X1/8", {2.825, 4.825, 0.116}},
            {"HSS5X2-1/2X1/4", {2.150, 4.650, 0.233}},
            {"HSS5X2-1/2X3/16", {2.240, 4.740, 0.174}},
            {"HSS5X2-1/2X1/8", {2.325, 4.825, 0.116}},
            {"HSS5X2X3/8", {1.476, 4.475, 0.349}},
            {"HSS5X2X5/16", {1.565, 4.565, 0.291}},
            {"HSS5X2X1/4", {1.650, 4.650, 0.233}},
            {"HSS5X2X3/16", {1.740, 4.740, 0.174}},
            {"HSS5X2X1/8", {1.825, 4.825, 0.116}},
            {"HSS4-1/2X4-1/2X1/2", {3.800, 3.800, 0.465}},
            {"HSS4-1/2X4-1/2X3/8", {3.975, 3.975, 0.349}},
            {"HSS4-1/2X4-1/2X5/16", {4.065, 4.065, 0.291}},
            {"HSS4-1/2X4-1/2X1/4", {4.150, 4.150, 0.233}},
            {"HSS4-1/2X4-1/2X3/16", {4.240, 4.240, 0.174}},
            {"HSS4-1/2X4-1/2X1/8", {4.325, 4.325, 0.116}},
            {"HSS4X4X1/2", {3.300, 3.300, 0.465}},
            {"HSS4X4X3/8", {3.475, 3.475, 0.349}},
            {"HSS4X4X5/16", {3.565, 3.565, 0.291}},
            {"HSS4X4X1/4", {3.650, 3.650, 0.233}},
            {"HSS4X4X3/16", {3.740, 3.740, 0.174}},
            {"HSS4X4X1/8", {3.825, 3.825, 0.116}},
            {"HSS4X3X3/8", {2.475, 3.475, 0.349}},
            {"HSS4X3X5/16", {2.565, 3.565, 0.291}},
            {"HSS4X3X1/4", {2.650, 3.650, 0.233}},
            {"HSS4X3X3/16", {2.740, 3.740, 0.174}},
            {"HSS4X3X1/8", {2.825, 3.825, 0.116}},
            {"HSS4X2-1/2X1/4", {2.150, 3.650, 0.233}},
            {"HSS4X2-1/2X3/16", {2.240, 3.740, 0.174}},
            {"HSS4X2-1/2X1/8", {2.325, 3.825, 0.116}},
            {"HSS4X2X3/8", {1.476, 3.475, 0.349}},
            {"HSS4X2X5/16", {1.565, 3.565, 0.291}},
            {"HSS4X2X1/4", {1.650, 3.650, 0.233}},
            {"HSS4X2X3/16", {1.740, 3.740, 0.174}},
            {"HSS4X2X1/8", {1.825, 3.825, 0.116}},
            {"HSS4X1-1/2X1/4", {1.151, 3.650, 0.233}},
            {"HSS4X1-1/2X3/16", {1.239, 3.740, 0.174}},
            {"HSS4X1-1/2X1/8", {1.325, 3.825, 0.116}},
            {"HSS3-1/2X3-1/2X3/8", {2.975, 2.975, 0.349}},
            {"HSS3-1/2X3-1/2X5/16", {3.065, 3.065, 0.291}},
            {"HSS3-1/2X3-1/2X1/4", {3.150, 3.150, 0.233}},
            {"HSS3-1/2X3-1/2X3/16", {3.240, 3.240, 0.174}},
            {"HSS3-1/2X3-1/2X1/8", {3.325, 3.325, 0.116}},
            {"HSS3-1/2X2-1/2X3/8", {1.975, 2.975, 0.349}},
            {"HSS3-1/2X2-1/2X5/16", {2.065, 3.065, 0.291}},
            {"HSS3-1/2X2-1/2X1/4", {2.150, 3.150, 0.233}},
            {"HSS3-1/2X2-1/2X3/16", {2.240, 3.240, 0.174}},
            {"HSS3-1/2X2-1/2X1/8", {2.325, 3.325, 0.116}},
            {"HSS3-1/2X2X1/4", {1.650, 3.150, 0.233}},
            {"HSS3-1/2X2X3/16", {1.740, 3.240, 0.174}},
            {"HSS3-1/2X2X1/8", {1.825, 3.325, 0.116}},
            {"HSS3-1/2X1-1/2X1/4", {1.151, 3.150, 0.233}},
            {"HSS3-1/2X1-1/2X3/16", {1.239, 3.240, 0.174}},
            {"HSS3-1/2X1-1/2X1/8", {1.325, 3.325, 0.116}},
            {"HSS3X3X3/8", {2.475, 2.475, 0.349}},
            {"HSS3X3X5/16", {2.565, 2.565, 0.291}},
            {"HSS3X3X1/4", {2.650, 2.650, 0.233}},
            {"HSS3X3X3/16", {2.740, 2.740, 0.174}},
            {"HSS3X3X1/8", {2.825, 2.825, 0.116}},
            {"HSS3X2-1/2X5/16", {2.065, 2.565, 0.291}},
            {"HSS3X2-1/2X1/4", {2.150, 2.650, 0.233}},
            {"HSS3X2-1/2X3/16", {2.240, 2.740, 0.174}},
            {"HSS3X2-1/2X1/8", {2.325, 2.825, 0.116}},
            {"HSS3X2X5/16", {1.565, 2.565, 0.291}},
            {"HSS3X2X1/4", {1.650, 2.650, 0.233}},
            {"HSS3X2X3/16", {1.740, 2.740, 0.174}},
            {"HSS3X2X1/8", {1.825, 2.825, 0.116}},
            {"HSS3X1-1/2X1/4", {1.151, 2.650, 0.233}},
            {"HSS3X1-1/2X3/16", {1.239, 2.740, 0.174}},
            {"HSS3X1-1/2X1/8", {1.325, 2.825, 0.116}},
            {"HSS3X1X3/16", {0.739, 2.740, 0.174}},
            {"HSS3X1X1/8", {0.826, 2.825, 0.116}},
            {"HSS2-1/2X2-1/2X5/16", {2.065, 2.065, 0.291}},
            {"HSS2-1/2X2-1/2X1/4", {2.150, 2.150, 0.233}},
            {"HSS2-1/2X2-1/2X3/16", {2.240, 2.240, 0.174}},
            {"HSS2-1/2X2-1/2X1/8", {2.325, 2.325, 0.116}},
            {"HSS2-1/2X2X1/4", {1.650, 2.150, 0.233}},
            {"HSS2-1/2X2X3/16", {1.740, 2.240, 0.174}},
            {"HSS2-1/2X2X1/8", {1.825, 2.325, 0.116}},
            {"HSS2-1/2X1-1/2X1/4", {1.151, 2.150, 0.233}},
            {"HSS2-1/2X1-1/2X3/16", {1.239, 2.240, 0.174}},
            {"HSS2-1/2X1-1/2X1/8", {1.325, 2.325, 0.116}},
            {"HSS2-1/2X1X3/16", {0.739, 2.240, 0.174}},
            {"HSS2-1/2X1X1/8", {0.826, 2.325, 0.116}},
            {"HSS2-1/4X2-1/4X1/4", {1.900, 1.900, 0.233}},
            {"HSS2-1/4X2-1/4X3/16", {1.990, 1.990, 0.174}},
            {"HSS2-1/4X2-1/4X1/8", {2.075, 2.075, 0.116}},
            {"HSS2X2X1/4", {1.650, 1.650, 0.233}},
            {"HSS2X2X3/16", {1.740, 1.740, 0.174}},
            {"HSS2X2X1/8", {1.825, 1.825, 0.116}},
            {"HSS2X1-1/2X3/16", {1.239, 1.740, 0.174}},
            {"HSS2X1-1/2X1/8", {1.325, 1.825, 0.116}},
            {"HSS2X1X3/16", {0.739, 1.740, 0.174}},
            {"HSS2X1X1/8", {0.826, 1.825, 0.116}},
            {"HSS1-1/2X1-1/2X1/4", {1.151, 1.151, 0.233}},
            {"HSS1-1/2X1-1/2X3/16", {1.239, 1.239, 0.174}},
            {"HSS1-1/2X1-1/2X1/8", {1.325, 1.325, 0.116}},
        };

        return registry.fetch(type);
    }

    // Adjusted registry: dimensions optimized against the AISC tabulated
    // section properties with section_optimizer.py.
    vec usrhssection_adjusted(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"HSS34X10X1", {9.0519, 32.4217, 0.9173}},
            {"HSS34X10X7/8", {9.1605, 32.6447, 0.8038}},
            {"HSS34X10X3/4", {9.2729, 32.8673, 0.6899}},
            {"HSS34X10X5/8", {9.4104, 33.0437, 0.576}},
            {"HSS30X10X1", {9.0223, 28.499, 0.9149}},
            {"HSS30X10X7/8", {9.1452, 28.6638, 0.8038}},
            {"HSS30X10X3/4", {9.2857, 28.8652, 0.6892}},
            {"HSS30X10X5/8", {9.3966, 29.0368, 0.5764}},
            {"HSS30X10X1/2", {9.5172, 29.2802, 0.4608}},
            {"HSS24X20X3/4", {19.2135, 23.1484, 0.6865}},
            {"HSS24X20X5/8", {19.3396, 23.2993, 0.5728}},
            {"HSS24X20X1/2", {19.4649, 23.4383, 0.4601}},
            {"HSS24X20X3/8", {19.584, 23.5704, 0.3462}},
            {"HSS24X20X5/16", {19.6582, 23.6487, 0.2893}},
            {"HSS24X18X3/4", {17.2197, 23.123, 0.686}},
            {"HSS24X18X5/8", {17.3557, 23.2604, 0.5733}},
            {"HSS24X18X1/2", {17.505, 23.4309, 0.4587}},
            {"HSS24X18X3/8", {17.6046, 23.5867, 0.3459}},
            {"HSS24X18X5/16", {17.652, 23.5992, 0.2896}},
            {"HSS24X16X3/4", {15.2491, 23.0929, 0.6854}},
            {"HSS24X16X5/8", {15.3567, 23.2346, 0.5733}},
            {"HSS24X16X1/2", {15.5069, 23.3698, 0.46}},
            {"HSS24X16X3/8", {15.5915, 23.5497, 0.346}},
            {"HSS24X16X5/16", {15.7229, 23.6218, 0.2887}},
            {"HSS24X14X3/4", {13.2681, 23.0483, 0.6851}},
            {"HSS24X14X5/8", {13.38, 23.2157, 0.5723}},
            {"HSS24X14X1/2", {13.5157, 23.3955, 0.4587}},
            {"HSS24X14X3/8", {13.6245, 23.5584, 0.3456}},
            {"HSS24X14X5/16", {13.6739, 23.6338, 0.2889}},
            {"HSS24X14X1/4", {13.7378, 23.6778, 0.2321}},
            {"HSS24X12X1", {11.0191, 22.6388, 0.9091}},
            {"HSS24X12X7/8", {11.1414, 22.8363, 0.7972}},
            {"HSS24X12X3/4", {11.2694, 22.991, 0.6862}},
            {"HSS24X12X5/8", {11.3716, 23.167, 0.5727}},
            {"HSS24X12X1/2", {11.5118, 23.3352, 0.4598}},
            {"HSS24X12X3/8", {11.6434, 23.5023, 0.3456}},
            {"HSS24X12X5/16", {11.6887, 23.565, 0.2893}},
            {"HSS24X12X1/4", {11.7712, 23.7245, 0.2308}},
            {"HSS24X8X1/2", {7.5227, 23.2553, 0.4599}},
            {"HSS24X8X3/8", {7.6462, 23.4412, 0.3459}},
            {"HSS24X8X5/16", {7.7012, 23.5588, 0.2889}},
            {"HSS24X8X1/4", {7.7499, 23.5655, 0.2329}},
            {"HSS22X22X1", {20.8904, 20.8904, 0.9109}},
            {"HSS22X22X7/8", {21.0406, 21.0406, 0.7987}},
            {"HSS22X22X3/4", {21.1662, 21.1662, 0.6873}},
            {"HSS22X22X5/8", {21.3372, 21.3372, 0.5719}},
            {"HSS22X22X1/2", {21.4557, 21.4557, 0.4596}},
            {"HSS22X20X3/4", {19.2023, 21.1684, 0.6856}},
            {"HSS22X20X5/8", {19.3142, 21.3, 0.5735}},
            {"HSS22X20X1/2", {19.5036, 21.4698, 0.4583}},
            {"HSS22X20X3/8", {19.6049, 21.5804, 0.346}},
            {"HSS22X20X5/16", {19.679, 21.6353, 0.2892}},
            {"HSS22X18X3/4", {17.2044, 21.1352, 0.6856}},
            {"HSS22X18X5/8", {17.3194, 21.2703, 0.5735}},
            {"HSS22X18X1/2", {17.4546, 21.4254, 0.4599}},
            {"HSS22X18X3/8", {17.5946, 21.5966, 0.3458}},
            {"HSS22X18X5/16", {17.6634, 21.6073, 0.2894}},
            {"HSS22X16X3/4", {15.2203, 21.0918, 0.6855}},
            {"HSS22X16X5/8", {15.3496, 21.2881, 0.5714}},
            {"HSS22X16X1/2", {15.522, 21.3997, 0.4583}},
            {"HSS22X16X3/8", {15.6254, 21.5334, 0.3462}},
            {"HSS22X16X5/16", {15.6739, 21.6453, 0.2886}},
            {"HSS22X16X1/4", {15.7192, 21.6363, 0.233}},
            {"HSS22X14X3/4", {13.2146, 21.0736, 0.6861}},
            {"HSS22X14X5/8", {13.3963, 21.2448, 0.5705}},
            {"HSS22X14X1/2", {13.4655, 21.3902, 0.4597}},
            {"HSS22X14X3/8", {13.6228, 21.5403, 0.3461}},
            {"HSS22X14X5/16", {13.6758, 21.6321, 0.2885}},
            {"HSS22X14X1/4", {13.7617, 21.663, 0.2313}},
            {"HSS22X10X5/8", {9.3897, 21.1527, 0.5722}},
            {"HSS22X10X1/2", {9.5245, 21.3505, 0.4579}},
            {"HSS22X10X3/8", {9.6476, 21.4794, 0.3453}},
            {"HSS22X10X5/16", {9.6839, 21.6081, 0.2891}},
            {"HSS22X10X1/4", {9.7535, 21.6312, 0.2322}},
            {"HSS20X20X1", {18.9081, 18.9081, 0.907}},
            {"HSS20X20X7/8", {19.0186, 19.0186, 0.7984}},
            {"HSS20X20X3/4", {19.1863, 19.1863, 0.6848}},
            {"HSS20X20X5/8", {19.2954, 19.2954, 0.5734}},
            {"HSS20X20X1/2", {19.4592, 19.4592, 0.4597}},
            {"HSS20X20X3/8", {19.6308, 19.6308, 0.3447}},
            {"HSS20X20X5/16", {19.6459, 19.6459, 0.2887}},
            {"HSS20X16X3/4", {15.2034, 19.1288, 0.6848}},
            {"HSS20X16X5/8", {15.326, 19.3049, 0.5711}},
            {"HSS20X16X1/2", {15.4565, 19.3937, 0.4601}},
            {"HSS20X16X3/8", {15.583, 19.6047, 0.3457}},
            {"HSS20X16X5/16", {15.6803, 19.6645, 0.2879}},
            {"HSS20X16X1/4", {15.7377, 19.7007, 0.2313}},
            {"HSS20X12X1", {10.9955, 18.7134, 0.9043}},
            {"HSS20X12X7/8", {11.1006, 18.9119, 0.7939}},
            {"HSS20X12X3/4", {11.2404, 19.0464, 0.6843}},
            {"HSS20X12X5/8", {11.3656, 19.2148, 0.5718}},
            {"HSS20X12X1/2", {11.5081, 19.4165, 0.4567}},
            {"HSS20X12X3/8", {11.6416, 19.5734, 0.3443}},
            {"HSS20X12X5/16", {11.6768, 19.583, 0.2895}},
            {"HSS20X8X1", {7.027, 18.4947, 0.9083}},
            {"HSS20X8X7/8", {7.1442, 18.7588, 0.7945}},
            {"HSS20X8X3/4", {7.2678, 18.9381, 0.6832}},
            {"HSS20X8X5/8", {7.391, 19.1218, 0.5711}},
            {"HSS20X8X1/2", {7.5129, 19.2773, 0.4592}},
            {"HSS20X8X3/8", {7.6393, 19.4881, 0.345}},
            {"HSS20X8X5/16", {7.6982, 19.5889, 0.2878}},
            {"HSS20X6X5/8", {5.4019, 19.0217, 0.5725}},
            {"HSS20X6X1/2", {5.5123, 19.2693, 0.4596}},
            {"HSS20X6X3/8", {5.6613, 19.4585, 0.3441}},
            {"HSS20X6X5/16", {5.7068, 19.5127, 0.2891}},
            {"HSS20X6X1/4", {5.7541, 19.6017, 0.2323}},
            {"HSS20X4X1/2", {3.5291, 19.1423, 0.4601}},
            {"HSS20X4X3/8", {3.6416, 19.3381, 0.3477}},
            {"HSS20X4X5/16", {3.7059, 19.4814, 0.289}},
            {"HSS20X4X1/4", {3.7691, 19.6181, 0.2308}},
            {"HSS18X18X1", {16.8899, 16.8899, 0.9055}},
            {"HSS18X18X7/8", {17.0353, 17.0353, 0.7956}},
            {"HSS18X18X3/4", {17.1899, 17.1899, 0.683}},
            {"HSS18X18X5/8", {17.3126, 17.3126, 0.5716}},
            {"HSS18X18X1/2", {17.4498, 17.4498, 0.4598}},
            {"HSS18X18X3/8", {17.566, 17.566, 0.3457}},
            {"HSS18X18X5/16", {17.6558, 17.6558, 0.2888}},
            {"HSS18X18X1/4", {17.7384, 17.7384, 0.231}},
            {"HSS18X10X5/8", {9.383, 17.231, 0.5685}},
            {"HSS18X10X1/2", {9.5091, 17.3959, 0.4563}},
            {"HSS18X10X3/8", {9.6413, 17.5386, 0.344}},
            {"HSS18X10X5/16", {9.6996, 17.6181, 0.2874}},
            {"HSS18X10X1/4", {9.7428, 17.6775, 0.2315}},
            {"HSS18X8X5/8", {7.3818, 17.1417, 0.5706}},
            {"HSS18X8X1/2", {7.5104, 17.2959, 0.4588}},
            {"HSS18X8X3/8", {7.6493, 17.5082, 0.3437}},
            {"HSS18X8X5/16", {7.6887, 17.5656, 0.2887}},
            {"HSS18X8X1/4", {7.7467, 17.6317, 0.2321}},
            {"HSS18X6X3/4", {5.2751, 16.84, 0.6845}},
            {"HSS18X6X5/8", {5.4001, 17.0398, 0.5714}},
            {"HSS18X6X1/2", {5.5195, 17.2468, 0.4589}},
            {"HSS18X6X3/8", {5.6225, 17.4096, 0.3471}},
            {"HSS18X6X5/16", {5.7053, 17.5544, 0.2878}},
            {"HSS18X6X1/4", {5.766, 17.6613, 0.2304}},
            {"HSS16X16X1", {14.8724, 14.8724, 0.9037}},
            {"HSS16X16X7/8", {15.034, 15.034, 0.7915}},
            {"HSS16X16X3/4", {15.1539, 15.1539, 0.6832}},
            {"HSS16X16X5/8", {15.3259, 15.3259, 0.5692}},
            {"HSS16X16X1/2", {15.4746, 15.4746, 0.4568}},
            {"HSS16X16X3/8", {15.6019, 15.6019, 0.3448}},
            {"HSS16X16X5/16", {15.6537, 15.6537, 0.2887}},
            {"HSS16X16X1/4", {15.7053, 15.7053, 0.2321}},
            {"HSS16X12X1", {10.9429, 14.7799, 0.901}},
            {"HSS16X12X7/8", {11.0874, 14.9899, 0.7876}},
            {"HSS16X12X3/4", {11.2128, 15.1318, 0.6801}},
            {"HSS16X12X5/8", {11.3492, 15.2577, 0.569}},
            {"HSS16X12X1/2", {11.4882, 15.417, 0.4567}},
            {"HSS16X12X3/8", {11.6226, 15.5807, 0.3437}},
            {"HSS16X12X5/16", {11.6959, 15.6528, 0.287}},
            {"HSS16X10X5/8", {9.3624, 15.2161, 0.5692}},
            {"HSS16X10X1/2", {9.4843, 15.3743, 0.4582}},
            {"HSS16X10X3/8", {9.6279, 15.555, 0.3435}},
            {"HSS16X10X5/16", {9.6848, 15.605, 0.2883}},
            {"HSS16X10X1/4", {9.7341, 15.6626, 0.232}},
            {"HSS16X8X7/8", {7.129, 14.7971, 0.7912}},
            {"HSS16X8X3/4", {7.261, 14.9822, 0.6798}},
            {"HSS16X8X5/8", {7.3726, 15.1522, 0.5699}},
            {"HSS16X8X1/2", {7.5033, 15.3269, 0.4571}},
            {"HSS16X8X3/8", {7.6136, 15.4896, 0.3458}},
            {"HSS16X8X5/16", {7.7016, 15.5923, 0.2874}},
            {"HSS16X8X1/4", {7.7665, 15.681, 0.2303}},
            {"HSS16X6X5/8", {5.402, 15.0926, 0.5681}},
            {"HSS16X6X1/2", {5.5311, 15.285, 0.4558}},
            {"HSS16X6X3/8", {5.6335, 15.4385, 0.3457}},
            {"HSS16X6X5/16", {5.7099, 15.5693, 0.2868}},
            {"HSS16X6X1/4", {5.7612, 15.6477, 0.2311}},
            {"HSS16X6X3/16", {5.8205, 15.7268, 0.1731}},
            {"HSS16X4X5/8", {3.4061, 14.9685, 0.5707}},
            {"HSS16X4X1/2", {3.5216, 15.1759, 0.4592}},
            {"HSS16X4X3/8", {3.6413, 15.3782, 0.3464}},
            {"HSS16X4X5/16", {3.7012, 15.5047, 0.2884}},
            {"HSS16X4X1/4", {3.7664, 15.6236, 0.2308}},
            {"HSS16X4X3/16", {3.824, 15.676, 0.1735}},
            {"HSS14X14X1", {12.8883, 12.8883, 0.8988}},
            {"HSS14X14X7/8", {13.0341, 13.0341, 0.789}},
            {"HSS14X14X3/4", {13.1704, 13.1704, 0.6806}},
            {"HSS14X14X5/8", {13.3201, 13.3201, 0.5679}},
            {"HSS14X14X1/2", {13.4564, 13.4564, 0.4567}},
            {"HSS14X14X3/8", {13.6021, 13.6021, 0.3437}},
            {"HSS14X14X5/16", {13.6828, 13.6828, 0.2867}},
            {"HSS14X14X1/4", {13.7302, 13.7302, 0.2311}},
            {"HSS14X12X5/8", {11.3402, 13.2795, 0.5674}},
            {"HSS14X12X1/2", {11.4656, 13.4275, 0.4564}},
            {"HSS14X12X3/8", {11.6222, 13.5881, 0.3432}},
            {"HSS14X12X5/16", {11.6639, 13.6394, 0.2882}},
            {"HSS14X12X1/4", {11.7053, 13.6898, 0.2321}},
            {"HSS14X10X7/8", {9.0883, 12.9302, 0.7865}},
            {"HSS14X10X3/4", {9.2165, 13.0803, 0.6793}},
            {"HSS14X10X5/8", {9.3415, 13.2371, 0.5684}},
            {"HSS14X10X1/2", {9.4918, 13.4052, 0.4558}},
            {"HSS14X10X3/8", {9.6071, 13.5342, 0.345}},
            {"HSS14X10X5/16", {9.677, 13.6307, 0.2874}},
            {"HSS14X10X1/4", {9.7618, 13.7074, 0.2301}},
            {"HSS14X8X5/8", {7.3784, 13.1961, 0.5653}},
            {"HSS14X8X1/2", {7.5081, 13.3649, 0.4548}},
            {"HSS14X8X3/8", {7.6195, 13.5016, 0.3451}},
            {"HSS14X8X5/16", {7.7171, 13.6136, 0.2861}},
            {"HSS14X8X1/4", {7.7571, 13.6762, 0.2307}},
            {"HSS14X8X3/16", {7.8149, 13.7662, 0.1726}},
            {"HSS14X6X5/8", {5.398, 13.1188, 0.5662}},
            {"HSS14X6X1/2", {5.4995, 13.2859, 0.4573}},
            {"HSS14X6X3/8", {5.6283, 13.4609, 0.3451}},
            {"HSS14X6X5/16", {5.6951, 13.5648, 0.2878}},
            {"HSS14X6X1/4", {5.7578, 13.6539, 0.2308}},
            {"HSS14X6X3/16", {5.8197, 13.7614, 0.1724}},
            {"HSS14X4X5/8", {3.3988, 12.9861, 0.5697}},
            {"HSS14X4X1/2", {3.5254, 13.2165, 0.4561}},
            {"HSS14X4X3/8", {3.6366, 13.4008, 0.3456}},
            {"HSS14X4X5/16", {3.7032, 13.5084, 0.2879}},
            {"HSS14X4X1/4", {3.7666, 13.6159, 0.2309}},
            {"HSS14X4X3/16", {3.8294, 13.7524, 0.172}},
            {"HSS12X12X1", {10.8719, 10.8719, 0.8931}},
            {"HSS12X12X7/8", {11.0273, 11.0273, 0.7847}},
            {"HSS12X12X3/4", {11.1567, 11.1567, 0.6792}},
            {"HSS12X12X5/8", {11.2961, 11.2961, 0.5689}},
            {"HSS12X12X1/2", {11.4476, 11.4476, 0.456}},
            {"HSS12X12X3/8", {11.5701, 11.5701, 0.3451}},
            {"HSS12X12X5/16", {11.6641, 11.6641, 0.2871}},
            {"HSS12X12X1/4", {11.7334, 11.7334, 0.2303}},
            {"HSS12X12X3/16", {11.7968, 11.7968, 0.1726}},
            {"HSS12X10X5/8", {9.3456, 11.2801, 0.5641}},
            {"HSS12X10X1/2", {9.4743, 11.4299, 0.4542}},
            {"HSS12X10X3/8", {9.5846, 11.5566, 0.3448}},
            {"HSS12X10X5/16", {9.6936, 11.6636, 0.2855}},
            {"HSS12X10X1/4", {9.7477, 11.7151, 0.2304}},
            {"HSS12X10X3/16", {9.8055, 11.7975, 0.1724}},
            {"HSS12X8X5/8", {7.3576, 11.2386, 0.5637}},
            {"HSS12X8X1/2", {7.4964, 11.3722, 0.4549}},
            {"HSS12X8X3/8", {7.5927, 11.5162, 0.345}},
            {"HSS12X8X5/16", {7.6713, 11.607, 0.2876}},
            {"HSS12X8X1/4", {7.7492, 11.7086, 0.2301}},
            {"HSS12X8X3/16", {7.8104, 11.7574, 0.1727}},
            {"HSS12X6X5/8", {5.391, 11.1511, 0.5638}},
            {"HSS12X6X1/2", {5.5076, 11.325, 0.4539}},
            {"HSS12X6X3/8", {5.6225, 11.4773, 0.3445}},
            {"HSS12X6X5/16", {5.6942, 11.5748, 0.2871}},
            {"HSS12X6X1/4", {5.7566, 11.65, 0.2306}},
            {"HSS12X6X3/16", {5.8191, 11.7556, 0.1723}},
            {"HSS12X4X5/8", {3.3946, 11.0201, 0.567}},
            {"HSS12X4X1/2", {3.5121, 11.2261, 0.4572}},
            {"HSS12X4X3/8", {3.6328, 11.4262, 0.3447}},
            {"HSS12X4X5/16", {3.7, 11.5199, 0.2874}},
            {"HSS12X4X1/4", {3.7568, 11.6202, 0.2306}},
            {"HSS12X4X3/16", {3.8203, 11.7251, 0.1727}},
            {"HSS12X3X5/16", {2.7029, 11.4757, 0.2881}},
            {"HSS12X3X1/4", {2.7673, 11.5995, 0.2305}},
            {"HSS12X3X3/16", {2.8238, 11.7013, 0.1727}},
            {"HSS12X2X5/16", {1.7039, 11.4202, 0.2891}},
            {"HSS12X2X1/4", {1.7651, 11.5545, 0.2312}},
            {"HSS12X2X3/16", {1.8236, 11.6772, 0.173}},
            {"HSS10X10X3/4", {9.16, 9.16, 0.673}},
            {"HSS10X10X5/8", {9.3051, 9.3051, 0.5636}},
            {"HSS10X10X1/2", {9.4438, 9.4438, 0.4544}},
            {"HSS10X10X3/8", {9.5821, 9.5821, 0.3435}},
            {"HSS10X10X5/16", {9.6389, 9.6389, 0.2878}},
            {"HSS10X10X1/4", {9.712, 9.712, 0.2309}},
            {"HSS10X10X3/16", {9.7902, 9.7902, 0.1726}},
            {"HSS10X8X5/8", {7.3326, 9.2685, 0.5619}},
            {"HSS10X8X1/2", {7.4728, 9.4268, 0.4522}},
            {"HSS10X8X3/8", {7.595, 9.5409, 0.3438}},
            {"HSS10X8X5/16", {7.6775, 9.6382, 0.2861}},
            {"HSS10X8X1/4", {7.7389, 9.702, 0.2302}},
            {"HSS10X8X3/16", {7.8108, 9.792, 0.172}},
            {"HSS10X6X5/8", {5.3582, 9.1814, 0.5627}},
            {"HSS10X6X1/2", {5.4853, 9.3351, 0.4545}},
            {"HSS10X6X3/8", {5.6147, 9.523, 0.3427}},
            {"HSS10X6X5/16", {5.6845, 9.6315, 0.2853}},
            {"HSS10X6X1/4", {5.745, 9.6858, 0.23}},
            {"HSS10X6X3/16", {5.8131, 9.7685, 0.1723}},
            {"HSS10X5X3/8", {4.6284, 9.4745, 0.3426}},
            {"HSS10X5X5/16", {4.6922, 9.5986, 0.2856}},
            {"HSS10X5X1/4", {4.7572, 9.6658, 0.2299}},
            {"HSS10X4X5/8", {3.3928, 9.0627, 0.5615}},
            {"HSS10X4X1/2", {3.514, 9.2594, 0.4533}},
            {"HSS10X4X3/8", {3.6383, 9.4527, 0.3418}},
            {"HSS10X4X5/16", {3.7009, 9.5505, 0.2861}},
            {"HSS10X4X1/4", {3.7542, 9.6366, 0.2305}},
            {"HSS10X4X3/16", {3.8255, 9.7476, 0.172}},
            {"HSS10X4X1/8", {3.8761, 9.8256, 0.1154}},
            {"HSS10X3-1/2X3/8", {3.1387, 9.4419, 0.3417}},
            {"HSS10X3-1/2X5/16", {3.1984, 9.535, 0.2862}},
            {"HSS10X3-1/2X1/4", {3.2548, 9.6318, 0.2302}},
            {"HSS10X3-1/2X3/16", {3.3191, 9.7296, 0.1723}},
            {"HSS10X3X3/8", {2.6382, 9.4123, 0.3429}},
            {"HSS10X3X5/16", {2.7063, 9.5093, 0.2865}},
            {"HSS10X3X1/4", {2.7623, 9.617, 0.2299}},
            {"HSS10X3X3/16", {2.8254, 9.724, 0.1718}},
            {"HSS10X3X1/8", {2.8808, 9.8139, 0.1153}},
            {"HSS10X2X3/8", {1.643, 9.3377, 0.3446}},
            {"HSS10X2X5/16", {1.7021, 9.4532, 0.288}},
            {"HSS10X2X1/4", {1.7625, 9.5681, 0.2308}},
            {"HSS10X2X3/16", {1.8247, 9.6842, 0.1726}},
            {"HSS10X2X1/8", {1.8838, 9.7889, 0.1155}},
            {"HSS9X9X5/8", {8.3104, 8.3104, 0.5613}},
            {"HSS9X9X1/2", {8.4637, 8.4637, 0.451}},
            {"HSS9X9X3/8", {8.5854, 8.5854, 0.3427}},
            {"HSS9X9X5/16", {8.6596, 8.6596, 0.2858}},
            {"HSS9X9X1/4", {8.731, 8.731, 0.2295}},
            {"HSS9X9X3/16", {8.7972, 8.7972, 0.1722}},
            {"HSS9X9X1/8", {8.8572, 8.8572, 0.1155}},
            {"HSS9X7X5/8", {6.3187, 8.2364, 0.5622}},
            {"HSS9X7X1/2", {6.4484, 8.4088, 0.4532}},
            {"HSS9X7X3/8", {6.5962, 8.5672, 0.3419}},
            {"HSS9X7X5/16", {6.6702, 8.6382, 0.2857}},
            {"HSS9X7X1/4", {6.7376, 8.7101, 0.2298}},
            {"HSS9X7X3/16", {6.8085, 8.7809, 0.1722}},
            {"HSS9X5X5/8", {4.3761, 8.1712, 0.5569}},
            {"HSS9X5X1/2", {4.4972, 8.3491, 0.4504}},
            {"HSS9X5X3/8", {4.6251, 8.5106, 0.341}},
            {"HSS9X5X5/16", {4.6918, 8.5989, 0.2849}},
            {"HSS9X5X1/4", {4.7458, 8.6811, 0.2294}},
            {"HSS9X5X3/16", {4.8205, 8.7606, 0.1717}},
            {"HSS9X3X1/2", {2.5127, 8.2154, 0.4528}},
            {"HSS9X3X3/8", {2.6396, 8.4197, 0.3421}},
            {"HSS9X3X5/16", {2.6992, 8.5178, 0.2861}},
            {"HSS9X3X1/4", {2.7576, 8.6167, 0.2302}},
            {"HSS9X3X3/16", {2.8202, 8.727, 0.1722}},
            {"HSS8X8X5/8", {7.2912, 7.2912, 0.5609}},
            {"HSS8X8X1/2", {7.4453, 7.4453, 0.4519}},
            {"HSS8X8X3/8", {7.5953, 7.5953, 0.3409}},
            {"HSS8X8X5/16", {7.6537, 7.6537, 0.2858}},
            {"HSS8X8X1/4", {7.7302, 7.7302, 0.2291}},
            {"HSS8X8X3/16", {7.7943, 7.7943, 0.1723}},
            {"HSS8X8X1/8", {7.8712, 7.8712, 0.1151}},
            {"HSS8X6X5/8", {5.3451, 7.2459, 0.5553}},
            {"HSS8X6X1/2", {5.4703, 7.399, 0.4497}},
            {"HSS8X6X3/8", {5.6078, 7.5487, 0.3406}},
            {"HSS8X6X5/16", {5.6753, 7.6306, 0.2848}},
            {"HSS8X6X1/4", {5.7428, 7.7034, 0.2291}},
            {"HSS8X6X3/16", {5.8083, 7.7745, 0.1722}},
            {"HSS8X4X5/8", {3.3743, 7.1236, 0.5558}},
            {"HSS8X4X1/2", {3.4975, 7.3079, 0.4501}},
            {"HSS8X4X3/8", {3.6253, 7.4895, 0.3407}},
            {"HSS8X4X5/16", {3.6952, 7.5767, 0.2846}},
            {"HSS8X4X1/4", {3.7483, 7.6598, 0.2294}},
            {"HSS8X4X3/16", {3.8163, 7.752, 0.1717}},
            {"HSS8X4X1/8", {3.8753, 7.8244, 0.1153}},
            {"HSS8X3X1/2", {2.5091, 7.238, 0.451}},
            {"HSS8X3X3/8", {2.642, 7.4416, 0.3404}},
            {"HSS8X3X5/16", {2.6948, 7.541, 0.2854}},
            {"HSS8X3X1/4", {2.7573, 7.6341, 0.2293}},
            {"HSS8X3X3/16", {2.8201, 7.7372, 0.1717}},
            {"HSS8X3X1/8", {2.8834, 7.8208, 0.1149}},
            {"HSS8X2X1/2", {1.52, 7.1475, 0.4525}},
            {"HSS8X2X3/8", {1.6403, 7.3595, 0.3427}},
            {"HSS8X2X5/16", {1.7015, 7.4796, 0.2858}},
            {"HSS8X2X1/4", {1.7613, 7.5959, 0.2295}},
            {"HSS8X2X3/16", {1.8239, 7.6966, 0.1721}},
            {"HSS8X2X1/8", {1.8812, 7.7934, 0.1153}},
            {"HSS7X7X5/8", {6.3042, 6.3042, 0.5545}},
            {"HSS7X7X1/2", {6.4415, 6.4415, 0.4489}},
            {"HSS7X7X3/8", {6.5878, 6.5878, 0.3398}},
            {"HSS7X7X5/16", {6.6556, 6.6556, 0.2847}},
            {"HSS7X7X1/4", {6.7239, 6.7239, 0.2289}},
            {"HSS7X7X3/16", {6.7995, 6.7995, 0.1716}},
            {"HSS7X7X1/8", {6.8613, 6.8613, 0.1151}},
            {"HSS7X5X1/2", {4.4734, 6.3796, 0.4481}},
            {"HSS7X5X3/8", {4.6099, 6.5452, 0.339}},
            {"HSS7X5X5/16", {4.6759, 6.6204, 0.2841}},
            {"HSS7X5X1/4", {4.7354, 6.7027, 0.2289}},
            {"HSS7X5X3/16", {4.801, 6.7789, 0.1717}},
            {"HSS7X5X1/8", {4.8764, 6.8439, 0.1151}},
            {"HSS7X4X1/2", {3.4915, 6.3374, 0.4464}},
            {"HSS7X4X3/8", {3.6264, 6.5112, 0.3381}},
            {"HSS7X4X5/16", {3.6877, 6.5993, 0.2839}},
            {"HSS7X4X1/4", {3.7542, 6.6759, 0.2284}},
            {"HSS7X4X3/16", {3.8067, 6.7578, 0.1717}},
            {"HSS7X4X1/8", {3.8806, 6.8539, 0.1145}},
            {"HSS7X3X1/2", {2.5051, 6.2576, 0.4484}},
            {"HSS7X3X3/8", {2.6315, 6.4634, 0.3389}},
            {"HSS7X3X5/16", {2.6944, 6.5574, 0.2838}},
            {"HSS7X3X1/4", {2.7576, 6.6528, 0.2282}},
            {"HSS7X3X3/16", {2.8182, 6.7478, 0.1713}},
            {"HSS7X3X1/8", {2.8789, 6.8273, 0.1148}},
            {"HSS7X2X1/4", {1.7609, 6.5925, 0.2294}},
            {"HSS7X2X3/16", {1.8249, 6.7056, 0.1716}},
            {"HSS7X2X1/8", {1.8842, 6.8141, 0.1148}},
            {"HSS6X6X5/8", {5.2914, 5.2914, 0.5519}},
            {"HSS6X6X1/2", {5.4399, 5.4399, 0.4464}},
            {"HSS6X6X3/8", {5.5853, 5.5853, 0.3383}},
            {"HSS6X6X5/16", {5.6523, 5.6523, 0.284}},
            {"HSS6X6X1/4", {5.7209, 5.7209, 0.2285}},
            {"HSS6X6X3/16", {5.7977, 5.7977, 0.1714}},
            {"HSS6X6X1/8", {5.8709, 5.8709, 0.1147}},
            {"HSS6X5X1/2", {4.4586, 5.402, 0.4456}},
            {"HSS6X5X3/8", {4.6023, 5.5616, 0.3373}},
            {"HSS6X5X5/16", {4.6707, 5.6412, 0.2826}},
            {"HSS6X5X1/4", {4.7379, 5.7043, 0.2281}},
            {"HSS6X5X3/16", {4.8003, 5.782, 0.1714}},
            {"HSS6X5X1/8", {4.8758, 5.8517, 0.1147}},
            {"HSS6X4X1/2", {3.4829, 5.3603, 0.4437}},
            {"HSS6X4X3/8", {3.6063, 5.5303, 0.3377}},
            {"HSS6X4X5/16", {3.6851, 5.6126, 0.2825}},
            {"HSS6X4X1/4", {3.7427, 5.6981, 0.2275}},
            {"HSS6X4X3/16", {3.8112, 5.7799, 0.1707}},
            {"HSS6X4X1/8", {3.8744, 5.8405, 0.1148}},
            {"HSS6X3X1/2", {2.4963, 5.2859, 0.4452}},
            {"HSS6X3X3/8", {2.6259, 5.4779, 0.3374}},
            {"HSS6X3X5/16", {2.6901, 5.5776, 0.2825}},
            {"HSS6X3X1/4", {2.7502, 5.6628, 0.2278}},
            {"HSS6X3X3/16", {2.8177, 5.7501, 0.1709}},
            {"HSS6X3X1/8", {2.878, 5.8361, 0.1146}},
            {"HSS6X2X3/8", {1.6371, 5.3963, 0.339}},
            {"HSS6X2X5/16", {1.695, 5.5056, 0.2842}},
            {"HSS6X2X1/4", {1.7584, 5.6122, 0.2283}},
            {"HSS6X2X3/16", {1.8196, 5.7343, 0.1707}},
            {"HSS6X2X1/8", {1.8795, 5.8103, 0.115}},
            {"HSS5-1/2X5-1/2X3/8", {5.0791, 5.0791, 0.3383}},
            {"HSS5-1/2X5-1/2X5/16", {5.1463, 5.1463, 0.2842}},
            {"HSS5-1/2X5-1/2X1/4", {5.2212, 5.2212, 0.2281}},
            {"HSS5-1/2X5-1/2X3/16", {5.3006, 5.3006, 0.1709}},
            {"HSS5-1/2X5-1/2X1/8", {5.3633, 5.3633, 0.1147}},
            {"HSS5X5X1/2", {4.4282, 4.4282, 0.4444}},
            {"HSS5X5X3/8", {4.5791, 4.5791, 0.3369}},
            {"HSS5X5X5/16", {4.648, 4.648, 0.2826}},
            {"HSS5X5X1/4", {4.72, 4.72, 0.2277}},
            {"HSS5X5X3/16", {4.7994, 4.7994, 0.1706}},
            {"HSS5X5X1/8", {4.8644, 4.8644, 0.1146}},
            {"HSS5X4X1/2", {3.4615, 4.3812, 0.4418}},
            {"HSS5X4X3/8", {3.5987, 4.5501, 0.3355}},
            {"HSS5X4X5/16", {3.6588, 4.6295, 0.2819}},
            {"HSS5X4X1/4", {3.735, 4.7092, 0.2269}},
            {"HSS5X4X3/16", {3.8046, 4.7956, 0.1702}},
            {"HSS5X4X1/8", {3.8675, 4.856, 0.1146}},
            {"HSS5X3X1/2", {2.4809, 4.3145, 0.4418}},
            {"HSS5X3X3/8", {2.6172, 4.5025, 0.3349}},
            {"HSS5X3X5/16", {2.6822, 4.5978, 0.2809}},
            {"HSS5X3X1/4", {2.7475, 4.6722, 0.2268}},
            {"HSS5X3X3/16", {2.8134, 4.767, 0.1701}},
            {"HSS5X3X1/8", {2.8755, 4.8372, 0.1146}},
            {"HSS5X2-1/2X1/4", {2.2501, 4.6584, 0.2267}},
            {"HSS5X2-1/2X3/16", {2.3145, 4.7494, 0.1704}},
            {"HSS5X2-1/2X1/8", {2.3765, 4.8352, 0.1144}},
            {"HSS5X2X3/8", {1.6289, 4.4374, 0.3356}},
            {"HSS5X2X5/16", {1.694, 4.5324, 0.2817}},
            {"HSS5X2X1/4", {1.7519, 4.6286, 0.2276}},
            {"HSS5X2X3/16", {1.8149, 4.7279, 0.1709}},
            {"HSS5X2X1/8", {1.8725, 4.8176, 0.115}},
            {"HSS4-1/2X4-1/2X1/2", {3.9255, 3.9255, 0.4418}},
            {"HSS4-1/2X4-1/2X3/8", {4.0818, 4.0818, 0.3347}},
            {"HSS4-1/2X4-1/2X5/16", {4.1526, 4.1526, 0.2812}},
            {"HSS4-1/2X4-1/2X1/4", {4.2157, 4.2157, 0.2274}},
            {"HSS4-1/2X4-1/2X3/16", {4.2958, 4.2958, 0.1703}},
            {"HSS4-1/2X4-1/2X1/8", {4.3637, 4.3637, 0.1145}},
            {"HSS4X4X1/2", {3.4082, 3.4082, 0.4418}},
            {"HSS4X4X3/8", {3.5835, 3.5835, 0.3322}},
            {"HSS4X4X5/16", {3.6488, 3.6488, 0.2802}},
            {"HSS4X4X1/4", {3.7214, 3.7214, 0.2259}},
            {"HSS4X4X3/16", {3.7972, 3.7972, 0.1697}},
            {"HSS4X4X1/8", {3.861, 3.861, 0.1145}},
            {"HSS4X3X3/8", {2.5975, 3.5322, 0.3325}},
            {"HSS4X3X5/16", {2.6665, 3.6173, 0.2792}},
            {"HSS4X3X1/4", {2.7332, 3.6962, 0.2256}},
            {"HSS4X3X3/16", {2.8052, 3.7737, 0.1698}},
            {"HSS4X3X1/8", {2.8697, 3.8466, 0.1145}},
            {"HSS4X2-1/2X1/4", {2.2443, 3.6797, 0.2251}},
            {"HSS4X2-1/2X3/16", {2.3114, 3.7664, 0.1693}},
            {"HSS4X2-1/2X1/8", {2.3715, 3.8477, 0.114}},
            {"HSS4X2X3/8", {1.6199, 3.4684, 0.3315}},
            {"HSS4X2X5/16", {1.687, 3.5542, 0.2794}},
            {"HSS4X2X1/4", {1.7503, 3.651, 0.2254}},
            {"HSS4X2X3/16", {1.8116, 3.7468, 0.1697}},
            {"HSS4X2X1/8", {1.878, 3.8409, 0.1136}},
            {"HSS4X1-1/2X1/4", {1.2542, 3.6093, 0.2268}},
            {"HSS4X1-1/2X3/16", {1.3195, 3.7244, 0.1694}},
            {"HSS4X1-1/2X1/8", {1.3783, 3.813, 0.1144}},
            {"HSS3-1/2X3-1/2X3/8", {3.0702, 3.0702, 0.3318}},
            {"HSS3-1/2X3-1/2X5/16", {3.1445, 3.1445, 0.2791}},
            {"HSS3-1/2X3-1/2X1/4", {3.2175, 3.2175, 0.2256}},
            {"HSS3-1/2X3-1/2X3/16", {3.2901, 3.2901, 0.17}},
            {"HSS3-1/2X3-1/2X1/8", {3.362, 3.362, 0.1142}},
            {"HSS3-1/2X2-1/2X3/8", {2.096, 3.0111, 0.3315}},
            {"HSS3-1/2X2-1/2X5/16", {2.1696, 3.1058, 0.2773}},
            {"HSS3-1/2X2-1/2X1/4", {2.2369, 3.189, 0.2243}},
            {"HSS3-1/2X2-1/2X3/16", {2.3008, 3.2752, 0.1691}},
            {"HSS3-1/2X2-1/2X1/8", {2.3762, 3.3543, 0.1133}},
            {"HSS3-1/2X2X1/4", {1.7466, 3.1636, 0.2243}},
            {"HSS3-1/2X2X3/16", {1.8159, 3.2651, 0.1679}},
            {"HSS3-1/2X2X1/8", {1.8742, 3.3361, 0.1139}},
            {"HSS3-1/2X1-1/2X1/4", {1.2524, 3.1323, 0.2242}},
            {"HSS3-1/2X1-1/2X3/16", {1.3161, 3.228, 0.1692}},
            {"HSS3-1/2X1-1/2X1/8", {1.379, 3.3292, 0.1136}},
            {"HSS3X3X3/8", {2.5579, 2.5579, 0.3315}},
            {"HSS3X3X5/16", {2.6411, 2.6411, 0.2773}},
            {"HSS3X3X1/4", {2.7174, 2.7174, 0.224}},
            {"HSS3X3X3/16", {2.7912, 2.7912, 0.1688}},
            {"HSS3X3X1/8", {2.8631, 2.8631, 0.1137}},
            {"HSS3X2-1/2X5/16", {2.1534, 2.6213, 0.2764}},
            {"HSS3X2-1/2X1/4", {2.2271, 2.6985, 0.2238}},
            {"HSS3X2-1/2X3/16", {2.3022, 2.7837, 0.1679}},
            {"HSS3X2-1/2X1/8", {2.3625, 2.8502, 0.1141}},
            {"HSS3X2X5/16", {1.6691, 2.5833, 0.2764}},
            {"HSS3X2X1/4", {1.7388, 2.6823, 0.2223}},
            {"HSS3X2X3/16", {1.8085, 2.7681, 0.168}},
            {"HSS3X2X1/8", {1.874, 2.8483, 0.1131}},
            {"HSS3X1-1/2X1/4", {1.2464, 2.6446, 0.2228}},
            {"HSS3X1-1/2X3/16", {1.3102, 2.7416, 0.1688}},
            {"HSS3X1-1/2X1/8", {1.3755, 2.8334, 0.1134}},
            {"HSS3X1X3/16", {0.8188, 2.7048, 0.1686}},
            {"HSS3X1X1/8", {0.8784, 2.8076, 0.1138}},
            {"HSS2-1/2X2-1/2X5/16", {2.13, 2.13, 0.2764}},
            {"HSS2-1/2X2-1/2X1/4", {2.2176, 2.2176, 0.2214}},
            {"HSS2-1/2X2-1/2X3/16", {2.288, 2.288, 0.168}},
            {"HSS2-1/2X2-1/2X1/8", {2.3636, 2.3636, 0.1131}},
            {"HSS2-1/2X2X1/4", {1.727, 2.1947, 0.2214}},
            {"HSS2-1/2X2X3/16", {1.7966, 2.2778, 0.1676}},
            {"HSS2-1/2X2X1/8", {1.8667, 2.3523, 0.1132}},
            {"HSS2-1/2X1-1/2X1/4", {1.2367, 2.1627, 0.2214}},
            {"HSS2-1/2X1-1/2X3/16", {1.3091, 2.2568, 0.1665}},
            {"HSS2-1/2X1-1/2X1/8", {1.3723, 2.3385, 0.1131}},
            {"HSS2-1/2X1X3/16", {0.813, 2.2114, 0.1684}},
            {"HSS2-1/2X1X1/8", {0.8765, 2.314, 0.1133}},
            {"HSS2-1/4X2-1/4X1/4", {1.9618, 1.9618, 0.2214}},
            {"HSS2-1/4X2-1/4X3/16", {2.0382, 2.0382, 0.1673}},
            {"HSS2-1/4X2-1/4X1/8", {2.111, 2.111, 0.1131}},
            {"HSS2X2X1/4", {1.7043, 1.7043, 0.2214}},
            {"HSS2X2X3/16", {1.7899, 1.7899, 0.1659}},
            {"HSS2X2X1/8", {1.8591, 1.8591, 0.1128}},
            {"HSS2X1-1/2X3/16", {1.2966, 1.7646, 0.1663}},
            {"HSS2X1-1/2X1/8", {1.3665, 1.8467, 0.1125}},
            {"HSS2X1X3/16", {0.8043, 1.7339, 0.1661}},
            {"HSS2X1X1/8", {0.8736, 1.8241, 0.1126}},
            {"HSS1-1/2X1-1/2X1/4", {1.2036, 1.2036, 0.2214}},
            {"HSS1-1/2X1-1/2X3/16", {1.2784, 1.2784, 0.1653}},
            {"HSS1-1/2X1-1/2X1/8", {1.3558, 1.3558, 0.112}},
        };

        return registry.fetch(type);
    }

    vec uschssection(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"HSS28.000X1.000", {14.000, 0.930}},
            {"HSS28.000X0.875", {14.000, 0.814}},
            {"HSS28.000X0.750", {14.000, 0.698}},
            {"HSS28.000X0.625", {14.000, 0.581}},
            {"HSS28.000X0.500", {14.000, 0.465}},
            {"HSS28.000X0.375", {14.000, 0.349}},
            {"HSS26.000X0.750", {13.000, 0.698}},
            {"HSS26.000X0.625", {13.000, 0.581}},
            {"HSS26.000X0.500", {13.000, 0.465}},
            {"HSS26.000X0.375", {13.000, 0.349}},
            {"HSS26.000X0.313", {13.000, 0.291}},
            {"HSS24.000X1.000", {12.000, 0.930}},
            {"HSS24.000X0.875", {12.000, 0.814}},
            {"HSS24.000X0.750", {12.000, 0.698}},
            {"HSS24.000X0.625", {12.000, 0.581}},
            {"HSS24.000X0.500", {12.000, 0.465}},
            {"HSS24.000X0.375", {12.000, 0.349}},
            {"HSS24.000X0.313", {12.000, 0.291}},
            {"HSS22.000X0.750", {11.000, 0.698}},
            {"HSS22.000X0.625", {11.000, 0.581}},
            {"HSS22.000X0.500", {11.000, 0.465}},
            {"HSS22.000X0.375", {11.000, 0.349}},
            {"HSS22.000X0.313", {11.000, 0.291}},
            {"HSS20.000X1.000", {10.000, 0.930}},
            {"HSS20.000X0.875", {10.000, 0.814}},
            {"HSS20.000X0.750", {10.000, 0.698}},
            {"HSS20.000X0.625", {10.000, 0.581}},
            {"HSS20.000X0.500", {10.000, 0.465}},
            {"HSS20.000X0.375", {10.000, 0.349}},
            {"HSS20.000X0.313", {10.000, 0.291}},
            {"HSS20.000X0.250", {10.000, 0.233}},
            {"HSS18.000X1.000", {9.000, 0.930}},
            {"HSS18.000X0.875", {9.000, 0.814}},
            {"HSS18.000X0.750", {9.000, 0.698}},
            {"HSS18.000X0.625", {9.000, 0.581}},
            {"HSS18.000X0.500", {9.000, 0.465}},
            {"HSS18.000X0.375", {9.000, 0.349}},
            {"HSS18.000X0.313", {9.000, 0.291}},
            {"HSS18.000X0.250", {9.000, 0.233}},
            {"HSS16.000X1.000", {8.000, 0.930}},
            {"HSS16.000X0.875", {8.000, 0.814}},
            {"HSS16.000X0.750", {8.000, 0.698}},
            {"HSS16.000X0.625", {8.000, 0.581}},
            {"HSS16.000X0.500", {8.000, 0.465}},
            {"HSS16.000X0.438", {8.000, 0.407}},
            {"HSS16.000X0.375", {8.000, 0.349}},
            {"HSS16.000X0.312", {8.000, 0.291}},
            {"HSS16.000X0.250", {8.000, 0.233}},
            {"HSS14.000X1.000", {7.000, 0.930}},
            {"HSS14.000X0.875", {7.000, 0.814}},
            {"HSS14.000X0.750", {7.000, 0.698}},
            {"HSS14.000X0.625", {7.000, 0.581}},
            {"HSS14.000X0.500", {7.000, 0.465}},
            {"HSS14.000X0.375", {7.000, 0.349}},
            {"HSS14.000X0.312", {7.000, 0.291}},
            {"HSS14.000X0.250", {7.000, 0.233}},
            {"HSS14.000X0.188", {7.000, 0.174}},
            {"HSS13.375X0.625", {6.700, 0.581}},
            {"HSS13.375X0.500", {6.700, 0.465}},
            {"HSS13.375X0.375", {6.700, 0.349}},
            {"HSS13.375X0.313", {6.700, 0.291}},
            {"HSS13.375X0.250", {6.700, 0.233}},
            {"HSS13.375X0.188", {6.700, 0.174}},
            {"HSS12.750X0.750", {6.400, 0.698}},
            {"HSS12.750X0.625", {6.400, 0.581}},
            {"HSS12.750X0.500", {6.400, 0.465}},
            {"HSS12.750X0.375", {6.400, 0.349}},
            {"HSS12.750X0.250", {6.400, 0.233}},
            {"HSS12.750X0.188", {6.400, 0.174}},
            {"HSS12.000X0.625", {6.000, 0.581}},
            {"HSS12.000X0.500", {6.000, 0.465}},
            {"HSS12.000X0.375", {6.000, 0.349}},
            {"HSS12.000X0.250", {6.000, 0.233}},
            {"HSS11.750X0.625", {5.900, 0.581}},
            {"HSS11.750X0.500", {5.900, 0.465}},
            {"HSS11.750X0.375", {5.900, 0.349}},
            {"HSS11.750X0.337", {5.900, 0.313}},
            {"HSS11.750X0.250", {5.900, 0.233}},
            {"HSS10.750X0.625", {5.400, 0.581}},
            {"HSS10.750X0.500", {5.400, 0.465}},
            {"HSS10.750X0.375", {5.400, 0.349}},
            {"HSS10.750X0.313", {5.400, 0.291}},
            {"HSS10.750X0.250", {5.400, 0.233}},
            {"HSS10.750X0.188", {5.400, 0.174}},
            {"HSS10.000X0.625", {5.000, 0.581}},
            {"HSS10.000X0.500", {5.000, 0.465}},
            {"HSS10.000X0.375", {5.000, 0.349}},
            {"HSS10.000X0.312", {5.000, 0.291}},
            {"HSS10.000X0.250", {5.000, 0.233}},
            {"HSS10.000X0.188", {5.000, 0.174}},
            {"HSS9.625X0.625", {4.815, 0.581}},
            {"HSS9.625X0.500", {4.815, 0.465}},
            {"HSS9.625X0.375", {4.815, 0.349}},
            {"HSS9.625X0.312", {4.815, 0.291}},
            {"HSS9.625X0.250", {4.815, 0.233}},
            {"HSS9.625X0.188", {4.815, 0.174}},
            {"HSS8.625X0.625", {4.315, 0.581}},
            {"HSS8.625X0.500", {4.315, 0.465}},
            {"HSS8.625X0.375", {4.315, 0.349}},
            {"HSS8.625X0.322", {4.315, 0.300}},
            {"HSS8.625X0.250", {4.315, 0.233}},
            {"HSS8.625X0.188", {4.315, 0.174}},
            {"HSS7.500X0.500", {3.750, 0.465}},
            {"HSS7.500X0.375", {3.750, 0.349}},
            {"HSS7.500X0.312", {3.750, 0.291}},
            {"HSS7.500X0.250", {3.750, 0.233}},
            {"HSS7.500X0.188", {3.750, 0.174}},
            {"HSS7.000X0.500", {3.500, 0.465}},
            {"HSS7.000X0.375", {3.500, 0.349}},
            {"HSS7.000X0.312", {3.500, 0.291}},
            {"HSS7.000X0.250", {3.500, 0.233}},
            {"HSS7.000X0.188", {3.500, 0.174}},
            {"HSS7.000X0.125", {3.500, 0.116}},
            {"HSS6.875X0.375", {3.440, 0.349}},
            {"HSS6.875X0.312", {3.440, 0.291}},
            {"HSS6.875X0.250", {3.440, 0.233}},
            {"HSS6.875X0.188", {3.440, 0.174}},
            {"HSS6.625X0.500", {3.315, 0.465}},
            {"HSS6.625X0.432", {3.315, 0.402}},
            {"HSS6.625X0.375", {3.315, 0.349}},
            {"HSS6.625X0.312", {3.315, 0.291}},
            {"HSS6.625X0.280", {3.315, 0.260}},
            {"HSS6.625X0.250", {3.315, 0.233}},
            {"HSS6.625X0.188", {3.315, 0.174}},
            {"HSS6.625X0.125", {3.315, 0.116}},
            {"HSS6.000X0.500", {3.000, 0.465}},
            {"HSS6.000X0.375", {3.000, 0.349}},
            {"HSS6.000X0.312", {3.000, 0.291}},
            {"HSS6.000X0.280", {3.000, 0.260}},
            {"HSS6.000X0.250", {3.000, 0.233}},
            {"HSS6.000X0.188", {3.000, 0.174}},
            {"HSS6.000X0.125", {3.000, 0.116}},
            {"HSS5.563X0.500", {2.780, 0.465}},
            {"HSS5.563X0.375", {2.780, 0.349}},
            {"HSS5.563X0.258", {2.780, 0.240}},
            {"HSS5.563X0.188", {2.780, 0.174}},
            {"HSS5.563X0.134", {2.780, 0.124}},
            {"HSS5.500X0.500", {2.750, 0.465}},
            {"HSS5.500X0.375", {2.750, 0.349}},
            {"HSS5.500X0.258", {2.750, 0.240}},
            {"HSS5.000X0.500", {2.500, 0.465}},
            {"HSS5.000X0.375", {2.500, 0.349}},
            {"HSS5.000X0.312", {2.500, 0.291}},
            {"HSS5.000X0.258", {2.500, 0.240}},
            {"HSS5.000X0.250", {2.500, 0.233}},
            {"HSS5.000X0.188", {2.500, 0.174}},
            {"HSS5.000X0.125", {2.500, 0.116}},
            {"HSS4.500X0.375", {2.250, 0.349}},
            {"HSS4.500X0.337", {2.250, 0.313}},
            {"HSS4.500X0.237", {2.250, 0.220}},
            {"HSS4.500X0.188", {2.250, 0.174}},
            {"HSS4.500X0.125", {2.250, 0.116}},
            {"HSS4.000X0.313", {2.000, 0.291}},
            {"HSS4.000X0.250", {2.000, 0.233}},
            {"HSS4.000X0.237", {2.000, 0.220}},
            {"HSS4.000X0.226", {2.000, 0.210}},
            {"HSS4.000X0.220", {2.000, 0.205}},
            {"HSS4.000X0.188", {2.000, 0.174}},
            {"HSS4.000X0.125", {2.000, 0.116}},
            {"HSS3.500X0.313", {1.750, 0.291}},
            {"HSS3.500X0.300", {1.750, 0.279}},
            {"HSS3.500X0.250", {1.750, 0.233}},
            {"HSS3.500X0.216", {1.750, 0.201}},
            {"HSS3.500X0.203", {1.750, 0.189}},
            {"HSS3.500X0.188", {1.750, 0.174}},
            {"HSS3.500X0.125", {1.750, 0.116}},
            {"HSS3.000X0.250", {1.500, 0.233}},
            {"HSS3.000X0.216", {1.500, 0.201}},
            {"HSS3.000X0.203", {1.500, 0.189}},
            {"HSS3.000X0.188", {1.500, 0.174}},
            {"HSS3.000X0.152", {1.500, 0.141}},
            {"HSS3.000X0.134", {1.500, 0.124}},
            {"HSS3.000X0.125", {1.500, 0.116}},
            {"HSS2.875X0.250", {1.440, 0.233}},
            {"HSS2.875X0.203", {1.440, 0.189}},
            {"HSS2.875X0.188", {1.440, 0.174}},
            {"HSS2.875X0.125", {1.440, 0.116}},
            {"HSS2.500X0.250", {1.250, 0.233}},
            {"HSS2.500X0.188", {1.250, 0.174}},
            {"HSS2.500X0.125", {1.250, 0.116}},
            {"HSS2.375X0.250", {1.190, 0.233}},
            {"HSS2.375X0.218", {1.190, 0.203}},
            {"HSS2.375X0.188", {1.190, 0.174}},
            {"HSS2.375X0.154", {1.190, 0.143}},
            {"HSS2.375X0.125", {1.190, 0.116}},
            {"HSS1.900X0.188", {0.950, 0.174}},
            {"HSS1.900X0.145", {0.950, 0.135}},
            {"HSS1.900X0.120", {0.950, 0.111}},
            {"HSS1.660X0.140", {0.830, 0.130}},
        };

        return registry.fetch(type);
    }

    // Adjusted registry: dimensions optimized against the AISC tabulated
    // section properties with section_optimizer.py.
    vec uschssection_adjusted(const std::string_view type) {
        static const TokenRegistry<vec> registry{
            {"HSS28.000X1.000", {13.9957, 0.9307}},
            {"HSS28.000X0.875", {14.0028, 0.8138}},
            {"HSS28.000X0.750", {14.005, 0.6969}},
            {"HSS28.000X0.625", {14.0012, 0.5813}},
            {"HSS28.000X0.500", {13.9948, 0.4654}},
            {"HSS28.000X0.375", {14.0077, 0.3487}},
            {"HSS26.000X0.750", {13.0022, 0.6973}},
            {"HSS26.000X0.625", {12.9988, 0.5815}},
            {"HSS26.000X0.500", {12.9985, 0.4649}},
            {"HSS26.000X0.375", {12.9983, 0.3484}},
            {"HSS26.000X0.313", {12.9962, 0.2908}},
            {"HSS24.000X1.000", {11.9988, 0.9299}},
            {"HSS24.000X0.875", {11.9998, 0.8143}},
            {"HSS24.000X0.750", {11.9979, 0.6981}},
            {"HSS24.000X0.625", {11.9877, 0.5826}},
            {"HSS24.000X0.500", {11.9912, 0.466}},
            {"HSS24.000X0.375", {11.9948, 0.3488}},
            {"HSS24.000X0.313", {12.0087, 0.2896}},
            {"HSS22.000X0.750", {10.9951, 0.6986}},
            {"HSS22.000X0.625", {10.9881, 0.5825}},
            {"HSS22.000X0.500", {10.9776, 0.4672}},
            {"HSS22.000X0.375", {11.0068, 0.3477}},
            {"HSS22.000X0.313", {11.0154, 0.2899}},
            {"HSS20.000X1.000", {10.005, 0.9284}},
            {"HSS20.000X0.875", {10.0004, 0.8136}},
            {"HSS20.000X0.750", {9.9931, 0.6983}},
            {"HSS20.000X0.625", {9.9877, 0.5824}},
            {"HSS20.000X0.500", {10.0004, 0.464}},
            {"HSS20.000X0.375", {10.0047, 0.3487}},
            {"HSS20.000X0.313", {9.9979, 0.2909}},
            {"HSS20.000X0.250", {10.0076, 0.232}},
            {"HSS18.000X1.000", {8.9944, 0.9308}},
            {"HSS18.000X0.875", {9.0127, 0.8123}},
            {"HSS18.000X0.750", {8.9972, 0.6978}},
            {"HSS18.000X0.625", {9.0124, 0.5791}},
            {"HSS18.000X0.500", {9.0011, 0.4647}},
            {"HSS18.000X0.375", {8.9887, 0.3505}},
            {"HSS18.000X0.313", {8.9921, 0.2914}},
            {"HSS18.000X0.250", {8.9918, 0.233}},
            {"HSS16.000X1.000", {7.9888, 0.9308}},
            {"HSS16.000X0.875", {7.9925, 0.8152}},
            {"HSS16.000X0.750", {8.007, 0.6951}},
            {"HSS16.000X0.625", {8.0081, 0.5791}},
            {"HSS16.000X0.500", {7.9999, 0.4645}},
            {"HSS16.000X0.438", {8.0031, 0.4064}},
            {"HSS16.000X0.375", {7.9953, 0.3497}},
            {"HSS16.000X0.312", {7.9902, 0.2919}},
            {"HSS16.000X0.250", {8.0156, 0.2318}},
            {"HSS14.000X1.000", {7.0016, 0.9294}},
            {"HSS14.000X0.875", {6.997, 0.8153}},
            {"HSS14.000X0.750", {7.0041, 0.6976}},
            {"HSS14.000X0.625", {6.9948, 0.5827}},
            {"HSS14.000X0.500", {6.9939, 0.4659}},
            {"HSS14.000X0.375", {6.9948, 0.35}},
            {"HSS14.000X0.312", {7.013, 0.2898}},
            {"HSS14.000X0.250", {6.9959, 0.2336}},
            {"HSS14.000X0.188", {7.0025, 0.1741}},
            {"HSS13.375X0.625", {6.6838, 0.5823}},
            {"HSS13.375X0.500", {6.6787, 0.4664}},
            {"HSS13.375X0.375", {6.6831, 0.3495}},
            {"HSS13.375X0.313", {6.7003, 0.2893}},
            {"HSS13.375X0.250", {6.6809, 0.233}},
            {"HSS13.375X0.188", {6.6984, 0.174}},
            {"HSS12.750X0.750", {6.3774, 0.6953}},
            {"HSS12.750X0.625", {6.3755, 0.5807}},
            {"HSS12.750X0.500", {6.3809, 0.4636}},
            {"HSS12.750X0.375", {6.3797, 0.3487}},
            {"HSS12.750X0.250", {6.3855, 0.2323}},
            {"HSS12.750X0.188", {6.3686, 0.1747}},
            {"HSS12.000X0.625", {5.9971, 0.5827}},
            {"HSS12.000X0.500", {5.9967, 0.4661}},
            {"HSS12.000X0.375", {5.9959, 0.3498}},
            {"HSS12.000X0.250", {6.0023, 0.2324}},
            {"HSS11.750X0.625", {5.8754, 0.5814}},
            {"HSS11.750X0.500", {5.874, 0.4656}},
            {"HSS11.750X0.375", {5.8717, 0.349}},
            {"HSS11.750X0.337", {5.8629, 0.315}},
            {"HSS11.750X0.250", {5.8871, 0.2315}},
            {"HSS10.750X0.625", {5.3739, 0.5822}},
            {"HSS10.750X0.500", {5.3775, 0.4641}},
            {"HSS10.750X0.375", {5.3682, 0.3498}},
            {"HSS10.750X0.313", {5.3817, 0.29}},
            {"HSS10.750X0.250", {5.3606, 0.2341}},
            {"HSS10.750X0.188", {5.3758, 0.1743}},
            {"HSS10.000X0.625", {4.9944, 0.5822}},
            {"HSS10.000X0.500", {5.0091, 0.463}},
            {"HSS10.000X0.375", {4.9895, 0.3504}},
            {"HSS10.000X0.312", {5.0081, 0.2901}},
            {"HSS10.000X0.250", {5.0007, 0.2328}},
            {"HSS10.000X0.188", {4.9986, 0.174}},
            {"HSS9.625X0.625", {4.8204, 0.5792}},
            {"HSS9.625X0.500", {4.8158, 0.4645}},
            {"HSS9.625X0.375", {4.8185, 0.3485}},
            {"HSS9.625X0.312", {4.8118, 0.2912}},
            {"HSS9.625X0.250", {4.8143, 0.233}},
            {"HSS9.625X0.188", {4.8127, 0.1738}},
            {"HSS8.625X0.625", {4.3039, 0.5838}},
            {"HSS8.625X0.500", {4.325, 0.4623}},
            {"HSS8.625X0.375", {4.3128, 0.3487}},
            {"HSS8.625X0.322", {4.313, 0.3}},
            {"HSS8.625X0.250", {4.3129, 0.2328}},
            {"HSS8.625X0.188", {4.316, 0.1736}},
            {"HSS7.500X0.500", {3.7491, 0.4652}},
            {"HSS7.500X0.375", {3.7479, 0.3497}},
            {"HSS7.500X0.312", {3.7517, 0.2904}},
            {"HSS7.500X0.250", {3.753, 0.2326}},
            {"HSS7.500X0.188", {3.7531, 0.1737}},
            {"HSS7.000X0.500", {3.4988, 0.4655}},
            {"HSS7.000X0.375", {3.4979, 0.3499}},
            {"HSS7.000X0.312", {3.5023, 0.2905}},
            {"HSS7.000X0.250", {3.5008, 0.2331}},
            {"HSS7.000X0.188", {3.4964, 0.1743}},
            {"HSS7.000X0.125", {3.5036, 0.1159}},
            {"HSS6.875X0.375", {3.4362, 0.3496}},
            {"HSS6.875X0.312", {3.4391, 0.2905}},
            {"HSS6.875X0.250", {3.4344, 0.2335}},
            {"HSS6.875X0.188", {3.4412, 0.1736}},
            {"HSS6.625X0.500", {3.311, 0.4657}},
            {"HSS6.625X0.432", {3.3121, 0.4023}},
            {"HSS6.625X0.375", {3.3125, 0.3495}},
            {"HSS6.625X0.312", {3.3121, 0.2913}},
            {"HSS6.625X0.280", {3.3154, 0.2591}},
            {"HSS6.625X0.250", {3.3104, 0.2332}},
            {"HSS6.625X0.188", {3.3155, 0.1738}},
            {"HSS6.625X0.125", {3.3181, 0.1157}},
            {"HSS6.000X0.500", {2.9999, 0.4655}},
            {"HSS6.000X0.375", {2.9967, 0.3505}},
            {"HSS6.000X0.312", {2.9986, 0.2912}},
            {"HSS6.000X0.280", {2.9958, 0.2606}},
            {"HSS6.000X0.250", {3.0024, 0.2326}},
            {"HSS6.000X0.188", {2.9987, 0.174}},
            {"HSS6.000X0.125", {3.0015, 0.1158}},
            {"HSS5.563X0.500", {2.7818, 0.4646}},
            {"HSS5.563X0.375", {2.7801, 0.3494}},
            {"HSS5.563X0.258", {2.7782, 0.2403}},
            {"HSS5.563X0.188", {2.7796, 0.1742}},
            {"HSS5.563X0.134", {2.781, 0.1241}},
            {"HSS5.500X0.500", {2.7495, 0.4649}},
            {"HSS5.500X0.375", {2.7485, 0.3493}},
            {"HSS5.500X0.258", {2.7446, 0.2408}},
            {"HSS5.000X0.500", {2.4998, 0.4649}},
            {"HSS5.000X0.375", {2.5029, 0.3482}},
            {"HSS5.000X0.312", {2.5033, 0.2901}},
            {"HSS5.000X0.258", {2.501, 0.2398}},
            {"HSS5.000X0.250", {2.5004, 0.233}},
            {"HSS5.000X0.188", {2.4996, 0.174}},
            {"HSS5.000X0.125", {2.4997, 0.1161}},
            {"HSS4.500X0.375", {2.2499, 0.349}},
            {"HSS4.500X0.337", {2.2491, 0.3134}},
            {"HSS4.500X0.237", {2.2495, 0.2201}},
            {"HSS4.500X0.188", {2.2514, 0.1737}},
            {"HSS4.500X0.125", {2.2484, 0.1162}},
            {"HSS4.000X0.313", {2.0008, 0.2907}},
            {"HSS4.000X0.250", {1.9994, 0.2332}},
            {"HSS4.000X0.237", {2.0001, 0.22}},
            {"HSS4.000X0.226", {1.9994, 0.2101}},
            {"HSS4.000X0.220", {2.0002, 0.205}},
            {"HSS4.000X0.188", {1.9991, 0.1741}},
            {"HSS4.000X0.125", {1.997, 0.1165}},
            {"HSS3.500X0.313", {1.7515, 0.2901}},
            {"HSS3.500X0.300", {1.751, 0.2784}},
            {"HSS3.500X0.250", {1.7513, 0.2326}},
            {"HSS3.500X0.216", {1.7496, 0.2009}},
            {"HSS3.500X0.203", {1.7482, 0.1893}},
            {"HSS3.500X0.188", {1.7487, 0.1745}},
            {"HSS3.500X0.125", {1.7529, 0.1156}},
            {"HSS3.000X0.250", {1.498, 0.234}},
            {"HSS3.000X0.216", {1.4992, 0.2015}},
            {"HSS3.000X0.203", {1.5011, 0.1891}},
            {"HSS3.000X0.188", {1.5027, 0.1731}},
            {"HSS3.000X0.152", {1.5007, 0.1409}},
            {"HSS3.000X0.134", {1.4992, 0.1243}},
            {"HSS3.000X0.125", {1.4974, 0.1162}},
            {"HSS2.875X0.250", {1.4383, 0.2325}},
            {"HSS2.875X0.203", {1.4408, 0.1882}},
            {"HSS2.875X0.188", {1.4354, 0.1746}},
            {"HSS2.875X0.125", {1.4347, 0.1167}},
            {"HSS2.500X0.250", {1.2514, 0.2324}},
            {"HSS2.500X0.188", {1.2507, 0.1737}},
            {"HSS2.500X0.125", {1.2502, 0.116}},
            {"HSS2.375X0.250", {1.1872, 0.2329}},
            {"HSS2.375X0.218", {1.1861, 0.2038}},
            {"HSS2.375X0.188", {1.1886, 0.1735}},
            {"HSS2.375X0.154", {1.1888, 0.1424}},
            {"HSS2.375X0.125", {1.1881, 0.1159}},
            {"HSS1.900X0.188", {0.9502, 0.1738}},
            {"HSS1.900X0.145", {0.9494, 0.1352}},
            {"HSS1.900X0.120", {0.9506, 0.1109}},
            {"HSS1.660X0.140", {0.8297, 0.1301}},
        };

        return registry.fetch(type);
    }

    template<bool adjusted> void new_eu2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string& type) {
        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto scale = 1.;
        if(!command.eof() && !get_input(command, scale)) {
            suanpan_error("A valid scale is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        const auto dim = adjusted ? euisection_adjusted(type) : euisection(type);

        if(dim.is_empty()) {
            suanpan_error("Cannot identify section type.\n");
            return;
        }

        return_obj = std::make_unique<ISection2D>(tag, scale * dim, material_id, int_pt, eccentricity);
    }

    template<bool adjusted> void new_eu3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string& type) {
        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto scale = 1.;
        if(!command.eof() && !get_input(command, scale)) {
            suanpan_error("A valid scale is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_y = 0.;
        if(!command.eof() && !get_input(command, eccentricity_y)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        auto eccentricity_z = 0.;
        if(!command.eof() && !get_input(command, eccentricity_z)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        const auto dim = adjusted ? euisection_adjusted(type) : euisection(type);

        if(dim.is_empty()) {
            suanpan_error("Cannot identify section type.\n");
            return;
        }

        return_obj = std::make_unique<ISection3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
    }

    void new_nz2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string& type) {
        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto scale = 1.;
        if(!command.eof() && !get_input(command, scale)) {
            suanpan_error("A valid scale is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        auto dim = nzisection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<ISection2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = nzchsection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<CircularHollow2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = nzrhsection(type);

        if(!dim.is_empty()) {
            dim[0] -= dim[2]; // account for centroid
            dim[1] -= dim[2]; // account for centroid
            return_obj = std::make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = nzshsection(type);

        if(!dim.is_empty()) {
            dim[0] -= dim[2]; // account for centroid
            dim[1] -= dim[2]; // account for centroid
            return_obj = std::make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        suanpan_error("Cannot identify section type.\n");
    }

    void new_nz3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string& type) {
        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto scale = 1.;
        if(!command.eof() && !get_input(command, scale)) {
            suanpan_error("A valid scale is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_y = 0.;
        if(!command.eof() && !get_input(command, eccentricity_y)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        auto eccentricity_z = 0.;
        if(!command.eof() && !get_input(command, eccentricity_z)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        auto dim = nzisection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<ISection3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = nzchsection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<CircularHollow3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = nzrhsection(type);

        if(!dim.is_empty()) {
            dim[0] -= dim[2]; // account for centroid
            dim[1] -= dim[2]; // account for centroid
            return_obj = std::make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = nzshsection(type);

        if(!dim.is_empty()) {
            dim[0] -= dim[2]; // account for centroid
            dim[1] -= dim[2]; // account for centroid
            return_obj = std::make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        suanpan_error("Cannot identify section type.\n");
    }

    template<bool recenter, bool adjusted> void new_us2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string& type) {
        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto scale = 1.;
        if(!command.eof() && !get_input(command, scale)) {
            suanpan_error("A valid scale is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity = 0.;
        if(!recenter && !command.eof() && !get_input(command, eccentricity)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        auto dim = adjusted ? usisection_adjusted(type) : usisection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<ISection2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = adjusted ? usrhssection_adjusted(type) : usrhssection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = adjusted ? uschssection_adjusted(type) : uschssection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<CircularHollow2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = adjusted ? ustsection_adjusted(type) : ustsection(type);

        if(!dim.is_empty()) {
            if constexpr(recenter) eccentricity = barycenter(dim *= scale);
            return_obj = std::make_unique<TSection2D>(tag, std::move(dim), material_id, int_pt, eccentricity);
            return;
        }

        suanpan_error("Cannot identify section type.\n");
    }

    template<bool recenter, bool adjusted> void new_us3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned tag, const std::string& type) {
        unsigned material_id;
        if(!get_input(command, material_id)) {
            suanpan_error("A valid material tag is required.\n");
            return;
        }

        auto scale = 1.;
        if(!command.eof() && !get_input(command, scale)) {
            suanpan_error("A valid scale is required.\n");
            return;
        }

        auto int_pt = 6u;
        if(!command.eof() && !get_input(command, int_pt)) {
            suanpan_error("A valid number of integration points is required.\n");
            return;
        }

        auto eccentricity_y = 0.;
        if(!recenter && !command.eof() && !get_input(command, eccentricity_y)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        auto eccentricity_z = 0.;
        if(!recenter && !command.eof() && !get_input(command, eccentricity_z)) {
            suanpan_error("A valid eccentricity is required.\n");
            return;
        }

        auto dim = adjusted ? usisection_adjusted(type) : usisection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<ISection3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = adjusted ? usrhssection_adjusted(type) : usrhssection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = adjusted ? uschssection_adjusted(type) : uschssection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<CircularHollow3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = adjusted ? ustsection_adjusted(type) : ustsection(type);

        if(!dim.is_empty()) {
            if constexpr(recenter) eccentricity_y = barycenter(dim *= scale);
            return_obj = std::make_unique<TSection3D>(tag, std::move(dim), material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        suanpan_error("Cannot identify section type.\n");
    }
} // namespace

int create_new_section(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    static const TokenRegistry<void (*)(unique_ptr<Section>&, std::istringstream&, unsigned, const std::string&)> registry{
        {"Box2D", new_box2d},
        {"Box3D", new_box3d},
        {"Cell2D", new_cell2d},
        {"Bar2D", new_cell2d},
        {"Cell3D", new_cell3d},
        {"Bar3D", new_cell3d},
        {"Cell3DOS", new_cell3dos},
        {"Circle1D", new_circle1d},
        {"Circle2D", new_circle2d},
        {"Circle3D", new_circle3d},
        {"CircularHollow2D", new_circularhollow2D},
        {"CircularHollow3D", new_circularhollow3D},
        {"Fibre1D", new_fibre1d},
        {"Fibre2D", new_fibre2d},
        {"Fibre3D", new_fibre3d<false>},
        {"Fibre3DOS", new_fibre3d<true>},
        {"HSection2D", new_hsection2d},
        {"ISection2D", new_isection2d<false>},
        {"ISection3D", new_isection3d<false>},
        {"ISection2DC", new_isection2d<true>},
        {"ISection3DC", new_isection3d<true>},
        {"Rectangle1D", new_rectangle1d},
        {"Rectangle2D", new_rectangle2d},
        {"Rectangle3D", new_rectangle3d},
        {"TrussSection", new_trusssection},
        {"TSection2D", new_tsection2d<false>},
        {"TSection3D", new_tsection3d<false>},
        {"TSection2DC", new_tsection2d<true>},
        {"TSection3DC", new_tsection3d<true>},
        {"NM2D1", new_nm2d<3>},
        {"NM2D2", new_nm2d<8>},
        {"NM2D3", new_nm2d<11>},
        {"NM3D1", new_nm3d<4>},
        {"NM3D2", new_nm3d<10>},
        {"NM3D3", new_nm3d<13>},
        {"NM2D3K", new_nmk<13>},
        {"NM3D3K", new_nmk<17>},
        {"EU2D", new_eu2d<false>},
        {"EU3D", new_eu3d<false>},
        {"NZ2D", new_nz2d},
        {"NZ3D", new_nz3d},
        // adjusted variants, optimized against the EN 10365 section properties
        {"EU2DA", new_eu2d<true>},
        {"EU3DA", new_eu3d<true>},
        {"US2D", new_us2d<false, false>},
        {"US3D", new_us3d<false, false>},
        {"US2DC", new_us2d<true, false>},
        {"US3DC", new_us3d<true, false>},
        // adjusted variants, optimized against the AISC section properties
        {"US2DA", new_us2d<false, true>},
        {"US3DA", new_us3d<false, true>},
        {"US2DAC", new_us2d<true, true>},
        {"US3DAC", new_us3d<true, true>},
    };

    std::string section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section type is required.\n");
        return SUANPAN_SUCCESS;
    }

    // code-based sections take an additional designation before the tag
    std::string designation;
    if(is_equal_any(section_id, "EU2D", "EU3D", "EU2DA", "EU3DA", "NZ2D", "NZ3D", "US2D", "US3D", "US2DC", "US3DC", "US2DA", "US3DA", "US2DAC", "US3DAC") && !get_input(command, designation)) {
        suanpan_error("A valid designation is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Section> new_section = nullptr;

    if(const auto factory = registry.fetch(section_id); nullptr != factory) factory(new_section, command, tag, designation);
    else external_module::object(new_section, domain, section_id, command, tag);

    if(new_section == nullptr || !domain->insert(std::move(new_section))) suanpan_error("Fail to create new section via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
