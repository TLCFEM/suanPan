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
#include <Toolbox/utility.h>

namespace {
    void new_cell2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_cell3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_box2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_box3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_cell3dos(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_circle1d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_circle2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_circle3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_circularhollow2D(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_circularhollow3D(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_fibre1d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Fibre1D>(tag, get_remaining<uword>(command));
    }

    void new_fibre2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Fibre2D>(tag, get_remaining<uword>(command));
    }

    template<bool if_os> void new_fibre3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        if constexpr(if_os) return_obj = std::make_unique<Fibre3DOS>(tag, get_remaining<uword>(command));
        else return_obj = std::make_unique<Fibre3D>(tag, get_remaining<uword>(command));
    }

    void new_hsection2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    template<bool recenter> void new_isection2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    template<bool recenter> void new_isection3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_rectangle1d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_rectangle2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_rectangle3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_trusssection(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    template<bool recenter> void new_tsection2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    template<bool recenter> void new_tsection3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    template<unsigned size> void new_nm2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    template<unsigned size> void new_nm3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    template<unsigned size> void new_nmk(unique_ptr<Section>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_eu2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        std::string type;
        if(!get_input(command, type)) {
            suanpan_error("A valid designation is required.\n");
            return;
        }

        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

        const auto dim = euisection(type);

        if(dim.is_empty()) {
            suanpan_error("Cannot identify section type.\n");
            return;
        }

        return_obj = std::make_unique<ISection2D>(tag, scale * dim, material_id, int_pt, eccentricity);
    }

    void new_eu3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        std::string type;
        if(!get_input(command, type)) {
            suanpan_error("A valid designation is required.\n");
            return;
        }

        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

        const auto dim = euisection(type);

        if(dim.is_empty()) {
            suanpan_error("Cannot identify section type.\n");
            return;
        }

        return_obj = std::make_unique<ISection3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
    }

    void new_nz2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        std::string type;
        if(!get_input(command, type)) {
            suanpan_error("A valid designation is required.\n");
            return;
        }

        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    void new_nz3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        std::string type;
        if(!get_input(command, type)) {
            suanpan_error("A valid designation is required.\n");
            return;
        }

        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

    template<bool recenter> void new_us2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        std::string type;
        if(!get_input(command, type)) {
            suanpan_error("A valid designation is required.\n");
            return;
        }

        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

        auto dim = usisection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<ISection2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = usrhssection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = uschssection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<CircularHollow2D>(tag, scale * dim, material_id, int_pt, eccentricity);
            return;
        }

        dim = ustsection(type);

        if(!dim.is_empty()) {
            if constexpr(recenter) eccentricity = barycenter(dim *= scale);
            return_obj = std::make_unique<TSection2D>(tag, std::move(dim), material_id, int_pt, eccentricity);
            return;
        }

        suanpan_error("Cannot identify section type.\n");
    }

    template<bool recenter> void new_us3d(unique_ptr<Section>& return_obj, std::istringstream& command) {
        std::string type;
        if(!get_input(command, type)) {
            suanpan_error("A valid designation is required.\n");
            return;
        }

        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

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

        auto dim = usisection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<ISection3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = usrhssection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = uschssection(type);

        if(!dim.is_empty()) {
            return_obj = std::make_unique<CircularHollow3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        dim = ustsection(type);

        if(!dim.is_empty()) {
            if constexpr(recenter) eccentricity_y = barycenter(dim *= scale);
            return_obj = std::make_unique<TSection3D>(tag, std::move(dim), material_id, int_pt, vec{eccentricity_y, eccentricity_z});
            return;
        }

        suanpan_error("Cannot identify section type.\n");
    }
} // namespace

int create_new_section(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    static const TokenRegistry<void (*)(unique_ptr<Section>&, std::istringstream&)> registry{
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
        {"EU2D", new_eu2d},
        {"EU3D", new_eu3d},
        {"NZ2D", new_nz2d},
        {"NZ3D", new_nz3d},
        {"US2D", new_us2d<false>},
        {"US3D", new_us3d<false>},
        {"US2DC", new_us2d<true>},
        {"US3DC", new_us3d<true>},
    };

    std::string section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Section> new_section = nullptr;

    if(const auto factory = registry.fetch(section_id); nullptr != factory) factory(new_section, command);
    else external_module::object(new_section, domain, section_id, command);

    if(new_section == nullptr || !domain->insert(std::move(new_section))) suanpan_error("Fail to create new section via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
