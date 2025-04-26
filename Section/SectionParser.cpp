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

#include "SectionParser.h"

#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Section/Section>
#include <Toolbox/utility.h>

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

    return_obj = make_unique<Cell2D>(tag, area, material_id, eccentricity);
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

    return_obj = make_unique<Cell3D>(tag, area, material_id, eccentricity_a, eccentricity_b);
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

    return_obj = make_unique<Box2D>(tag, width, height, thickness, material_id, int_pt, eccentricity);
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

    return_obj = make_unique<Box3D>(tag, width, height, thickness, material_id, int_pt, eccentricity_a, eccentricity_b);
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

    return_obj = make_unique<Cell3DOS>(tag, area, omega, py, pz, material_id, eccentricity_a, eccentricity_b);
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

    return_obj = make_unique<Circle1D>(tag, radius, material_id);
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

    return_obj = make_unique<Circle2D>(tag, radius, material_id, int_pt, eccentricity);
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

    return_obj = make_unique<Circle3D>(tag, radius, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
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

    return_obj = make_unique<CircularHollow2D>(tag, radius, thickness, material_id, int_pt, eccentricity);
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

    return_obj = make_unique<CircularHollow3D>(tag, radius, thickness, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
}

void new_fibre1d(unique_ptr<Section>& return_obj, std::istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    return_obj = make_unique<Fibre1D>(tag, get_remaining<uword>(command));
}

void new_fibre2d(unique_ptr<Section>& return_obj, std::istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    return_obj = make_unique<Fibre2D>(tag, get_remaining<uword>(command));
}

void new_fibre3d(unique_ptr<Section>& return_obj, std::istringstream& command, const bool if_os) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    if(if_os) return_obj = make_unique<Fibre3DOS>(tag, get_remaining<uword>(command));
    else return_obj = make_unique<Fibre3D>(tag, get_remaining<uword>(command));
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

    return_obj = make_unique<HSection2D>(tag, dim(0), dim(1), dim(2), dim(3), dim(4), dim(5), material_id, int_pt, eccentricity);
}

double barycenter(const vec& dim) {
    if(4llu == dim.n_elem) {
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

void new_isection2d(unique_ptr<Section>& return_obj, std::istringstream& command, const bool recenter) {
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
    if(recenter) eccentricity = barycenter(dim);
    else if(!command.eof() && !get_input(command, eccentricity)) {
        suanpan_error("A valid eccentricity is required.\n");
        return;
    }

    return_obj = make_unique<ISection2D>(tag, dim(0), dim(1), dim(2), dim(3), dim(4), dim(5), material_id, int_pt, eccentricity);
}

void new_isection3d(unique_ptr<Section>& return_obj, std::istringstream& command, const bool recenter) {
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
    if(recenter) eccentricity_y = barycenter(dim);
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

    return_obj = make_unique<ISection3D>(tag, std::move(dim), material_id, int_pt, vec{eccentricity_y, eccentricity_z});
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

    return_obj = make_unique<Rectangle1D>(tag, width, height, material_id);
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

    return_obj = make_unique<Rectangle2D>(tag, width, height, material_id, int_pt, eccentricity);
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

    return_obj = make_unique<Rectangle3D>(tag, width, height, material_id, int_pt, eccentricity_a, eccentricity_b);
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

    return_obj = make_unique<TrussSection>(tag, area, material_id);
}

void new_tsection2d(unique_ptr<Section>& return_obj, std::istringstream& command, const bool recenter) {
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
    if(recenter) eccentricity = barycenter(dim);
    else if(!command.eof() && !get_input(command, eccentricity)) {
        suanpan_error("A valid eccentricity is required.\n");
        return;
    }

    return_obj = make_unique<TSection2D>(tag, dim(0), dim(1), dim(2), dim(3), material_id, int_pt, eccentricity);
}

void new_tsection3d(unique_ptr<Section>& return_obj, std::istringstream& command, const bool recenter) {
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
    if(recenter) eccentricity_y = barycenter(dim);
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

    return_obj = make_unique<TSection3D>(tag, std::move(dim), material_id, int_pt, vec{eccentricity_y, eccentricity_z});
}

void new_nm2d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned size) {
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

    if(3u == size) {
        return_obj = make_unique<NM2D1>(tag, P(0), P(1), P(2));
        return;
    }

    const auto para_set = get_remaining<double>(command);

    if(para_set.size() % 3 != 0) {
        suanpan_error("A valid parameter set is required.\n");
        return;
    }

    mat poly_set(para_set);
    poly_set.reshape(3, poly_set.n_elem / 3);
    inplace_trans(poly_set);

    if(8 == size) return_obj = make_unique<NM2D2>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), std::move(poly_set));
    else if(11 == size) return_obj = make_unique<NM2D3>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), vec{P(8), P(8)}, vec{P(9), P(9)}, P(10), std::move(poly_set));
}

void new_nm3d(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned size) {
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

    if(4u == size) {
        return_obj = make_unique<NM3D1>(tag, P(0), P(1), P(2), P(3));
        return;
    }

    const auto para_set = get_remaining<double>(command);

    if(para_set.size() % 4 != 0) {
        suanpan_error("A valid parameter set is required.\n");
        return;
    }

    mat poly_set(para_set);
    poly_set.reshape(4, poly_set.n_elem / 4);
    inplace_trans(poly_set);

    if(10 == size) return_obj = make_unique<NM3D2>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), std::move(poly_set));
    else if(13 == size) return_obj = make_unique<NM3D3>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), vec{P(10), P(10), P(10)}, vec{P(11), P(11), P(11)}, P(12), std::move(poly_set));
}

void new_nmk(unique_ptr<Section>& return_obj, std::istringstream& command, const unsigned size) {
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

    const auto para_set = get_remaining<double>(command);

    const auto p_size = 13u == size ? 3 : 4;

    if(para_set.size() % p_size != 0) {
        suanpan_error("A valid parameter set is required.\n");
        return;
    }

    mat poly_set(para_set);
    poly_set.reshape(p_size, poly_set.n_elem / p_size);
    inplace_trans(poly_set);

    if(13 == size) return_obj = make_unique<NM2D3>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), vec{P(8), P(9)}, vec{P(10), P(11)}, P(12), std::move(poly_set));
    else if(17 == size) return_obj = make_unique<NM3D3>(tag, P(0), P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), vec{P(10), P(11), P(12)}, vec{P(13), P(14), P(15)}, P(16), std::move(poly_set));
}

vec euisection(const std::string_view type) {
    if(is_equal(type, "HEA100")) return {100.0, 8.0, 100.0, 8.0, 80.0, 5.0};
    if(is_equal(type, "HEA120")) return {120.0, 8.0, 120.0, 8.0, 98.0, 5.0};
    if(is_equal(type, "HEA140")) return {140.0, 8.5, 140.0, 8.5, 116.0, 5.5};
    if(is_equal(type, "HEA160")) return {160.0, 9.0, 160.0, 9.0, 134.0, 6.0};
    if(is_equal(type, "HEA180")) return {180.0, 9.5, 180.0, 9.5, 152.0, 6.0};
    if(is_equal(type, "HEA200")) return {200.0, 10.0, 200.0, 10.0, 170.0, 6.5};
    if(is_equal(type, "HEA220")) return {220.0, 11.0, 220.0, 11.0, 188.0, 7.0};
    if(is_equal(type, "HEA240")) return {240.0, 12.0, 240.0, 12.0, 206.0, 7.5};
    if(is_equal(type, "HEA260")) return {260.0, 12.5, 260.0, 12.5, 225.0, 7.5};
    if(is_equal(type, "HEA280")) return {280.0, 13.0, 280.0, 13.0, 244.0, 8.0};
    if(is_equal(type, "HEA300")) return {300.0, 14.0, 300.0, 14.0, 262.0, 8.5};
    if(is_equal(type, "HEA320")) return {300.0, 15.5, 300.0, 15.5, 279.0, 9.0};
    if(is_equal(type, "HEA340")) return {300.0, 16.5, 300.0, 16.5, 297.0, 9.5};
    if(is_equal(type, "HEA360")) return {300.0, 17.5, 300.0, 17.5, 315.0, 10.0};
    if(is_equal(type, "HEA400")) return {300.0, 19.0, 300.0, 19.0, 352.0, 11.0};
    if(is_equal(type, "HEA450")) return {300.0, 21.0, 300.0, 21.0, 398.0, 11.5};
    if(is_equal(type, "HEA500")) return {300.0, 23.0, 300.0, 23.0, 444.0, 12.0};
    if(is_equal(type, "HEA550")) return {300.0, 24.0, 300.0, 24.0, 492.0, 12.5};
    if(is_equal(type, "HEA600")) return {300.0, 25.0, 300.0, 25.0, 540.0, 13.0};
    if(is_equal(type, "HEA650")) return {300.0, 26.0, 300.0, 26.0, 588.0, 13.5};
    if(is_equal(type, "HEA700")) return {300.0, 27.0, 300.0, 27.0, 636.0, 14.5};
    if(is_equal(type, "HEA800")) return {300.0, 28.0, 300.0, 28.0, 734.0, 15.0};
    if(is_equal(type, "HEA900")) return {300.0, 30.0, 300.0, 30.0, 830.0, 16.0};
    if(is_equal(type, "HEA1000")) return {300.0, 31.0, 300.0, 31.0, 928.0, 16.5};
    if(is_equal(type, "HEB100")) return {100.0, 10.0, 100.0, 10.0, 80.0, 6.0};
    if(is_equal(type, "HEB120")) return {120.0, 11.0, 120.0, 11.0, 98.0, 6.5};
    if(is_equal(type, "HEB140")) return {140.0, 12.0, 140.0, 12.0, 116.0, 7.0};
    if(is_equal(type, "HEB160")) return {160.0, 13.0, 160.0, 13.0, 134.0, 8.0};
    if(is_equal(type, "HEB180")) return {180.0, 14.0, 180.0, 14.0, 152.0, 8.5};
    if(is_equal(type, "HEB200")) return {200.0, 15.0, 200.0, 15.0, 170.0, 9.0};
    if(is_equal(type, "HEB220")) return {220.0, 16.0, 220.0, 16.0, 188.0, 9.5};
    if(is_equal(type, "HEB240")) return {240.0, 17.0, 240.0, 17.0, 206.0, 10.0};
    if(is_equal(type, "HEB260")) return {260.0, 17.5, 260.0, 17.5, 225.0, 10.0};
    if(is_equal(type, "HEB280")) return {280.0, 18.0, 280.0, 18.0, 244.0, 10.5};
    if(is_equal(type, "HEB300")) return {300.0, 19.0, 300.0, 19.0, 262.0, 11.0};
    if(is_equal(type, "HEB320")) return {300.0, 20.5, 300.0, 20.5, 279.0, 11.5};
    if(is_equal(type, "HEB340")) return {300.0, 21.5, 300.0, 21.5, 297.0, 12.0};
    if(is_equal(type, "HEB360")) return {300.0, 22.5, 300.0, 22.5, 315.0, 12.5};
    if(is_equal(type, "HEB400")) return {300.0, 24.0, 300.0, 24.0, 352.0, 13.5};
    if(is_equal(type, "HEB450")) return {300.0, 26.0, 300.0, 26.0, 398.0, 14.0};
    if(is_equal(type, "HEB500")) return {300.0, 28.0, 300.0, 28.0, 444.0, 14.5};
    if(is_equal(type, "HEB550")) return {300.0, 29.0, 300.0, 29.0, 492.0, 15.0};
    if(is_equal(type, "HEB600")) return {300.0, 30.0, 300.0, 30.0, 540.0, 15.5};
    if(is_equal(type, "HEB650")) return {300.0, 31.0, 300.0, 31.0, 588.0, 16.0};
    if(is_equal(type, "HEB700")) return {300.0, 32.0, 300.0, 32.0, 636.0, 17.0};
    if(is_equal(type, "HEB800")) return {300.0, 33.0, 300.0, 33.0, 734.0, 17.5};
    if(is_equal(type, "HEB900")) return {300.0, 35.0, 300.0, 35.0, 830.0, 18.5};
    if(is_equal(type, "HEB1000")) return {300.0, 36.0, 300.0, 36.0, 928.0, 19.0};
    if(is_equal(type, "IPE80")) return {46.0, 5.2, 46.0, 5.2, 69.6, 3.8};
    if(is_equal(type, "IPE100")) return {55.0, 5.7, 55.0, 5.7, 88.6, 4.1};
    if(is_equal(type, "IPE120")) return {64.0, 6.3, 64.0, 6.3, 107.4, 4.4};
    if(is_equal(type, "IPE140")) return {73.0, 6.9, 73.0, 6.9, 126.2, 4.7};
    if(is_equal(type, "IPE160")) return {82.0, 7.4, 82.0, 7.4, 145.2, 5.0};
    if(is_equal(type, "IPE180")) return {91.0, 8.0, 91.0, 8.0, 164.0, 5.3};
    if(is_equal(type, "IPE200")) return {100.0, 8.5, 100.0, 8.5, 183.0, 5.6};
    if(is_equal(type, "IPE220")) return {110.0, 9.2, 110.0, 9.2, 201.6, 5.9};
    if(is_equal(type, "IPE240")) return {120.0, 9.8, 120.0, 9.8, 220.4, 6.2};
    if(is_equal(type, "IPE270")) return {135.0, 10.2, 135.0, 10.2, 249.6, 6.6};
    if(is_equal(type, "IPE300")) return {150.0, 10.7, 150.0, 10.7, 278.6, 7.1};
    if(is_equal(type, "IPE330")) return {160.0, 11.5, 160.0, 11.5, 307.0, 7.5};
    if(is_equal(type, "IPE360")) return {170.0, 12.7, 170.0, 12.7, 334.6, 8.0};
    if(is_equal(type, "IPE400")) return {180.0, 13.5, 180.0, 13.5, 373.0, 8.6};
    if(is_equal(type, "IPE450")) return {190.0, 14.6, 190.0, 14.6, 420.8, 9.4};
    if(is_equal(type, "IPE500")) return {200.0, 16.0, 200.0, 16.0, 468.0, 10.2};
    if(is_equal(type, "IPE550")) return {210.0, 17.2, 210.0, 17.2, 515.6, 11.1};
    if(is_equal(type, "IPE600")) return {220.0, 19.0, 220.0, 19.0, 562.0, 12.0};
    if(is_equal(type, "HEM100")) return {106.0, 20.0, 106.0, 20.0, 80.0, 12.0};
    if(is_equal(type, "HEM120")) return {126.0, 21.0, 126.0, 21.0, 98.0, 12.5};
    if(is_equal(type, "HEM140")) return {146.0, 22.0, 146.0, 22.0, 116.0, 13.0};
    if(is_equal(type, "HEM160")) return {166.0, 23.0, 166.0, 23.0, 134.0, 14.0};
    if(is_equal(type, "HEM180")) return {186.0, 24.0, 186.0, 24.0, 152.0, 14.5};
    if(is_equal(type, "HEM200")) return {206.0, 25.0, 206.0, 25.0, 170.0, 15.0};
    if(is_equal(type, "HEM220")) return {226.0, 26.0, 226.0, 26.0, 188.0, 15.5};
    if(is_equal(type, "HEM240")) return {248.0, 32.0, 248.0, 32.0, 206.0, 18.0};
    if(is_equal(type, "HEM260")) return {268.0, 32.5, 268.0, 32.5, 225.0, 18.0};
    if(is_equal(type, "HEM280")) return {288.0, 33.0, 288.0, 33.0, 244.0, 18.5};
    if(is_equal(type, "HEM300")) return {310.0, 39.0, 310.0, 39.0, 262.0, 21.0};
    if(is_equal(type, "HEM320")) return {309.0, 40.0, 309.0, 40.0, 279.0, 21.0};
    if(is_equal(type, "HEM340")) return {309.0, 40.0, 309.0, 40.0, 297.0, 21.0};
    if(is_equal(type, "HEM360")) return {308.0, 40.0, 308.0, 40.0, 315.0, 21.0};
    if(is_equal(type, "HEM400")) return {307.0, 40.0, 307.0, 40.0, 352.0, 21.0};
    if(is_equal(type, "HEM450")) return {307.0, 40.0, 307.0, 40.0, 398.0, 21.0};
    if(is_equal(type, "HEM500")) return {306.0, 40.0, 306.0, 40.0, 444.0, 21.0};
    if(is_equal(type, "HEM550")) return {306.0, 40.0, 306.0, 40.0, 492.0, 21.0};
    if(is_equal(type, "HEM600")) return {305.0, 40.0, 305.0, 40.0, 540.0, 21.0};
    if(is_equal(type, "HEM650")) return {305.0, 40.0, 305.0, 40.0, 588.0, 21.0};
    if(is_equal(type, "HEM700")) return {304.0, 40.0, 304.0, 40.0, 636.0, 21.0};
    if(is_equal(type, "HEM800")) return {303.0, 40.0, 303.0, 40.0, 734.0, 21.0};
    if(is_equal(type, "HEM900")) return {302.0, 40.0, 302.0, 40.0, 830.0, 21.0};
    if(is_equal(type, "HEM1000")) return {302.0, 40.0, 302.0, 40.0, 928.0, 21.0};
    return {};
}

vec nzchsection(const std::string_view type) {
    if(is_equal(type, "610.0X12.7CHS")) return {305.00, 12.7};
    if(is_equal(type, "610.0X9.5CHS")) return {305.00, 9.5};
    if(is_equal(type, "610.0X6.4CHS")) return {305.00, 6.4};
    if(is_equal(type, "508.0X12.7CHS")) return {254.00, 12.7};
    if(is_equal(type, "508.0X9.5CHS")) return {254.00, 9.5};
    if(is_equal(type, "508.0X6.4CHS")) return {254.00, 6.4};
    if(is_equal(type, "165.1X5.4CHS")) return {82.55, 5.4};
    if(is_equal(type, "165.1X5.0CHS")) return {82.55, 5.0};
    if(is_equal(type, "139.7X5.4CHS")) return {69.85, 5.4};
    if(is_equal(type, "139.7X5.0CHS")) return {69.85, 5.0};
    if(is_equal(type, "114.3X5.4CHS")) return {57.15, 5.4};
    if(is_equal(type, "114.3X4.5CHS")) return {57.15, 4.5};
    if(is_equal(type, "101.6X5.0CHS")) return {50.80, 5.0};
    if(is_equal(type, "101.6X4.0CHS")) return {50.80, 4.0};
    if(is_equal(type, "88.9X5.9CHS")) return {44.45, 5.9};
    if(is_equal(type, "88.9X5.0CHS")) return {44.45, 5.0};
    if(is_equal(type, "88.9X4.0CHS")) return {44.45, 4.0};
    if(is_equal(type, "76.1X5.9CHS")) return {38.05, 5.9};
    if(is_equal(type, "76.1X4.5CHS")) return {38.05, 4.5};
    if(is_equal(type, "76.1X3.6CHS")) return {38.05, 3.6};
    if(is_equal(type, "60.3X5.4CHS")) return {30.15, 5.4};
    if(is_equal(type, "60.3X4.5CHS")) return {30.15, 4.5};
    if(is_equal(type, "60.3X3.6CHS")) return {30.15, 3.6};
    if(is_equal(type, "48.3X5.4CHS")) return {24.15, 5.4};
    if(is_equal(type, "48.3X4.0CHS")) return {24.15, 4.0};
    if(is_equal(type, "48.3X3.2CHS")) return {24.15, 3.2};
    if(is_equal(type, "42.4X4.9CHS")) return {21.20, 4.9};
    if(is_equal(type, "42.4X4.0CHS")) return {21.20, 4.0};
    if(is_equal(type, "42.4X3.2CHS")) return {21.20, 3.2};
    if(is_equal(type, "457.0X12.7CHS")) return {228.50, 12.7};
    if(is_equal(type, "457.0X9.5CHS")) return {228.50, 9.5};
    if(is_equal(type, "457.0X6.4CHS")) return {228.50, 6.4};
    if(is_equal(type, "406.4X12.7CHS")) return {203.20, 12.7};
    if(is_equal(type, "406.4X9.5CHS")) return {203.20, 9.5};
    if(is_equal(type, "406.4X6.4CHS")) return {203.20, 6.4};
    if(is_equal(type, "355.6X12.7CHS")) return {177.80, 12.7};
    if(is_equal(type, "355.6X9.5CHS")) return {177.80, 9.5};
    if(is_equal(type, "355.6X6.4CHS")) return {177.80, 6.4};
    if(is_equal(type, "323.9X12.7CHS")) return {161.95, 12.7};
    if(is_equal(type, "323.9X9.5CHS")) return {161.95, 9.5};
    if(is_equal(type, "323.9X6.4CHS")) return {161.95, 6.4};
    if(is_equal(type, "273.1X9.3CHS")) return {136.55, 9.3};
    if(is_equal(type, "273.1X6.4CHS")) return {136.55, 6.4};
    if(is_equal(type, "273.1X4.8CHS")) return {136.55, 4.8};
    if(is_equal(type, "219.1X8.2CHS")) return {109.55, 8.2};
    if(is_equal(type, "219.1X6.4CHS")) return {109.55, 6.4};
    if(is_equal(type, "219.1X4.8CHS")) return {109.55, 4.8};
    if(is_equal(type, "168.3X7.1CHS")) return {84.15, 7.1};
    if(is_equal(type, "168.3X6.4CHS")) return {84.15, 6.4};
    if(is_equal(type, "168.3X4.8CHS")) return {84.15, 4.8};
    if(is_equal(type, "165.1X3.5CHS")) return {82.55, 3.5};
    if(is_equal(type, "165.1X3.0CHS")) return {82.55, 3.0};
    if(is_equal(type, "139.7X3.5CHS")) return {69.85, 3.5};
    if(is_equal(type, "139.7X3.0CHS")) return {69.85, 3.0};
    if(is_equal(type, "114.3X6.0CHS")) return {57.15, 6.0};
    if(is_equal(type, "114.3X4.8CHS")) return {57.15, 4.8};
    if(is_equal(type, "114.3X3.6CHS")) return {57.15, 3.6};
    if(is_equal(type, "114.3X3.2CHS")) return {57.15, 3.2};
    if(is_equal(type, "101.6X3.2CHS")) return {50.80, 3.2};
    if(is_equal(type, "101.6X2.6CHS")) return {50.80, 2.6};
    if(is_equal(type, "88.9X5.5CHS")) return {44.45, 5.5};
    if(is_equal(type, "88.9X4.8CHS")) return {44.45, 4.8};
    if(is_equal(type, "88.9X3.2CHS")) return {44.45, 3.2};
    if(is_equal(type, "88.9X2.6CHS")) return {44.45, 2.6};
    if(is_equal(type, "76.1X3.2CHS")) return {38.05, 3.2};
    if(is_equal(type, "76.1X2.3CHS")) return {38.05, 2.3};
    return {};
}

vec nzisection(const std::string_view type) {
    if(is_equal(type, "1200WB455")) return {500.0, 40.0, 500.0, 40.0, 1120.0, 16.0};
    if(is_equal(type, "1200WB423")) return {500.0, 36.0, 500.0, 36.0, 1120.0, 16.0};
    if(is_equal(type, "1200WB392")) return {500.0, 32.0, 500.0, 32.0, 1120.0, 16.0};
    if(is_equal(type, "1200WB342")) return {400.0, 32.0, 400.0, 32.0, 1120.0, 16.0};
    if(is_equal(type, "1200WB317")) return {400.0, 28.0, 400.0, 28.0, 1120.0, 16.0};
    if(is_equal(type, "1200WB278")) return {350.0, 25.0, 350.0, 25.0, 1120.0, 16.0};
    if(is_equal(type, "1200WB249")) return {275.0, 25.0, 275.0, 25.0, 1120.0, 16.0};
    if(is_equal(type, "1000WB322")) return {400.0, 32.0, 400.0, 32.0, 960.0, 16.0};
    if(is_equal(type, "1000WB296")) return {400.0, 28.0, 400.0, 28.0, 960.0, 16.0};
    if(is_equal(type, "1000WB258")) return {350.0, 25.0, 350.0, 25.0, 960.0, 16.0};
    if(is_equal(type, "1000WB215")) return {300.0, 20.0, 300.0, 20.0, 960.0, 16.0};
    if(is_equal(type, "900WB282")) return {400.0, 32.0, 400.0, 32.0, 860.0, 12.0};
    if(is_equal(type, "900WB257")) return {400.0, 28.0, 400.0, 28.0, 860.0, 12.0};
    if(is_equal(type, "900WB218")) return {350.0, 25.0, 350.0, 25.0, 860.0, 12.0};
    if(is_equal(type, "900WB175")) return {300.0, 20.0, 300.0, 20.0, 860.0, 12.0};
    if(is_equal(type, "800WB192")) return {300.0, 28.0, 300.0, 28.0, 760.0, 10.0};
    if(is_equal(type, "800WB168")) return {275.0, 25.0, 275.0, 25.0, 760.0, 10.0};
    if(is_equal(type, "800WB146")) return {275.0, 20.0, 275.0, 20.0, 760.0, 10.0};
    if(is_equal(type, "800WB122")) return {250.0, 16.0, 250.0, 16.0, 760.0, 10.0};
    if(is_equal(type, "700WB173")) return {275.0, 28.0, 275.0, 28.0, 660.0, 10.0};
    if(is_equal(type, "700WB150")) return {250.0, 25.0, 250.0, 25.0, 660.0, 10.0};
    if(is_equal(type, "700WB130")) return {250.0, 20.0, 250.0, 20.0, 660.0, 10.0};
    if(is_equal(type, "700WB115")) return {250.0, 16.0, 250.0, 16.0, 660.0, 10.0};
    if(is_equal(type, "500WC440")) return {500.0, 40.0, 500.0, 40.0, 400.0, 40.0};
    if(is_equal(type, "500WC414")) return {500.0, 40.0, 500.0, 40.0, 400.0, 32.0};
    if(is_equal(type, "500WC383")) return {500.0, 36.0, 500.0, 36.0, 400.0, 32.0};
    if(is_equal(type, "500WC340")) return {500.0, 32.0, 500.0, 32.0, 450.0, 25.0};
    if(is_equal(type, "500WC290")) return {500.0, 28.0, 500.0, 28.0, 450.0, 20.0};
    if(is_equal(type, "500WC267")) return {500.0, 25.0, 500.0, 25.0, 450.0, 20.0};
    if(is_equal(type, "500WC228")) return {500.0, 20.0, 500.0, 20.0, 450.0, 20.0};
    if(is_equal(type, "400WC361")) return {400.0, 40.0, 400.0, 40.0, 350.0, 40.0};
    if(is_equal(type, "400WC328")) return {400.0, 40.0, 400.0, 40.0, 350.0, 28.0};
    if(is_equal(type, "400WC303")) return {400.0, 36.0, 400.0, 36.0, 350.0, 28.0};
    if(is_equal(type, "400WC270")) return {400.0, 32.0, 400.0, 32.0, 350.0, 25.0};
    if(is_equal(type, "400WC212")) return {400.0, 25.0, 400.0, 25.0, 350.0, 20.0};
    if(is_equal(type, "400WC181")) return {400.0, 20.0, 400.0, 20.0, 350.0, 20.0};
    if(is_equal(type, "400WC144")) return {400.0, 16.0, 400.0, 16.0, 350.0, 16.0};
    if(is_equal(type, "350WC280")) return {350.0, 40.0, 350.0, 40.0, 275.0, 28.0};
    if(is_equal(type, "350WC258")) return {350.0, 36.0, 350.0, 36.0, 275.0, 28.0};
    if(is_equal(type, "350WC230")) return {350.0, 32.0, 350.0, 32.0, 275.0, 25.0};
    if(is_equal(type, "350WC197")) return {350.0, 28.0, 350.0, 28.0, 275.0, 20.0};
    if(is_equal(type, "610UB125")) return {229.0, 19.6, 229.0, 19.6, 572.0, 11.9};
    if(is_equal(type, "610UB113")) return {228.0, 17.3, 228.0, 17.3, 572.0, 11.2};
    if(is_equal(type, "610UB101")) return {228.0, 14.8, 228.0, 14.8, 572.0, 10.6};
    if(is_equal(type, "530UB92.4")) return {209.0, 15.6, 209.0, 15.6, 502.0, 10.2};
    if(is_equal(type, "530UB82.0")) return {209.0, 13.2, 209.0, 13.2, 502.0, 9.6};
    if(is_equal(type, "460UB82.1")) return {191.0, 16.0, 191.0, 16.0, 428.0, 9.9};
    if(is_equal(type, "460UB74.6")) return {190.0, 14.5, 190.0, 14.5, 428.0, 9.1};
    if(is_equal(type, "460UB67.1")) return {190.0, 12.7, 190.0, 12.7, 428.0, 8.5};
    if(is_equal(type, "410UB59.7")) return {178.0, 12.8, 178.0, 12.8, 381.0, 7.8};
    if(is_equal(type, "410UB53.7")) return {178.0, 10.9, 178.0, 10.9, 381.0, 7.6};
    if(is_equal(type, "360UB56.7")) return {172.0, 13.0, 172.0, 13.0, 333.0, 8.0};
    if(is_equal(type, "360UB50.7")) return {171.0, 11.5, 171.0, 11.5, 333.0, 7.3};
    if(is_equal(type, "360UB44.7")) return {171.0, 9.7, 171.0, 9.7, 333.0, 6.9};
    if(is_equal(type, "310UB46.2")) return {166.0, 11.8, 166.0, 11.8, 284.0, 6.7};
    if(is_equal(type, "310UB40.4")) return {165.0, 10.2, 165.0, 10.2, 284.0, 6.1};
    if(is_equal(type, "310UB32.0")) return {149.0, 8.0, 149.0, 8.0, 282.0, 5.5};
    if(is_equal(type, "250UB37.3")) return {146.0, 10.9, 146.0, 10.9, 234.0, 6.4};
    if(is_equal(type, "250UB31.4")) return {146.0, 8.6, 146.0, 8.6, 234.0, 6.1};
    if(is_equal(type, "250UB25.7")) return {124.0, 8.0, 124.0, 8.0, 232.0, 5.0};
    if(is_equal(type, "200UB29.8")) return {134.0, 9.6, 134.0, 9.6, 188.0, 6.3};
    if(is_equal(type, "200UB25.4")) return {133.0, 7.8, 133.0, 7.8, 188.0, 5.8};
    if(is_equal(type, "200UB22.3")) return {133.0, 7.0, 133.0, 7.0, 188.0, 5.0};
    if(is_equal(type, "200UB18.2")) return {99.0, 7.0, 99.0, 7.0, 184.0, 4.5};
    if(is_equal(type, "180UB22.2")) return {90.0, 10.0, 90.0, 10.0, 159.0, 6.0};
    if(is_equal(type, "180UB18.1")) return {90.0, 8.0, 90.0, 8.0, 159.0, 5.0};
    if(is_equal(type, "180UB16.1")) return {90.0, 7.0, 90.0, 7.0, 159.0, 4.5};
    if(is_equal(type, "150UB18.0")) return {75.0, 9.5, 75.0, 9.5, 136.0, 6.0};
    if(is_equal(type, "150UB14.0")) return {75.0, 7.0, 75.0, 7.0, 136.0, 5.0};
    if(is_equal(type, "320UC158")) return {311.0, 25.0, 311.0, 25.0, 277.0, 15.7};
    if(is_equal(type, "310UC137")) return {309.0, 21.7, 309.0, 21.7, 277.0, 13.8};
    if(is_equal(type, "310UC118")) return {307.0, 18.7, 307.0, 18.7, 277.0, 11.9};
    if(is_equal(type, "310UC96.8")) return {305.0, 15.4, 305.0, 15.4, 277.0, 9.9};
    if(is_equal(type, "250UC89.5")) return {256.0, 17.3, 256.0, 17.3, 225.0, 10.5};
    if(is_equal(type, "250UC72.9")) return {254.0, 14.2, 254.0, 14.2, 225.0, 8.6};
    if(is_equal(type, "200UC59.5")) return {205.0, 14.2, 205.0, 14.2, 181.0, 9.3};
    if(is_equal(type, "200UC52.2")) return {204.0, 12.5, 204.0, 12.5, 181.0, 8.0};
    if(is_equal(type, "200UC46.2")) return {203.0, 11.0, 203.0, 11.0, 181.0, 7.3};
    if(is_equal(type, "150UC37.2")) return {154.0, 11.5, 154.0, 11.5, 139.0, 8.1};
    if(is_equal(type, "150UC30.0")) return {153.0, 9.4, 153.0, 9.4, 139.0, 6.6};
    if(is_equal(type, "150UC23.4")) return {152.0, 6.8, 152.0, 6.8, 139.0, 6.1};
    if(is_equal(type, "100UC14.8")) return {99.0, 7.0, 99.0, 7.0, 83.0, 5.0};
    return {};
}

vec nzrhsection(const std::string_view type) {
    if(is_equal(type, "250X150X9.0RHS")) return {150.0, 250.0, 9.0};
    if(is_equal(type, "250X150X6.0RHS")) return {150.0, 250.0, 6.0};
    if(is_equal(type, "250X150X5.0RHS")) return {150.0, 250.0, 5.0};
    if(is_equal(type, "200X100X9.0RHS")) return {100.0, 200.0, 9.0};
    if(is_equal(type, "200X100X6.0RHS")) return {100.0, 200.0, 6.0};
    if(is_equal(type, "200X100X5.0RHS")) return {100.0, 200.0, 5.0};
    if(is_equal(type, "200X100X4.0RHS")) return {100.0, 200.0, 4.0};
    if(is_equal(type, "150X100X6.0RHS")) return {100.0, 150.0, 6.0};
    if(is_equal(type, "150X100X5.0RHS")) return {100.0, 150.0, 5.0};
    if(is_equal(type, "150X100X4.0RHS")) return {100.0, 150.0, 4.0};
    if(is_equal(type, "150X50X5.0RHS")) return {50.0, 150.0, 5.0};
    if(is_equal(type, "150X50X4.0RHS")) return {50.0, 150.0, 4.0};
    if(is_equal(type, "150X50X3.0RHS")) return {50.0, 150.0, 3.0};
    if(is_equal(type, "125X75X5.0RHS")) return {75.0, 125.0, 5.0};
    if(is_equal(type, "125X75X4.0RHS")) return {75.0, 125.0, 4.0};
    if(is_equal(type, "125X75X3.0RHS")) return {75.0, 125.0, 3.0};
    if(is_equal(type, "100X50X6.0RHS")) return {50.0, 100.0, 6.0};
    if(is_equal(type, "100X50X5.0RHS")) return {50.0, 100.0, 5.0};
    if(is_equal(type, "100X50X4.0RHS")) return {50.0, 100.0, 4.0};
    if(is_equal(type, "100X50X3.5RHS")) return {50.0, 100.0, 3.5};
    if(is_equal(type, "100X50X3.0RHS")) return {50.0, 100.0, 3.0};
    if(is_equal(type, "100X50X2.5RHS")) return {50.0, 100.0, 2.5};
    if(is_equal(type, "100X50X2.0RHS")) return {50.0, 100.0, 2.0};
    if(is_equal(type, "75X50X4.0RHS")) return {50.0, 75.0, 4.0};
    if(is_equal(type, "75X50X3.0RHS")) return {50.0, 75.0, 3.0};
    if(is_equal(type, "75X50X2.5RHS")) return {50.0, 75.0, 2.5};
    if(is_equal(type, "75X50X2.0RHS")) return {50.0, 75.0, 2.0};
    if(is_equal(type, "75X25X2.5RHS")) return {25.0, 75.0, 2.5};
    if(is_equal(type, "75X25X2.0RHS")) return {25.0, 75.0, 2.0};
    if(is_equal(type, "75X25X1.6RHS")) return {25.0, 75.0, 1.6};
    if(is_equal(type, "65X35X3.0RHS")) return {35.0, 65.0, 3.0};
    if(is_equal(type, "65X35X2.5RHS")) return {35.0, 65.0, 2.5};
    if(is_equal(type, "65X35X2.0RHS")) return {35.0, 65.0, 2.0};
    if(is_equal(type, "50X25X3.0RHS")) return {25.0, 50.0, 3.0};
    if(is_equal(type, "50X25X2.5RHS")) return {25.0, 50.0, 2.5};
    if(is_equal(type, "50X25X2.0RHS")) return {25.0, 50.0, 2.0};
    if(is_equal(type, "50X25X1.6RHS")) return {25.0, 50.0, 1.6};
    if(is_equal(type, "50X20X3.0RHS")) return {20.0, 50.0, 3.0};
    if(is_equal(type, "50X20X2.5RHS")) return {20.0, 50.0, 2.5};
    if(is_equal(type, "50X20X2.0RHS")) return {20.0, 50.0, 2.0};
    if(is_equal(type, "50X20X1.6RHS")) return {20.0, 50.0, 1.6};
    return {};
}

vec nzshsection(const std::string_view type) {
    if(is_equal(type, "250X9.0SHS")) return {250.0, 250.0, 9.0};
    if(is_equal(type, "250X6.0SHS")) return {250.0, 250.0, 6.0};
    if(is_equal(type, "200X9.0SHS")) return {200.0, 200.0, 9.0};
    if(is_equal(type, "200X6.0SHS")) return {200.0, 200.0, 6.0};
    if(is_equal(type, "200X5.0SHS")) return {200.0, 200.0, 5.0};
    if(is_equal(type, "150X9.0SHS")) return {150.0, 150.0, 9.0};
    if(is_equal(type, "150X6.0SHS")) return {150.0, 150.0, 6.0};
    if(is_equal(type, "150X5.0SHS")) return {150.0, 150.0, 5.0};
    if(is_equal(type, "125X9.0SHS")) return {125.0, 125.0, 9.0};
    if(is_equal(type, "125X6.0SHS")) return {125.0, 125.0, 6.0};
    if(is_equal(type, "125X5.0SHS")) return {125.0, 125.0, 5.0};
    if(is_equal(type, "125X4.0SHS")) return {125.0, 125.0, 4.0};
    if(is_equal(type, "100X9.0SHS")) return {100.0, 100.0, 9.0};
    if(is_equal(type, "100X6.0SHS")) return {100.0, 100.0, 6.0};
    if(is_equal(type, "100X5.0SHS")) return {100.0, 100.0, 5.0};
    if(is_equal(type, "100X4.0SHS")) return {100.0, 100.0, 4.0};
    if(is_equal(type, "100X3.0SHS")) return {100.0, 100.0, 3.0};
    if(is_equal(type, "89X6.0SHS")) return {89.0, 89.0, 6.0};
    if(is_equal(type, "89X5.0SHS")) return {89.0, 89.0, 5.0};
    if(is_equal(type, "89X3.5SHS")) return {89.0, 89.0, 3.5};
    if(is_equal(type, "75X6.0SHS")) return {75.0, 75.0, 6.0};
    if(is_equal(type, "75X5.0SHS")) return {75.0, 75.0, 5.0};
    if(is_equal(type, "75X4.0SHS")) return {75.0, 75.0, 4.0};
    if(is_equal(type, "75X3.5SHS")) return {75.0, 75.0, 3.5};
    if(is_equal(type, "75X3.0SHS")) return {75.0, 75.0, 3.0};
    if(is_equal(type, "75X2.5SHS")) return {75.0, 75.0, 2.5};
    if(is_equal(type, "65X3.0SHS")) return {65.0, 65.0, 3.0};
    if(is_equal(type, "65X2.5SHS")) return {65.0, 65.0, 2.5};
    if(is_equal(type, "65X2.0SHS")) return {65.0, 65.0, 2.0};
    if(is_equal(type, "50X4.0SHS")) return {50.0, 50.0, 4.0};
    if(is_equal(type, "50X3.0SHS")) return {50.0, 50.0, 3.0};
    if(is_equal(type, "50X2.5SHS")) return {50.0, 50.0, 2.5};
    if(is_equal(type, "50X2.0SHS")) return {50.0, 50.0, 2.0};
    if(is_equal(type, "50X1.6SHS")) return {50.0, 50.0, 1.6};
    if(is_equal(type, "40X4.0SHS")) return {40.0, 40.0, 4.0};
    if(is_equal(type, "40X2.5SHS")) return {40.0, 40.0, 2.5};
    if(is_equal(type, "40X2.0SHS")) return {40.0, 40.0, 2.0};
    if(is_equal(type, "40X1.6SHS")) return {40.0, 40.0, 1.6};
    if(is_equal(type, "35X3.0SHS")) return {35.0, 35.0, 3.0};
    if(is_equal(type, "35X2.5SHS")) return {35.0, 35.0, 2.5};
    if(is_equal(type, "35X2.0SHS")) return {35.0, 35.0, 2.0};
    if(is_equal(type, "35X1.6SHS")) return {35.0, 35.0, 1.6};
    if(is_equal(type, "30X2.0SHS")) return {30.0, 30.0, 2.0};
    if(is_equal(type, "30X1.6SHS")) return {30.0, 30.0, 1.6};
    if(is_equal(type, "25X3.0SHS")) return {25.0, 25.0, 3.0};
    if(is_equal(type, "25X2.5SHS")) return {25.0, 25.0, 2.5};
    if(is_equal(type, "25X2.0SHS")) return {25.0, 25.0, 2.0};
    if(is_equal(type, "25X1.6SHS")) return {25.0, 25.0, 1.6};
    if(is_equal(type, "20X1.6SHS")) return {20.0, 20.0, 1.6};
    return {};
}

vec usisection(const std::string_view type) {
    if(is_equal(type, "W44X408")) return {16.100, 2.170, 16.100, 2.170, 40.460, 1.220};
    if(is_equal(type, "W44X368")) return {16.000, 1.970, 16.000, 1.970, 40.460, 1.100};
    if(is_equal(type, "W44X335")) return {15.900, 1.770, 15.900, 1.770, 40.460, 1.030};
    if(is_equal(type, "W44X290")) return {15.800, 1.580, 15.800, 1.580, 40.440, 0.865};
    if(is_equal(type, "W44X262")) return {15.800, 1.420, 15.800, 1.420, 40.460, 0.785};
    if(is_equal(type, "W44X230")) return {15.800, 1.220, 15.800, 1.220, 40.460, 0.710};
    if(is_equal(type, "W40X655")) return {16.900, 3.540, 16.900, 3.540, 36.520, 1.970};
    if(is_equal(type, "W40X593")) return {16.700, 3.230, 16.700, 3.230, 36.540, 1.790};
    if(is_equal(type, "W40X503")) return {16.400, 2.760, 16.400, 2.760, 36.580, 1.540};
    if(is_equal(type, "W40X431")) return {16.200, 2.360, 16.200, 2.360, 36.580, 1.340};
    if(is_equal(type, "W40X397")) return {16.100, 2.200, 16.100, 2.200, 36.600, 1.220};
    if(is_equal(type, "W40X372")) return {16.100, 2.050, 16.100, 2.050, 36.500, 1.160};
    if(is_equal(type, "W40X362")) return {16.000, 2.010, 16.000, 2.010, 36.580, 1.120};
    if(is_equal(type, "W40X324")) return {15.900, 1.810, 15.900, 1.810, 36.580, 1.000};
    if(is_equal(type, "W40X297")) return {15.800, 1.650, 15.800, 1.650, 36.500, 0.930};
    if(is_equal(type, "W40X277")) return {15.800, 1.580, 15.800, 1.580, 36.540, 0.830};
    if(is_equal(type, "W40X249")) return {15.800, 1.420, 15.800, 1.420, 36.560, 0.750};
    if(is_equal(type, "W40X215")) return {15.800, 1.220, 15.800, 1.220, 36.560, 0.650};
    if(is_equal(type, "W40X199")) return {15.800, 1.070, 15.800, 1.070, 36.560, 0.650};
    if(is_equal(type, "W40X392")) return {12.400, 2.520, 12.400, 2.520, 36.560, 1.420};
    if(is_equal(type, "W40X331")) return {12.200, 2.130, 12.200, 2.130, 36.540, 1.220};
    if(is_equal(type, "W40X327")) return {12.100, 2.130, 12.100, 2.130, 36.540, 1.180};
    if(is_equal(type, "W40X294")) return {12.000, 1.930, 12.000, 1.930, 36.540, 1.060};
    if(is_equal(type, "W40X278")) return {12.000, 1.810, 12.000, 1.810, 36.580, 1.030};
    if(is_equal(type, "W40X264")) return {11.900, 1.730, 11.900, 1.730, 36.540, 0.960};
    if(is_equal(type, "W40X235")) return {11.900, 1.580, 11.900, 1.580, 36.540, 0.830};
    if(is_equal(type, "W40X211")) return {11.800, 1.420, 11.800, 1.420, 36.560, 0.750};
    if(is_equal(type, "W40X183")) return {11.800, 1.200, 11.800, 1.200, 36.600, 0.650};
    if(is_equal(type, "W40X167")) return {11.800, 1.030, 11.800, 1.030, 36.540, 0.650};
    if(is_equal(type, "W40X149")) return {11.800, 0.830, 11.800, 0.830, 36.540, 0.630};
    if(is_equal(type, "W36X925")) return {18.600, 4.530, 18.600, 4.530, 34.040, 3.020};
    if(is_equal(type, "W36X853")) return {18.200, 4.530, 18.200, 4.530, 34.040, 2.520};
    if(is_equal(type, "W36X802")) return {18.000, 4.290, 18.000, 4.290, 34.020, 2.380};
    if(is_equal(type, "W36X723")) return {17.800, 3.900, 17.800, 3.900, 34.000, 2.170};
    if(is_equal(type, "W36X652")) return {17.600, 3.540, 17.600, 3.540, 34.020, 1.970};
    if(is_equal(type, "W36X529")) return {17.200, 2.910, 17.200, 2.910, 33.980, 1.610};
    if(is_equal(type, "W36X487")) return {17.100, 2.680, 17.100, 2.680, 33.940, 1.500};
    if(is_equal(type, "W36X441")) return {17.000, 2.440, 17.000, 2.440, 34.020, 1.360};
    if(is_equal(type, "W36X395")) return {16.800, 2.200, 16.800, 2.200, 34.000, 1.220};
    if(is_equal(type, "W36X361")) return {16.700, 2.010, 16.700, 2.010, 33.980, 1.120};
    if(is_equal(type, "W36X330")) return {16.600, 1.850, 16.600, 1.850, 34.000, 1.020};
    if(is_equal(type, "W36X302")) return {16.700, 1.680, 16.700, 1.680, 33.940, 0.945};
    if(is_equal(type, "W36X282")) return {16.600, 1.570, 16.600, 1.570, 33.960, 0.885};
    if(is_equal(type, "W36X262")) return {16.600, 1.440, 16.600, 1.440, 34.020, 0.840};
    if(is_equal(type, "W36X247")) return {16.500, 1.350, 16.500, 1.350, 34.000, 0.800};
    if(is_equal(type, "W36X231")) return {16.500, 1.260, 16.500, 1.260, 33.980, 0.760};
    if(is_equal(type, "W36X387")) return {12.700, 2.560, 12.700, 2.560, 33.980, 1.420};
    if(is_equal(type, "W36X350")) return {12.600, 2.320, 12.600, 2.320, 33.960, 1.300};
    if(is_equal(type, "W36X318")) return {12.400, 2.130, 12.400, 2.130, 33.940, 1.180};
    if(is_equal(type, "W36X286")) return {12.300, 1.930, 12.300, 1.930, 33.940, 1.060};
    if(is_equal(type, "W36X256")) return {12.200, 1.730, 12.200, 1.730, 33.940, 0.960};
    if(is_equal(type, "W36X232")) return {12.100, 1.570, 12.100, 1.570, 33.960, 0.870};
    if(is_equal(type, "W36X210")) return {12.200, 1.360, 12.200, 1.360, 33.980, 0.830};
    if(is_equal(type, "W36X194")) return {12.100, 1.260, 12.100, 1.260, 33.980, 0.765};
    if(is_equal(type, "W36X182")) return {12.100, 1.180, 12.100, 1.180, 33.940, 0.725};
    if(is_equal(type, "W36X170")) return {12.000, 1.100, 12.000, 1.100, 34.000, 0.680};
    if(is_equal(type, "W36X160")) return {12.000, 1.020, 12.000, 1.020, 33.960, 0.650};
    if(is_equal(type, "W36X150")) return {12.000, 0.940, 12.000, 0.940, 34.020, 0.625};
    if(is_equal(type, "W36X135")) return {12.000, 0.790, 12.000, 0.790, 34.020, 0.600};
    if(is_equal(type, "W33X387")) return {16.200, 2.280, 16.200, 2.280, 31.440, 1.260};
    if(is_equal(type, "W33X354")) return {16.100, 2.090, 16.100, 2.090, 31.420, 1.160};
    if(is_equal(type, "W33X318")) return {16.000, 1.890, 16.000, 1.890, 31.420, 1.040};
    if(is_equal(type, "W33X291")) return {15.900, 1.730, 15.900, 1.730, 31.340, 0.960};
    if(is_equal(type, "W33X263")) return {15.800, 1.570, 15.800, 1.570, 31.360, 0.870};
    if(is_equal(type, "W33X241")) return {15.900, 1.400, 15.900, 1.400, 31.400, 0.830};
    if(is_equal(type, "W33X221")) return {15.800, 1.280, 15.800, 1.280, 31.340, 0.775};
    if(is_equal(type, "W33X201")) return {15.700, 1.150, 15.700, 1.150, 31.400, 0.715};
    if(is_equal(type, "W33X169")) return {11.500, 1.220, 11.500, 1.220, 31.360, 0.670};
    if(is_equal(type, "W33X152")) return {11.600, 1.060, 11.600, 1.060, 31.380, 0.635};
    if(is_equal(type, "W33X141")) return {11.500, 0.960, 11.500, 0.960, 31.380, 0.605};
    if(is_equal(type, "W33X130")) return {11.500, 0.855, 11.500, 0.855, 31.390, 0.580};
    if(is_equal(type, "W33X118")) return {11.500, 0.740, 11.500, 0.740, 31.420, 0.550};
    if(is_equal(type, "W30X391")) return {15.600, 2.440, 15.600, 2.440, 28.320, 1.360};
    if(is_equal(type, "W30X357")) return {15.500, 2.240, 15.500, 2.240, 28.320, 1.240};
    if(is_equal(type, "W30X326")) return {15.400, 2.050, 15.400, 2.050, 28.300, 1.140};
    if(is_equal(type, "W30X292")) return {15.300, 1.850, 15.300, 1.850, 28.300, 1.020};
    if(is_equal(type, "W30X261")) return {15.200, 1.650, 15.200, 1.650, 28.300, 0.930};
    if(is_equal(type, "W30X235")) return {15.100, 1.500, 15.100, 1.500, 28.300, 0.830};
    if(is_equal(type, "W30X211")) return {15.100, 1.320, 15.100, 1.320, 28.260, 0.775};
    if(is_equal(type, "W30X191")) return {15.000, 1.190, 15.000, 1.190, 28.320, 0.710};
    if(is_equal(type, "W30X173")) return {15.000, 1.070, 15.000, 1.070, 28.260, 0.655};
    if(is_equal(type, "W30X148")) return {10.500, 1.180, 10.500, 1.180, 28.340, 0.650};
    if(is_equal(type, "W30X132")) return {10.500, 1.000, 10.500, 1.000, 28.300, 0.615};
    if(is_equal(type, "W30X124")) return {10.500, 0.930, 10.500, 0.930, 28.340, 0.585};
    if(is_equal(type, "W30X116")) return {10.500, 0.850, 10.500, 0.850, 28.300, 0.565};
    if(is_equal(type, "W30X108")) return {10.500, 0.760, 10.500, 0.760, 28.280, 0.545};
    if(is_equal(type, "W30X99")) return {10.500, 0.670, 10.500, 0.670, 28.360, 0.520};
    if(is_equal(type, "W30X90")) return {10.400, 0.610, 10.400, 0.610, 28.280, 0.470};
    if(is_equal(type, "W27X539")) return {15.300, 3.540, 15.300, 3.540, 25.420, 1.970};
    if(is_equal(type, "W27X368")) return {14.700, 2.480, 14.700, 2.480, 25.440, 1.380};
    if(is_equal(type, "W27X336")) return {14.600, 2.280, 14.600, 2.280, 25.440, 1.260};
    if(is_equal(type, "W27X307")) return {14.400, 2.090, 14.400, 2.090, 25.420, 1.160};
    if(is_equal(type, "W27X281")) return {14.400, 1.930, 14.400, 1.930, 25.440, 1.060};
    if(is_equal(type, "W27X258")) return {14.300, 1.770, 14.300, 1.770, 25.460, 0.980};
    if(is_equal(type, "W27X235")) return {14.200, 1.610, 14.200, 1.610, 25.480, 0.910};
    if(is_equal(type, "W27X217")) return {14.100, 1.500, 14.100, 1.500, 25.400, 0.830};
    if(is_equal(type, "W27X194")) return {14.000, 1.340, 14.000, 1.340, 25.420, 0.750};
    if(is_equal(type, "W27X178")) return {14.100, 1.190, 14.100, 1.190, 25.420, 0.725};
    if(is_equal(type, "W27X161")) return {14.000, 1.080, 14.000, 1.080, 25.440, 0.660};
    if(is_equal(type, "W27X146")) return {14.000, 0.975, 14.000, 0.975, 25.450, 0.605};
    if(is_equal(type, "W27X129")) return {10.000, 1.100, 10.000, 1.100, 25.400, 0.610};
    if(is_equal(type, "W27X114")) return {10.100, 0.930, 10.100, 0.930, 25.440, 0.570};
    if(is_equal(type, "W27X102")) return {10.000, 0.830, 10.000, 0.830, 25.440, 0.515};
    if(is_equal(type, "W27X94")) return {10.000, 0.745, 10.000, 0.745, 25.410, 0.490};
    if(is_equal(type, "W27X84")) return {10.000, 0.640, 10.000, 0.640, 25.420, 0.460};
    if(is_equal(type, "W24X370")) return {13.700, 2.720, 13.700, 2.720, 22.560, 1.520};
    if(is_equal(type, "W24X335")) return {13.500, 2.480, 13.500, 2.480, 22.540, 1.380};
    if(is_equal(type, "W24X306")) return {13.400, 2.280, 13.400, 2.280, 22.540, 1.260};
    if(is_equal(type, "W24X279")) return {13.300, 2.090, 13.300, 2.090, 22.520, 1.160};
    if(is_equal(type, "W24X250")) return {13.200, 1.890, 13.200, 1.890, 22.520, 1.040};
    if(is_equal(type, "W24X229")) return {13.100, 1.730, 13.100, 1.730, 22.540, 0.960};
    if(is_equal(type, "W24X207")) return {13.000, 1.570, 13.000, 1.570, 22.560, 0.870};
    if(is_equal(type, "W24X192")) return {13.000, 1.460, 13.000, 1.460, 22.580, 0.810};
    if(is_equal(type, "W24X176")) return {12.900, 1.340, 12.900, 1.340, 22.520, 0.750};
    if(is_equal(type, "W24X162")) return {13.000, 1.220, 13.000, 1.220, 22.560, 0.705};
    if(is_equal(type, "W24X146")) return {12.900, 1.090, 12.900, 1.090, 22.520, 0.650};
    if(is_equal(type, "W24X131")) return {12.900, 0.960, 12.900, 0.960, 22.580, 0.605};
    if(is_equal(type, "W24X117")) return {12.800, 0.850, 12.800, 0.850, 22.600, 0.550};
    if(is_equal(type, "W24X104")) return {12.800, 0.750, 12.800, 0.750, 22.600, 0.500};
    if(is_equal(type, "W24X103")) return {9.000, 0.980, 9.000, 0.980, 22.540, 0.550};
    if(is_equal(type, "W24X94")) return {9.070, 0.875, 9.070, 0.875, 22.550, 0.515};
    if(is_equal(type, "W24X84")) return {9.020, 0.770, 9.020, 0.770, 22.560, 0.470};
    if(is_equal(type, "W24X76")) return {8.990, 0.680, 8.990, 0.680, 22.540, 0.440};
    if(is_equal(type, "W24X68")) return {8.970, 0.585, 8.970, 0.585, 22.530, 0.415};
    if(is_equal(type, "W24X62")) return {7.040, 0.590, 7.040, 0.590, 22.520, 0.430};
    if(is_equal(type, "W24X55")) return {7.010, 0.505, 7.010, 0.505, 22.590, 0.395};
    if(is_equal(type, "W21X275")) return {12.900, 2.190, 12.900, 2.190, 19.720, 1.220};
    if(is_equal(type, "W21X248")) return {12.800, 1.990, 12.800, 1.990, 19.720, 1.100};
    if(is_equal(type, "W21X223")) return {12.700, 1.790, 12.700, 1.790, 19.820, 1.000};
    if(is_equal(type, "W21X201")) return {12.600, 1.630, 12.600, 1.630, 19.740, 0.910};
    if(is_equal(type, "W21X182")) return {12.500, 1.480, 12.500, 1.480, 19.740, 0.830};
    if(is_equal(type, "W21X166")) return {12.400, 1.360, 12.400, 1.360, 19.780, 0.750};
    if(is_equal(type, "W21X147")) return {12.500, 1.150, 12.500, 1.150, 19.800, 0.720};
    if(is_equal(type, "W21X132")) return {12.400, 1.040, 12.400, 1.040, 19.720, 0.650};
    if(is_equal(type, "W21X122")) return {12.400, 0.960, 12.400, 0.960, 19.780, 0.600};
    if(is_equal(type, "W21X111")) return {12.300, 0.875, 12.300, 0.875, 19.750, 0.550};
    if(is_equal(type, "W21X101")) return {12.300, 0.800, 12.300, 0.800, 19.800, 0.500};
    if(is_equal(type, "W21X93")) return {8.420, 0.930, 8.420, 0.930, 19.740, 0.580};
    if(is_equal(type, "W21X83")) return {8.360, 0.835, 8.360, 0.835, 19.730, 0.515};
    if(is_equal(type, "W21X73")) return {8.300, 0.740, 8.300, 0.740, 19.720, 0.455};
    if(is_equal(type, "W21X68")) return {8.270, 0.685, 8.270, 0.685, 19.730, 0.430};
    if(is_equal(type, "W21X62")) return {8.240, 0.615, 8.240, 0.615, 19.770, 0.400};
    if(is_equal(type, "W21X55")) return {8.220, 0.522, 8.220, 0.522, 19.756, 0.375};
    if(is_equal(type, "W21X48")) return {8.140, 0.430, 8.140, 0.430, 19.740, 0.350};
    if(is_equal(type, "W21X57")) return {6.560, 0.650, 6.560, 0.650, 19.800, 0.405};
    if(is_equal(type, "W21X50")) return {6.530, 0.535, 6.530, 0.535, 19.730, 0.380};
    if(is_equal(type, "W21X44")) return {6.500, 0.450, 6.500, 0.450, 19.800, 0.350};
    if(is_equal(type, "W18X311")) return {12.000, 2.740, 12.000, 2.740, 16.820, 1.520};
    if(is_equal(type, "W18X283")) return {11.900, 2.500, 11.900, 2.500, 16.900, 1.400};
    if(is_equal(type, "W18X258")) return {11.800, 2.300, 11.800, 2.300, 16.900, 1.280};
    if(is_equal(type, "W18X234")) return {11.700, 2.110, 11.700, 2.110, 16.880, 1.160};
    if(is_equal(type, "W18X211")) return {11.600, 1.910, 11.600, 1.910, 16.880, 1.060};
    if(is_equal(type, "W18X192")) return {11.500, 1.750, 11.500, 1.750, 16.900, 0.960};
    if(is_equal(type, "W18X175")) return {11.400, 1.590, 11.400, 1.590, 16.820, 0.890};
    if(is_equal(type, "W18X158")) return {11.300, 1.440, 11.300, 1.440, 16.820, 0.810};
    if(is_equal(type, "W18X143")) return {11.200, 1.320, 11.200, 1.320, 16.860, 0.730};
    if(is_equal(type, "W18X130")) return {11.200, 1.200, 11.200, 1.200, 16.900, 0.670};
    if(is_equal(type, "W18X119")) return {11.300, 1.060, 11.300, 1.060, 16.880, 0.655};
    if(is_equal(type, "W18X106")) return {11.200, 0.940, 11.200, 0.940, 16.820, 0.590};
    if(is_equal(type, "W18X97")) return {11.100, 0.870, 11.100, 0.870, 16.860, 0.535};
    if(is_equal(type, "W18X86")) return {11.100, 0.770, 11.100, 0.770, 16.860, 0.480};
    if(is_equal(type, "W18X76")) return {11.000, 0.680, 11.000, 0.680, 16.840, 0.425};
    if(is_equal(type, "W18X71")) return {7.640, 0.810, 7.640, 0.810, 16.880, 0.495};
    if(is_equal(type, "W18X65")) return {7.590, 0.750, 7.590, 0.750, 16.900, 0.450};
    if(is_equal(type, "W18X60")) return {7.560, 0.695, 7.560, 0.695, 16.810, 0.415};
    if(is_equal(type, "W18X55")) return {7.530, 0.630, 7.530, 0.630, 16.840, 0.390};
    if(is_equal(type, "W18X50")) return {7.500, 0.570, 7.500, 0.570, 16.860, 0.355};
    if(is_equal(type, "W18X46")) return {6.060, 0.605, 6.060, 0.605, 16.890, 0.360};
    if(is_equal(type, "W18X40")) return {6.020, 0.525, 6.020, 0.525, 16.850, 0.315};
    if(is_equal(type, "W18X35")) return {6.000, 0.425, 6.000, 0.425, 16.850, 0.300};
    if(is_equal(type, "W16X100")) return {10.400, 0.985, 10.400, 0.985, 15.030, 0.585};
    if(is_equal(type, "W16X89")) return {10.400, 0.875, 10.400, 0.875, 15.050, 0.525};
    if(is_equal(type, "W16X77")) return {10.300, 0.760, 10.300, 0.760, 14.980, 0.455};
    if(is_equal(type, "W16X67")) return {10.200, 0.665, 10.200, 0.665, 14.970, 0.395};
    if(is_equal(type, "W16X57")) return {7.120, 0.715, 7.120, 0.715, 14.970, 0.430};
    if(is_equal(type, "W16X50")) return {7.070, 0.630, 7.070, 0.630, 15.040, 0.380};
    if(is_equal(type, "W16X45")) return {7.040, 0.565, 7.040, 0.565, 14.970, 0.345};
    if(is_equal(type, "W16X40")) return {7.000, 0.505, 7.000, 0.505, 14.990, 0.305};
    if(is_equal(type, "W16X36")) return {6.990, 0.430, 6.990, 0.430, 15.040, 0.295};
    if(is_equal(type, "W16X31")) return {5.530, 0.440, 5.530, 0.440, 15.020, 0.275};
    if(is_equal(type, "W16X26")) return {5.500, 0.345, 5.500, 0.345, 15.010, 0.250};
    if(is_equal(type, "W14X873")) return {18.800, 5.510, 18.800, 5.510, 12.580, 3.940};
    if(is_equal(type, "W14X808")) return {18.600, 5.120, 18.600, 5.120, 12.560, 3.740};
    if(is_equal(type, "W14X730")) return {17.900, 4.910, 17.900, 4.910, 12.580, 3.070};
    if(is_equal(type, "W14X665")) return {17.700, 4.520, 17.700, 4.520, 12.560, 2.830};
    if(is_equal(type, "W14X605")) return {17.400, 4.160, 17.400, 4.160, 12.580, 2.600};
    if(is_equal(type, "W14X550")) return {17.200, 3.820, 17.200, 3.820, 12.560, 2.380};
    if(is_equal(type, "W14X500")) return {17.000, 3.500, 17.000, 3.500, 12.600, 2.190};
    if(is_equal(type, "W14X455")) return {16.800, 3.210, 16.800, 3.210, 12.580, 2.020};
    if(is_equal(type, "W14X426")) return {16.700, 3.040, 16.700, 3.040, 12.620, 1.880};
    if(is_equal(type, "W14X398")) return {16.600, 2.850, 16.600, 2.850, 12.600, 1.770};
    if(is_equal(type, "W14X370")) return {16.500, 2.660, 16.500, 2.660, 12.580, 1.660};
    if(is_equal(type, "W14X342")) return {16.400, 2.470, 16.400, 2.470, 12.560, 1.540};
    if(is_equal(type, "W14X311")) return {16.200, 2.260, 16.200, 2.260, 12.580, 1.410};
    if(is_equal(type, "W14X283")) return {16.100, 2.070, 16.100, 2.070, 12.560, 1.290};
    if(is_equal(type, "W14X257")) return {16.000, 1.890, 16.000, 1.890, 12.620, 1.180};
    if(is_equal(type, "W14X233")) return {15.900, 1.720, 15.900, 1.720, 12.560, 1.070};
    if(is_equal(type, "W14X211")) return {15.800, 1.560, 15.800, 1.560, 12.580, 0.980};
    if(is_equal(type, "W14X193")) return {15.700, 1.440, 15.700, 1.440, 12.620, 0.890};
    if(is_equal(type, "W14X176")) return {15.700, 1.310, 15.700, 1.310, 12.580, 0.830};
    if(is_equal(type, "W14X159")) return {15.600, 1.190, 15.600, 1.190, 12.620, 0.745};
    if(is_equal(type, "W14X145")) return {15.500, 1.090, 15.500, 1.090, 12.620, 0.680};
    if(is_equal(type, "W14X132")) return {14.700, 1.030, 14.700, 1.030, 12.640, 0.645};
    if(is_equal(type, "W14X120")) return {14.700, 0.940, 14.700, 0.940, 12.620, 0.590};
    if(is_equal(type, "W14X109")) return {14.600, 0.860, 14.600, 0.860, 12.580, 0.525};
    if(is_equal(type, "W14X99")) return {14.600, 0.780, 14.600, 0.780, 12.640, 0.485};
    if(is_equal(type, "W14X90")) return {14.500, 0.710, 14.500, 0.710, 12.580, 0.440};
    if(is_equal(type, "W14X82")) return {10.100, 0.855, 10.100, 0.855, 12.590, 0.510};
    if(is_equal(type, "W14X74")) return {10.100, 0.785, 10.100, 0.785, 12.630, 0.450};
    if(is_equal(type, "W14X68")) return {10.000, 0.720, 10.000, 0.720, 12.560, 0.415};
    if(is_equal(type, "W14X61")) return {10.000, 0.645, 10.000, 0.645, 12.610, 0.375};
    if(is_equal(type, "W14X53")) return {8.060, 0.660, 8.060, 0.660, 12.580, 0.370};
    if(is_equal(type, "W14X48")) return {8.030, 0.595, 8.030, 0.595, 12.610, 0.340};
    if(is_equal(type, "W14X43")) return {8.000, 0.530, 8.000, 0.530, 12.640, 0.305};
    if(is_equal(type, "W14X38")) return {6.770, 0.515, 6.770, 0.515, 13.070, 0.310};
    if(is_equal(type, "W14X34")) return {6.750, 0.455, 6.750, 0.455, 13.090, 0.285};
    if(is_equal(type, "W14X30")) return {6.730, 0.385, 6.730, 0.385, 13.030, 0.270};
    if(is_equal(type, "W14X26")) return {5.030, 0.420, 5.030, 0.420, 13.060, 0.255};
    if(is_equal(type, "W14X22")) return {5.000, 0.335, 5.000, 0.335, 13.030, 0.230};
    if(is_equal(type, "W12X336")) return {13.400, 2.960, 13.400, 2.960, 10.880, 1.780};
    if(is_equal(type, "W12X305")) return {13.200, 2.710, 13.200, 2.710, 10.880, 1.630};
    if(is_equal(type, "W12X279")) return {13.100, 2.470, 13.100, 2.470, 10.960, 1.530};
    if(is_equal(type, "W12X252")) return {13.000, 2.250, 13.000, 2.250, 10.900, 1.400};
    if(is_equal(type, "W12X230")) return {12.900, 2.070, 12.900, 2.070, 10.960, 1.290};
    if(is_equal(type, "W12X210")) return {12.800, 1.900, 12.800, 1.900, 10.900, 1.180};
    if(is_equal(type, "W12X190")) return {12.700, 1.740, 12.700, 1.740, 10.920, 1.060};
    if(is_equal(type, "W12X170")) return {12.600, 1.560, 12.600, 1.560, 10.880, 0.960};
    if(is_equal(type, "W12X152")) return {12.500, 1.400, 12.500, 1.400, 10.900, 0.870};
    if(is_equal(type, "W12X136")) return {12.400, 1.250, 12.400, 1.250, 10.900, 0.790};
    if(is_equal(type, "W12X120")) return {12.300, 1.110, 12.300, 1.110, 10.880, 0.710};
    if(is_equal(type, "W12X106")) return {12.200, 0.990, 12.200, 0.990, 10.920, 0.610};
    if(is_equal(type, "W12X96")) return {12.200, 0.900, 12.200, 0.900, 10.900, 0.550};
    if(is_equal(type, "W12X87")) return {12.100, 0.810, 12.100, 0.810, 10.880, 0.515};
    if(is_equal(type, "W12X79")) return {12.100, 0.735, 12.100, 0.735, 10.930, 0.470};
    if(is_equal(type, "W12X72")) return {12.000, 0.670, 12.000, 0.670, 10.960, 0.430};
    if(is_equal(type, "W12X65")) return {12.000, 0.605, 12.000, 0.605, 10.890, 0.390};
    if(is_equal(type, "W12X58")) return {10.000, 0.640, 10.000, 0.640, 10.920, 0.360};
    if(is_equal(type, "W12X53")) return {10.000, 0.575, 10.000, 0.575, 10.950, 0.345};
    if(is_equal(type, "W12X50")) return {8.080, 0.640, 8.080, 0.640, 10.920, 0.370};
    if(is_equal(type, "W12X45")) return {8.050, 0.575, 8.050, 0.575, 10.950, 0.335};
    if(is_equal(type, "W12X40")) return {8.010, 0.515, 8.010, 0.515, 10.870, 0.295};
    if(is_equal(type, "W12X35")) return {6.560, 0.520, 6.560, 0.520, 11.460, 0.300};
    if(is_equal(type, "W12X30")) return {6.520, 0.440, 6.520, 0.440, 11.420, 0.260};
    if(is_equal(type, "W12X26")) return {6.490, 0.380, 6.490, 0.380, 11.440, 0.230};
    if(is_equal(type, "W12X22")) return {4.030, 0.425, 4.030, 0.425, 11.450, 0.260};
    if(is_equal(type, "W12X19")) return {4.010, 0.350, 4.010, 0.350, 11.500, 0.235};
    if(is_equal(type, "W12X16")) return {3.990, 0.265, 3.990, 0.265, 11.470, 0.220};
    if(is_equal(type, "W12X14")) return {3.970, 0.225, 3.970, 0.225, 11.450, 0.200};
    if(is_equal(type, "W10X112")) return {10.400, 1.250, 10.400, 1.250, 8.900, 0.755};
    if(is_equal(type, "W10X100")) return {10.300, 1.120, 10.300, 1.120, 8.860, 0.680};
    if(is_equal(type, "W10X88")) return {10.300, 0.990, 10.300, 0.990, 8.820, 0.605};
    if(is_equal(type, "W10X77")) return {10.200, 0.870, 10.200, 0.870, 8.860, 0.530};
    if(is_equal(type, "W10X68")) return {10.100, 0.770, 10.100, 0.770, 8.860, 0.470};
    if(is_equal(type, "W10X60")) return {10.100, 0.680, 10.100, 0.680, 8.840, 0.420};
    if(is_equal(type, "W10X54")) return {10.000, 0.615, 10.000, 0.615, 8.870, 0.370};
    if(is_equal(type, "W10X49")) return {10.000, 0.560, 10.000, 0.560, 8.880, 0.340};
    if(is_equal(type, "W10X45")) return {8.020, 0.620, 8.020, 0.620, 8.860, 0.350};
    if(is_equal(type, "W10X39")) return {7.990, 0.530, 7.990, 0.530, 8.860, 0.315};
    if(is_equal(type, "W10X33")) return {7.960, 0.435, 7.960, 0.435, 8.860, 0.290};
    if(is_equal(type, "W10X30")) return {5.810, 0.510, 5.810, 0.510, 9.480, 0.300};
    if(is_equal(type, "W10X26")) return {5.770, 0.440, 5.770, 0.440, 9.420, 0.260};
    if(is_equal(type, "W10X22")) return {5.750, 0.360, 5.750, 0.360, 9.480, 0.240};
    if(is_equal(type, "W10X19")) return {4.020, 0.395, 4.020, 0.395, 9.410, 0.250};
    if(is_equal(type, "W10X17")) return {4.010, 0.330, 4.010, 0.330, 9.440, 0.240};
    if(is_equal(type, "W10X15")) return {4.000, 0.270, 4.000, 0.270, 9.450, 0.230};
    if(is_equal(type, "W10X12")) return {3.960, 0.210, 3.960, 0.210, 9.450, 0.190};
    if(is_equal(type, "W8X67")) return {8.280, 0.935, 8.280, 0.935, 7.130, 0.570};
    if(is_equal(type, "W8X58")) return {8.220, 0.810, 8.220, 0.810, 7.130, 0.510};
    if(is_equal(type, "W8X48")) return {8.110, 0.685, 8.110, 0.685, 7.130, 0.400};
    if(is_equal(type, "W8X40")) return {8.070, 0.560, 8.070, 0.560, 7.130, 0.360};
    if(is_equal(type, "W8X35")) return {8.020, 0.495, 8.020, 0.495, 7.130, 0.310};
    if(is_equal(type, "W8X31")) return {8.000, 0.435, 8.000, 0.435, 7.130, 0.285};
    if(is_equal(type, "W8X28")) return {6.540, 0.465, 6.540, 0.465, 7.130, 0.285};
    if(is_equal(type, "W8X24")) return {6.500, 0.400, 6.500, 0.400, 7.130, 0.245};
    if(is_equal(type, "W8X21")) return {5.270, 0.400, 5.270, 0.400, 7.480, 0.250};
    if(is_equal(type, "W8X18")) return {5.250, 0.330, 5.250, 0.330, 7.480, 0.230};
    if(is_equal(type, "W8X15")) return {4.020, 0.315, 4.020, 0.315, 7.480, 0.245};
    if(is_equal(type, "W8X13")) return {4.000, 0.255, 4.000, 0.255, 7.480, 0.230};
    if(is_equal(type, "W8X10")) return {3.940, 0.205, 3.940, 0.205, 7.480, 0.170};
    if(is_equal(type, "W6X25")) return {6.080, 0.455, 6.080, 0.455, 5.470, 0.320};
    if(is_equal(type, "W6X20")) return {6.020, 0.365, 6.020, 0.365, 5.470, 0.260};
    if(is_equal(type, "W6X15")) return {5.990, 0.260, 5.990, 0.260, 5.470, 0.230};
    if(is_equal(type, "W6X16")) return {4.030, 0.405, 4.030, 0.405, 5.470, 0.260};
    if(is_equal(type, "W6X12")) return {4.000, 0.280, 4.000, 0.280, 5.470, 0.230};
    if(is_equal(type, "W6X9")) return {3.940, 0.215, 3.940, 0.215, 5.470, 0.170};
    if(is_equal(type, "W6X8.5")) return {3.940, 0.195, 3.940, 0.195, 5.440, 0.170};
    if(is_equal(type, "W5X19")) return {5.030, 0.430, 5.030, 0.430, 4.290, 0.270};
    if(is_equal(type, "W5X16")) return {5.000, 0.360, 5.000, 0.360, 4.290, 0.240};
    if(is_equal(type, "W4X13")) return {4.060, 0.345, 4.060, 0.345, 3.470, 0.280};
    if(is_equal(type, "M12.5X12.4")) return {3.750, 0.228, 3.750, 0.228, 12.044, 0.155};
    if(is_equal(type, "M12.5X11.6")) return {3.500, 0.211, 3.500, 0.211, 12.078, 0.155};
    if(is_equal(type, "M12X11.8")) return {3.070, 0.225, 3.070, 0.225, 11.550, 0.177};
    if(is_equal(type, "M12X10.8")) return {3.070, 0.210, 3.070, 0.210, 11.580, 0.160};
    if(is_equal(type, "M12X10")) return {3.250, 0.180, 3.250, 0.180, 11.640, 0.149};
    if(is_equal(type, "M10X9")) return {2.690, 0.206, 2.690, 0.206, 9.588, 0.157};
    if(is_equal(type, "M10X8")) return {2.690, 0.182, 2.690, 0.182, 9.586, 0.141};
    if(is_equal(type, "M10X7.5")) return {2.690, 0.173, 2.690, 0.173, 9.644, 0.130};
    if(is_equal(type, "M8X6.5")) return {2.280, 0.189, 2.280, 0.189, 7.622, 0.135};
    if(is_equal(type, "M8X6.2")) return {2.280, 0.177, 2.280, 0.177, 7.646, 0.129};
    if(is_equal(type, "M6X4.4")) return {1.840, 0.171, 1.840, 0.171, 5.658, 0.114};
    if(is_equal(type, "M6X3.7")) return {2.000, 0.129, 2.000, 0.129, 5.662, 0.098};
    if(is_equal(type, "M5X18.9")) return {5.000, 0.416, 5.000, 0.416, 4.168, 0.316};
    if(is_equal(type, "M4X6")) return {3.800, 0.160, 3.800, 0.160, 3.480, 0.130};
    if(is_equal(type, "M4X4.08")) return {2.250, 0.170, 2.250, 0.170, 3.660, 0.115};
    if(is_equal(type, "M3X2.9")) return {2.250, 0.130, 2.250, 0.130, 2.740, 0.090};
    if(is_equal(type, "S24X121")) return {8.050, 1.090, 8.050, 1.090, 22.320, 0.800};
    if(is_equal(type, "S24X106")) return {7.870, 1.090, 7.870, 1.090, 22.320, 0.620};
    if(is_equal(type, "S24X100")) return {7.250, 0.870, 7.250, 0.870, 22.260, 0.745};
    if(is_equal(type, "S24X90")) return {7.130, 0.870, 7.130, 0.870, 22.260, 0.625};
    if(is_equal(type, "S24X80")) return {7.000, 0.870, 7.000, 0.870, 22.260, 0.500};
    if(is_equal(type, "S20X96")) return {7.200, 0.920, 7.200, 0.920, 18.460, 0.800};
    if(is_equal(type, "S20X86")) return {7.060, 0.920, 7.060, 0.920, 18.460, 0.660};
    if(is_equal(type, "S20X75")) return {6.390, 0.795, 6.390, 0.795, 18.410, 0.635};
    if(is_equal(type, "S20X66")) return {6.260, 0.795, 6.260, 0.795, 18.410, 0.505};
    if(is_equal(type, "S18X70")) return {6.250, 0.691, 6.250, 0.691, 16.618, 0.711};
    if(is_equal(type, "S18X54.7")) return {6.000, 0.691, 6.000, 0.691, 16.618, 0.461};
    if(is_equal(type, "S15X50")) return {5.640, 0.622, 5.640, 0.622, 13.756, 0.550};
    if(is_equal(type, "S15X42.9")) return {5.500, 0.622, 5.500, 0.622, 13.756, 0.411};
    if(is_equal(type, "S12X50")) return {5.480, 0.659, 5.480, 0.659, 10.682, 0.687};
    if(is_equal(type, "S12X40.8")) return {5.250, 0.659, 5.250, 0.659, 10.682, 0.462};
    if(is_equal(type, "S12X35")) return {5.080, 0.544, 5.080, 0.544, 10.912, 0.428};
    if(is_equal(type, "S12X31.8")) return {5.000, 0.544, 5.000, 0.544, 10.912, 0.350};
    if(is_equal(type, "S10X35")) return {4.940, 0.491, 4.940, 0.491, 9.018, 0.594};
    if(is_equal(type, "S10X25.4")) return {4.660, 0.491, 4.660, 0.491, 9.018, 0.311};
    if(is_equal(type, "S8X23")) return {4.170, 0.425, 4.170, 0.425, 7.150, 0.441};
    if(is_equal(type, "S8X18.4")) return {4.000, 0.425, 4.000, 0.425, 7.150, 0.271};
    if(is_equal(type, "S6X17.25")) return {3.570, 0.359, 3.570, 0.359, 5.282, 0.465};
    if(is_equal(type, "S6X12.5")) return {3.330, 0.359, 3.330, 0.359, 5.282, 0.232};
    if(is_equal(type, "S5X10")) return {3.000, 0.326, 3.000, 0.326, 4.348, 0.214};
    if(is_equal(type, "S4X9.5")) return {2.800, 0.293, 2.800, 0.293, 3.414, 0.326};
    if(is_equal(type, "S4X7.7")) return {2.660, 0.293, 2.660, 0.293, 3.414, 0.193};
    if(is_equal(type, "S3X7.5")) return {2.510, 0.260, 2.510, 0.260, 2.480, 0.349};
    if(is_equal(type, "S3X5.7")) return {2.330, 0.260, 2.330, 0.260, 2.480, 0.170};
    if(is_equal(type, "HP18X204")) return {18.100, 1.130, 18.100, 1.130, 16.040, 1.130};
    if(is_equal(type, "HP18X181")) return {18.000, 1.000, 18.000, 1.000, 16.000, 1.000};
    if(is_equal(type, "HP18X157")) return {17.900, 0.870, 17.900, 0.870, 15.960, 0.870};
    if(is_equal(type, "HP18X135")) return {17.800, 0.750, 17.800, 0.750, 16.000, 0.750};
    if(is_equal(type, "HP16X183")) return {16.300, 1.130, 16.300, 1.130, 14.240, 1.130};
    if(is_equal(type, "HP16X162")) return {16.100, 1.000, 16.100, 1.000, 14.300, 1.000};
    if(is_equal(type, "HP16X141")) return {16.000, 0.875, 16.000, 0.875, 14.250, 0.875};
    if(is_equal(type, "HP16X121")) return {15.900, 0.750, 15.900, 0.750, 14.300, 0.750};
    if(is_equal(type, "HP16X101")) return {15.800, 0.625, 15.800, 0.625, 14.250, 0.625};
    if(is_equal(type, "HP16X88")) return {15.700, 0.540, 15.700, 0.540, 14.220, 0.540};
    if(is_equal(type, "HP14X117")) return {14.900, 0.805, 14.900, 0.805, 12.590, 0.805};
    if(is_equal(type, "HP14X102")) return {14.800, 0.705, 14.800, 0.705, 12.590, 0.705};
    if(is_equal(type, "HP14X89")) return {14.700, 0.615, 14.700, 0.615, 12.570, 0.615};
    if(is_equal(type, "HP14X73")) return {14.600, 0.505, 14.600, 0.505, 12.590, 0.505};
    if(is_equal(type, "HP12X89")) return {12.300, 0.720, 12.300, 0.720, 10.960, 0.720};
    if(is_equal(type, "HP12X84")) return {12.300, 0.685, 12.300, 0.685, 10.930, 0.685};
    if(is_equal(type, "HP12X74")) return {12.200, 0.610, 12.200, 0.610, 10.880, 0.605};
    if(is_equal(type, "HP12X63")) return {12.100, 0.515, 12.100, 0.515, 10.870, 0.515};
    if(is_equal(type, "HP12X53")) return {12.000, 0.435, 12.000, 0.435, 10.930, 0.435};
    if(is_equal(type, "HP10X57")) return {10.200, 0.565, 10.200, 0.565, 8.860, 0.565};
    if(is_equal(type, "HP10X42")) return {10.100, 0.420, 10.100, 0.420, 8.860, 0.415};
    if(is_equal(type, "HP8X36")) return {8.160, 0.445, 8.160, 0.445, 7.130, 0.445};
    return {};
}

vec ustsection(const std::string_view type) {
    if(is_equal(type, "WT22X167.5")) return {15.900, 1.770, 20.230, 1.030};
    if(is_equal(type, "WT22X145")) return {15.800, 1.580, 20.220, 0.865};
    if(is_equal(type, "WT22X131")) return {15.800, 1.420, 20.280, 0.785};
    if(is_equal(type, "WT22X115")) return {15.800, 1.220, 20.280, 0.710};
    if(is_equal(type, "WT20X327.5")) return {16.900, 3.540, 18.260, 1.970};
    if(is_equal(type, "WT20X296.5")) return {16.700, 3.230, 18.270, 1.790};
    if(is_equal(type, "WT20X251.5")) return {16.400, 2.760, 18.240, 1.540};
    if(is_equal(type, "WT20X215.5")) return {16.200, 2.360, 18.240, 1.340};
    if(is_equal(type, "WT20X198.5")) return {16.100, 2.200, 18.300, 1.220};
    if(is_equal(type, "WT20X186")) return {16.100, 2.050, 18.250, 1.160};
    if(is_equal(type, "WT20X181")) return {16.000, 2.010, 18.290, 1.120};
    if(is_equal(type, "WT20X162")) return {15.900, 1.810, 18.290, 1.000};
    if(is_equal(type, "WT20X148.5")) return {15.800, 1.650, 18.250, 0.930};
    if(is_equal(type, "WT20X138.5")) return {15.800, 1.580, 18.220, 0.830};
    if(is_equal(type, "WT20X124.5")) return {15.800, 1.420, 18.280, 0.750};
    if(is_equal(type, "WT20X107.5")) return {15.800, 1.220, 18.280, 0.650};
    if(is_equal(type, "WT20X99.5")) return {15.800, 1.070, 18.230, 0.650};
    if(is_equal(type, "WT20X196")) return {12.400, 2.520, 18.280, 1.420};
    if(is_equal(type, "WT20X165.5")) return {12.200, 2.130, 18.270, 1.220};
    if(is_equal(type, "WT20X163.5")) return {12.100, 2.130, 18.270, 1.180};
    if(is_equal(type, "WT20X147")) return {12.000, 1.930, 18.270, 1.060};
    if(is_equal(type, "WT20X139")) return {12.000, 1.810, 18.290, 1.030};
    if(is_equal(type, "WT20X132")) return {11.900, 1.730, 18.270, 0.960};
    if(is_equal(type, "WT20X117.5")) return {11.900, 1.580, 18.220, 0.830};
    if(is_equal(type, "WT20X105.5")) return {11.800, 1.420, 18.280, 0.750};
    if(is_equal(type, "WT20X91.5")) return {11.800, 1.200, 18.300, 0.650};
    if(is_equal(type, "WT20X83.5")) return {11.800, 1.030, 18.270, 0.650};
    if(is_equal(type, "WT20X74.5")) return {11.800, 0.830, 18.270, 0.630};
    if(is_equal(type, "WT18X462.5")) return {18.600, 4.530, 17.070, 3.020};
    if(is_equal(type, "WT18X426.5")) return {18.200, 4.530, 17.070, 2.520};
    if(is_equal(type, "WT18X401")) return {18.000, 4.290, 17.010, 2.380};
    if(is_equal(type, "WT18X361.5")) return {17.800, 3.900, 17.000, 2.170};
    if(is_equal(type, "WT18X326")) return {17.600, 3.540, 16.960, 1.970};
    if(is_equal(type, "WT18X264.5")) return {17.200, 2.910, 16.990, 1.610};
    if(is_equal(type, "WT18X243.5")) return {17.100, 2.680, 17.020, 1.500};
    if(is_equal(type, "WT18X220.5")) return {17.000, 2.440, 16.960, 1.360};
    if(is_equal(type, "WT18X197.5")) return {16.800, 2.200, 17.000, 1.220};
    if(is_equal(type, "WT18X180.5")) return {16.700, 2.010, 16.990, 1.120};
    if(is_equal(type, "WT18X165")) return {16.600, 1.850, 16.950, 1.020};
    if(is_equal(type, "WT18X151")) return {16.700, 1.680, 17.020, 0.945};
    if(is_equal(type, "WT18X141")) return {16.600, 1.570, 17.030, 0.885};
    if(is_equal(type, "WT18X131")) return {16.600, 1.440, 16.960, 0.840};
    if(is_equal(type, "WT18X123.5")) return {16.500, 1.350, 16.950, 0.800};
    if(is_equal(type, "WT18X115.5")) return {16.500, 1.260, 16.940, 0.760};
    if(is_equal(type, "WT18X128")) return {12.200, 1.730, 16.970, 0.960};
    if(is_equal(type, "WT18X116")) return {12.100, 1.570, 17.030, 0.870};
    if(is_equal(type, "WT18X105")) return {12.200, 1.360, 16.940, 0.830};
    if(is_equal(type, "WT18X97")) return {12.100, 1.260, 16.940, 0.765};
    if(is_equal(type, "WT18X91")) return {12.100, 1.180, 17.020, 0.725};
    if(is_equal(type, "WT18X85")) return {12.000, 1.100, 17.000, 0.680};
    if(is_equal(type, "WT18X80")) return {12.000, 1.020, 16.980, 0.650};
    if(is_equal(type, "WT18X75")) return {12.000, 0.940, 16.960, 0.625};
    if(is_equal(type, "WT18X67.5")) return {12.000, 0.790, 17.010, 0.600};
    if(is_equal(type, "WT16.5X193.5")) return {16.200, 2.280, 15.720, 1.260};
    if(is_equal(type, "WT16.5X177")) return {16.100, 2.090, 15.710, 1.160};
    if(is_equal(type, "WT16.5X159")) return {16.000, 1.890, 15.710, 1.040};
    if(is_equal(type, "WT16.5X145.5")) return {15.900, 1.730, 15.670, 0.960};
    if(is_equal(type, "WT16.5X131.5")) return {15.800, 1.570, 15.730, 0.870};
    if(is_equal(type, "WT16.5X120.5")) return {15.900, 1.400, 15.700, 0.830};
    if(is_equal(type, "WT16.5X110.5")) return {15.800, 1.280, 15.720, 0.775};
    if(is_equal(type, "WT16.5X100.5")) return {15.700, 1.150, 15.650, 0.715};
    if(is_equal(type, "WT16.5X84.5")) return {11.500, 1.220, 15.680, 0.670};
    if(is_equal(type, "WT16.5X76")) return {11.600, 1.060, 15.640, 0.635};
    if(is_equal(type, "WT16.5X70.5")) return {11.500, 0.960, 15.740, 0.605};
    if(is_equal(type, "WT16.5X65")) return {11.500, 0.855, 15.645, 0.580};
    if(is_equal(type, "WT16.5X59")) return {11.500, 0.740, 15.660, 0.550};
    if(is_equal(type, "WT15X195.5")) return {15.600, 2.440, 14.160, 1.360};
    if(is_equal(type, "WT15X178.5")) return {15.500, 2.240, 14.160, 1.240};
    if(is_equal(type, "WT15X163")) return {15.400, 2.050, 14.150, 1.140};
    if(is_equal(type, "WT15X146")) return {15.300, 1.850, 14.150, 1.020};
    if(is_equal(type, "WT15X130.5")) return {15.200, 1.650, 14.150, 0.930};
    if(is_equal(type, "WT15X117.5")) return {15.100, 1.500, 14.200, 0.830};
    if(is_equal(type, "WT15X105.5")) return {15.100, 1.320, 14.180, 0.775};
    if(is_equal(type, "WT15X95.5")) return {15.000, 1.190, 14.110, 0.710};
    if(is_equal(type, "WT15X86.5")) return {15.000, 1.070, 14.130, 0.655};
    if(is_equal(type, "WT15X74")) return {10.500, 1.180, 14.120, 0.650};
    if(is_equal(type, "WT15X66")) return {10.500, 1.000, 14.200, 0.615};
    if(is_equal(type, "WT15X62")) return {10.500, 0.930, 14.170, 0.585};
    if(is_equal(type, "WT15X58")) return {10.500, 0.850, 14.150, 0.565};
    if(is_equal(type, "WT15X54")) return {10.500, 0.760, 14.140, 0.545};
    if(is_equal(type, "WT15X49.5")) return {10.500, 0.670, 14.130, 0.520};
    if(is_equal(type, "WT15X45")) return {10.400, 0.610, 14.190, 0.470};
    if(is_equal(type, "WT13.5X269.5")) return {15.300, 3.540, 12.760, 1.970};
    if(is_equal(type, "WT13.5X184")) return {14.700, 2.480, 12.720, 1.380};
    if(is_equal(type, "WT13.5X168")) return {14.600, 2.280, 12.720, 1.260};
    if(is_equal(type, "WT13.5X153.5")) return {14.400, 2.090, 12.710, 1.160};
    if(is_equal(type, "WT13.5X140.5")) return {14.400, 1.930, 12.670, 1.060};
    if(is_equal(type, "WT13.5X129")) return {14.300, 1.770, 12.730, 0.980};
    if(is_equal(type, "WT13.5X117.5")) return {14.200, 1.610, 12.690, 0.910};
    if(is_equal(type, "WT13.5X108.5")) return {14.100, 1.500, 12.700, 0.830};
    if(is_equal(type, "WT13.5X97")) return {14.000, 1.340, 12.760, 0.750};
    if(is_equal(type, "WT13.5X89")) return {14.100, 1.190, 12.710, 0.725};
    if(is_equal(type, "WT13.5X80.5")) return {14.000, 1.080, 12.720, 0.660};
    if(is_equal(type, "WT13.5X73")) return {14.000, 0.975, 12.725, 0.605};
    if(is_equal(type, "WT13.5X64.5")) return {10.000, 1.100, 12.700, 0.610};
    if(is_equal(type, "WT13.5X57")) return {10.100, 0.930, 12.670, 0.570};
    if(is_equal(type, "WT13.5X51")) return {10.000, 0.830, 12.670, 0.515};
    if(is_equal(type, "WT13.5X47")) return {10.000, 0.745, 12.755, 0.490};
    if(is_equal(type, "WT13.5X42")) return {10.000, 0.640, 12.760, 0.460};
    if(is_equal(type, "WT12X185")) return {13.700, 2.720, 11.280, 1.520};
    if(is_equal(type, "WT12X167.5")) return {13.500, 2.480, 11.320, 1.380};
    if(is_equal(type, "WT12X153")) return {13.400, 2.280, 11.320, 1.260};
    if(is_equal(type, "WT12X139.5")) return {13.300, 2.090, 11.310, 1.160};
    if(is_equal(type, "WT12X125")) return {13.200, 1.890, 11.310, 1.040};
    if(is_equal(type, "WT12X114.5")) return {13.100, 1.730, 11.270, 0.960};
    if(is_equal(type, "WT12X103.5")) return {13.000, 1.570, 11.330, 0.870};
    if(is_equal(type, "WT12X96")) return {13.000, 1.460, 11.240, 0.810};
    if(is_equal(type, "WT12X88")) return {12.900, 1.340, 11.260, 0.750};
    if(is_equal(type, "WT12X81")) return {13.000, 1.220, 11.280, 0.705};
    if(is_equal(type, "WT12X73")) return {12.900, 1.090, 11.310, 0.650};
    if(is_equal(type, "WT12X65.5")) return {12.900, 0.960, 11.240, 0.605};
    if(is_equal(type, "WT12X58.5")) return {12.800, 0.850, 11.250, 0.550};
    if(is_equal(type, "WT12X52")) return {12.800, 0.750, 11.250, 0.500};
    if(is_equal(type, "WT12X51.5")) return {9.000, 0.980, 11.320, 0.550};
    if(is_equal(type, "WT12X47")) return {9.070, 0.875, 11.325, 0.515};
    if(is_equal(type, "WT12X42")) return {9.020, 0.770, 11.330, 0.470};
    if(is_equal(type, "WT12X38")) return {8.990, 0.680, 11.320, 0.440};
    if(is_equal(type, "WT12X34")) return {8.970, 0.585, 11.315, 0.415};
    if(is_equal(type, "WT12X31")) return {7.040, 0.590, 11.310, 0.430};
    if(is_equal(type, "WT12X27.5")) return {7.010, 0.505, 11.295, 0.395};
    if(is_equal(type, "WT10.5X137.5")) return {12.900, 2.190, 9.910, 1.220};
    if(is_equal(type, "WT10.5X124")) return {12.800, 1.990, 9.910, 1.100};
    if(is_equal(type, "WT10.5X111.5")) return {12.700, 1.790, 9.910, 1.000};
    if(is_equal(type, "WT10.5X100.5")) return {12.600, 1.630, 9.870, 0.910};
    if(is_equal(type, "WT10.5X91")) return {12.500, 1.480, 9.920, 0.830};
    if(is_equal(type, "WT10.5X83")) return {12.400, 1.360, 9.840, 0.750};
    if(is_equal(type, "WT10.5X73.5")) return {12.500, 1.150, 9.850, 0.720};
    if(is_equal(type, "WT10.5X66")) return {12.400, 1.040, 9.860, 0.650};
    if(is_equal(type, "WT10.5X61")) return {12.400, 0.960, 9.840, 0.600};
    if(is_equal(type, "WT10.5X55.5")) return {12.300, 0.875, 9.925, 0.550};
    if(is_equal(type, "WT10.5X50.5")) return {12.300, 0.800, 9.900, 0.500};
    if(is_equal(type, "WT10.5X46.5")) return {8.420, 0.930, 9.870, 0.580};
    if(is_equal(type, "WT10.5X41.5")) return {8.360, 0.835, 9.865, 0.515};
    if(is_equal(type, "WT10.5X36.5")) return {8.300, 0.740, 9.860, 0.455};
    if(is_equal(type, "WT10.5X34")) return {8.270, 0.685, 9.915, 0.430};
    if(is_equal(type, "WT10.5X31")) return {8.240, 0.615, 9.885, 0.400};
    if(is_equal(type, "WT10.5X27.5")) return {8.220, 0.522, 9.878, 0.375};
    if(is_equal(type, "WT10.5X24")) return {8.140, 0.430, 9.870, 0.350};
    if(is_equal(type, "WT10.5X28.5")) return {6.560, 0.650, 9.850, 0.405};
    if(is_equal(type, "WT10.5X25")) return {6.530, 0.535, 9.865, 0.380};
    if(is_equal(type, "WT10.5X22")) return {6.500, 0.450, 9.850, 0.350};
    if(is_equal(type, "WT9X155.5")) return {12.000, 2.740, 8.460, 1.520};
    if(is_equal(type, "WT9X141.5")) return {11.900, 2.500, 8.400, 1.400};
    if(is_equal(type, "WT9X129")) return {11.800, 2.300, 8.400, 1.280};
    if(is_equal(type, "WT9X117")) return {11.700, 2.110, 8.390, 1.160};
    if(is_equal(type, "WT9X105.5")) return {11.600, 1.910, 8.390, 1.060};
    if(is_equal(type, "WT9X96")) return {11.500, 1.750, 8.450, 0.960};
    if(is_equal(type, "WT9X87.5")) return {11.400, 1.590, 8.410, 0.890};
    if(is_equal(type, "WT9X79")) return {11.300, 1.440, 8.420, 0.810};
    if(is_equal(type, "WT9X71.5")) return {11.200, 1.320, 8.430, 0.730};
    if(is_equal(type, "WT9X65")) return {11.200, 1.200, 8.430, 0.670};
    if(is_equal(type, "WT9X59.5")) return {11.300, 1.060, 8.430, 0.655};
    if(is_equal(type, "WT9X53")) return {11.200, 0.940, 8.430, 0.590};
    if(is_equal(type, "WT9X48.5")) return {11.100, 0.870, 8.430, 0.535};
    if(is_equal(type, "WT9X43")) return {11.100, 0.770, 8.430, 0.480};
    if(is_equal(type, "WT9X38")) return {11.000, 0.680, 8.430, 0.425};
    if(is_equal(type, "WT9X35.5")) return {7.640, 0.810, 8.430, 0.495};
    if(is_equal(type, "WT9X32.5")) return {7.590, 0.750, 8.430, 0.450};
    if(is_equal(type, "WT9X30")) return {7.560, 0.695, 8.425, 0.415};
    if(is_equal(type, "WT9X27.5")) return {7.530, 0.630, 8.430, 0.390};
    if(is_equal(type, "WT9X25")) return {7.500, 0.570, 8.430, 0.355};
    if(is_equal(type, "WT9X23")) return {6.060, 0.605, 8.425, 0.360};
    if(is_equal(type, "WT9X20")) return {6.020, 0.525, 8.425, 0.315};
    if(is_equal(type, "WT9X17.5")) return {6.000, 0.425, 8.425, 0.300};
    if(is_equal(type, "WT8X50")) return {10.400, 0.985, 7.505, 0.585};
    if(is_equal(type, "WT8X44.5")) return {10.400, 0.875, 7.505, 0.525};
    if(is_equal(type, "WT8X38.5")) return {10.300, 0.760, 7.500, 0.455};
    if(is_equal(type, "WT8X33.5")) return {10.200, 0.665, 7.505, 0.395};
    if(is_equal(type, "WT8X28.5")) return {7.120, 0.715, 7.505, 0.430};
    if(is_equal(type, "WT8X25")) return {7.070, 0.630, 7.500, 0.380};
    if(is_equal(type, "WT8X22.5")) return {7.040, 0.565, 7.505, 0.345};
    if(is_equal(type, "WT8X20")) return {7.000, 0.505, 7.505, 0.305};
    if(is_equal(type, "WT8X18")) return {6.990, 0.430, 7.500, 0.295};
    if(is_equal(type, "WT8X15.5")) return {5.530, 0.440, 7.500, 0.275};
    if(is_equal(type, "WT8X13")) return {5.500, 0.345, 7.505, 0.250};
    if(is_equal(type, "WT7X436.5")) return {18.800, 5.510, 6.290, 3.940};
    if(is_equal(type, "WT7X404")) return {18.600, 5.120, 6.280, 3.740};
    if(is_equal(type, "WT7X365")) return {17.900, 4.910, 6.290, 3.070};
    if(is_equal(type, "WT7X332.5")) return {17.700, 4.520, 6.280, 2.830};
    if(is_equal(type, "WT7X302.5")) return {17.400, 4.160, 6.340, 2.600};
    if(is_equal(type, "WT7X275")) return {17.200, 3.820, 6.280, 2.380};
    if(is_equal(type, "WT7X250")) return {17.000, 3.500, 6.300, 2.190};
    if(is_equal(type, "WT7X227.5")) return {16.800, 3.210, 6.300, 2.020};
    if(is_equal(type, "WT7X213")) return {16.700, 3.040, 6.300, 1.880};
    if(is_equal(type, "WT7X199")) return {16.600, 2.850, 6.300, 1.770};
    if(is_equal(type, "WT7X185")) return {16.500, 2.660, 6.300, 1.660};
    if(is_equal(type, "WT7X171")) return {16.400, 2.470, 6.300, 1.540};
    if(is_equal(type, "WT7X155.5")) return {16.200, 2.260, 6.300, 1.410};
    if(is_equal(type, "WT7X141.5")) return {16.100, 2.070, 6.300, 1.290};
    if(is_equal(type, "WT7X128.5")) return {16.000, 1.890, 6.300, 1.180};
    if(is_equal(type, "WT7X116.5")) return {15.900, 1.720, 6.300, 1.070};
    if(is_equal(type, "WT7X105.5")) return {15.800, 1.560, 6.300, 0.980};
    if(is_equal(type, "WT7X96.5")) return {15.700, 1.440, 6.300, 0.890};
    if(is_equal(type, "WT7X88")) return {15.700, 1.310, 6.300, 0.830};
    if(is_equal(type, "WT7X79.5")) return {15.600, 1.190, 6.300, 0.745};
    if(is_equal(type, "WT7X72.5")) return {15.500, 1.090, 6.300, 0.680};
    if(is_equal(type, "WT7X66")) return {14.700, 1.030, 6.300, 0.645};
    if(is_equal(type, "WT7X60")) return {14.700, 0.940, 6.300, 0.590};
    if(is_equal(type, "WT7X54.5")) return {14.600, 0.860, 6.300, 0.525};
    if(is_equal(type, "WT7X49.5")) return {14.600, 0.780, 6.300, 0.485};
    if(is_equal(type, "WT7X45")) return {14.500, 0.710, 6.300, 0.440};
    if(is_equal(type, "WT7X41")) return {10.100, 0.855, 6.305, 0.510};
    if(is_equal(type, "WT7X37")) return {10.100, 0.785, 6.305, 0.450};
    if(is_equal(type, "WT7X34")) return {10.000, 0.720, 6.300, 0.415};
    if(is_equal(type, "WT7X30.5")) return {10.000, 0.645, 6.305, 0.375};
    if(is_equal(type, "WT7X26.5")) return {8.060, 0.660, 6.300, 0.370};
    if(is_equal(type, "WT7X24")) return {8.030, 0.595, 6.305, 0.340};
    if(is_equal(type, "WT7X21.5")) return {8.000, 0.530, 6.300, 0.305};
    if(is_equal(type, "WT7X19")) return {6.770, 0.515, 6.535, 0.310};
    if(is_equal(type, "WT7X17")) return {6.750, 0.455, 6.535, 0.285};
    if(is_equal(type, "WT7X15")) return {6.730, 0.385, 6.535, 0.270};
    if(is_equal(type, "WT7X13")) return {5.030, 0.420, 6.540, 0.255};
    if(is_equal(type, "WT7X11")) return {5.000, 0.335, 6.535, 0.230};
    if(is_equal(type, "WT6X168")) return {13.400, 2.960, 5.450, 1.780};
    if(is_equal(type, "WT6X152.5")) return {13.200, 2.710, 5.450, 1.630};
    if(is_equal(type, "WT6X139.5")) return {13.100, 2.470, 5.460, 1.530};
    if(is_equal(type, "WT6X126")) return {13.000, 2.250, 5.460, 1.400};
    if(is_equal(type, "WT6X115")) return {12.900, 2.070, 5.460, 1.290};
    if(is_equal(type, "WT6X105")) return {12.800, 1.900, 5.460, 1.180};
    if(is_equal(type, "WT6X95")) return {12.700, 1.740, 5.450, 1.060};
    if(is_equal(type, "WT6X85")) return {12.600, 1.560, 5.460, 0.960};
    if(is_equal(type, "WT6X76")) return {12.500, 1.400, 5.460, 0.870};
    if(is_equal(type, "WT6X68")) return {12.400, 1.250, 5.460, 0.790};
    if(is_equal(type, "WT6X60")) return {12.300, 1.110, 5.450, 0.710};
    if(is_equal(type, "WT6X53")) return {12.200, 0.990, 5.460, 0.610};
    if(is_equal(type, "WT6X48")) return {12.200, 0.900, 5.460, 0.550};
    if(is_equal(type, "WT6X43.5")) return {12.100, 0.810, 5.460, 0.515};
    if(is_equal(type, "WT6X39.5")) return {12.100, 0.735, 5.455, 0.470};
    if(is_equal(type, "WT6X36")) return {12.000, 0.670, 5.460, 0.430};
    if(is_equal(type, "WT6X32.5")) return {12.000, 0.605, 5.455, 0.390};
    if(is_equal(type, "WT6X29")) return {10.000, 0.640, 5.460, 0.360};
    if(is_equal(type, "WT6X26.5")) return {10.000, 0.575, 5.455, 0.345};
    if(is_equal(type, "WT6X25")) return {8.080, 0.640, 5.460, 0.370};
    if(is_equal(type, "WT6X22.5")) return {8.050, 0.575, 5.455, 0.335};
    if(is_equal(type, "WT6X20")) return {8.010, 0.515, 5.455, 0.295};
    if(is_equal(type, "WT6X17.5")) return {6.560, 0.520, 5.730, 0.300};
    if(is_equal(type, "WT6X15")) return {6.520, 0.440, 5.730, 0.260};
    if(is_equal(type, "WT6X13")) return {6.490, 0.380, 5.730, 0.230};
    if(is_equal(type, "WT6X11")) return {4.030, 0.425, 5.735, 0.260};
    if(is_equal(type, "WT6X9.5")) return {4.010, 0.350, 5.730, 0.235};
    if(is_equal(type, "WT6X8")) return {3.990, 0.265, 5.735, 0.220};
    if(is_equal(type, "WT6X7")) return {3.970, 0.225, 5.735, 0.200};
    if(is_equal(type, "WT5X56")) return {10.400, 1.250, 4.430, 0.755};
    if(is_equal(type, "WT5X50")) return {10.300, 1.120, 4.430, 0.680};
    if(is_equal(type, "WT5X44")) return {10.300, 0.990, 4.430, 0.605};
    if(is_equal(type, "WT5X38.5")) return {10.200, 0.870, 4.430, 0.530};
    if(is_equal(type, "WT5X34")) return {10.100, 0.770, 4.430, 0.470};
    if(is_equal(type, "WT5X30")) return {10.100, 0.680, 4.430, 0.420};
    if(is_equal(type, "WT5X27")) return {10.000, 0.615, 4.435, 0.370};
    if(is_equal(type, "WT5X24.5")) return {10.000, 0.560, 4.430, 0.340};
    if(is_equal(type, "WT5X22.5")) return {8.020, 0.620, 4.430, 0.350};
    if(is_equal(type, "WT5X19.5")) return {7.990, 0.530, 4.430, 0.315};
    if(is_equal(type, "WT5X16.5")) return {7.960, 0.435, 4.435, 0.290};
    if(is_equal(type, "WT5X15")) return {5.810, 0.510, 4.730, 0.300};
    if(is_equal(type, "WT5X13")) return {5.770, 0.440, 4.730, 0.260};
    if(is_equal(type, "WT5X11")) return {5.750, 0.360, 4.730, 0.240};
    if(is_equal(type, "WT5X9.5")) return {4.020, 0.395, 4.725, 0.250};
    if(is_equal(type, "WT5X8.5")) return {4.010, 0.330, 4.730, 0.240};
    if(is_equal(type, "WT5X7.5")) return {4.000, 0.270, 4.730, 0.230};
    if(is_equal(type, "WT5X6")) return {3.960, 0.210, 4.730, 0.190};
    if(is_equal(type, "WT4X33.5")) return {8.280, 0.935, 3.565, 0.570};
    if(is_equal(type, "WT4X29")) return {8.220, 0.810, 3.570, 0.510};
    if(is_equal(type, "WT4X24")) return {8.110, 0.685, 3.565, 0.400};
    if(is_equal(type, "WT4X20")) return {8.070, 0.560, 3.570, 0.360};
    if(is_equal(type, "WT4X17.5")) return {8.020, 0.495, 3.565, 0.310};
    if(is_equal(type, "WT4X15.5")) return {8.000, 0.435, 3.565, 0.285};
    if(is_equal(type, "WT4X14")) return {6.540, 0.465, 3.565, 0.285};
    if(is_equal(type, "WT4X12")) return {6.500, 0.400, 3.570, 0.245};
    if(is_equal(type, "WT4X10.5")) return {5.270, 0.400, 3.740, 0.250};
    if(is_equal(type, "WT4X9")) return {5.250, 0.330, 3.740, 0.230};
    if(is_equal(type, "WT4X7.5")) return {4.020, 0.315, 3.745, 0.245};
    if(is_equal(type, "WT4X6.5")) return {4.000, 0.255, 3.745, 0.230};
    if(is_equal(type, "WT4X5")) return {3.940, 0.205, 3.745, 0.170};
    if(is_equal(type, "WT3X12.5")) return {6.080, 0.455, 2.735, 0.320};
    if(is_equal(type, "WT3X10")) return {6.020, 0.365, 2.735, 0.260};
    if(is_equal(type, "WT3X7.5")) return {5.990, 0.260, 2.740, 0.230};
    if(is_equal(type, "WT3X8")) return {4.030, 0.405, 2.735, 0.260};
    if(is_equal(type, "WT3X6")) return {4.000, 0.280, 2.740, 0.230};
    if(is_equal(type, "WT3X4.5")) return {3.940, 0.215, 2.735, 0.170};
    if(is_equal(type, "WT3X4.25")) return {3.940, 0.195, 2.725, 0.170};
    if(is_equal(type, "WT2.5X9.5")) return {5.030, 0.430, 2.150, 0.270};
    if(is_equal(type, "WT2.5X8")) return {5.000, 0.360, 2.150, 0.240};
    if(is_equal(type, "WT2X6.5")) return {4.060, 0.345, 1.735, 0.280};
    if(is_equal(type, "MT6.25X6.2")) return {3.750, 0.228, 6.042, 0.155};
    if(is_equal(type, "MT6.25X5.8")) return {3.500, 0.211, 6.039, 0.155};
    if(is_equal(type, "MT6X5.9")) return {3.070, 0.225, 5.775, 0.177};
    if(is_equal(type, "MT6X5.4")) return {3.070, 0.210, 5.780, 0.160};
    if(is_equal(type, "MT6X5")) return {3.250, 0.180, 5.810, 0.149};
    if(is_equal(type, "MT5X4.5")) return {2.690, 0.206, 4.794, 0.157};
    if(is_equal(type, "MT5X4")) return {2.690, 0.182, 4.798, 0.141};
    if(is_equal(type, "MT5X3.75")) return {2.690, 0.173, 4.827, 0.130};
    if(is_equal(type, "MT4X3.25")) return {2.280, 0.189, 3.811, 0.135};
    if(is_equal(type, "MT4X3.1")) return {2.280, 0.177, 3.823, 0.129};
    if(is_equal(type, "MT3X2.2")) return {1.840, 0.171, 2.829, 0.114};
    if(is_equal(type, "MT3X1.85")) return {2.000, 0.129, 2.831, 0.098};
    if(is_equal(type, "MT2.5X9.45")) return {5.000, 0.416, 2.084, 0.316};
    if(is_equal(type, "MT2X3")) return {3.800, 0.160, 1.740, 0.130};
    if(is_equal(type, "ST12X60.5")) return {8.050, 1.090, 11.210, 0.800};
    if(is_equal(type, "ST12X53")) return {7.870, 1.090, 11.210, 0.620};
    if(is_equal(type, "ST12X50")) return {7.250, 0.870, 11.130, 0.745};
    if(is_equal(type, "ST12X45")) return {7.130, 0.870, 11.130, 0.625};
    if(is_equal(type, "ST12X40")) return {7.000, 0.870, 11.130, 0.500};
    if(is_equal(type, "ST10X48")) return {7.200, 0.920, 9.280, 0.800};
    if(is_equal(type, "ST10X43")) return {7.060, 0.920, 9.280, 0.660};
    if(is_equal(type, "ST10X37.5")) return {6.390, 0.795, 9.205, 0.635};
    if(is_equal(type, "ST10X33")) return {6.260, 0.795, 9.205, 0.505};
    if(is_equal(type, "ST9X35")) return {6.250, 0.691, 8.309, 0.711};
    if(is_equal(type, "ST9X27.35")) return {6.000, 0.691, 8.309, 0.461};
    if(is_equal(type, "ST7.5X25")) return {5.640, 0.622, 6.878, 0.550};
    if(is_equal(type, "ST7.5X21.45")) return {5.500, 0.622, 6.878, 0.411};
    if(is_equal(type, "ST6X25")) return {5.480, 0.659, 5.341, 0.687};
    if(is_equal(type, "ST6X20.4")) return {5.250, 0.659, 5.341, 0.462};
    if(is_equal(type, "ST6X17.5")) return {5.080, 0.544, 5.456, 0.428};
    if(is_equal(type, "ST6X15.9")) return {5.000, 0.544, 5.456, 0.350};
    if(is_equal(type, "ST5X17.5")) return {4.940, 0.491, 4.509, 0.594};
    if(is_equal(type, "ST5X12.7")) return {4.660, 0.491, 4.509, 0.311};
    if(is_equal(type, "ST4X11.5")) return {4.170, 0.425, 3.575, 0.441};
    if(is_equal(type, "ST4X9.2")) return {4.000, 0.425, 3.575, 0.271};
    if(is_equal(type, "ST3X8.6")) return {3.570, 0.359, 2.641, 0.465};
    if(is_equal(type, "ST3X6.25")) return {3.330, 0.359, 2.641, 0.232};
    if(is_equal(type, "ST2.5X5")) return {3.000, 0.326, 2.174, 0.214};
    if(is_equal(type, "ST2X4.75")) return {2.800, 0.293, 1.707, 0.326};
    if(is_equal(type, "ST2X3.85")) return {2.660, 0.293, 1.707, 0.193};
    if(is_equal(type, "ST1.5X3.75")) return {2.510, 0.260, 1.240, 0.349};
    if(is_equal(type, "ST1.5X2.85")) return {2.330, 0.260, 1.240, 0.170};
    return {};
}

vec usrhssection(const std::string_view type) {
    if(is_equal(type, "HSS34X10X1")) return {8.605, 32.600, 0.930};
    if(is_equal(type, "HSS34X10X7/8")) return {8.780, 32.800, 0.814};
    if(is_equal(type, "HSS34X10X3/4")) return {8.955, 32.950, 0.698};
    if(is_equal(type, "HSS34X10X5/8")) return {9.130, 33.150, 0.581};
    if(is_equal(type, "HSS30X10X1")) return {8.605, 28.600, 0.930};
    if(is_equal(type, "HSS30X10X7/8")) return {8.780, 28.800, 0.814};
    if(is_equal(type, "HSS30X10X3/4")) return {8.955, 28.950, 0.698};
    if(is_equal(type, "HSS30X10X5/8")) return {9.130, 29.150, 0.581};
    if(is_equal(type, "HSS30X10X1/2")) return {9.305, 29.300, 0.465};
    if(is_equal(type, "HSS24X20X3/4")) return {18.950, 22.950, 0.698};
    if(is_equal(type, "HSS24X20X5/8")) return {19.150, 23.150, 0.581};
    if(is_equal(type, "HSS24X20X1/2")) return {19.300, 23.300, 0.465};
    if(is_equal(type, "HSS24X20X3/8")) return {19.500, 23.500, 0.349};
    if(is_equal(type, "HSS24X20X5/16")) return {19.550, 23.550, 0.291};
    if(is_equal(type, "HSS24X18X3/4")) return {16.950, 22.950, 0.698};
    if(is_equal(type, "HSS24X18X5/8")) return {17.150, 23.150, 0.581};
    if(is_equal(type, "HSS24X18X1/2")) return {17.300, 23.300, 0.465};
    if(is_equal(type, "HSS24X18X3/8")) return {17.500, 23.500, 0.349};
    if(is_equal(type, "HSS24X18X5/16")) return {17.550, 23.550, 0.291};
    if(is_equal(type, "HSS24X16X3/4")) return {14.950, 22.950, 0.698};
    if(is_equal(type, "HSS24X16X5/8")) return {15.150, 23.150, 0.581};
    if(is_equal(type, "HSS24X16X1/2")) return {15.300, 23.300, 0.465};
    if(is_equal(type, "HSS24X16X3/8")) return {15.500, 23.500, 0.349};
    if(is_equal(type, "HSS24X16X5/16")) return {15.550, 23.550, 0.291};
    if(is_equal(type, "HSS24X14X3/4")) return {12.950, 22.950, 0.698};
    if(is_equal(type, "HSS24X14X5/8")) return {13.150, 23.150, 0.581};
    if(is_equal(type, "HSS24X14X1/2")) return {13.300, 23.300, 0.465};
    if(is_equal(type, "HSS24X14X3/8")) return {13.500, 23.500, 0.349};
    if(is_equal(type, "HSS24X14X5/16")) return {13.550, 23.550, 0.291};
    if(is_equal(type, "HSS24X14X1/4")) return {13.650, 23.650, 0.233};
    if(is_equal(type, "HSS24X12X1")) return {10.605, 22.600, 0.930};
    if(is_equal(type, "HSS24X12X7/8")) return {10.780, 22.800, 0.814};
    if(is_equal(type, "HSS24X12X3/4")) return {10.955, 22.950, 0.698};
    if(is_equal(type, "HSS24X12X5/8")) return {11.150, 23.150, 0.581};
    if(is_equal(type, "HSS24X12X1/2")) return {11.300, 23.300, 0.465};
    if(is_equal(type, "HSS24X12X3/8")) return {11.500, 23.500, 0.349};
    if(is_equal(type, "HSS24X12X5/16")) return {11.550, 23.550, 0.291};
    if(is_equal(type, "HSS24X12X1/4")) return {11.650, 23.650, 0.233};
    if(is_equal(type, "HSS24X8X1/2")) return {7.305, 23.300, 0.465};
    if(is_equal(type, "HSS24X8X3/8")) return {7.475, 23.500, 0.349};
    if(is_equal(type, "HSS24X8X5/16")) return {7.565, 23.550, 0.291};
    if(is_equal(type, "HSS24X8X1/4")) return {7.650, 23.650, 0.233};
    if(is_equal(type, "HSS22X22X1")) return {20.600, 20.600, 0.930};
    if(is_equal(type, "HSS22X22X7/8")) return {20.800, 20.800, 0.814};
    if(is_equal(type, "HSS22X22X3/4")) return {20.950, 20.950, 0.698};
    if(is_equal(type, "HSS22X22X5/8")) return {21.150, 21.150, 0.581};
    if(is_equal(type, "HSS22X22X1/2")) return {21.300, 21.300, 0.465};
    if(is_equal(type, "HSS22X20X3/4")) return {18.950, 20.950, 0.698};
    if(is_equal(type, "HSS22X20X5/8")) return {19.150, 21.150, 0.581};
    if(is_equal(type, "HSS22X20X1/2")) return {19.300, 21.300, 0.465};
    if(is_equal(type, "HSS22X20X3/8")) return {19.500, 21.500, 0.349};
    if(is_equal(type, "HSS22X20X5/16")) return {19.550, 21.550, 0.291};
    if(is_equal(type, "HSS22X18X3/4")) return {16.950, 20.950, 0.698};
    if(is_equal(type, "HSS22X18X5/8")) return {17.150, 21.150, 0.581};
    if(is_equal(type, "HSS22X18X1/2")) return {17.300, 21.300, 0.465};
    if(is_equal(type, "HSS22X18X3/8")) return {17.500, 21.500, 0.349};
    if(is_equal(type, "HSS22X18X5/16")) return {17.550, 21.550, 0.291};
    if(is_equal(type, "HSS22X16X3/4")) return {14.950, 20.950, 0.698};
    if(is_equal(type, "HSS22X16X5/8")) return {15.150, 21.150, 0.581};
    if(is_equal(type, "HSS22X16X1/2")) return {15.300, 21.300, 0.465};
    if(is_equal(type, "HSS22X16X3/8")) return {15.500, 21.500, 0.349};
    if(is_equal(type, "HSS22X16X5/16")) return {15.550, 21.550, 0.291};
    if(is_equal(type, "HSS22X16X1/4")) return {15.650, 21.650, 0.233};
    if(is_equal(type, "HSS22X14X3/4")) return {12.950, 20.950, 0.698};
    if(is_equal(type, "HSS22X14X5/8")) return {13.150, 21.150, 0.581};
    if(is_equal(type, "HSS22X14X1/2")) return {13.300, 21.300, 0.465};
    if(is_equal(type, "HSS22X14X3/8")) return {13.500, 21.500, 0.349};
    if(is_equal(type, "HSS22X14X5/16")) return {13.550, 21.550, 0.291};
    if(is_equal(type, "HSS22X14X1/4")) return {13.650, 21.650, 0.233};
    if(is_equal(type, "HSS22X10X5/8")) return {9.130, 21.150, 0.581};
    if(is_equal(type, "HSS22X10X1/2")) return {9.305, 21.300, 0.465};
    if(is_equal(type, "HSS22X10X3/8")) return {9.475, 21.500, 0.349};
    if(is_equal(type, "HSS22X10X5/16")) return {9.565, 21.550, 0.291};
    if(is_equal(type, "HSS22X10X1/4")) return {9.650, 21.650, 0.233};
    if(is_equal(type, "HSS20X20X1")) return {18.600, 18.600, 0.930};
    if(is_equal(type, "HSS20X20X7/8")) return {18.800, 18.800, 0.814};
    if(is_equal(type, "HSS20X20X3/4")) return {18.950, 18.950, 0.698};
    if(is_equal(type, "HSS20X20X5/8")) return {19.150, 19.150, 0.581};
    if(is_equal(type, "HSS20X20X1/2")) return {19.300, 19.300, 0.465};
    if(is_equal(type, "HSS20X20X3/8")) return {19.500, 19.500, 0.349};
    if(is_equal(type, "HSS20X20X5/16")) return {19.550, 19.550, 0.291};
    if(is_equal(type, "HSS20X16X3/4")) return {14.950, 18.950, 0.698};
    if(is_equal(type, "HSS20X16X5/8")) return {15.150, 19.150, 0.581};
    if(is_equal(type, "HSS20X16X1/2")) return {15.300, 19.300, 0.465};
    if(is_equal(type, "HSS20X16X3/8")) return {15.500, 19.500, 0.349};
    if(is_equal(type, "HSS20X16X5/16")) return {15.550, 19.550, 0.291};
    if(is_equal(type, "HSS20X16X1/4")) return {15.650, 19.650, 0.233};
    if(is_equal(type, "HSS20X12X1")) return {10.605, 18.600, 0.930};
    if(is_equal(type, "HSS20X12X7/8")) return {10.780, 18.800, 0.814};
    if(is_equal(type, "HSS20X12X3/4")) return {10.955, 18.950, 0.698};
    if(is_equal(type, "HSS20X12X5/8")) return {11.150, 19.150, 0.581};
    if(is_equal(type, "HSS20X12X1/2")) return {11.300, 19.300, 0.465};
    if(is_equal(type, "HSS20X12X3/8")) return {11.500, 19.500, 0.349};
    if(is_equal(type, "HSS20X12X5/16")) return {11.550, 19.550, 0.291};
    if(is_equal(type, "HSS20X8X1")) return {6.605, 18.600, 0.930};
    if(is_equal(type, "HSS20X8X7/8")) return {6.780, 18.800, 0.814};
    if(is_equal(type, "HSS20X8X3/4")) return {6.955, 18.950, 0.698};
    if(is_equal(type, "HSS20X8X5/8")) return {7.130, 19.150, 0.581};
    if(is_equal(type, "HSS20X8X1/2")) return {7.300, 19.300, 0.465};
    if(is_equal(type, "HSS20X8X3/8")) return {7.475, 19.500, 0.349};
    if(is_equal(type, "HSS20X8X5/16")) return {7.565, 19.550, 0.291};
    if(is_equal(type, "HSS20X6X5/8")) return {5.130, 19.150, 0.581};
    if(is_equal(type, "HSS20X6X1/2")) return {5.305, 19.300, 0.465};
    if(is_equal(type, "HSS20X6X3/8")) return {5.475, 19.500, 0.349};
    if(is_equal(type, "HSS20X6X5/16")) return {5.565, 19.550, 0.291};
    if(is_equal(type, "HSS20X6X1/4")) return {5.650, 19.650, 0.233};
    if(is_equal(type, "HSS20X4X1/2")) return {3.300, 19.300, 0.465};
    if(is_equal(type, "HSS20X4X3/8")) return {3.475, 19.500, 0.349};
    if(is_equal(type, "HSS20X4X5/16")) return {3.565, 19.550, 0.291};
    if(is_equal(type, "HSS20X4X1/4")) return {3.650, 19.650, 0.233};
    if(is_equal(type, "HSS18X18X1")) return {16.600, 16.600, 0.930};
    if(is_equal(type, "HSS18X18X7/8")) return {16.800, 16.800, 0.814};
    if(is_equal(type, "HSS18X18X3/4")) return {16.950, 16.950, 0.698};
    if(is_equal(type, "HSS18X18X5/8")) return {17.150, 17.150, 0.581};
    if(is_equal(type, "HSS18X18X1/2")) return {17.300, 17.300, 0.465};
    if(is_equal(type, "HSS18X18X3/8")) return {17.500, 17.500, 0.349};
    if(is_equal(type, "HSS18X18X5/16")) return {17.550, 17.550, 0.291};
    if(is_equal(type, "HSS18X18X1/4")) return {17.650, 17.650, 0.233};
    if(is_equal(type, "HSS18X10X5/8")) return {9.130, 17.150, 0.581};
    if(is_equal(type, "HSS18X10X1/2")) return {9.305, 17.300, 0.465};
    if(is_equal(type, "HSS18X10X3/8")) return {9.475, 17.500, 0.349};
    if(is_equal(type, "HSS18X10X5/16")) return {9.565, 17.550, 0.291};
    if(is_equal(type, "HSS18X10X1/4")) return {9.650, 17.650, 0.233};
    if(is_equal(type, "HSS18X8X5/8")) return {7.130, 17.150, 0.581};
    if(is_equal(type, "HSS18X8X1/2")) return {7.305, 17.300, 0.465};
    if(is_equal(type, "HSS18X8X3/8")) return {7.475, 17.500, 0.349};
    if(is_equal(type, "HSS18X8X5/16")) return {7.565, 17.550, 0.291};
    if(is_equal(type, "HSS18X8X1/4")) return {7.650, 17.650, 0.233};
    if(is_equal(type, "HSS18X6X3/4")) return {4.955, 16.950, 0.698};
    if(is_equal(type, "HSS18X6X5/8")) return {5.130, 17.150, 0.581};
    if(is_equal(type, "HSS18X6X1/2")) return {5.305, 17.300, 0.465};
    if(is_equal(type, "HSS18X6X3/8")) return {5.475, 17.500, 0.349};
    if(is_equal(type, "HSS18X6X5/16")) return {5.565, 17.550, 0.291};
    if(is_equal(type, "HSS18X6X1/4")) return {5.650, 17.650, 0.233};
    if(is_equal(type, "HSS16X16X1")) return {14.600, 14.600, 0.930};
    if(is_equal(type, "HSS16X16X7/8")) return {14.800, 14.800, 0.814};
    if(is_equal(type, "HSS16X16X3/4")) return {14.950, 14.950, 0.698};
    if(is_equal(type, "HSS16X16X5/8")) return {15.150, 15.150, 0.581};
    if(is_equal(type, "HSS16X16X1/2")) return {15.300, 15.300, 0.465};
    if(is_equal(type, "HSS16X16X3/8")) return {15.500, 15.500, 0.349};
    if(is_equal(type, "HSS16X16X5/16")) return {15.550, 15.550, 0.291};
    if(is_equal(type, "HSS16X16X1/4")) return {15.650, 15.650, 0.233};
    if(is_equal(type, "HSS16X12X1")) return {10.605, 14.600, 0.930};
    if(is_equal(type, "HSS16X12X7/8")) return {10.780, 14.800, 0.814};
    if(is_equal(type, "HSS16X12X3/4")) return {10.955, 14.950, 0.698};
    if(is_equal(type, "HSS16X12X5/8")) return {11.150, 15.150, 0.581};
    if(is_equal(type, "HSS16X12X1/2")) return {11.300, 15.300, 0.465};
    if(is_equal(type, "HSS16X12X3/8")) return {11.500, 15.500, 0.349};
    if(is_equal(type, "HSS16X12X5/16")) return {11.550, 15.550, 0.291};
    if(is_equal(type, "HSS16X10X5/8")) return {9.130, 15.150, 0.581};
    if(is_equal(type, "HSS16X10X1/2")) return {9.305, 15.300, 0.465};
    if(is_equal(type, "HSS16X10X3/8")) return {9.475, 15.500, 0.349};
    if(is_equal(type, "HSS16X10X5/16")) return {9.565, 15.550, 0.291};
    if(is_equal(type, "HSS16X10X1/4")) return {9.650, 15.650, 0.233};
    if(is_equal(type, "HSS16X8X7/8")) return {6.780, 14.800, 0.814};
    if(is_equal(type, "HSS16X8X3/4")) return {6.955, 14.950, 0.698};
    if(is_equal(type, "HSS16X8X5/8")) return {7.130, 15.150, 0.581};
    if(is_equal(type, "HSS16X8X1/2")) return {7.300, 15.300, 0.465};
    if(is_equal(type, "HSS16X8X3/8")) return {7.475, 15.500, 0.349};
    if(is_equal(type, "HSS16X8X5/16")) return {7.565, 15.550, 0.291};
    if(is_equal(type, "HSS16X8X1/4")) return {7.650, 15.650, 0.233};
    if(is_equal(type, "HSS16X6X5/8")) return {5.130, 15.150, 0.581};
    if(is_equal(type, "HSS16X6X1/2")) return {5.305, 15.300, 0.465};
    if(is_equal(type, "HSS16X6X3/8")) return {5.475, 15.500, 0.349};
    if(is_equal(type, "HSS16X6X5/16")) return {5.565, 15.550, 0.291};
    if(is_equal(type, "HSS16X6X1/4")) return {5.650, 15.650, 0.233};
    if(is_equal(type, "HSS16X6X3/16")) return {5.740, 15.750, 0.174};
    if(is_equal(type, "HSS16X4X5/8")) return {3.130, 15.150, 0.581};
    if(is_equal(type, "HSS16X4X1/2")) return {3.300, 15.300, 0.465};
    if(is_equal(type, "HSS16X4X3/8")) return {3.475, 15.500, 0.349};
    if(is_equal(type, "HSS16X4X5/16")) return {3.565, 15.550, 0.291};
    if(is_equal(type, "HSS16X4X1/4")) return {3.650, 15.650, 0.233};
    if(is_equal(type, "HSS16X4X3/16")) return {3.740, 15.750, 0.174};
    if(is_equal(type, "HSS14X14X1")) return {12.600, 12.600, 0.930};
    if(is_equal(type, "HSS14X14X7/8")) return {12.800, 12.800, 0.814};
    if(is_equal(type, "HSS14X14X3/4")) return {12.950, 12.950, 0.698};
    if(is_equal(type, "HSS14X14X5/8")) return {13.150, 13.150, 0.581};
    if(is_equal(type, "HSS14X14X1/2")) return {13.300, 13.300, 0.465};
    if(is_equal(type, "HSS14X14X3/8")) return {13.500, 13.500, 0.349};
    if(is_equal(type, "HSS14X14X5/16")) return {13.550, 13.550, 0.291};
    if(is_equal(type, "HSS14X14X1/4")) return {13.650, 13.650, 0.233};
    if(is_equal(type, "HSS14X12X5/8")) return {11.150, 13.150, 0.581};
    if(is_equal(type, "HSS14X12X1/2")) return {11.300, 13.300, 0.465};
    if(is_equal(type, "HSS14X12X3/8")) return {11.500, 13.500, 0.349};
    if(is_equal(type, "HSS14X12X5/16")) return {11.550, 13.550, 0.291};
    if(is_equal(type, "HSS14X12X1/4")) return {11.650, 13.650, 0.233};
    if(is_equal(type, "HSS14X10X7/8")) return {8.780, 12.800, 0.814};
    if(is_equal(type, "HSS14X10X3/4")) return {8.955, 12.950, 0.698};
    if(is_equal(type, "HSS14X10X5/8")) return {9.130, 13.150, 0.581};
    if(is_equal(type, "HSS14X10X1/2")) return {9.300, 13.300, 0.465};
    if(is_equal(type, "HSS14X10X3/8")) return {9.475, 13.500, 0.349};
    if(is_equal(type, "HSS14X10X5/16")) return {9.565, 13.550, 0.291};
    if(is_equal(type, "HSS14X10X1/4")) return {9.650, 13.650, 0.233};
    if(is_equal(type, "HSS14X8X5/8")) return {7.130, 13.150, 0.581};
    if(is_equal(type, "HSS14X8X1/2")) return {7.305, 13.300, 0.465};
    if(is_equal(type, "HSS14X8X3/8")) return {7.475, 13.500, 0.349};
    if(is_equal(type, "HSS14X8X5/16")) return {7.565, 13.550, 0.291};
    if(is_equal(type, "HSS14X8X1/4")) return {7.650, 13.650, 0.233};
    if(is_equal(type, "HSS14X8X3/16")) return {7.740, 13.750, 0.174};
    if(is_equal(type, "HSS14X6X5/8")) return {5.130, 13.150, 0.581};
    if(is_equal(type, "HSS14X6X1/2")) return {5.305, 13.300, 0.465};
    if(is_equal(type, "HSS14X6X3/8")) return {5.475, 13.500, 0.349};
    if(is_equal(type, "HSS14X6X5/16")) return {5.565, 13.550, 0.291};
    if(is_equal(type, "HSS14X6X1/4")) return {5.650, 13.650, 0.233};
    if(is_equal(type, "HSS14X6X3/16")) return {5.740, 13.750, 0.174};
    if(is_equal(type, "HSS14X4X5/8")) return {3.130, 13.150, 0.581};
    if(is_equal(type, "HSS14X4X1/2")) return {3.300, 13.300, 0.465};
    if(is_equal(type, "HSS14X4X3/8")) return {3.475, 13.500, 0.349};
    if(is_equal(type, "HSS14X4X5/16")) return {3.565, 13.550, 0.291};
    if(is_equal(type, "HSS14X4X1/4")) return {3.650, 13.650, 0.233};
    if(is_equal(type, "HSS14X4X3/16")) return {3.740, 13.750, 0.174};
    if(is_equal(type, "HSS12X12X1")) return {10.605, 10.605, 0.930};
    if(is_equal(type, "HSS12X12X7/8")) return {10.780, 10.780, 0.814};
    if(is_equal(type, "HSS12X12X3/4")) return {10.955, 10.955, 0.698};
    if(is_equal(type, "HSS12X12X5/8")) return {11.150, 11.150, 0.581};
    if(is_equal(type, "HSS12X12X1/2")) return {11.300, 11.300, 0.465};
    if(is_equal(type, "HSS12X12X3/8")) return {11.500, 11.500, 0.349};
    if(is_equal(type, "HSS12X12X5/16")) return {11.550, 11.550, 0.291};
    if(is_equal(type, "HSS12X12X1/4")) return {11.650, 11.650, 0.233};
    if(is_equal(type, "HSS12X12X3/16")) return {11.750, 11.750, 0.174};
    if(is_equal(type, "HSS12X10X5/8")) return {9.130, 11.150, 0.581};
    if(is_equal(type, "HSS12X10X1/2")) return {9.300, 11.300, 0.465};
    if(is_equal(type, "HSS12X10X3/8")) return {9.475, 11.500, 0.349};
    if(is_equal(type, "HSS12X10X5/16")) return {9.565, 11.550, 0.291};
    if(is_equal(type, "HSS12X10X1/4")) return {9.650, 11.650, 0.233};
    if(is_equal(type, "HSS12X10X3/16")) return {9.740, 11.750, 0.174};
    if(is_equal(type, "HSS12X8X5/8")) return {7.130, 11.150, 0.581};
    if(is_equal(type, "HSS12X8X1/2")) return {7.300, 11.300, 0.465};
    if(is_equal(type, "HSS12X8X3/8")) return {7.475, 11.500, 0.349};
    if(is_equal(type, "HSS12X8X5/16")) return {7.565, 11.550, 0.291};
    if(is_equal(type, "HSS12X8X1/4")) return {7.650, 11.650, 0.233};
    if(is_equal(type, "HSS12X8X3/16")) return {7.740, 11.750, 0.174};
    if(is_equal(type, "HSS12X6X5/8")) return {5.130, 11.150, 0.581};
    if(is_equal(type, "HSS12X6X1/2")) return {5.305, 11.300, 0.465};
    if(is_equal(type, "HSS12X6X3/8")) return {5.475, 11.500, 0.349};
    if(is_equal(type, "HSS12X6X5/16")) return {5.565, 11.550, 0.291};
    if(is_equal(type, "HSS12X6X1/4")) return {5.650, 11.650, 0.233};
    if(is_equal(type, "HSS12X6X3/16")) return {5.740, 11.750, 0.174};
    if(is_equal(type, "HSS12X4X5/8")) return {3.130, 11.150, 0.581};
    if(is_equal(type, "HSS12X4X1/2")) return {3.300, 11.300, 0.465};
    if(is_equal(type, "HSS12X4X3/8")) return {3.475, 11.500, 0.349};
    if(is_equal(type, "HSS12X4X5/16")) return {3.565, 11.550, 0.291};
    if(is_equal(type, "HSS12X4X1/4")) return {3.650, 11.650, 0.233};
    if(is_equal(type, "HSS12X4X3/16")) return {3.740, 11.750, 0.174};
    if(is_equal(type, "HSS12X3X5/16")) return {2.565, 11.550, 0.291};
    if(is_equal(type, "HSS12X3X1/4")) return {2.650, 11.650, 0.233};
    if(is_equal(type, "HSS12X3X3/16")) return {2.740, 11.750, 0.174};
    if(is_equal(type, "HSS12X2X5/16")) return {1.565, 11.550, 0.291};
    if(is_equal(type, "HSS12X2X1/4")) return {1.650, 11.650, 0.233};
    if(is_equal(type, "HSS12X2X3/16")) return {1.740, 11.750, 0.174};
    if(is_equal(type, "HSS10X10X3/4")) return {8.955, 8.955, 0.698};
    if(is_equal(type, "HSS10X10X5/8")) return {9.130, 9.130, 0.581};
    if(is_equal(type, "HSS10X10X1/2")) return {9.300, 9.300, 0.465};
    if(is_equal(type, "HSS10X10X3/8")) return {9.475, 9.475, 0.349};
    if(is_equal(type, "HSS10X10X5/16")) return {9.565, 9.565, 0.291};
    if(is_equal(type, "HSS10X10X1/4")) return {9.650, 9.650, 0.233};
    if(is_equal(type, "HSS10X10X3/16")) return {9.740, 9.740, 0.174};
    if(is_equal(type, "HSS10X8X5/8")) return {7.130, 9.130, 0.581};
    if(is_equal(type, "HSS10X8X1/2")) return {7.300, 9.300, 0.465};
    if(is_equal(type, "HSS10X8X3/8")) return {7.475, 9.475, 0.349};
    if(is_equal(type, "HSS10X8X5/16")) return {7.565, 9.565, 0.291};
    if(is_equal(type, "HSS10X8X1/4")) return {7.650, 9.650, 0.233};
    if(is_equal(type, "HSS10X8X3/16")) return {7.740, 9.740, 0.174};
    if(is_equal(type, "HSS10X6X5/8")) return {5.130, 9.130, 0.581};
    if(is_equal(type, "HSS10X6X1/2")) return {5.305, 9.300, 0.465};
    if(is_equal(type, "HSS10X6X3/8")) return {5.475, 9.475, 0.349};
    if(is_equal(type, "HSS10X6X5/16")) return {5.565, 9.565, 0.291};
    if(is_equal(type, "HSS10X6X1/4")) return {5.650, 9.650, 0.233};
    if(is_equal(type, "HSS10X6X3/16")) return {5.740, 9.740, 0.174};
    if(is_equal(type, "HSS10X5X3/8")) return {4.475, 9.475, 0.349};
    if(is_equal(type, "HSS10X5X5/16")) return {4.565, 9.565, 0.291};
    if(is_equal(type, "HSS10X5X1/4")) return {4.650, 9.650, 0.233};
    if(is_equal(type, "HSS10X4X5/8")) return {3.130, 9.130, 0.581};
    if(is_equal(type, "HSS10X4X1/2")) return {3.300, 9.300, 0.465};
    if(is_equal(type, "HSS10X4X3/8")) return {3.475, 9.475, 0.349};
    if(is_equal(type, "HSS10X4X5/16")) return {3.565, 9.565, 0.291};
    if(is_equal(type, "HSS10X4X1/4")) return {3.650, 9.650, 0.233};
    if(is_equal(type, "HSS10X4X3/16")) return {3.740, 9.740, 0.174};
    if(is_equal(type, "HSS10X4X1/8")) return {3.825, 9.825, 0.116};
    if(is_equal(type, "HSS10X3-1/2X3/8")) return {2.975, 9.475, 0.349};
    if(is_equal(type, "HSS10X3-1/2X5/16")) return {3.065, 9.565, 0.291};
    if(is_equal(type, "HSS10X3-1/2X1/4")) return {3.150, 9.650, 0.233};
    if(is_equal(type, "HSS10X3-1/2X3/16")) return {3.240, 9.740, 0.174};
    if(is_equal(type, "HSS10X3X3/8")) return {2.475, 9.475, 0.349};
    if(is_equal(type, "HSS10X3X5/16")) return {2.565, 9.565, 0.291};
    if(is_equal(type, "HSS10X3X1/4")) return {2.650, 9.650, 0.233};
    if(is_equal(type, "HSS10X3X3/16")) return {2.740, 9.740, 0.174};
    if(is_equal(type, "HSS10X3X1/8")) return {2.825, 9.825, 0.116};
    if(is_equal(type, "HSS10X2X3/8")) return {1.476, 9.475, 0.349};
    if(is_equal(type, "HSS10X2X5/16")) return {1.565, 9.565, 0.291};
    if(is_equal(type, "HSS10X2X1/4")) return {1.650, 9.650, 0.233};
    if(is_equal(type, "HSS10X2X3/16")) return {1.740, 9.740, 0.174};
    if(is_equal(type, "HSS10X2X1/8")) return {1.825, 9.825, 0.116};
    if(is_equal(type, "HSS9X9X5/8")) return {8.130, 8.130, 0.581};
    if(is_equal(type, "HSS9X9X1/2")) return {8.300, 8.300, 0.465};
    if(is_equal(type, "HSS9X9X3/8")) return {8.475, 8.475, 0.349};
    if(is_equal(type, "HSS9X9X5/16")) return {8.565, 8.565, 0.291};
    if(is_equal(type, "HSS9X9X1/4")) return {8.650, 8.650, 0.233};
    if(is_equal(type, "HSS9X9X3/16")) return {8.740, 8.740, 0.174};
    if(is_equal(type, "HSS9X9X1/8")) return {8.825, 8.825, 0.116};
    if(is_equal(type, "HSS9X7X5/8")) return {6.130, 8.130, 0.581};
    if(is_equal(type, "HSS9X7X1/2")) return {6.300, 8.300, 0.465};
    if(is_equal(type, "HSS9X7X3/8")) return {6.475, 8.475, 0.349};
    if(is_equal(type, "HSS9X7X5/16")) return {6.565, 8.565, 0.291};
    if(is_equal(type, "HSS9X7X1/4")) return {6.650, 8.650, 0.233};
    if(is_equal(type, "HSS9X7X3/16")) return {6.740, 8.740, 0.174};
    if(is_equal(type, "HSS9X5X5/8")) return {4.130, 8.130, 0.581};
    if(is_equal(type, "HSS9X5X1/2")) return {4.300, 8.300, 0.465};
    if(is_equal(type, "HSS9X5X3/8")) return {4.475, 8.475, 0.349};
    if(is_equal(type, "HSS9X5X5/16")) return {4.565, 8.565, 0.291};
    if(is_equal(type, "HSS9X5X1/4")) return {4.650, 8.650, 0.233};
    if(is_equal(type, "HSS9X5X3/16")) return {4.740, 8.740, 0.174};
    if(is_equal(type, "HSS9X3X1/2")) return {2.300, 8.300, 0.465};
    if(is_equal(type, "HSS9X3X3/8")) return {2.475, 8.475, 0.349};
    if(is_equal(type, "HSS9X3X5/16")) return {2.565, 8.565, 0.291};
    if(is_equal(type, "HSS9X3X1/4")) return {2.650, 8.650, 0.233};
    if(is_equal(type, "HSS9X3X3/16")) return {2.740, 8.740, 0.174};
    if(is_equal(type, "HSS8X8X5/8")) return {7.130, 7.130, 0.581};
    if(is_equal(type, "HSS8X8X1/2")) return {7.300, 7.300, 0.465};
    if(is_equal(type, "HSS8X8X3/8")) return {7.475, 7.475, 0.349};
    if(is_equal(type, "HSS8X8X5/16")) return {7.565, 7.565, 0.291};
    if(is_equal(type, "HSS8X8X1/4")) return {7.650, 7.650, 0.233};
    if(is_equal(type, "HSS8X8X3/16")) return {7.740, 7.740, 0.174};
    if(is_equal(type, "HSS8X8X1/8")) return {7.825, 7.825, 0.116};
    if(is_equal(type, "HSS8X6X5/8")) return {5.130, 7.130, 0.581};
    if(is_equal(type, "HSS8X6X1/2")) return {5.305, 7.300, 0.465};
    if(is_equal(type, "HSS8X6X3/8")) return {5.475, 7.475, 0.349};
    if(is_equal(type, "HSS8X6X5/16")) return {5.565, 7.565, 0.291};
    if(is_equal(type, "HSS8X6X1/4")) return {5.650, 7.650, 0.233};
    if(is_equal(type, "HSS8X6X3/16")) return {5.740, 7.740, 0.174};
    if(is_equal(type, "HSS8X4X5/8")) return {3.130, 7.130, 0.581};
    if(is_equal(type, "HSS8X4X1/2")) return {3.300, 7.300, 0.465};
    if(is_equal(type, "HSS8X4X3/8")) return {3.475, 7.475, 0.349};
    if(is_equal(type, "HSS8X4X5/16")) return {3.565, 7.565, 0.291};
    if(is_equal(type, "HSS8X4X1/4")) return {3.650, 7.650, 0.233};
    if(is_equal(type, "HSS8X4X3/16")) return {3.740, 7.740, 0.174};
    if(is_equal(type, "HSS8X4X1/8")) return {3.825, 7.825, 0.116};
    if(is_equal(type, "HSS8X3X1/2")) return {2.300, 7.300, 0.465};
    if(is_equal(type, "HSS8X3X3/8")) return {2.475, 7.475, 0.349};
    if(is_equal(type, "HSS8X3X5/16")) return {2.565, 7.565, 0.291};
    if(is_equal(type, "HSS8X3X1/4")) return {2.650, 7.650, 0.233};
    if(is_equal(type, "HSS8X3X3/16")) return {2.740, 7.740, 0.174};
    if(is_equal(type, "HSS8X3X1/8")) return {2.825, 7.825, 0.116};
    if(is_equal(type, "HSS8X2X1/2")) return {1.302, 7.305, 0.465};
    if(is_equal(type, "HSS8X2X3/8")) return {1.476, 7.475, 0.349};
    if(is_equal(type, "HSS8X2X5/16")) return {1.565, 7.565, 0.291};
    if(is_equal(type, "HSS8X2X1/4")) return {1.650, 7.650, 0.233};
    if(is_equal(type, "HSS8X2X3/16")) return {1.740, 7.740, 0.174};
    if(is_equal(type, "HSS8X2X1/8")) return {1.825, 7.825, 0.116};
    if(is_equal(type, "HSS7X7X5/8")) return {6.130, 6.130, 0.581};
    if(is_equal(type, "HSS7X7X1/2")) return {6.300, 6.300, 0.465};
    if(is_equal(type, "HSS7X7X3/8")) return {6.475, 6.475, 0.349};
    if(is_equal(type, "HSS7X7X5/16")) return {6.565, 6.565, 0.291};
    if(is_equal(type, "HSS7X7X1/4")) return {6.650, 6.650, 0.233};
    if(is_equal(type, "HSS7X7X3/16")) return {6.740, 6.740, 0.174};
    if(is_equal(type, "HSS7X7X1/8")) return {6.825, 6.825, 0.116};
    if(is_equal(type, "HSS7X5X1/2")) return {4.300, 6.300, 0.465};
    if(is_equal(type, "HSS7X5X3/8")) return {4.475, 6.475, 0.349};
    if(is_equal(type, "HSS7X5X5/16")) return {4.565, 6.565, 0.291};
    if(is_equal(type, "HSS7X5X1/4")) return {4.650, 6.650, 0.233};
    if(is_equal(type, "HSS7X5X3/16")) return {4.740, 6.740, 0.174};
    if(is_equal(type, "HSS7X5X1/8")) return {4.825, 6.825, 0.116};
    if(is_equal(type, "HSS7X4X1/2")) return {3.300, 6.300, 0.465};
    if(is_equal(type, "HSS7X4X3/8")) return {3.475, 6.475, 0.349};
    if(is_equal(type, "HSS7X4X5/16")) return {3.565, 6.565, 0.291};
    if(is_equal(type, "HSS7X4X1/4")) return {3.650, 6.650, 0.233};
    if(is_equal(type, "HSS7X4X3/16")) return {3.740, 6.740, 0.174};
    if(is_equal(type, "HSS7X4X1/8")) return {3.825, 6.825, 0.116};
    if(is_equal(type, "HSS7X3X1/2")) return {2.300, 6.300, 0.465};
    if(is_equal(type, "HSS7X3X3/8")) return {2.475, 6.475, 0.349};
    if(is_equal(type, "HSS7X3X5/16")) return {2.565, 6.565, 0.291};
    if(is_equal(type, "HSS7X3X1/4")) return {2.650, 6.650, 0.233};
    if(is_equal(type, "HSS7X3X3/16")) return {2.740, 6.740, 0.174};
    if(is_equal(type, "HSS7X3X1/8")) return {2.825, 6.825, 0.116};
    if(is_equal(type, "HSS7X2X1/4")) return {1.650, 6.650, 0.233};
    if(is_equal(type, "HSS7X2X3/16")) return {1.740, 6.740, 0.174};
    if(is_equal(type, "HSS7X2X1/8")) return {1.825, 6.825, 0.116};
    if(is_equal(type, "HSS6X6X5/8")) return {5.130, 5.130, 0.581};
    if(is_equal(type, "HSS6X6X1/2")) return {5.305, 5.305, 0.465};
    if(is_equal(type, "HSS6X6X3/8")) return {5.475, 5.475, 0.349};
    if(is_equal(type, "HSS6X6X5/16")) return {5.565, 5.565, 0.291};
    if(is_equal(type, "HSS6X6X1/4")) return {5.650, 5.650, 0.233};
    if(is_equal(type, "HSS6X6X3/16")) return {5.740, 5.740, 0.174};
    if(is_equal(type, "HSS6X6X1/8")) return {5.825, 5.825, 0.116};
    if(is_equal(type, "HSS6X5X1/2")) return {4.300, 5.305, 0.465};
    if(is_equal(type, "HSS6X5X3/8")) return {4.475, 5.475, 0.349};
    if(is_equal(type, "HSS6X5X5/16")) return {4.565, 5.565, 0.291};
    if(is_equal(type, "HSS6X5X1/4")) return {4.650, 5.650, 0.233};
    if(is_equal(type, "HSS6X5X3/16")) return {4.740, 5.740, 0.174};
    if(is_equal(type, "HSS6X5X1/8")) return {4.825, 5.825, 0.116};
    if(is_equal(type, "HSS6X4X1/2")) return {3.300, 5.305, 0.465};
    if(is_equal(type, "HSS6X4X3/8")) return {3.475, 5.475, 0.349};
    if(is_equal(type, "HSS6X4X5/16")) return {3.565, 5.565, 0.291};
    if(is_equal(type, "HSS6X4X1/4")) return {3.650, 5.650, 0.233};
    if(is_equal(type, "HSS6X4X3/16")) return {3.740, 5.740, 0.174};
    if(is_equal(type, "HSS6X4X1/8")) return {3.825, 5.825, 0.116};
    if(is_equal(type, "HSS6X3X1/2")) return {2.300, 5.305, 0.465};
    if(is_equal(type, "HSS6X3X3/8")) return {2.475, 5.475, 0.349};
    if(is_equal(type, "HSS6X3X5/16")) return {2.565, 5.565, 0.291};
    if(is_equal(type, "HSS6X3X1/4")) return {2.650, 5.650, 0.233};
    if(is_equal(type, "HSS6X3X3/16")) return {2.740, 5.740, 0.174};
    if(is_equal(type, "HSS6X3X1/8")) return {2.825, 5.825, 0.116};
    if(is_equal(type, "HSS6X2X3/8")) return {1.476, 5.475, 0.349};
    if(is_equal(type, "HSS6X2X5/16")) return {1.565, 5.565, 0.291};
    if(is_equal(type, "HSS6X2X1/4")) return {1.650, 5.650, 0.233};
    if(is_equal(type, "HSS6X2X3/16")) return {1.740, 5.740, 0.174};
    if(is_equal(type, "HSS6X2X1/8")) return {1.825, 5.825, 0.116};
    if(is_equal(type, "HSS5-1/2X5-1/2X3/8")) return {4.975, 4.975, 0.349};
    if(is_equal(type, "HSS5-1/2X5-1/2X5/16")) return {5.065, 5.065, 0.291};
    if(is_equal(type, "HSS5-1/2X5-1/2X1/4")) return {5.150, 5.150, 0.233};
    if(is_equal(type, "HSS5-1/2X5-1/2X3/16")) return {5.240, 5.240, 0.174};
    if(is_equal(type, "HSS5-1/2X5-1/2X1/8")) return {5.325, 5.325, 0.116};
    if(is_equal(type, "HSS5X5X1/2")) return {4.300, 4.300, 0.465};
    if(is_equal(type, "HSS5X5X3/8")) return {4.475, 4.475, 0.349};
    if(is_equal(type, "HSS5X5X5/16")) return {4.565, 4.565, 0.291};
    if(is_equal(type, "HSS5X5X1/4")) return {4.650, 4.650, 0.233};
    if(is_equal(type, "HSS5X5X3/16")) return {4.740, 4.740, 0.174};
    if(is_equal(type, "HSS5X5X1/8")) return {4.825, 4.825, 0.116};
    if(is_equal(type, "HSS5X4X1/2")) return {3.300, 4.300, 0.465};
    if(is_equal(type, "HSS5X4X3/8")) return {3.475, 4.475, 0.349};
    if(is_equal(type, "HSS5X4X5/16")) return {3.565, 4.565, 0.291};
    if(is_equal(type, "HSS5X4X1/4")) return {3.650, 4.650, 0.233};
    if(is_equal(type, "HSS5X4X3/16")) return {3.740, 4.740, 0.174};
    if(is_equal(type, "HSS5X4X1/8")) return {3.825, 4.825, 0.116};
    if(is_equal(type, "HSS5X3X1/2")) return {2.300, 4.300, 0.465};
    if(is_equal(type, "HSS5X3X3/8")) return {2.475, 4.475, 0.349};
    if(is_equal(type, "HSS5X3X5/16")) return {2.565, 4.565, 0.291};
    if(is_equal(type, "HSS5X3X1/4")) return {2.650, 4.650, 0.233};
    if(is_equal(type, "HSS5X3X3/16")) return {2.740, 4.740, 0.174};
    if(is_equal(type, "HSS5X3X1/8")) return {2.825, 4.825, 0.116};
    if(is_equal(type, "HSS5X2-1/2X1/4")) return {2.150, 4.650, 0.233};
    if(is_equal(type, "HSS5X2-1/2X3/16")) return {2.240, 4.740, 0.174};
    if(is_equal(type, "HSS5X2-1/2X1/8")) return {2.325, 4.825, 0.116};
    if(is_equal(type, "HSS5X2X3/8")) return {1.476, 4.475, 0.349};
    if(is_equal(type, "HSS5X2X5/16")) return {1.565, 4.565, 0.291};
    if(is_equal(type, "HSS5X2X1/4")) return {1.650, 4.650, 0.233};
    if(is_equal(type, "HSS5X2X3/16")) return {1.740, 4.740, 0.174};
    if(is_equal(type, "HSS5X2X1/8")) return {1.825, 4.825, 0.116};
    if(is_equal(type, "HSS4-1/2X4-1/2X1/2")) return {3.800, 3.800, 0.465};
    if(is_equal(type, "HSS4-1/2X4-1/2X3/8")) return {3.975, 3.975, 0.349};
    if(is_equal(type, "HSS4-1/2X4-1/2X5/16")) return {4.065, 4.065, 0.291};
    if(is_equal(type, "HSS4-1/2X4-1/2X1/4")) return {4.150, 4.150, 0.233};
    if(is_equal(type, "HSS4-1/2X4-1/2X3/16")) return {4.240, 4.240, 0.174};
    if(is_equal(type, "HSS4-1/2X4-1/2X1/8")) return {4.325, 4.325, 0.116};
    if(is_equal(type, "HSS4X4X1/2")) return {3.300, 3.300, 0.465};
    if(is_equal(type, "HSS4X4X3/8")) return {3.475, 3.475, 0.349};
    if(is_equal(type, "HSS4X4X5/16")) return {3.565, 3.565, 0.291};
    if(is_equal(type, "HSS4X4X1/4")) return {3.650, 3.650, 0.233};
    if(is_equal(type, "HSS4X4X3/16")) return {3.740, 3.740, 0.174};
    if(is_equal(type, "HSS4X4X1/8")) return {3.825, 3.825, 0.116};
    if(is_equal(type, "HSS4X3X3/8")) return {2.475, 3.475, 0.349};
    if(is_equal(type, "HSS4X3X5/16")) return {2.565, 3.565, 0.291};
    if(is_equal(type, "HSS4X3X1/4")) return {2.650, 3.650, 0.233};
    if(is_equal(type, "HSS4X3X3/16")) return {2.740, 3.740, 0.174};
    if(is_equal(type, "HSS4X3X1/8")) return {2.825, 3.825, 0.116};
    if(is_equal(type, "HSS4X2-1/2X1/4")) return {2.150, 3.650, 0.233};
    if(is_equal(type, "HSS4X2-1/2X3/16")) return {2.240, 3.740, 0.174};
    if(is_equal(type, "HSS4X2-1/2X1/8")) return {2.325, 3.825, 0.116};
    if(is_equal(type, "HSS4X2X3/8")) return {1.476, 3.475, 0.349};
    if(is_equal(type, "HSS4X2X5/16")) return {1.565, 3.565, 0.291};
    if(is_equal(type, "HSS4X2X1/4")) return {1.650, 3.650, 0.233};
    if(is_equal(type, "HSS4X2X3/16")) return {1.740, 3.740, 0.174};
    if(is_equal(type, "HSS4X2X1/8")) return {1.825, 3.825, 0.116};
    if(is_equal(type, "HSS4X1-1/2X1/4")) return {1.151, 3.650, 0.233};
    if(is_equal(type, "HSS4X1-1/2X3/16")) return {1.239, 3.740, 0.174};
    if(is_equal(type, "HSS4X1-1/2X1/8")) return {1.325, 3.825, 0.116};
    if(is_equal(type, "HSS3-1/2X3-1/2X3/8")) return {2.975, 2.975, 0.349};
    if(is_equal(type, "HSS3-1/2X3-1/2X5/16")) return {3.065, 3.065, 0.291};
    if(is_equal(type, "HSS3-1/2X3-1/2X1/4")) return {3.150, 3.150, 0.233};
    if(is_equal(type, "HSS3-1/2X3-1/2X3/16")) return {3.240, 3.240, 0.174};
    if(is_equal(type, "HSS3-1/2X3-1/2X1/8")) return {3.325, 3.325, 0.116};
    if(is_equal(type, "HSS3-1/2X2-1/2X3/8")) return {1.975, 2.975, 0.349};
    if(is_equal(type, "HSS3-1/2X2-1/2X5/16")) return {2.065, 3.065, 0.291};
    if(is_equal(type, "HSS3-1/2X2-1/2X1/4")) return {2.150, 3.150, 0.233};
    if(is_equal(type, "HSS3-1/2X2-1/2X3/16")) return {2.240, 3.240, 0.174};
    if(is_equal(type, "HSS3-1/2X2-1/2X1/8")) return {2.325, 3.325, 0.116};
    if(is_equal(type, "HSS3-1/2X2X1/4")) return {1.650, 3.150, 0.233};
    if(is_equal(type, "HSS3-1/2X2X3/16")) return {1.740, 3.240, 0.174};
    if(is_equal(type, "HSS3-1/2X2X1/8")) return {1.825, 3.325, 0.116};
    if(is_equal(type, "HSS3-1/2X1-1/2X1/4")) return {1.151, 3.150, 0.233};
    if(is_equal(type, "HSS3-1/2X1-1/2X3/16")) return {1.239, 3.240, 0.174};
    if(is_equal(type, "HSS3-1/2X1-1/2X1/8")) return {1.325, 3.325, 0.116};
    if(is_equal(type, "HSS3X3X3/8")) return {2.475, 2.475, 0.349};
    if(is_equal(type, "HSS3X3X5/16")) return {2.565, 2.565, 0.291};
    if(is_equal(type, "HSS3X3X1/4")) return {2.650, 2.650, 0.233};
    if(is_equal(type, "HSS3X3X3/16")) return {2.740, 2.740, 0.174};
    if(is_equal(type, "HSS3X3X1/8")) return {2.825, 2.825, 0.116};
    if(is_equal(type, "HSS3X2-1/2X5/16")) return {2.065, 2.565, 0.291};
    if(is_equal(type, "HSS3X2-1/2X1/4")) return {2.150, 2.650, 0.233};
    if(is_equal(type, "HSS3X2-1/2X3/16")) return {2.240, 2.740, 0.174};
    if(is_equal(type, "HSS3X2-1/2X1/8")) return {2.325, 2.825, 0.116};
    if(is_equal(type, "HSS3X2X5/16")) return {1.565, 2.565, 0.291};
    if(is_equal(type, "HSS3X2X1/4")) return {1.650, 2.650, 0.233};
    if(is_equal(type, "HSS3X2X3/16")) return {1.740, 2.740, 0.174};
    if(is_equal(type, "HSS3X2X1/8")) return {1.825, 2.825, 0.116};
    if(is_equal(type, "HSS3X1-1/2X1/4")) return {1.151, 2.650, 0.233};
    if(is_equal(type, "HSS3X1-1/2X3/16")) return {1.239, 2.740, 0.174};
    if(is_equal(type, "HSS3X1-1/2X1/8")) return {1.325, 2.825, 0.116};
    if(is_equal(type, "HSS3X1X3/16")) return {0.739, 2.740, 0.174};
    if(is_equal(type, "HSS3X1X1/8")) return {0.826, 2.825, 0.116};
    if(is_equal(type, "HSS2-1/2X2-1/2X5/16")) return {2.065, 2.065, 0.291};
    if(is_equal(type, "HSS2-1/2X2-1/2X1/4")) return {2.150, 2.150, 0.233};
    if(is_equal(type, "HSS2-1/2X2-1/2X3/16")) return {2.240, 2.240, 0.174};
    if(is_equal(type, "HSS2-1/2X2-1/2X1/8")) return {2.325, 2.325, 0.116};
    if(is_equal(type, "HSS2-1/2X2X1/4")) return {1.650, 2.150, 0.233};
    if(is_equal(type, "HSS2-1/2X2X3/16")) return {1.740, 2.240, 0.174};
    if(is_equal(type, "HSS2-1/2X2X1/8")) return {1.825, 2.325, 0.116};
    if(is_equal(type, "HSS2-1/2X1-1/2X1/4")) return {1.151, 2.150, 0.233};
    if(is_equal(type, "HSS2-1/2X1-1/2X3/16")) return {1.239, 2.240, 0.174};
    if(is_equal(type, "HSS2-1/2X1-1/2X1/8")) return {1.325, 2.325, 0.116};
    if(is_equal(type, "HSS2-1/2X1X3/16")) return {0.739, 2.240, 0.174};
    if(is_equal(type, "HSS2-1/2X1X1/8")) return {0.826, 2.325, 0.116};
    if(is_equal(type, "HSS2-1/4X2-1/4X1/4")) return {1.900, 1.900, 0.233};
    if(is_equal(type, "HSS2-1/4X2-1/4X3/16")) return {1.990, 1.990, 0.174};
    if(is_equal(type, "HSS2-1/4X2-1/4X1/8")) return {2.075, 2.075, 0.116};
    if(is_equal(type, "HSS2X2X1/4")) return {1.650, 1.650, 0.233};
    if(is_equal(type, "HSS2X2X3/16")) return {1.740, 1.740, 0.174};
    if(is_equal(type, "HSS2X2X1/8")) return {1.825, 1.825, 0.116};
    if(is_equal(type, "HSS2X1-1/2X3/16")) return {1.239, 1.740, 0.174};
    if(is_equal(type, "HSS2X1-1/2X1/8")) return {1.325, 1.825, 0.116};
    if(is_equal(type, "HSS2X1X3/16")) return {0.739, 1.740, 0.174};
    if(is_equal(type, "HSS2X1X1/8")) return {0.826, 1.825, 0.116};
    if(is_equal(type, "HSS1-1/2X1-1/2X1/4")) return {1.151, 1.151, 0.233};
    if(is_equal(type, "HSS1-1/2X1-1/2X3/16")) return {1.239, 1.239, 0.174};
    if(is_equal(type, "HSS1-1/2X1-1/2X1/8")) return {1.325, 1.325, 0.116};
    return {};
}

vec uschssection(const std::string_view type) {
    if(is_equal(type, "HSS28.000X1.000")) return {14.000, 0.930};
    if(is_equal(type, "HSS28.000X0.875")) return {14.000, 0.814};
    if(is_equal(type, "HSS28.000X0.750")) return {14.000, 0.698};
    if(is_equal(type, "HSS28.000X0.625")) return {14.000, 0.581};
    if(is_equal(type, "HSS28.000X0.500")) return {14.000, 0.465};
    if(is_equal(type, "HSS28.000X0.375")) return {14.000, 0.349};
    if(is_equal(type, "HSS26.000X0.750")) return {13.000, 0.698};
    if(is_equal(type, "HSS26.000X0.625")) return {13.000, 0.581};
    if(is_equal(type, "HSS26.000X0.500")) return {13.000, 0.465};
    if(is_equal(type, "HSS26.000X0.375")) return {13.000, 0.349};
    if(is_equal(type, "HSS26.000X0.313")) return {13.000, 0.291};
    if(is_equal(type, "HSS24.000X1.000")) return {12.000, 0.930};
    if(is_equal(type, "HSS24.000X0.875")) return {12.000, 0.814};
    if(is_equal(type, "HSS24.000X0.750")) return {12.000, 0.698};
    if(is_equal(type, "HSS24.000X0.625")) return {12.000, 0.581};
    if(is_equal(type, "HSS24.000X0.500")) return {12.000, 0.465};
    if(is_equal(type, "HSS24.000X0.375")) return {12.000, 0.349};
    if(is_equal(type, "HSS24.000X0.313")) return {12.000, 0.291};
    if(is_equal(type, "HSS22.000X0.750")) return {11.000, 0.698};
    if(is_equal(type, "HSS22.000X0.625")) return {11.000, 0.581};
    if(is_equal(type, "HSS22.000X0.500")) return {11.000, 0.465};
    if(is_equal(type, "HSS22.000X0.375")) return {11.000, 0.349};
    if(is_equal(type, "HSS22.000X0.313")) return {11.000, 0.291};
    if(is_equal(type, "HSS20.000X1.000")) return {10.000, 0.930};
    if(is_equal(type, "HSS20.000X0.875")) return {10.000, 0.814};
    if(is_equal(type, "HSS20.000X0.750")) return {10.000, 0.698};
    if(is_equal(type, "HSS20.000X0.625")) return {10.000, 0.581};
    if(is_equal(type, "HSS20.000X0.500")) return {10.000, 0.465};
    if(is_equal(type, "HSS20.000X0.375")) return {10.000, 0.349};
    if(is_equal(type, "HSS20.000X0.313")) return {10.000, 0.291};
    if(is_equal(type, "HSS20.000X0.250")) return {10.000, 0.233};
    if(is_equal(type, "HSS18.000X1.000")) return {9.000, 0.930};
    if(is_equal(type, "HSS18.000X0.875")) return {9.000, 0.814};
    if(is_equal(type, "HSS18.000X0.750")) return {9.000, 0.698};
    if(is_equal(type, "HSS18.000X0.625")) return {9.000, 0.581};
    if(is_equal(type, "HSS18.000X0.500")) return {9.000, 0.465};
    if(is_equal(type, "HSS18.000X0.375")) return {9.000, 0.349};
    if(is_equal(type, "HSS18.000X0.313")) return {9.000, 0.291};
    if(is_equal(type, "HSS18.000X0.250")) return {9.000, 0.233};
    if(is_equal(type, "HSS16.000X1.000")) return {8.000, 0.930};
    if(is_equal(type, "HSS16.000X0.875")) return {8.000, 0.814};
    if(is_equal(type, "HSS16.000X0.750")) return {8.000, 0.698};
    if(is_equal(type, "HSS16.000X0.625")) return {8.000, 0.581};
    if(is_equal(type, "HSS16.000X0.500")) return {8.000, 0.465};
    if(is_equal(type, "HSS16.000X0.438")) return {8.000, 0.407};
    if(is_equal(type, "HSS16.000X0.375")) return {8.000, 0.349};
    if(is_equal(type, "HSS16.000X0.312")) return {8.000, 0.291};
    if(is_equal(type, "HSS16.000X0.250")) return {8.000, 0.233};
    if(is_equal(type, "HSS14.000X1.000")) return {7.000, 0.930};
    if(is_equal(type, "HSS14.000X0.875")) return {7.000, 0.814};
    if(is_equal(type, "HSS14.000X0.750")) return {7.000, 0.698};
    if(is_equal(type, "HSS14.000X0.625")) return {7.000, 0.581};
    if(is_equal(type, "HSS14.000X0.500")) return {7.000, 0.465};
    if(is_equal(type, "HSS14.000X0.375")) return {7.000, 0.349};
    if(is_equal(type, "HSS14.000X0.312")) return {7.000, 0.291};
    if(is_equal(type, "HSS14.000X0.250")) return {7.000, 0.233};
    if(is_equal(type, "HSS14.000X0.188")) return {7.000, 0.174};
    if(is_equal(type, "HSS13.375X0.625")) return {6.700, 0.581};
    if(is_equal(type, "HSS13.375X0.500")) return {6.700, 0.465};
    if(is_equal(type, "HSS13.375X0.375")) return {6.700, 0.349};
    if(is_equal(type, "HSS13.375X0.313")) return {6.700, 0.291};
    if(is_equal(type, "HSS13.375X0.250")) return {6.700, 0.233};
    if(is_equal(type, "HSS13.375X0.188")) return {6.700, 0.174};
    if(is_equal(type, "HSS12.750X0.750")) return {6.400, 0.698};
    if(is_equal(type, "HSS12.750X0.625")) return {6.400, 0.581};
    if(is_equal(type, "HSS12.750X0.500")) return {6.400, 0.465};
    if(is_equal(type, "HSS12.750X0.375")) return {6.400, 0.349};
    if(is_equal(type, "HSS12.750X0.250")) return {6.400, 0.233};
    if(is_equal(type, "HSS12.750X0.188")) return {6.400, 0.174};
    if(is_equal(type, "HSS12.000X0.625")) return {6.000, 0.581};
    if(is_equal(type, "HSS12.000X0.500")) return {6.000, 0.465};
    if(is_equal(type, "HSS12.000X0.375")) return {6.000, 0.349};
    if(is_equal(type, "HSS12.000X0.250")) return {6.000, 0.233};
    if(is_equal(type, "HSS11.750X0.625")) return {5.900, 0.581};
    if(is_equal(type, "HSS11.750X0.500")) return {5.900, 0.465};
    if(is_equal(type, "HSS11.750X0.375")) return {5.900, 0.349};
    if(is_equal(type, "HSS11.750X0.337")) return {5.900, 0.313};
    if(is_equal(type, "HSS11.750X0.250")) return {5.900, 0.233};
    if(is_equal(type, "HSS10.750X0.625")) return {5.400, 0.581};
    if(is_equal(type, "HSS10.750X0.500")) return {5.400, 0.465};
    if(is_equal(type, "HSS10.750X0.375")) return {5.400, 0.349};
    if(is_equal(type, "HSS10.750X0.313")) return {5.400, 0.291};
    if(is_equal(type, "HSS10.750X0.250")) return {5.400, 0.233};
    if(is_equal(type, "HSS10.750X0.188")) return {5.400, 0.174};
    if(is_equal(type, "HSS10.000X0.625")) return {5.000, 0.581};
    if(is_equal(type, "HSS10.000X0.500")) return {5.000, 0.465};
    if(is_equal(type, "HSS10.000X0.375")) return {5.000, 0.349};
    if(is_equal(type, "HSS10.000X0.312")) return {5.000, 0.291};
    if(is_equal(type, "HSS10.000X0.250")) return {5.000, 0.233};
    if(is_equal(type, "HSS10.000X0.188")) return {5.000, 0.174};
    if(is_equal(type, "HSS9.625X0.625")) return {4.815, 0.581};
    if(is_equal(type, "HSS9.625X0.500")) return {4.815, 0.465};
    if(is_equal(type, "HSS9.625X0.375")) return {4.815, 0.349};
    if(is_equal(type, "HSS9.625X0.312")) return {4.815, 0.291};
    if(is_equal(type, "HSS9.625X0.250")) return {4.815, 0.233};
    if(is_equal(type, "HSS9.625X0.188")) return {4.815, 0.174};
    if(is_equal(type, "HSS8.625X0.625")) return {4.315, 0.581};
    if(is_equal(type, "HSS8.625X0.500")) return {4.315, 0.465};
    if(is_equal(type, "HSS8.625X0.375")) return {4.315, 0.349};
    if(is_equal(type, "HSS8.625X0.322")) return {4.315, 0.300};
    if(is_equal(type, "HSS8.625X0.250")) return {4.315, 0.233};
    if(is_equal(type, "HSS8.625X0.188")) return {4.315, 0.174};
    if(is_equal(type, "HSS7.500X0.500")) return {3.750, 0.465};
    if(is_equal(type, "HSS7.500X0.375")) return {3.750, 0.349};
    if(is_equal(type, "HSS7.500X0.312")) return {3.750, 0.291};
    if(is_equal(type, "HSS7.500X0.250")) return {3.750, 0.233};
    if(is_equal(type, "HSS7.500X0.188")) return {3.750, 0.174};
    if(is_equal(type, "HSS7.000X0.500")) return {3.500, 0.465};
    if(is_equal(type, "HSS7.000X0.375")) return {3.500, 0.349};
    if(is_equal(type, "HSS7.000X0.312")) return {3.500, 0.291};
    if(is_equal(type, "HSS7.000X0.250")) return {3.500, 0.233};
    if(is_equal(type, "HSS7.000X0.188")) return {3.500, 0.174};
    if(is_equal(type, "HSS7.000X0.125")) return {3.500, 0.116};
    if(is_equal(type, "HSS6.875X0.375")) return {3.440, 0.349};
    if(is_equal(type, "HSS6.875X0.312")) return {3.440, 0.291};
    if(is_equal(type, "HSS6.875X0.250")) return {3.440, 0.233};
    if(is_equal(type, "HSS6.875X0.188")) return {3.440, 0.174};
    if(is_equal(type, "HSS6.625X0.500")) return {3.315, 0.465};
    if(is_equal(type, "HSS6.625X0.432")) return {3.315, 0.402};
    if(is_equal(type, "HSS6.625X0.375")) return {3.315, 0.349};
    if(is_equal(type, "HSS6.625X0.312")) return {3.315, 0.291};
    if(is_equal(type, "HSS6.625X0.280")) return {3.315, 0.260};
    if(is_equal(type, "HSS6.625X0.250")) return {3.315, 0.233};
    if(is_equal(type, "HSS6.625X0.188")) return {3.315, 0.174};
    if(is_equal(type, "HSS6.625X0.125")) return {3.315, 0.116};
    if(is_equal(type, "HSS6.000X0.500")) return {3.000, 0.465};
    if(is_equal(type, "HSS6.000X0.375")) return {3.000, 0.349};
    if(is_equal(type, "HSS6.000X0.312")) return {3.000, 0.291};
    if(is_equal(type, "HSS6.000X0.280")) return {3.000, 0.260};
    if(is_equal(type, "HSS6.000X0.250")) return {3.000, 0.233};
    if(is_equal(type, "HSS6.000X0.188")) return {3.000, 0.174};
    if(is_equal(type, "HSS6.000X0.125")) return {3.000, 0.116};
    if(is_equal(type, "HSS5.563X0.500")) return {2.780, 0.465};
    if(is_equal(type, "HSS5.563X0.375")) return {2.780, 0.349};
    if(is_equal(type, "HSS5.563X0.258")) return {2.780, 0.240};
    if(is_equal(type, "HSS5.563X0.188")) return {2.780, 0.174};
    if(is_equal(type, "HSS5.563X0.134")) return {2.780, 0.124};
    if(is_equal(type, "HSS5.500X0.500")) return {2.750, 0.465};
    if(is_equal(type, "HSS5.500X0.375")) return {2.750, 0.349};
    if(is_equal(type, "HSS5.500X0.258")) return {2.750, 0.240};
    if(is_equal(type, "HSS5.000X0.500")) return {2.500, 0.465};
    if(is_equal(type, "HSS5.000X0.375")) return {2.500, 0.349};
    if(is_equal(type, "HSS5.000X0.312")) return {2.500, 0.291};
    if(is_equal(type, "HSS5.000X0.258")) return {2.500, 0.240};
    if(is_equal(type, "HSS5.000X0.250")) return {2.500, 0.233};
    if(is_equal(type, "HSS5.000X0.188")) return {2.500, 0.174};
    if(is_equal(type, "HSS5.000X0.125")) return {2.500, 0.116};
    if(is_equal(type, "HSS4.500X0.375")) return {2.250, 0.349};
    if(is_equal(type, "HSS4.500X0.337")) return {2.250, 0.313};
    if(is_equal(type, "HSS4.500X0.237")) return {2.250, 0.220};
    if(is_equal(type, "HSS4.500X0.188")) return {2.250, 0.174};
    if(is_equal(type, "HSS4.500X0.125")) return {2.250, 0.116};
    if(is_equal(type, "HSS4.000X0.313")) return {2.000, 0.291};
    if(is_equal(type, "HSS4.000X0.250")) return {2.000, 0.233};
    if(is_equal(type, "HSS4.000X0.237")) return {2.000, 0.220};
    if(is_equal(type, "HSS4.000X0.226")) return {2.000, 0.210};
    if(is_equal(type, "HSS4.000X0.220")) return {2.000, 0.205};
    if(is_equal(type, "HSS4.000X0.188")) return {2.000, 0.174};
    if(is_equal(type, "HSS4.000X0.125")) return {2.000, 0.116};
    if(is_equal(type, "HSS3.500X0.313")) return {1.750, 0.291};
    if(is_equal(type, "HSS3.500X0.300")) return {1.750, 0.279};
    if(is_equal(type, "HSS3.500X0.250")) return {1.750, 0.233};
    if(is_equal(type, "HSS3.500X0.216")) return {1.750, 0.201};
    if(is_equal(type, "HSS3.500X0.203")) return {1.750, 0.189};
    if(is_equal(type, "HSS3.500X0.188")) return {1.750, 0.174};
    if(is_equal(type, "HSS3.500X0.125")) return {1.750, 0.116};
    if(is_equal(type, "HSS3.000X0.250")) return {1.500, 0.233};
    if(is_equal(type, "HSS3.000X0.216")) return {1.500, 0.201};
    if(is_equal(type, "HSS3.000X0.203")) return {1.500, 0.189};
    if(is_equal(type, "HSS3.000X0.188")) return {1.500, 0.174};
    if(is_equal(type, "HSS3.000X0.152")) return {1.500, 0.141};
    if(is_equal(type, "HSS3.000X0.134")) return {1.500, 0.124};
    if(is_equal(type, "HSS3.000X0.125")) return {1.500, 0.116};
    if(is_equal(type, "HSS2.875X0.250")) return {1.440, 0.233};
    if(is_equal(type, "HSS2.875X0.203")) return {1.440, 0.189};
    if(is_equal(type, "HSS2.875X0.188")) return {1.440, 0.174};
    if(is_equal(type, "HSS2.875X0.125")) return {1.440, 0.116};
    if(is_equal(type, "HSS2.500X0.250")) return {1.250, 0.233};
    if(is_equal(type, "HSS2.500X0.188")) return {1.250, 0.174};
    if(is_equal(type, "HSS2.500X0.125")) return {1.250, 0.116};
    if(is_equal(type, "HSS2.375X0.250")) return {1.190, 0.233};
    if(is_equal(type, "HSS2.375X0.218")) return {1.190, 0.203};
    if(is_equal(type, "HSS2.375X0.188")) return {1.190, 0.174};
    if(is_equal(type, "HSS2.375X0.154")) return {1.190, 0.143};
    if(is_equal(type, "HSS2.375X0.125")) return {1.190, 0.116};
    if(is_equal(type, "HSS1.900X0.188")) return {0.950, 0.174};
    if(is_equal(type, "HSS1.900X0.145")) return {0.950, 0.135};
    if(is_equal(type, "HSS1.900X0.120")) return {0.950, 0.111};
    if(is_equal(type, "HSS1.660X0.140")) return {0.830, 0.130};
    return {};
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

    return_obj = make_unique<ISection2D>(tag, scale * dim, material_id, int_pt, eccentricity);
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

    return_obj = make_unique<ISection3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
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
        return_obj = make_unique<ISection2D>(tag, scale * dim, material_id, int_pt, eccentricity);
        return;
    }

    dim = nzchsection(type);

    if(!dim.is_empty()) {
        return_obj = make_unique<CircularHollow2D>(tag, scale * dim, material_id, int_pt, eccentricity);
        return;
    }

    dim = nzrhsection(type);

    if(!dim.is_empty()) {
        dim[0] -= dim[2]; // account for centroid
        dim[1] -= dim[2]; // account for centroid
        return_obj = make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
        return;
    }

    dim = nzshsection(type);

    if(!dim.is_empty()) {
        dim[0] -= dim[2]; // account for centroid
        dim[1] -= dim[2]; // account for centroid
        return_obj = make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
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
        return_obj = make_unique<ISection3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    dim = nzchsection(type);

    if(!dim.is_empty()) {
        return_obj = make_unique<CircularHollow3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    dim = nzrhsection(type);

    if(!dim.is_empty()) {
        dim[0] -= dim[2]; // account for centroid
        dim[1] -= dim[2]; // account for centroid
        return_obj = make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    dim = nzshsection(type);

    if(!dim.is_empty()) {
        dim[0] -= dim[2]; // account for centroid
        dim[1] -= dim[2]; // account for centroid
        return_obj = make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    suanpan_error("Cannot identify section type.\n");
}

void new_us2d(unique_ptr<Section>& return_obj, std::istringstream& command, const bool recenter) {
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
        return_obj = make_unique<ISection2D>(tag, scale * dim, material_id, int_pt, eccentricity);
        return;
    }

    dim = usrhssection(type);

    if(!dim.is_empty()) {
        return_obj = make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
        return;
    }

    dim = uschssection(type);

    if(!dim.is_empty()) {
        return_obj = make_unique<CircularHollow2D>(tag, scale * dim, material_id, int_pt, eccentricity);
        return;
    }

    dim = ustsection(type);

    if(!dim.is_empty()) {
        if(recenter) eccentricity = barycenter(dim *= scale);
        return_obj = make_unique<TSection2D>(tag, std::move(dim), material_id, int_pt, eccentricity);
        return;
    }

    suanpan_error("Cannot identify section type.\n");
}

void new_us3d(unique_ptr<Section>& return_obj, std::istringstream& command, const bool recenter) {
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
        return_obj = make_unique<ISection3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    dim = usrhssection(type);

    if(!dim.is_empty()) {
        return_obj = make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    dim = uschssection(type);

    if(!dim.is_empty()) {
        return_obj = make_unique<CircularHollow3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    dim = ustsection(type);

    if(!dim.is_empty()) {
        if(recenter) eccentricity_y = barycenter(dim *= scale);
        return_obj = make_unique<TSection3D>(tag, std::move(dim), material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    suanpan_error("Cannot identify section type.\n");
}

int create_new_section(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section type is required.\n");
        return 0;
    }

    unique_ptr<Section> new_section = nullptr;

    if(is_equal(section_id, "Box2D")) new_box2d(new_section, command);
    else if(is_equal(section_id, "Box3D")) new_box3d(new_section, command);
    else if(is_equal(section_id, "Cell2D") || is_equal(section_id, "Bar2D")) new_cell2d(new_section, command);
    else if(is_equal(section_id, "Cell3D") || is_equal(section_id, "Bar3D")) new_cell3d(new_section, command);
    else if(is_equal(section_id, "Cell3DOS")) new_cell3dos(new_section, command);
    else if(is_equal(section_id, "Circle1D")) new_circle1d(new_section, command);
    else if(is_equal(section_id, "Circle2D")) new_circle2d(new_section, command);
    else if(is_equal(section_id, "Circle3D")) new_circle3d(new_section, command);
    else if(is_equal(section_id, "CircularHollow2D")) new_circularhollow2D(new_section, command);
    else if(is_equal(section_id, "CircularHollow3D")) new_circularhollow3D(new_section, command);
    else if(is_equal(section_id, "Fibre1D")) new_fibre1d(new_section, command);
    else if(is_equal(section_id, "Fibre2D")) new_fibre2d(new_section, command);
    else if(is_equal(section_id, "Fibre3D")) new_fibre3d(new_section, command, false);
    else if(is_equal(section_id, "Fibre3DOS")) new_fibre3d(new_section, command, true);
    else if(is_equal(section_id, "HSection2D")) new_hsection2d(new_section, command);
    else if(is_equal(section_id, "ISection2D")) new_isection2d(new_section, command, false);
    else if(is_equal(section_id, "ISection3D")) new_isection3d(new_section, command, false);
    else if(is_equal(section_id, "ISection2DC")) new_isection2d(new_section, command, true);
    else if(is_equal(section_id, "ISection3DC")) new_isection3d(new_section, command, true);
    else if(is_equal(section_id, "Rectangle1D")) new_rectangle1d(new_section, command);
    else if(is_equal(section_id, "Rectangle2D")) new_rectangle2d(new_section, command);
    else if(is_equal(section_id, "Rectangle3D")) new_rectangle3d(new_section, command);
    else if(is_equal(section_id, "TrussSection")) new_trusssection(new_section, command);
    else if(is_equal(section_id, "TSection2D")) new_tsection2d(new_section, command, false);
    else if(is_equal(section_id, "TSection3D")) new_tsection3d(new_section, command, false);
    else if(is_equal(section_id, "TSection2DC")) new_tsection2d(new_section, command, true);
    else if(is_equal(section_id, "TSection3DC")) new_tsection3d(new_section, command, true);
    else if(is_equal(section_id, "NM2D1")) new_nm2d(new_section, command, 3);
    else if(is_equal(section_id, "NM2D2")) new_nm2d(new_section, command, 8);
    else if(is_equal(section_id, "NM2D3")) new_nm2d(new_section, command, 11);
    else if(is_equal(section_id, "NM3D1")) new_nm3d(new_section, command, 4);
    else if(is_equal(section_id, "NM3D2")) new_nm3d(new_section, command, 10);
    else if(is_equal(section_id, "NM3D3")) new_nm3d(new_section, command, 13);
    else if(is_equal(section_id, "NM2D3K")) new_nmk(new_section, command, 13);
    else if(is_equal(section_id, "NM3D3K")) new_nmk(new_section, command, 17);
    else if(is_equal(section_id, "EU2D")) new_eu2d(new_section, command);
    else if(is_equal(section_id, "EU3D")) new_eu3d(new_section, command);
    else if(is_equal(section_id, "NZ2D")) new_nz2d(new_section, command);
    else if(is_equal(section_id, "NZ3D")) new_nz3d(new_section, command);
    else if(is_equal(section_id, "US2D")) new_us2d(new_section, command, false);
    else if(is_equal(section_id, "US3D")) new_us3d(new_section, command, false);
    else if(is_equal(section_id, "US2DC")) new_us2d(new_section, command, true);
    else if(is_equal(section_id, "US3DC")) new_us3d(new_section, command, true);
    else load::object(new_section, domain, section_id, command);

    if(new_section == nullptr || !domain->insert(std::move(new_section)))
        suanpan_error("Fail to create new section via \"{}\".\n", command.str());

    return 0;
}
