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

#include "SectionParser.h"
#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Section/Section>
#include <Toolbox/utility.h>

void new_bar2d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    return_obj = make_unique<Bar2D>(tag, area, material_id, eccentricity);
}

void new_bar3d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    return_obj = make_unique<Bar3D>(tag, area, material_id, eccentricity_a, eccentricity_b);
}

void new_box2d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 6;
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

void new_box3d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 3;
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

void new_circle1d(unique_ptr<Section>& return_obj, istringstream& command) {
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

void new_circle2d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 6;
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

void new_circle3d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 6;
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

void new_circularhollow2D(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 10;
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

void new_circularhollow3D(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 10;
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

void new_fibre1d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<uword> tag_vector;
    while(!command.eof())
        if(uword section_tag; get_input(command, section_tag)) tag_vector.emplace_back(section_tag);
        else {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

    return_obj = make_unique<Fibre1D>(tag, std::move(tag_vector));
}

void new_fibre2d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<uword> tag_vector;
    while(!command.eof())
        if(uword section_tag; get_input(command, section_tag)) tag_vector.emplace_back(section_tag);
        else {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

    return_obj = make_unique<Fibre2D>(tag, std::move(tag_vector));
}

void new_fibre3d(unique_ptr<Section>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    vector<uword> tag_vector;
    while(!command.eof())
        if(uword section_tag; get_input(command, section_tag)) tag_vector.emplace_back(section_tag);
        else {
            suanpan_error("A valid parameter is required.\n");
            return;
        }

    return_obj = make_unique<Fibre3D>(tag, std::move(tag_vector));
}

void new_hsection2d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 6;
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
        return .5 * flange_area * (dim(1) + dim(2)) / (flange_area + dim(2) * dim(3));
    }

    // dim(0): top flange width
    // dim(1): top flange thickness
    // dim(2): bottom flange width
    // dim(3): bottom flange thickness
    // dim(4): web height
    // dim(5): web thickness
    const auto top_flange_area = dim(0) * dim(1);
    const auto bottom_flange_area = dim(2) * dim(3);
    return .5 * (top_flange_area * (dim(1) + dim(4)) + bottom_flange_area * (dim(3) + dim(4))) / (top_flange_area + bottom_flange_area + dim(4) * dim(5));
}

void new_isection2d(unique_ptr<Section>& return_obj, istringstream& command, const bool recenter) {
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

    unsigned int_pt = 6;
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

void new_isection3d(unique_ptr<Section>& return_obj, istringstream& command, const bool recenter) {
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

    unsigned int_pt = 6;
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

void new_rectangle1d(unique_ptr<Section>& return_obj, istringstream& command) {
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

void new_rectangle2d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 6;
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

void new_rectangle3d(unique_ptr<Section>& return_obj, istringstream& command) {
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

    unsigned int_pt = 3;
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

void new_trusssection(unique_ptr<Section>& return_obj, istringstream& command) {
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

void new_tsection2d(unique_ptr<Section>& return_obj, istringstream& command, const bool recenter) {
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

    unsigned int_pt = 4;
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

void new_tsection3d(unique_ptr<Section>& return_obj, istringstream& command, const bool recenter) {
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

    unsigned int_pt = 3;
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

void new_nm2d(unique_ptr<Section>& return_obj, istringstream& command, const unsigned size) {
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

    if(3 == size) {
        return_obj = make_unique<NM2D1>(tag, P(0), P(1), P(2));
        return;
    }

    vector<double> para_set;
    double para;
    while(!command.eof() && get_input(command, para)) para_set.emplace_back(para);

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

void new_nm3d(unique_ptr<Section>& return_obj, istringstream& command, const unsigned size) {
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

    if(4 == size) {
        return_obj = make_unique<NM3D1>(tag, P(0), P(1), P(2), P(3));
        return;
    }

    vector<double> para_set;
    double para;
    while(!command.eof() && get_input(command, para)) para_set.emplace_back(para);

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

void new_nmk(unique_ptr<Section>& return_obj, istringstream& command, const unsigned size) {
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

    vector<double> para_set;
    double para;
    while(!command.eof() && get_input(command, para)) para_set.emplace_back(para);

    const auto p_size = 13 == size ? 3 : 4;

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

vec euisection(const string& type) {
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

vec nzchsection(const string& type) {
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

vec nzisection(const string& type) {
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

vec nzrhsection(const string& type) {
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

vec nzshsection(const string& type) {
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

vec usisection(const string& type) {
    if(is_equal(type, "W44X335")) return {15.9, 1.77, 15.9, 1.77, 40.46, 1.03};
    if(is_equal(type, "W44X290")) return {15.8, 1.58, 15.8, 1.58, 40.44, 0.865};
    if(is_equal(type, "W44X262")) return {15.8, 1.42, 15.8, 1.42, 40.46, 0.785};
    if(is_equal(type, "W44X230")) return {15.8, 1.22, 15.8, 1.22, 40.46, 0.71};
    if(is_equal(type, "W40X655")) return {16.9, 3.54, 16.9, 3.54, 36.52, 1.97};
    if(is_equal(type, "W40X593")) return {16.7, 3.23, 16.7, 3.23, 36.54, 1.79};
    if(is_equal(type, "W40X503")) return {16.4, 2.76, 16.4, 2.76, 36.58, 1.54};
    if(is_equal(type, "W40X431")) return {16.2, 2.36, 16.2, 2.36, 36.58, 1.34};
    if(is_equal(type, "W40X397")) return {16.1, 2.2, 16.1, 2.2, 36.6, 1.22};
    if(is_equal(type, "W40X372")) return {16.1, 2.05, 16.1, 2.05, 36.5, 1.16};
    if(is_equal(type, "W40X362")) return {16, 2.01, 16, 2.01, 36.58, 1.12};
    if(is_equal(type, "W40X324")) return {15.9, 1.81, 15.9, 1.81, 36.58, 1};
    if(is_equal(type, "W40X297")) return {15.8, 1.65, 15.8, 1.65, 36.5, 0.93};
    if(is_equal(type, "W40X277")) return {15.8, 1.58, 15.8, 1.58, 36.54, 0.83};
    if(is_equal(type, "W40X249")) return {15.8, 1.42, 15.8, 1.42, 36.56, 0.75};
    if(is_equal(type, "W40X215")) return {15.8, 1.22, 15.8, 1.22, 36.56, 0.65};
    if(is_equal(type, "W40X199")) return {15.8, 1.07, 15.8, 1.07, 36.56, 0.65};
    if(is_equal(type, "W40X392")) return {12.4, 2.52, 12.4, 2.52, 36.56, 1.42};
    if(is_equal(type, "W40X331")) return {12.2, 2.13, 12.2, 2.13, 36.54, 1.22};
    if(is_equal(type, "W40X327")) return {12.1, 2.13, 12.1, 2.13, 36.54, 1.18};
    if(is_equal(type, "W40X294")) return {12, 1.93, 12, 1.93, 36.54, 1.06};
    if(is_equal(type, "W40X278")) return {12, 1.81, 12, 1.81, 36.58, 1.03};
    if(is_equal(type, "W40X264")) return {11.9, 1.73, 11.9, 1.73, 36.54, 0.96};
    if(is_equal(type, "W40X235")) return {11.9, 1.58, 11.9, 1.58, 36.54, 0.83};
    if(is_equal(type, "W40X211")) return {11.8, 1.42, 11.8, 1.42, 36.56, 0.75};
    if(is_equal(type, "W40X183")) return {11.8, 1.2, 11.8, 1.2, 36.6, 0.65};
    if(is_equal(type, "W40X167")) return {11.8, 1.03, 11.8, 1.03, 36.54, 0.65};
    if(is_equal(type, "W40X149")) return {11.8, 0.83, 11.8, 0.83, 36.54, 0.63};
    if(is_equal(type, "W36X925")) return {18.6, 4.53, 18.6, 4.53, 34.04, 3.02};
    if(is_equal(type, "W36X853")) return {18.2, 4.53, 18.2, 4.53, 34.04, 2.52};
    if(is_equal(type, "W36X802")) return {18, 4.29, 18, 4.29, 34.02, 2.38};
    if(is_equal(type, "W36X723")) return {17.8, 3.9, 17.8, 3.9, 34, 2.17};
    if(is_equal(type, "W36X652")) return {17.6, 3.54, 17.6, 3.54, 34.02, 1.97};
    if(is_equal(type, "W36X529")) return {17.2, 2.91, 17.2, 2.91, 33.98, 1.61};
    if(is_equal(type, "W36X487")) return {17.1, 2.68, 17.1, 2.68, 33.94, 1.5};
    if(is_equal(type, "W36X441")) return {17, 2.44, 17, 2.44, 34.02, 1.36};
    if(is_equal(type, "W36X395")) return {16.8, 2.2, 16.8, 2.2, 34, 1.22};
    if(is_equal(type, "W36X361")) return {16.7, 2.01, 16.7, 2.01, 33.98, 1.12};
    if(is_equal(type, "W36X330")) return {16.6, 1.85, 16.6, 1.85, 34, 1.02};
    if(is_equal(type, "W36X302")) return {16.7, 1.68, 16.7, 1.68, 33.94, 0.945};
    if(is_equal(type, "W36X282")) return {16.6, 1.57, 16.6, 1.57, 33.96, 0.885};
    if(is_equal(type, "W36X262")) return {16.6, 1.44, 16.6, 1.44, 34.02, 0.84};
    if(is_equal(type, "W36X247")) return {16.5, 1.35, 16.5, 1.35, 34, 0.8};
    if(is_equal(type, "W36X231")) return {16.5, 1.26, 16.5, 1.26, 33.98, 0.76};
    if(is_equal(type, "W36X256")) return {12.2, 1.73, 12.2, 1.73, 33.94, 0.96};
    if(is_equal(type, "W36X232")) return {12.1, 1.57, 12.1, 1.57, 33.96, 0.87};
    if(is_equal(type, "W36X210")) return {12.2, 1.36, 12.2, 1.36, 33.98, 0.83};
    if(is_equal(type, "W36X194")) return {12.1, 1.26, 12.1, 1.26, 33.98, 0.765};
    if(is_equal(type, "W36X182")) return {12.1, 1.18, 12.1, 1.18, 33.94, 0.725};
    if(is_equal(type, "W36X170")) return {12, 1.1, 12, 1.1, 34, 0.68};
    if(is_equal(type, "W36X160")) return {12, 1.02, 12, 1.02, 33.96, 0.65};
    if(is_equal(type, "W36X150")) return {12, 0.94, 12, 0.94, 34.02, 0.625};
    if(is_equal(type, "W36X135")) return {12, 0.79, 12, 0.79, 34.02, 0.6};
    if(is_equal(type, "W33X387")) return {16.2, 2.28, 16.2, 2.28, 31.44, 1.26};
    if(is_equal(type, "W33X354")) return {16.1, 2.09, 16.1, 2.09, 31.42, 1.16};
    if(is_equal(type, "W33X318")) return {16, 1.89, 16, 1.89, 31.42, 1.04};
    if(is_equal(type, "W33X291")) return {15.9, 1.73, 15.9, 1.73, 31.34, 0.96};
    if(is_equal(type, "W33X263")) return {15.8, 1.57, 15.8, 1.57, 31.36, 0.87};
    if(is_equal(type, "W33X241")) return {15.9, 1.4, 15.9, 1.4, 31.4, 0.83};
    if(is_equal(type, "W33X221")) return {15.8, 1.28, 15.8, 1.28, 31.34, 0.775};
    if(is_equal(type, "W33X201")) return {15.7, 1.15, 15.7, 1.15, 31.4, 0.715};
    if(is_equal(type, "W33X169")) return {11.5, 1.22, 11.5, 1.22, 31.36, 0.67};
    if(is_equal(type, "W33X152")) return {11.6, 1.06, 11.6, 1.06, 31.38, 0.635};
    if(is_equal(type, "W33X141")) return {11.5, 0.96, 11.5, 0.96, 31.38, 0.605};
    if(is_equal(type, "W33X130")) return {11.5, 0.855, 11.5, 0.855, 31.39, 0.58};
    if(is_equal(type, "W33X118")) return {11.5, 0.74, 11.5, 0.74, 31.42, 0.55};
    if(is_equal(type, "W30X391")) return {15.6, 2.44, 15.6, 2.44, 28.32, 1.36};
    if(is_equal(type, "W30X357")) return {15.5, 2.24, 15.5, 2.24, 28.32, 1.24};
    if(is_equal(type, "W30X326")) return {15.4, 2.05, 15.4, 2.05, 28.3, 1.14};
    if(is_equal(type, "W30X292")) return {15.3, 1.85, 15.3, 1.85, 28.3, 1.02};
    if(is_equal(type, "W30X261")) return {15.2, 1.65, 15.2, 1.65, 28.3, 0.93};
    if(is_equal(type, "W30X235")) return {15.1, 1.5, 15.1, 1.5, 28.3, 0.83};
    if(is_equal(type, "W30X211")) return {15.1, 1.32, 15.1, 1.32, 28.26, 0.775};
    if(is_equal(type, "W30X191")) return {15, 1.19, 15, 1.19, 28.32, 0.71};
    if(is_equal(type, "W30X173")) return {15, 1.07, 15, 1.07, 28.26, 0.655};
    if(is_equal(type, "W30X148")) return {10.5, 1.18, 10.5, 1.18, 28.34, 0.65};
    if(is_equal(type, "W30X132")) return {10.5, 1, 10.5, 1, 28.3, 0.615};
    if(is_equal(type, "W30X124")) return {10.5, 0.93, 10.5, 0.93, 28.34, 0.585};
    if(is_equal(type, "W30X116")) return {10.5, 0.85, 10.5, 0.85, 28.3, 0.565};
    if(is_equal(type, "W30X108")) return {10.5, 0.76, 10.5, 0.76, 28.28, 0.545};
    if(is_equal(type, "W30X99")) return {10.5, 0.67, 10.5, 0.67, 28.36, 0.52};
    if(is_equal(type, "W30X90")) return {10.4, 0.61, 10.4, 0.61, 28.28, 0.47};
    if(is_equal(type, "W27X539")) return {15.3, 3.54, 15.3, 3.54, 25.42, 1.97};
    if(is_equal(type, "W27X368")) return {14.7, 2.48, 14.7, 2.48, 25.44, 1.38};
    if(is_equal(type, "W27X336")) return {14.6, 2.28, 14.6, 2.28, 25.44, 1.26};
    if(is_equal(type, "W27X307")) return {14.4, 2.09, 14.4, 2.09, 25.42, 1.16};
    if(is_equal(type, "W27X281")) return {14.4, 1.93, 14.4, 1.93, 25.44, 1.06};
    if(is_equal(type, "W27X258")) return {14.3, 1.77, 14.3, 1.77, 25.46, 0.98};
    if(is_equal(type, "W27X235")) return {14.2, 1.61, 14.2, 1.61, 25.48, 0.91};
    if(is_equal(type, "W27X217")) return {14.1, 1.5, 14.1, 1.5, 25.4, 0.83};
    if(is_equal(type, "W27X194")) return {14, 1.34, 14, 1.34, 25.42, 0.75};
    if(is_equal(type, "W27X178")) return {14.1, 1.19, 14.1, 1.19, 25.42, 0.725};
    if(is_equal(type, "W27X161")) return {14, 1.08, 14, 1.08, 25.44, 0.66};
    if(is_equal(type, "W27X146")) return {14, 0.975, 14, 0.975, 25.45, 0.605};
    if(is_equal(type, "W27X129")) return {10, 1.1, 10, 1.1, 25.4, 0.61};
    if(is_equal(type, "W27X114")) return {10.1, 0.93, 10.1, 0.93, 25.44, 0.57};
    if(is_equal(type, "W27X102")) return {10, 0.83, 10, 0.83, 25.44, 0.515};
    if(is_equal(type, "W27X94")) return {10, 0.745, 10, 0.745, 25.41, 0.49};
    if(is_equal(type, "W27X84")) return {10, 0.64, 10, 0.64, 25.42, 0.46};
    if(is_equal(type, "W24X370")) return {13.7, 2.72, 13.7, 2.72, 22.56, 1.52};
    if(is_equal(type, "W24X335")) return {13.5, 2.48, 13.5, 2.48, 22.54, 1.38};
    if(is_equal(type, "W24X306")) return {13.4, 2.28, 13.4, 2.28, 22.54, 1.26};
    if(is_equal(type, "W24X279")) return {13.3, 2.09, 13.3, 2.09, 22.52, 1.16};
    if(is_equal(type, "W24X250")) return {13.2, 1.89, 13.2, 1.89, 22.52, 1.04};
    if(is_equal(type, "W24X229")) return {13.1, 1.73, 13.1, 1.73, 22.54, 0.96};
    if(is_equal(type, "W24X207")) return {13, 1.57, 13, 1.57, 22.56, 0.87};
    if(is_equal(type, "W24X192")) return {13, 1.46, 13, 1.46, 22.58, 0.81};
    if(is_equal(type, "W24X176")) return {12.9, 1.34, 12.9, 1.34, 22.52, 0.75};
    if(is_equal(type, "W24X162")) return {13, 1.22, 13, 1.22, 22.56, 0.705};
    if(is_equal(type, "W24X146")) return {12.9, 1.09, 12.9, 1.09, 22.52, 0.65};
    if(is_equal(type, "W24X131")) return {12.9, 0.96, 12.9, 0.96, 22.58, 0.605};
    if(is_equal(type, "W24X117")) return {12.8, 0.85, 12.8, 0.85, 22.6, 0.55};
    if(is_equal(type, "W24X104")) return {12.8, 0.75, 12.8, 0.75, 22.6, 0.5};
    if(is_equal(type, "W24X103")) return {9, 0.98, 9, 0.98, 22.54, 0.55};
    if(is_equal(type, "W24X94")) return {9.07, 0.875, 9.07, 0.875, 22.55, 0.515};
    if(is_equal(type, "W24X84")) return {9.02, 0.77, 9.02, 0.77, 22.56, 0.47};
    if(is_equal(type, "W24X76")) return {8.99, 0.68, 8.99, 0.68, 22.54, 0.44};
    if(is_equal(type, "W24X68")) return {8.97, 0.585, 8.97, 0.585, 22.53, 0.415};
    if(is_equal(type, "W24X62")) return {7.04, 0.59, 7.04, 0.59, 22.52, 0.43};
    if(is_equal(type, "W24X55")) return {7.01, 0.505, 7.01, 0.505, 22.59, 0.395};
    if(is_equal(type, "W21X275")) return {12.9, 2.19, 12.9, 2.19, 19.72, 1.22};
    if(is_equal(type, "W21X248")) return {12.8, 1.99, 12.8, 1.99, 19.72, 1.1};
    if(is_equal(type, "W21X223")) return {12.7, 1.79, 12.7, 1.79, 19.82, 1};
    if(is_equal(type, "W21X201")) return {12.6, 1.63, 12.6, 1.63, 19.74, 0.91};
    if(is_equal(type, "W21X182")) return {12.5, 1.48, 12.5, 1.48, 19.74, 0.83};
    if(is_equal(type, "W21X166")) return {12.4, 1.36, 12.4, 1.36, 19.78, 0.75};
    if(is_equal(type, "W21X147")) return {12.5, 1.15, 12.5, 1.15, 19.8, 0.72};
    if(is_equal(type, "W21X132")) return {12.4, 1.04, 12.4, 1.04, 19.72, 0.65};
    if(is_equal(type, "W21X122")) return {12.4, 0.96, 12.4, 0.96, 19.78, 0.6};
    if(is_equal(type, "W21X111")) return {12.3, 0.875, 12.3, 0.875, 19.75, 0.55};
    if(is_equal(type, "W21X101")) return {12.3, 0.8, 12.3, 0.8, 19.8, 0.5};
    if(is_equal(type, "W21X93")) return {8.42, 0.93, 8.42, 0.93, 19.74, 0.58};
    if(is_equal(type, "W21X83")) return {8.36, 0.835, 8.36, 0.835, 19.73, 0.515};
    if(is_equal(type, "W21X73")) return {8.3, 0.74, 8.3, 0.74, 19.72, 0.455};
    if(is_equal(type, "W21X68")) return {8.27, 0.685, 8.27, 0.685, 19.73, 0.43};
    if(is_equal(type, "W21X62")) return {8.24, 0.615, 8.24, 0.615, 19.77, 0.4};
    if(is_equal(type, "W21X55")) return {8.22, 0.522, 8.22, 0.522, 19.756, 0.375};
    if(is_equal(type, "W21X48")) return {8.14, 0.43, 8.14, 0.43, 19.74, 0.35};
    if(is_equal(type, "W21X57")) return {6.56, 0.65, 6.56, 0.65, 19.8, 0.405};
    if(is_equal(type, "W21X50")) return {6.53, 0.535, 6.53, 0.535, 19.73, 0.38};
    if(is_equal(type, "W21X44")) return {6.5, 0.45, 6.5, 0.45, 19.8, 0.35};
    if(is_equal(type, "W18X311")) return {12, 2.74, 12, 2.74, 16.82, 1.52};
    if(is_equal(type, "W18X283")) return {11.9, 2.5, 11.9, 2.5, 16.9, 1.4};
    if(is_equal(type, "W18X258")) return {11.8, 2.3, 11.8, 2.3, 16.9, 1.28};
    if(is_equal(type, "W18X234")) return {11.7, 2.11, 11.7, 2.11, 16.88, 1.16};
    if(is_equal(type, "W18X211")) return {11.6, 1.91, 11.6, 1.91, 16.88, 1.06};
    if(is_equal(type, "W18X192")) return {11.5, 1.75, 11.5, 1.75, 16.9, 0.96};
    if(is_equal(type, "W18X175")) return {11.4, 1.59, 11.4, 1.59, 16.82, 0.89};
    if(is_equal(type, "W18X158")) return {11.3, 1.44, 11.3, 1.44, 16.82, 0.81};
    if(is_equal(type, "W18X143")) return {11.2, 1.32, 11.2, 1.32, 16.86, 0.73};
    if(is_equal(type, "W18X130")) return {11.2, 1.2, 11.2, 1.2, 16.9, 0.67};
    if(is_equal(type, "W18X119")) return {11.3, 1.06, 11.3, 1.06, 16.88, 0.655};
    if(is_equal(type, "W18X106")) return {11.2, 0.94, 11.2, 0.94, 16.82, 0.59};
    if(is_equal(type, "W18X97")) return {11.1, 0.87, 11.1, 0.87, 16.86, 0.535};
    if(is_equal(type, "W18X86")) return {11.1, 0.77, 11.1, 0.77, 16.86, 0.48};
    if(is_equal(type, "W18X76")) return {11, 0.68, 11, 0.68, 16.84, 0.425};
    if(is_equal(type, "W18X71")) return {7.64, 0.81, 7.64, 0.81, 16.88, 0.495};
    if(is_equal(type, "W18X65")) return {7.59, 0.75, 7.59, 0.75, 16.9, 0.45};
    if(is_equal(type, "W18X60")) return {7.56, 0.695, 7.56, 0.695, 16.81, 0.415};
    if(is_equal(type, "W18X55")) return {7.53, 0.63, 7.53, 0.63, 16.84, 0.39};
    if(is_equal(type, "W18X50")) return {7.5, 0.57, 7.5, 0.57, 16.86, 0.355};
    if(is_equal(type, "W18X46")) return {6.06, 0.605, 6.06, 0.605, 16.89, 0.36};
    if(is_equal(type, "W18X40")) return {6.02, 0.525, 6.02, 0.525, 16.85, 0.315};
    if(is_equal(type, "W18X35")) return {6, 0.425, 6, 0.425, 16.85, 0.3};
    if(is_equal(type, "W16X100")) return {10.4, 0.985, 10.4, 0.985, 15.03, 0.585};
    if(is_equal(type, "W16X89")) return {10.4, 0.875, 10.4, 0.875, 15.05, 0.525};
    if(is_equal(type, "W16X77")) return {10.3, 0.76, 10.3, 0.76, 14.98, 0.455};
    if(is_equal(type, "W16X67")) return {10.2, 0.665, 10.2, 0.665, 14.97, 0.395};
    if(is_equal(type, "W16X57")) return {7.12, 0.715, 7.12, 0.715, 14.97, 0.43};
    if(is_equal(type, "W16X50")) return {7.07, 0.63, 7.07, 0.63, 15.04, 0.38};
    if(is_equal(type, "W16X45")) return {7.04, 0.565, 7.04, 0.565, 14.97, 0.345};
    if(is_equal(type, "W16X40")) return {7, 0.505, 7, 0.505, 14.99, 0.305};
    if(is_equal(type, "W16X36")) return {6.99, 0.43, 6.99, 0.43, 15.04, 0.295};
    if(is_equal(type, "W16X31")) return {5.53, 0.44, 5.53, 0.44, 15.02, 0.275};
    if(is_equal(type, "W16X26")) return {5.5, 0.345, 5.5, 0.345, 15.01, 0.25};
    if(is_equal(type, "W14X873")) return {18.8, 5.51, 18.8, 5.51, 12.58, 3.94};
    if(is_equal(type, "W14X808")) return {18.6, 5.12, 18.6, 5.12, 12.56, 3.74};
    if(is_equal(type, "W14X730")) return {17.9, 4.91, 17.9, 4.91, 12.58, 3.07};
    if(is_equal(type, "W14X665")) return {17.7, 4.52, 17.7, 4.52, 12.56, 2.83};
    if(is_equal(type, "W14X605")) return {17.4, 4.16, 17.4, 4.16, 12.58, 2.6};
    if(is_equal(type, "W14X550")) return {17.2, 3.82, 17.2, 3.82, 12.56, 2.38};
    if(is_equal(type, "W14X500")) return {17, 3.5, 17, 3.5, 12.6, 2.19};
    if(is_equal(type, "W14X455")) return {16.8, 3.21, 16.8, 3.21, 12.58, 2.02};
    if(is_equal(type, "W14X426")) return {16.7, 3.04, 16.7, 3.04, 12.62, 1.88};
    if(is_equal(type, "W14X398")) return {16.6, 2.85, 16.6, 2.85, 12.6, 1.77};
    if(is_equal(type, "W14X370")) return {16.5, 2.66, 16.5, 2.66, 12.58, 1.66};
    if(is_equal(type, "W14X342")) return {16.4, 2.47, 16.4, 2.47, 12.56, 1.54};
    if(is_equal(type, "W14X311")) return {16.2, 2.26, 16.2, 2.26, 12.58, 1.41};
    if(is_equal(type, "W14X283")) return {16.1, 2.07, 16.1, 2.07, 12.56, 1.29};
    if(is_equal(type, "W14X257")) return {16, 1.89, 16, 1.89, 12.62, 1.18};
    if(is_equal(type, "W14X233")) return {15.9, 1.72, 15.9, 1.72, 12.56, 1.07};
    if(is_equal(type, "W14X211")) return {15.8, 1.56, 15.8, 1.56, 12.58, 0.98};
    if(is_equal(type, "W14X193")) return {15.7, 1.44, 15.7, 1.44, 12.62, 0.89};
    if(is_equal(type, "W14X176")) return {15.7, 1.31, 15.7, 1.31, 12.58, 0.83};
    if(is_equal(type, "W14X159")) return {15.6, 1.19, 15.6, 1.19, 12.62, 0.745};
    if(is_equal(type, "W14X145")) return {15.5, 1.09, 15.5, 1.09, 12.62, 0.68};
    if(is_equal(type, "W14X132")) return {14.7, 1.03, 14.7, 1.03, 12.64, 0.645};
    if(is_equal(type, "W14X120")) return {14.7, 0.94, 14.7, 0.94, 12.62, 0.59};
    if(is_equal(type, "W14X109")) return {14.6, 0.86, 14.6, 0.86, 12.58, 0.525};
    if(is_equal(type, "W14X99")) return {14.6, 0.78, 14.6, 0.78, 12.64, 0.485};
    if(is_equal(type, "W14X90")) return {14.5, 0.71, 14.5, 0.71, 12.58, 0.44};
    if(is_equal(type, "W14X82")) return {10.1, 0.855, 10.1, 0.855, 12.59, 0.51};
    if(is_equal(type, "W14X74")) return {10.1, 0.785, 10.1, 0.785, 12.63, 0.45};
    if(is_equal(type, "W14X68")) return {10, 0.72, 10, 0.72, 12.56, 0.415};
    if(is_equal(type, "W14X61")) return {10, 0.645, 10, 0.645, 12.61, 0.375};
    if(is_equal(type, "W14X53")) return {8.06, 0.66, 8.06, 0.66, 12.58, 0.37};
    if(is_equal(type, "W14X48")) return {8.03, 0.595, 8.03, 0.595, 12.61, 0.34};
    if(is_equal(type, "W14X43")) return {8, 0.53, 8, 0.53, 12.64, 0.305};
    if(is_equal(type, "W14X38")) return {6.77, 0.515, 6.77, 0.515, 13.07, 0.31};
    if(is_equal(type, "W14X34")) return {6.75, 0.455, 6.75, 0.455, 13.09, 0.285};
    if(is_equal(type, "W14X30")) return {6.73, 0.385, 6.73, 0.385, 13.03, 0.27};
    if(is_equal(type, "W14X26")) return {5.03, 0.42, 5.03, 0.42, 13.06, 0.255};
    if(is_equal(type, "W14X22")) return {5, 0.335, 5, 0.335, 13.03, 0.23};
    if(is_equal(type, "W12X336")) return {13.4, 2.96, 13.4, 2.96, 10.88, 1.78};
    if(is_equal(type, "W12X305")) return {13.2, 2.71, 13.2, 2.71, 10.88, 1.63};
    if(is_equal(type, "W12X279")) return {13.1, 2.47, 13.1, 2.47, 10.96, 1.53};
    if(is_equal(type, "W12X252")) return {13, 2.25, 13, 2.25, 10.9, 1.4};
    if(is_equal(type, "W12X230")) return {12.9, 2.07, 12.9, 2.07, 10.96, 1.29};
    if(is_equal(type, "W12X210")) return {12.8, 1.9, 12.8, 1.9, 10.9, 1.18};
    if(is_equal(type, "W12X190")) return {12.7, 1.74, 12.7, 1.74, 10.92, 1.06};
    if(is_equal(type, "W12X170")) return {12.6, 1.56, 12.6, 1.56, 10.88, 0.96};
    if(is_equal(type, "W12X152")) return {12.5, 1.4, 12.5, 1.4, 10.9, 0.87};
    if(is_equal(type, "W12X136")) return {12.4, 1.25, 12.4, 1.25, 10.9, 0.79};
    if(is_equal(type, "W12X120")) return {12.3, 1.11, 12.3, 1.11, 10.88, 0.71};
    if(is_equal(type, "W12X106")) return {12.2, 0.99, 12.2, 0.99, 10.92, 0.61};
    if(is_equal(type, "W12X96")) return {12.2, 0.9, 12.2, 0.9, 10.9, 0.55};
    if(is_equal(type, "W12X87")) return {12.1, 0.81, 12.1, 0.81, 10.88, 0.515};
    if(is_equal(type, "W12X79")) return {12.1, 0.735, 12.1, 0.735, 10.93, 0.47};
    if(is_equal(type, "W12X72")) return {12, 0.67, 12, 0.67, 10.96, 0.43};
    if(is_equal(type, "W12X65")) return {12, 0.605, 12, 0.605, 10.89, 0.39};
    if(is_equal(type, "W12X58")) return {10, 0.64, 10, 0.64, 10.92, 0.36};
    if(is_equal(type, "W12X53")) return {10, 0.575, 10, 0.575, 10.95, 0.345};
    if(is_equal(type, "W12X50")) return {8.08, 0.64, 8.08, 0.64, 10.92, 0.37};
    if(is_equal(type, "W12X45")) return {8.05, 0.575, 8.05, 0.575, 10.95, 0.335};
    if(is_equal(type, "W12X40")) return {8.01, 0.515, 8.01, 0.515, 10.87, 0.295};
    if(is_equal(type, "W12X35")) return {6.56, 0.52, 6.56, 0.52, 11.46, 0.3};
    if(is_equal(type, "W12X30")) return {6.52, 0.44, 6.52, 0.44, 11.42, 0.26};
    if(is_equal(type, "W12X26")) return {6.49, 0.38, 6.49, 0.38, 11.44, 0.23};
    if(is_equal(type, "W12X22")) return {4.03, 0.425, 4.03, 0.425, 11.45, 0.26};
    if(is_equal(type, "W12X19")) return {4.01, 0.35, 4.01, 0.35, 11.5, 0.235};
    if(is_equal(type, "W12X16")) return {3.99, 0.265, 3.99, 0.265, 11.47, 0.22};
    if(is_equal(type, "W12X14")) return {3.97, 0.225, 3.97, 0.225, 11.45, 0.2};
    if(is_equal(type, "W10X112")) return {10.4, 1.25, 10.4, 1.25, 8.9, 0.755};
    if(is_equal(type, "W10X100")) return {10.3, 1.12, 10.3, 1.12, 8.86, 0.68};
    if(is_equal(type, "W10X88")) return {10.3, 0.99, 10.3, 0.99, 8.82, 0.605};
    if(is_equal(type, "W10X77")) return {10.2, 0.87, 10.2, 0.87, 8.86, 0.53};
    if(is_equal(type, "W10X68")) return {10.1, 0.77, 10.1, 0.77, 8.86, 0.47};
    if(is_equal(type, "W10X60")) return {10.1, 0.68, 10.1, 0.68, 8.84, 0.42};
    if(is_equal(type, "W10X54")) return {10, 0.615, 10, 0.615, 8.87, 0.37};
    if(is_equal(type, "W10X49")) return {10, 0.56, 10, 0.56, 8.88, 0.34};
    if(is_equal(type, "W10X45")) return {8.02, 0.62, 8.02, 0.62, 8.86, 0.35};
    if(is_equal(type, "W10X39")) return {7.99, 0.53, 7.99, 0.53, 8.86, 0.315};
    if(is_equal(type, "W10X33")) return {7.96, 0.435, 7.96, 0.435, 8.86, 0.29};
    if(is_equal(type, "W10X30")) return {5.81, 0.51, 5.81, 0.51, 9.48, 0.3};
    if(is_equal(type, "W10X26")) return {5.77, 0.44, 5.77, 0.44, 9.42, 0.26};
    if(is_equal(type, "W10X22")) return {5.75, 0.36, 5.75, 0.36, 9.48, 0.24};
    if(is_equal(type, "W10X19")) return {4.02, 0.395, 4.02, 0.395, 9.41, 0.25};
    if(is_equal(type, "W10X17")) return {4.01, 0.33, 4.01, 0.33, 9.44, 0.24};
    if(is_equal(type, "W10X15")) return {4, 0.27, 4, 0.27, 9.45, 0.23};
    if(is_equal(type, "W10X12")) return {3.96, 0.21, 3.96, 0.21, 9.45, 0.19};
    if(is_equal(type, "W8X67")) return {8.28, 0.935, 8.28, 0.935, 7.13, 0.57};
    if(is_equal(type, "W8X58")) return {8.22, 0.81, 8.22, 0.81, 7.13, 0.51};
    if(is_equal(type, "W8X48")) return {8.11, 0.685, 8.11, 0.685, 7.13, 0.4};
    if(is_equal(type, "W8X40")) return {8.07, 0.56, 8.07, 0.56, 7.13, 0.36};
    if(is_equal(type, "W8X35")) return {8.02, 0.495, 8.02, 0.495, 7.13, 0.31};
    if(is_equal(type, "W8X31")) return {8, 0.435, 8, 0.435, 7.13, 0.285};
    if(is_equal(type, "W8X28")) return {6.54, 0.465, 6.54, 0.465, 7.13, 0.285};
    if(is_equal(type, "W8X24")) return {6.5, 0.4, 6.5, 0.4, 7.13, 0.245};
    if(is_equal(type, "W8X21")) return {5.27, 0.4, 5.27, 0.4, 7.48, 0.25};
    if(is_equal(type, "W8X18")) return {5.25, 0.33, 5.25, 0.33, 7.48, 0.23};
    if(is_equal(type, "W8X15")) return {4.02, 0.315, 4.02, 0.315, 7.48, 0.245};
    if(is_equal(type, "W8X13")) return {4, 0.255, 4, 0.255, 7.48, 0.23};
    if(is_equal(type, "W8X10")) return {3.94, 0.205, 3.94, 0.205, 7.48, 0.17};
    if(is_equal(type, "W6X25")) return {6.08, 0.455, 6.08, 0.455, 5.47, 0.32};
    if(is_equal(type, "W6X20")) return {6.02, 0.365, 6.02, 0.365, 5.47, 0.26};
    if(is_equal(type, "W6X15")) return {5.99, 0.26, 5.99, 0.26, 5.47, 0.23};
    if(is_equal(type, "W6X16")) return {4.03, 0.405, 4.03, 0.405, 5.47, 0.26};
    if(is_equal(type, "W6X12")) return {4, 0.28, 4, 0.28, 5.47, 0.23};
    if(is_equal(type, "W6X9")) return {3.94, 0.215, 3.94, 0.215, 5.47, 0.17};
    if(is_equal(type, "W6X8.5")) return {3.94, 0.195, 3.94, 0.195, 5.44, 0.17};
    if(is_equal(type, "W5X19")) return {5.03, 0.43, 5.03, 0.43, 4.29, 0.27};
    if(is_equal(type, "W5X16")) return {5, 0.36, 5, 0.36, 4.29, 0.24};
    if(is_equal(type, "W4X13")) return {4.06, 0.345, 4.06, 0.345, 3.47, 0.28};
    if(is_equal(type, "M12.5X12.4")) return {3.75, 0.228, 3.75, 0.228, 12.044, 0.155};
    if(is_equal(type, "M12.5X11.6")) return {3.50, 0.211, 3.50, 0.211, 12.078, 0.155};
    if(is_equal(type, "M12X11.8")) return {3.07, 0.225, 3.07, 0.225, 11.550, 0.177};
    if(is_equal(type, "M12X10.8")) return {3.07, 0.210, 3.07, 0.210, 11.580, 0.160};
    if(is_equal(type, "M12X10")) return {3.25, 0.180, 3.25, 0.180, 11.640, 0.149};
    if(is_equal(type, "M10X9")) return {2.69, 0.206, 2.69, 0.206, 9.588, 0.157};
    if(is_equal(type, "M10X8")) return {2.69, 0.182, 2.69, 0.182, 9.586, 0.141};
    if(is_equal(type, "M10X7.5")) return {2.69, 0.173, 2.69, 0.173, 9.644, 0.130};
    if(is_equal(type, "M8X6.5")) return {2.28, 0.189, 2.28, 0.189, 7.622, 0.135};
    if(is_equal(type, "M8X6.2")) return {2.28, 0.177, 2.28, 0.177, 7.646, 0.129};
    if(is_equal(type, "M6X4.4")) return {1.84, 0.171, 1.84, 0.171, 5.658, 0.114};
    if(is_equal(type, "M6X3.7")) return {2.00, 0.129, 2.00, 0.129, 5.662, 0.0980};
    if(is_equal(type, "M5X18.9")) return {5.00, 0.416, 5.00, 0.416, 4.168, 0.316};
    if(is_equal(type, "M4X6")) return {3.80, 0.160, 3.80, 0.160, 3.480, 0.130};
    if(is_equal(type, "M4X4.08")) return {2.25, 0.170, 2.25, 0.170, 3.660, 0.115};
    if(is_equal(type, "M4X3.45")) return {2.25, 0.130, 2.25, 0.130, 3.740, 0.0920};
    if(is_equal(type, "M4X3.2")) return {2.25, 0.130, 2.25, 0.130, 3.740, 0.0920};
    if(is_equal(type, "M3X2.9")) return {2.25, 0.130, 2.25, 0.130, 2.740, 0.0900};
    if(is_equal(type, "S24X121")) return {8.05, 1.09, 8.05, 1.09, 22.320, 0.800};
    if(is_equal(type, "S24X106")) return {7.87, 1.09, 7.87, 1.09, 22.320, 0.620};
    if(is_equal(type, "S24X100")) return {7.25, 0.870, 7.25, 0.870, 22.260, 0.745};
    if(is_equal(type, "S24X90")) return {7.13, 0.870, 7.13, 0.870, 22.260, 0.625};
    if(is_equal(type, "S24X80")) return {7.00, 0.870, 7.00, 0.870, 22.260, 0.500};
    if(is_equal(type, "S20X96")) return {7.20, 0.920, 7.20, 0.920, 18.460, 0.800};
    if(is_equal(type, "S20X86")) return {7.06, 0.920, 7.06, 0.920, 18.460, 0.660};
    if(is_equal(type, "S20X75")) return {6.39, 0.795, 6.39, 0.795, 18.410, 0.635};
    if(is_equal(type, "S20X66")) return {6.26, 0.795, 6.26, 0.795, 18.410, 0.505};
    if(is_equal(type, "S18X70")) return {6.25, 0.691, 6.25, 0.691, 16.618, 0.711};
    if(is_equal(type, "S18X54.7")) return {6.00, 0.691, 6.00, 0.691, 16.618, 0.461};
    if(is_equal(type, "S15X50")) return {5.64, 0.622, 5.64, 0.622, 13.756, 0.550};
    if(is_equal(type, "S15X42.9")) return {5.50, 0.622, 5.50, 0.622, 13.756, 0.411};
    if(is_equal(type, "S12X50")) return {5.48, 0.659, 5.48, 0.659, 10.682, 0.687};
    if(is_equal(type, "S12X40.8")) return {5.25, 0.659, 5.25, 0.659, 10.682, 0.462};
    if(is_equal(type, "S12X35")) return {5.08, 0.544, 5.08, 0.544, 10.912, 0.428};
    if(is_equal(type, "S12X31.8")) return {5.00, 0.544, 5.00, 0.544, 10.912, 0.350};
    if(is_equal(type, "S10X35")) return {4.94, 0.491, 4.94, 0.491, 9.018, 0.594};
    if(is_equal(type, "S10X25.4")) return {4.66, 0.491, 4.66, 0.491, 9.018, 0.311};
    if(is_equal(type, "S8X23")) return {4.17, 0.425, 4.17, 0.425, 7.150, 0.441};
    if(is_equal(type, "S8X18.4")) return {4.00, 0.425, 4.00, 0.425, 7.150, 0.271};
    if(is_equal(type, "S6X17.25")) return {3.57, 0.359, 3.57, 0.359, 5.282, 0.465};
    if(is_equal(type, "S6X12.5")) return {3.33, 0.359, 3.33, 0.359, 5.282, 0.232};
    if(is_equal(type, "S5X10")) return {3.00, 0.326, 3.00, 0.326, 4.348, 0.214};
    if(is_equal(type, "S4X9.5")) return {2.80, 0.293, 2.80, 0.293, 3.414, 0.326};
    if(is_equal(type, "S4X7.7")) return {2.66, 0.293, 2.66, 0.293, 3.414, 0.193};
    if(is_equal(type, "S3X7.5")) return {2.51, 0.260, 2.51, 0.260, 2.480, 0.349};
    if(is_equal(type, "S3X5.7")) return {2.33, 0.260, 2.33, 0.260, 2.480, 0.170};
    if(is_equal(type, "HP18X204")) return {18.1, 1.13, 18.1, 1.13, 16.040, 1.13};
    if(is_equal(type, "HP18X181")) return {18.0, 1.00, 18.0, 1.00, 16.000, 1.00};
    if(is_equal(type, "HP18X157")) return {17.9, 0.870, 17.9, 0.870, 15.960, 0.870};
    if(is_equal(type, "HP18X135")) return {17.8, 0.750, 17.8, 0.750, 16.000, 0.750};
    if(is_equal(type, "HP16X183")) return {16.3, 1.13, 16.3, 1.13, 14.240, 1.13};
    if(is_equal(type, "HP16X162")) return {16.1, 1.00, 16.1, 1.00, 14.300, 1.00};
    if(is_equal(type, "HP16X141")) return {16.0, 0.875, 16.0, 0.875, 14.250, 0.875};
    if(is_equal(type, "HP16X121")) return {15.9, 0.750, 15.9, 0.750, 14.300, 0.750};
    if(is_equal(type, "HP16X101")) return {15.8, 0.625, 15.8, 0.625, 14.250, 0.625};
    if(is_equal(type, "HP16X88")) return {15.7, 0.540, 15.7, 0.540, 14.220, 0.540};
    if(is_equal(type, "HP14X117")) return {14.9, 0.805, 14.9, 0.805, 12.590, 0.805};
    if(is_equal(type, "HP14X102")) return {14.8, 0.705, 14.8, 0.705, 12.590, 0.705};
    if(is_equal(type, "HP14X89")) return {14.7, 0.615, 14.7, 0.615, 12.570, 0.615};
    if(is_equal(type, "HP14X73")) return {14.6, 0.505, 14.6, 0.505, 12.590, 0.505};
    if(is_equal(type, "HP12X89")) return {12.3, 0.720, 12.3, 0.720, 10.960, 0.720};
    if(is_equal(type, "HP12X84")) return {12.3, 0.685, 12.3, 0.685, 10.930, 0.685};
    if(is_equal(type, "HP12X74")) return {12.2, 0.610, 12.2, 0.610, 10.880, 0.605};
    if(is_equal(type, "HP12X63")) return {12.1, 0.515, 12.1, 0.515, 10.870, 0.515};
    if(is_equal(type, "HP12X53")) return {12.0, 0.435, 12.0, 0.435, 10.930, 0.435};
    if(is_equal(type, "HP10X57")) return {10.2, 0.565, 10.2, 0.565, 8.860, 0.565};
    if(is_equal(type, "HP10X42")) return {10.1, 0.420, 10.1, 0.420, 8.860, 0.415};
    if(is_equal(type, "HP8X36")) return {8.16, 0.445, 8.16, 0.445, 7.130, 0.445};
    return {};
}

vec ustsection(const string& type) {
    if(is_equal(type, "WT22X167.5")) return {15.9, 1.77, 20.23, 1.03};
    if(is_equal(type, "WT22X145")) return {15.8, 1.58, 20.22, 0.865};
    if(is_equal(type, "WT22X131")) return {15.8, 1.42, 20.28, 0.785};
    if(is_equal(type, "WT22X115")) return {15.8, 1.22, 20.28, 0.710};
    if(is_equal(type, "WT20X327.5")) return {16.9, 3.54, 18.26, 1.97};
    if(is_equal(type, "WT20X296.5")) return {16.7, 3.23, 18.27, 1.79};
    if(is_equal(type, "WT20X251.5")) return {16.4, 2.76, 18.24, 1.54};
    if(is_equal(type, "WT20X215.5")) return {16.2, 2.36, 18.24, 1.34};
    if(is_equal(type, "WT20X198.5")) return {16.1, 2.20, 18.30, 1.22};
    if(is_equal(type, "WT20X186")) return {16.1, 2.05, 18.25, 1.16};
    if(is_equal(type, "WT20X181")) return {16.0, 2.01, 18.29, 1.12};
    if(is_equal(type, "WT20X162")) return {15.9, 1.81, 18.29, 1.00};
    if(is_equal(type, "WT20X148.5")) return {15.8, 1.65, 18.25, 0.930};
    if(is_equal(type, "WT20X138.5")) return {15.8, 1.58, 18.22, 0.830};
    if(is_equal(type, "WT20X124.5")) return {15.8, 1.42, 18.28, 0.750};
    if(is_equal(type, "WT20X107.5")) return {15.8, 1.22, 18.28, 0.650};
    if(is_equal(type, "WT20X99.5")) return {15.8, 1.07, 18.23, 0.650};
    if(is_equal(type, "WT20X196")) return {12.4, 2.52, 18.28, 1.42};
    if(is_equal(type, "WT20X165.5")) return {12.2, 2.13, 18.27, 1.22};
    if(is_equal(type, "WT20X163.5")) return {12.1, 2.13, 18.27, 1.18};
    if(is_equal(type, "WT20X147")) return {12.0, 1.93, 18.27, 1.06};
    if(is_equal(type, "WT20X139")) return {12.0, 1.81, 18.29, 1.03};
    if(is_equal(type, "WT20X132")) return {11.9, 1.73, 18.27, 0.960};
    if(is_equal(type, "WT20X117.5")) return {11.9, 1.58, 18.22, 0.830};
    if(is_equal(type, "WT20X105.5")) return {11.8, 1.42, 18.28, 0.750};
    if(is_equal(type, "WT20X91.5")) return {11.8, 1.20, 18.30, 0.650};
    if(is_equal(type, "WT20X83.5")) return {11.8, 1.03, 18.27, 0.650};
    if(is_equal(type, "WT20X74.5")) return {11.8, 0.830, 18.27, 0.630};
    if(is_equal(type, "WT18X462.5")) return {18.6, 4.53, 17.07, 3.02};
    if(is_equal(type, "WT18X426.5")) return {18.2, 4.53, 17.07, 2.52};
    if(is_equal(type, "WT18X401")) return {18.0, 4.29, 17.01, 2.38};
    if(is_equal(type, "WT18X361.5")) return {17.8, 3.90, 17.00, 2.17};
    if(is_equal(type, "WT18X326")) return {17.6, 3.54, 16.96, 1.97};
    if(is_equal(type, "WT18X264.5")) return {17.2, 2.91, 16.99, 1.61};
    if(is_equal(type, "WT18X243.5")) return {17.1, 2.68, 17.02, 1.50};
    if(is_equal(type, "WT18X220.5")) return {17.0, 2.44, 16.96, 1.36};
    if(is_equal(type, "WT18X197.5")) return {16.8, 2.20, 17.00, 1.22};
    if(is_equal(type, "WT18X180.5")) return {16.7, 2.01, 16.99, 1.12};
    if(is_equal(type, "WT18X165")) return {16.6, 1.85, 16.95, 1.02};
    if(is_equal(type, "WT18X151")) return {16.7, 1.68, 17.02, 0.945};
    if(is_equal(type, "WT18X141")) return {16.6, 1.57, 17.03, 0.885};
    if(is_equal(type, "WT18X131")) return {16.6, 1.44, 16.96, 0.840};
    if(is_equal(type, "WT18X123.5")) return {16.5, 1.35, 16.95, 0.800};
    if(is_equal(type, "WT18X115.5")) return {16.5, 1.26, 16.94, 0.760};
    if(is_equal(type, "WT18X128")) return {12.2, 1.73, 16.97, 0.960};
    if(is_equal(type, "WT18X116")) return {12.1, 1.57, 17.03, 0.870};
    if(is_equal(type, "WT18X105")) return {12.2, 1.36, 16.94, 0.830};
    if(is_equal(type, "WT18X97")) return {12.1, 1.26, 16.94, 0.765};
    if(is_equal(type, "WT18X91")) return {12.1, 1.18, 17.02, 0.725};
    if(is_equal(type, "WT18X85")) return {12.0, 1.10, 17.00, 0.680};
    if(is_equal(type, "WT18X80")) return {12.0, 1.02, 16.98, 0.650};
    if(is_equal(type, "WT18X75")) return {12.0, 0.940, 16.96, 0.625};
    if(is_equal(type, "WT18X67.5")) return {12.0, 0.790, 17.01, 0.600};
    if(is_equal(type, "WT16.5X193.5")) return {16.2, 2.28, 15.72, 1.26};
    if(is_equal(type, "WT16.5X177")) return {16.1, 2.09, 15.71, 1.16};
    if(is_equal(type, "WT16.5X159")) return {16.0, 1.89, 15.71, 1.04};
    if(is_equal(type, "WT16.5X145.5")) return {15.9, 1.73, 15.67, 0.960};
    if(is_equal(type, "WT16.5X131.5")) return {15.8, 1.57, 15.73, 0.870};
    if(is_equal(type, "WT16.5X120.5")) return {15.9, 1.40, 15.70, 0.830};
    if(is_equal(type, "WT16.5X110.5")) return {15.8, 1.28, 15.72, 0.775};
    if(is_equal(type, "WT16.5X100.5")) return {15.7, 1.15, 15.65, 0.715};
    if(is_equal(type, "WT16.5X84.5")) return {11.5, 1.22, 15.68, 0.670};
    if(is_equal(type, "WT16.5X76")) return {11.6, 1.06, 15.64, 0.635};
    if(is_equal(type, "WT16.5X70.5")) return {11.5, 0.960, 15.74, 0.605};
    if(is_equal(type, "WT16.5X65")) return {11.5, 0.855, 15.65, 0.580};
    if(is_equal(type, "WT16.5X59")) return {11.5, 0.740, 15.66, 0.550};
    if(is_equal(type, "WT15X195.5")) return {15.6, 2.44, 14.16, 1.36};
    if(is_equal(type, "WT15X178.5")) return {15.5, 2.24, 14.16, 1.24};
    if(is_equal(type, "WT15X163")) return {15.4, 2.05, 14.15, 1.14};
    if(is_equal(type, "WT15X146")) return {15.3, 1.85, 14.15, 1.02};
    if(is_equal(type, "WT15X130.5")) return {15.2, 1.65, 14.15, 0.930};
    if(is_equal(type, "WT15X117.5")) return {15.1, 1.50, 14.20, 0.830};
    if(is_equal(type, "WT15X105.5")) return {15.1, 1.32, 14.18, 0.775};
    if(is_equal(type, "WT15X95.5")) return {15.0, 1.19, 14.11, 0.710};
    if(is_equal(type, "WT15X86.5")) return {15.0, 1.07, 14.13, 0.655};
    if(is_equal(type, "WT15X74")) return {10.5, 1.18, 14.12, 0.650};
    if(is_equal(type, "WT15X66")) return {10.5, 1.00, 14.20, 0.615};
    if(is_equal(type, "WT15X62")) return {10.5, 0.930, 14.17, 0.585};
    if(is_equal(type, "WT15X58")) return {10.5, 0.850, 14.15, 0.565};
    if(is_equal(type, "WT15X54")) return {10.5, 0.760, 14.14, 0.545};
    if(is_equal(type, "WT15X49.5")) return {10.5, 0.670, 14.13, 0.520};
    if(is_equal(type, "WT15X45")) return {10.4, 0.610, 14.19, 0.470};
    if(is_equal(type, "WT13.5X269.5")) return {15.3, 3.54, 12.76, 1.97};
    if(is_equal(type, "WT13.5X184")) return {14.7, 2.48, 12.72, 1.38};
    if(is_equal(type, "WT13.5X168")) return {14.6, 2.28, 12.72, 1.26};
    if(is_equal(type, "WT13.5X153.5")) return {14.4, 2.09, 12.71, 1.16};
    if(is_equal(type, "WT13.5X140.5")) return {14.4, 1.93, 12.67, 1.06};
    if(is_equal(type, "WT13.5X129")) return {14.3, 1.77, 12.73, 0.980};
    if(is_equal(type, "WT13.5X117.5")) return {14.2, 1.61, 12.69, 0.910};
    if(is_equal(type, "WT13.5X108.5")) return {14.1, 1.50, 12.70, 0.830};
    if(is_equal(type, "WT13.5X97")) return {14.0, 1.34, 12.76, 0.750};
    if(is_equal(type, "WT13.5X89")) return {14.1, 1.19, 12.71, 0.725};
    if(is_equal(type, "WT13.5X80.5")) return {14.0, 1.08, 12.72, 0.660};
    if(is_equal(type, "WT13.5X73")) return {14.0, 0.975, 12.73, 0.605};
    if(is_equal(type, "WT13.5X64.5")) return {10.0, 1.10, 12.70, 0.610};
    if(is_equal(type, "WT13.5X57")) return {10.1, 0.930, 12.67, 0.570};
    if(is_equal(type, "WT13.5X51")) return {10.0, 0.830, 12.67, 0.515};
    if(is_equal(type, "WT13.5X47")) return {10.0, 0.745, 12.76, 0.490};
    if(is_equal(type, "WT13.5X42")) return {10.0, 0.640, 12.76, 0.460};
    if(is_equal(type, "WT12X185")) return {13.7, 2.72, 11.28, 1.52};
    if(is_equal(type, "WT12X167.5")) return {13.5, 2.48, 11.32, 1.38};
    if(is_equal(type, "WT12X153")) return {13.4, 2.28, 11.32, 1.26};
    if(is_equal(type, "WT12X139.5")) return {13.3, 2.09, 11.31, 1.16};
    if(is_equal(type, "WT12X125")) return {13.2, 1.89, 11.31, 1.04};
    if(is_equal(type, "WT12X114.5")) return {13.1, 1.73, 11.27, 0.960};
    if(is_equal(type, "WT12X103.5")) return {13.0, 1.57, 11.33, 0.870};
    if(is_equal(type, "WT12X96")) return {13.0, 1.46, 11.24, 0.810};
    if(is_equal(type, "WT12X88")) return {12.9, 1.34, 11.26, 0.750};
    if(is_equal(type, "WT12X81")) return {13.0, 1.22, 11.28, 0.705};
    if(is_equal(type, "WT12X73")) return {12.9, 1.09, 11.31, 0.650};
    if(is_equal(type, "WT12X65.5")) return {12.9, 0.960, 11.24, 0.605};
    if(is_equal(type, "WT12X58.5")) return {12.8, 0.850, 11.25, 0.550};
    if(is_equal(type, "WT12X52")) return {12.8, 0.750, 11.25, 0.500};
    if(is_equal(type, "WT12X51.5")) return {9.00, 0.980, 11.32, 0.550};
    if(is_equal(type, "WT12X47")) return {9.07, 0.875, 11.33, 0.515};
    if(is_equal(type, "WT12X42")) return {9.02, 0.770, 11.33, 0.470};
    if(is_equal(type, "WT12X38")) return {8.99, 0.680, 11.32, 0.440};
    if(is_equal(type, "WT12X34")) return {8.97, 0.585, 11.32, 0.415};
    if(is_equal(type, "WT12X31")) return {7.04, 0.590, 11.31, 0.430};
    if(is_equal(type, "WT12X27.5")) return {7.01, 0.505, 11.30, 0.395};
    if(is_equal(type, "WT10.5X137.5")) return {12.9, 2.19, 9.91, 1.22};
    if(is_equal(type, "WT10.5X124")) return {12.8, 1.99, 9.91, 1.10};
    if(is_equal(type, "WT10.5X111.5")) return {12.7, 1.79, 9.91, 1.00};
    if(is_equal(type, "WT10.5X100.5")) return {12.6, 1.63, 9.87, 0.910};
    if(is_equal(type, "WT10.5X91")) return {12.5, 1.48, 9.92, 0.830};
    if(is_equal(type, "WT10.5X83")) return {12.4, 1.36, 9.84, 0.750};
    if(is_equal(type, "WT10.5X73.5")) return {12.5, 1.15, 9.85, 0.720};
    if(is_equal(type, "WT10.5X66")) return {12.4, 1.04, 9.86, 0.650};
    if(is_equal(type, "WT10.5X61")) return {12.4, 0.960, 9.84, 0.600};
    if(is_equal(type, "WT10.5X55.5")) return {12.3, 0.875, 9.93, 0.550};
    if(is_equal(type, "WT10.5X50.5")) return {12.3, 0.800, 9.90, 0.500};
    if(is_equal(type, "WT10.5X46.5")) return {8.42, 0.930, 9.87, 0.580};
    if(is_equal(type, "WT10.5X41.5")) return {8.36, 0.835, 9.87, 0.515};
    if(is_equal(type, "WT10.5X36.5")) return {8.30, 0.740, 9.86, 0.455};
    if(is_equal(type, "WT10.5X34")) return {8.27, 0.685, 9.92, 0.430};
    if(is_equal(type, "WT10.5X31")) return {8.24, 0.615, 9.89, 0.400};
    if(is_equal(type, "WT10.5X27.5")) return {8.22, 0.522, 9.88, 0.375};
    if(is_equal(type, "WT10.5X24")) return {8.14, 0.430, 9.87, 0.350};
    if(is_equal(type, "WT10.5X28.5")) return {6.56, 0.650, 9.85, 0.405};
    if(is_equal(type, "WT10.5X25")) return {6.53, 0.535, 9.87, 0.380};
    if(is_equal(type, "WT10.5X22")) return {6.50, 0.450, 9.85, 0.350};
    if(is_equal(type, "WT9X155.5")) return {12.0, 2.74, 8.46, 1.52};
    if(is_equal(type, "WT9X141.5")) return {11.9, 2.50, 8.40, 1.40};
    if(is_equal(type, "WT9X129")) return {11.8, 2.30, 8.40, 1.28};
    if(is_equal(type, "WT9X117")) return {11.7, 2.11, 8.39, 1.16};
    if(is_equal(type, "WT9X105.5")) return {11.6, 1.91, 8.39, 1.06};
    if(is_equal(type, "WT9X96")) return {11.5, 1.75, 8.45, 0.960};
    if(is_equal(type, "WT9X87.5")) return {11.4, 1.59, 8.41, 0.890};
    if(is_equal(type, "WT9X79")) return {11.3, 1.44, 8.42, 0.810};
    if(is_equal(type, "WT9X71.5")) return {11.2, 1.32, 8.43, 0.730};
    if(is_equal(type, "WT9X65")) return {11.2, 1.20, 8.43, 0.670};
    if(is_equal(type, "WT9X59.5")) return {11.3, 1.06, 8.43, 0.655};
    if(is_equal(type, "WT9X53")) return {11.2, 0.940, 8.43, 0.590};
    if(is_equal(type, "WT9X48.5")) return {11.1, 0.870, 8.43, 0.535};
    if(is_equal(type, "WT9X43")) return {11.1, 0.770, 8.43, 0.480};
    if(is_equal(type, "WT9X38")) return {11.0, 0.680, 8.43, 0.425};
    if(is_equal(type, "WT9X35.5")) return {7.64, 0.810, 8.43, 0.495};
    if(is_equal(type, "WT9X32.5")) return {7.59, 0.750, 8.43, 0.450};
    if(is_equal(type, "WT9X30")) return {7.56, 0.695, 8.43, 0.415};
    if(is_equal(type, "WT9X27.5")) return {7.53, 0.630, 8.43, 0.390};
    if(is_equal(type, "WT9X25")) return {7.50, 0.570, 8.43, 0.355};
    if(is_equal(type, "WT9X23")) return {6.06, 0.605, 8.43, 0.360};
    if(is_equal(type, "WT9X20")) return {6.02, 0.525, 8.43, 0.315};
    if(is_equal(type, "WT9X17.5")) return {6.00, 0.425, 8.43, 0.300};
    if(is_equal(type, "WT8X50")) return {10.4, 0.985, 7.51, 0.585};
    if(is_equal(type, "WT8X44.5")) return {10.4, 0.875, 7.51, 0.525};
    if(is_equal(type, "WT8X38.5")) return {10.3, 0.760, 7.50, 0.455};
    if(is_equal(type, "WT8X33.5")) return {10.2, 0.665, 7.51, 0.395};
    if(is_equal(type, "WT8X28.5")) return {7.12, 0.715, 7.51, 0.430};
    if(is_equal(type, "WT8X25")) return {7.07, 0.630, 7.50, 0.380};
    if(is_equal(type, "WT8X22.5")) return {7.04, 0.565, 7.51, 0.345};
    if(is_equal(type, "WT8X20")) return {7.00, 0.505, 7.51, 0.305};
    if(is_equal(type, "WT8X18")) return {6.99, 0.430, 7.50, 0.295};
    if(is_equal(type, "WT8X15.5")) return {5.53, 0.440, 7.50, 0.275};
    if(is_equal(type, "WT8X13")) return {5.50, 0.345, 7.51, 0.250};
    if(is_equal(type, "WT7X436.5")) return {18.8, 5.51, 6.29, 3.94};
    if(is_equal(type, "WT7X404")) return {18.6, 5.12, 6.28, 3.74};
    if(is_equal(type, "WT7X365")) return {17.9, 4.91, 6.29, 3.07};
    if(is_equal(type, "WT7X332.5")) return {17.7, 4.52, 6.28, 2.83};
    if(is_equal(type, "WT7X302.5")) return {17.4, 4.16, 6.34, 2.60};
    if(is_equal(type, "WT7X275")) return {17.2, 3.82, 6.28, 2.38};
    if(is_equal(type, "WT7X250")) return {17.0, 3.50, 6.30, 2.19};
    if(is_equal(type, "WT7X227.5")) return {16.8, 3.21, 6.30, 2.02};
    if(is_equal(type, "WT7X213")) return {16.7, 3.04, 6.30, 1.88};
    if(is_equal(type, "WT7X199")) return {16.6, 2.85, 6.30, 1.77};
    if(is_equal(type, "WT7X185")) return {16.5, 2.66, 6.30, 1.66};
    if(is_equal(type, "WT7X171")) return {16.4, 2.47, 6.30, 1.54};
    if(is_equal(type, "WT7X155.5")) return {16.2, 2.26, 6.30, 1.41};
    if(is_equal(type, "WT7X141.5")) return {16.1, 2.07, 6.30, 1.29};
    if(is_equal(type, "WT7X128.5")) return {16.0, 1.89, 6.30, 1.18};
    if(is_equal(type, "WT7X116.5")) return {15.9, 1.72, 6.30, 1.07};
    if(is_equal(type, "WT7X105.5")) return {15.8, 1.56, 6.30, 0.980};
    if(is_equal(type, "WT7X96.5")) return {15.7, 1.44, 6.30, 0.890};
    if(is_equal(type, "WT7X88")) return {15.7, 1.31, 6.30, 0.830};
    if(is_equal(type, "WT7X79.5")) return {15.6, 1.19, 6.30, 0.745};
    if(is_equal(type, "WT7X72.5")) return {15.5, 1.09, 6.30, 0.680};
    if(is_equal(type, "WT7X66")) return {14.7, 1.03, 6.30, 0.645};
    if(is_equal(type, "WT7X60")) return {14.7, 0.940, 6.30, 0.590};
    if(is_equal(type, "WT7X54.5")) return {14.6, 0.860, 6.30, 0.525};
    if(is_equal(type, "WT7X49.5")) return {14.6, 0.780, 6.30, 0.485};
    if(is_equal(type, "WT7X45")) return {14.5, 0.710, 6.30, 0.440};
    if(is_equal(type, "WT7X41")) return {10.1, 0.855, 6.31, 0.510};
    if(is_equal(type, "WT7X37")) return {10.1, 0.785, 6.31, 0.450};
    if(is_equal(type, "WT7X34")) return {10.0, 0.720, 6.30, 0.415};
    if(is_equal(type, "WT7X30.5")) return {10.0, 0.645, 6.31, 0.375};
    if(is_equal(type, "WT7X26.5")) return {8.06, 0.660, 6.30, 0.370};
    if(is_equal(type, "WT7X24")) return {8.03, 0.595, 6.31, 0.340};
    if(is_equal(type, "WT7X21.5")) return {8.00, 0.530, 6.30, 0.305};
    if(is_equal(type, "WT7X19")) return {6.77, 0.515, 6.54, 0.310};
    if(is_equal(type, "WT7X17")) return {6.75, 0.455, 6.54, 0.285};
    if(is_equal(type, "WT7X15")) return {6.73, 0.385, 6.54, 0.270};
    if(is_equal(type, "WT7X13")) return {5.03, 0.420, 6.54, 0.255};
    if(is_equal(type, "WT7X11")) return {5.00, 0.335, 6.54, 0.230};
    if(is_equal(type, "WT6X168")) return {13.4, 2.96, 5.45, 1.78};
    if(is_equal(type, "WT6X152.5")) return {13.2, 2.71, 5.45, 1.63};
    if(is_equal(type, "WT6X139.5")) return {13.1, 2.47, 5.46, 1.53};
    if(is_equal(type, "WT6X126")) return {13.0, 2.25, 5.46, 1.40};
    if(is_equal(type, "WT6X115")) return {12.9, 2.07, 5.46, 1.29};
    if(is_equal(type, "WT6X105")) return {12.8, 1.90, 5.46, 1.18};
    if(is_equal(type, "WT6X95")) return {12.7, 1.74, 5.45, 1.06};
    if(is_equal(type, "WT6X85")) return {12.6, 1.56, 5.46, 0.960};
    if(is_equal(type, "WT6X76")) return {12.5, 1.40, 5.46, 0.870};
    if(is_equal(type, "WT6X68")) return {12.4, 1.25, 5.46, 0.790};
    if(is_equal(type, "WT6X60")) return {12.3, 1.11, 5.45, 0.710};
    if(is_equal(type, "WT6X53")) return {12.2, 0.990, 5.46, 0.610};
    if(is_equal(type, "WT6X48")) return {12.2, 0.900, 5.46, 0.550};
    if(is_equal(type, "WT6X43.5")) return {12.1, 0.810, 5.46, 0.515};
    if(is_equal(type, "WT6X39.5")) return {12.1, 0.735, 5.46, 0.470};
    if(is_equal(type, "WT6X36")) return {12.0, 0.670, 5.46, 0.430};
    if(is_equal(type, "WT6X32.5")) return {12.0, 0.605, 5.46, 0.390};
    if(is_equal(type, "WT6X29")) return {10.0, 0.640, 5.46, 0.360};
    if(is_equal(type, "WT6X26.5")) return {10.0, 0.575, 5.46, 0.345};
    if(is_equal(type, "WT6X25")) return {8.08, 0.640, 5.46, 0.370};
    if(is_equal(type, "WT6X22.5")) return {8.05, 0.575, 5.46, 0.335};
    if(is_equal(type, "WT6X20")) return {8.01, 0.515, 5.46, 0.295};
    if(is_equal(type, "WT6X17.5")) return {6.56, 0.520, 5.73, 0.300};
    if(is_equal(type, "WT6X15")) return {6.52, 0.440, 5.73, 0.260};
    if(is_equal(type, "WT6X13")) return {6.49, 0.380, 5.73, 0.230};
    if(is_equal(type, "WT6X11")) return {4.03, 0.425, 5.74, 0.260};
    if(is_equal(type, "WT6X9.5")) return {4.01, 0.350, 5.73, 0.235};
    if(is_equal(type, "WT6X8")) return {3.99, 0.265, 5.74, 0.220};
    if(is_equal(type, "WT6X7")) return {3.97, 0.225, 5.74, 0.200};
    if(is_equal(type, "WT5X56")) return {10.4, 1.25, 4.43, 0.755};
    if(is_equal(type, "WT5X50")) return {10.3, 1.12, 4.43, 0.680};
    if(is_equal(type, "WT5X44")) return {10.3, 0.990, 4.43, 0.605};
    if(is_equal(type, "WT5X38.5")) return {10.2, 0.870, 4.43, 0.530};
    if(is_equal(type, "WT5X34")) return {10.1, 0.770, 4.43, 0.470};
    if(is_equal(type, "WT5X30")) return {10.1, 0.680, 4.43, 0.420};
    if(is_equal(type, "WT5X27")) return {10.0, 0.615, 4.44, 0.370};
    if(is_equal(type, "WT5X24.5")) return {10.0, 0.560, 4.43, 0.340};
    if(is_equal(type, "WT5X22.5")) return {8.02, 0.620, 4.43, 0.350};
    if(is_equal(type, "WT5X19.5")) return {7.99, 0.530, 4.43, 0.315};
    if(is_equal(type, "WT5X16.5")) return {7.96, 0.435, 4.44, 0.290};
    if(is_equal(type, "WT5X15")) return {5.81, 0.510, 4.73, 0.300};
    if(is_equal(type, "WT5X13")) return {5.77, 0.440, 4.73, 0.260};
    if(is_equal(type, "WT5X11")) return {5.75, 0.360, 4.73, 0.240};
    if(is_equal(type, "WT5X9.5")) return {4.02, 0.395, 4.73, 0.250};
    if(is_equal(type, "WT5X8.5")) return {4.01, 0.330, 4.73, 0.240};
    if(is_equal(type, "WT5X7.5")) return {4.00, 0.270, 4.73, 0.230};
    if(is_equal(type, "WT5X6")) return {3.96, 0.210, 4.73, 0.190};
    if(is_equal(type, "WT4X33.5")) return {8.28, 0.935, 3.57, 0.570};
    if(is_equal(type, "WT4X29")) return {8.22, 0.810, 3.57, 0.510};
    if(is_equal(type, "WT4X24")) return {8.11, 0.685, 3.57, 0.400};
    if(is_equal(type, "WT4X20")) return {8.07, 0.560, 3.57, 0.360};
    if(is_equal(type, "WT4X17.5")) return {8.02, 0.495, 3.57, 0.310};
    if(is_equal(type, "WT4X15.5")) return {8.00, 0.435, 3.57, 0.285};
    if(is_equal(type, "WT4X14")) return {6.54, 0.465, 3.57, 0.285};
    if(is_equal(type, "WT4X12")) return {6.50, 0.400, 3.57, 0.245};
    if(is_equal(type, "WT4X10.5")) return {5.27, 0.400, 3.74, 0.250};
    if(is_equal(type, "WT4X9")) return {5.25, 0.330, 3.74, 0.230};
    if(is_equal(type, "WT4X7.5")) return {4.02, 0.315, 3.75, 0.245};
    if(is_equal(type, "WT4X6.5")) return {4.00, 0.255, 3.75, 0.230};
    if(is_equal(type, "WT4X5")) return {3.94, 0.205, 3.75, 0.170};
    if(is_equal(type, "WT3X12.5")) return {6.08, 0.455, 2.74, 0.320};
    if(is_equal(type, "WT3X10")) return {6.02, 0.365, 2.74, 0.260};
    if(is_equal(type, "WT3X7.5")) return {5.99, 0.260, 2.74, 0.230};
    if(is_equal(type, "WT3X8")) return {4.03, 0.405, 2.74, 0.260};
    if(is_equal(type, "WT3X6")) return {4.00, 0.280, 2.74, 0.230};
    if(is_equal(type, "WT3X4.5")) return {3.94, 0.215, 2.74, 0.170};
    if(is_equal(type, "WT3X4.25")) return {3.94, 0.195, 2.73, 0.170};
    if(is_equal(type, "WT2.5X9.5")) return {5.03, 0.430, 2.15, 0.270};
    if(is_equal(type, "WT2.5X8")) return {5.00, 0.360, 2.15, 0.240};
    if(is_equal(type, "WT2X6.5")) return {4.06, 0.345, 1.74, 0.280};
    if(is_equal(type, "MT6.25X6.2")) return {3.75, 0.228, 6.04, 0.155};
    if(is_equal(type, "MT6.25X5.8")) return {3.50, 0.211, 6.04, 0.155};
    if(is_equal(type, "MT6X5.9")) return {3.07, 0.225, 5.78, 0.177};
    if(is_equal(type, "MT6X5.4")) return {3.07, 0.210, 5.78, 0.160};
    if(is_equal(type, "MT6X5")) return {3.25, 0.180, 5.81, 0.149};
    if(is_equal(type, "MT5X4.5")) return {2.69, 0.206, 4.79, 0.157};
    if(is_equal(type, "MT5X4")) return {2.69, 0.182, 4.80, 0.141};
    if(is_equal(type, "MT5X3.75")) return {2.69, 0.173, 4.83, 0.130};
    if(is_equal(type, "MT4X3.25")) return {2.28, 0.189, 3.81, 0.135};
    if(is_equal(type, "MT4X3.1")) return {2.28, 0.177, 3.82, 0.129};
    if(is_equal(type, "MT3X2.2")) return {1.84, 0.171, 2.83, 0.114};
    if(is_equal(type, "MT3X1.85")) return {2.00, 0.129, 2.83, 0.0980};
    if(is_equal(type, "MT2.5X9.45")) return {5.00, 0.416, 2.08, 0.316};
    if(is_equal(type, "MT2X3")) return {3.80, 0.160, 1.74, 0.130};
    if(is_equal(type, "ST12X60.5")) return {8.05, 1.09, 11.21, 0.800};
    if(is_equal(type, "ST12X53")) return {7.87, 1.09, 11.21, 0.620};
    if(is_equal(type, "ST12X50")) return {7.25, 0.870, 11.13, 0.745};
    if(is_equal(type, "ST12X45")) return {7.13, 0.870, 11.13, 0.625};
    if(is_equal(type, "ST12X40")) return {7.00, 0.870, 11.13, 0.500};
    if(is_equal(type, "ST10X48")) return {7.20, 0.920, 9.28, 0.800};
    if(is_equal(type, "ST10X43")) return {7.06, 0.920, 9.28, 0.660};
    if(is_equal(type, "ST10X37.5")) return {6.39, 0.795, 9.21, 0.635};
    if(is_equal(type, "ST10X33")) return {6.26, 0.795, 9.21, 0.505};
    if(is_equal(type, "ST9X35")) return {6.25, 0.691, 8.31, 0.711};
    if(is_equal(type, "ST9X27.35")) return {6.00, 0.691, 8.31, 0.461};
    if(is_equal(type, "ST7.5X25")) return {5.64, 0.622, 6.88, 0.550};
    if(is_equal(type, "ST7.5X21.45")) return {5.50, 0.622, 6.88, 0.411};
    if(is_equal(type, "ST6X25")) return {5.48, 0.659, 5.34, 0.687};
    if(is_equal(type, "ST6X20.4")) return {5.25, 0.659, 5.34, 0.462};
    if(is_equal(type, "ST6X17.5")) return {5.08, 0.544, 5.46, 0.428};
    if(is_equal(type, "ST6X15.9")) return {5.00, 0.544, 5.46, 0.350};
    if(is_equal(type, "ST5X17.5")) return {4.94, 0.491, 4.51, 0.594};
    if(is_equal(type, "ST5X12.7")) return {4.66, 0.491, 4.51, 0.311};
    if(is_equal(type, "ST4X11.5")) return {4.17, 0.425, 3.58, 0.441};
    if(is_equal(type, "ST4X9.2")) return {4.00, 0.425, 3.58, 0.271};
    if(is_equal(type, "ST3X8.6")) return {3.57, 0.359, 2.64, 0.465};
    if(is_equal(type, "ST3X6.25")) return {3.33, 0.359, 2.64, 0.232};
    if(is_equal(type, "ST2.5X5")) return {3.00, 0.326, 2.17, 0.214};
    if(is_equal(type, "ST2X4.75")) return {2.80, 0.293, 1.71, 0.326};
    if(is_equal(type, "ST2X3.85")) return {2.66, 0.293, 1.71, 0.193};
    if(is_equal(type, "ST1.5X3.75")) return {2.51, 0.260, 1.24, 0.349};
    if(is_equal(type, "ST1.5X2.85")) return {2.33, 0.260, 1.24, 0.170};
    return {};
}

void new_eu2d(unique_ptr<Section>& return_obj, istringstream& command) {
    string type;
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

    unsigned int_pt = 6;
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

void new_eu3d(unique_ptr<Section>& return_obj, istringstream& command) {
    string type;
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

    unsigned int_pt = 6;
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

void new_nz2d(unique_ptr<Section>& return_obj, istringstream& command) {
    string type;
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

    unsigned int_pt = 6;
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
        return_obj = make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
        return;
    }

    dim = nzshsection(type);

    if(!dim.is_empty()) {
        return_obj = make_unique<Box2D>(tag, scale * dim, material_id, int_pt, eccentricity);
        return;
    }

    suanpan_error("Cannot identify section type.\n");
}

void new_nz3d(unique_ptr<Section>& return_obj, istringstream& command) {
    string type;
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

    unsigned int_pt = 6;
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
        return_obj = make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    dim = nzshsection(type);

    if(!dim.is_empty()) {
        return_obj = make_unique<Box3D>(tag, scale * dim, material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    suanpan_error("Cannot identify section type.\n");
}

void new_us2d(unique_ptr<Section>& return_obj, istringstream& command, const bool recenter) {
    string type;
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

    unsigned int_pt = 6;
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

    dim = ustsection(type);

    if(!dim.is_empty()) {
        if(recenter) eccentricity = barycenter(dim *= scale);
        return_obj = make_unique<TSection2D>(tag, std::move(dim), material_id, int_pt, eccentricity);
        return;
    }

    suanpan_error("Cannot identify section type.\n");
}

void new_us3d(unique_ptr<Section>& return_obj, istringstream& command, const bool recenter) {
    string type;
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

    unsigned int_pt = 6;
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

    dim = ustsection(type);

    if(!dim.is_empty()) {
        if(recenter) eccentricity_y = barycenter(dim *= scale);
        return_obj = make_unique<TSection3D>(tag, std::move(dim), material_id, int_pt, vec{eccentricity_y, eccentricity_z});
        return;
    }

    suanpan_error("Cannot identify section type.\n");
}

int create_new_section(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string section_id;
    if(!get_input(command, section_id)) {
        suanpan_error("A valid section type is required.\n");
        return 0;
    }

    unique_ptr<Section> new_section = nullptr;

    if(is_equal(section_id, "Bar2D")) new_bar2d(new_section, command);
    else if(is_equal(section_id, "Bar3D")) new_bar3d(new_section, command);
    else if(is_equal(section_id, "Box2D")) new_box2d(new_section, command);
    else if(is_equal(section_id, "Box3D")) new_box3d(new_section, command);
    else if(is_equal(section_id, "Circle1D")) new_circle1d(new_section, command);
    else if(is_equal(section_id, "Circle2D")) new_circle2d(new_section, command);
    else if(is_equal(section_id, "Circle3D")) new_circle3d(new_section, command);
    else if(is_equal(section_id, "CircularHollow2D")) new_circularhollow2D(new_section, command);
    else if(is_equal(section_id, "CircularHollow3D")) new_circularhollow3D(new_section, command);
    else if(is_equal(section_id, "Fibre1D")) new_fibre1d(new_section, command);
    else if(is_equal(section_id, "Fibre2D")) new_fibre2d(new_section, command);
    else if(is_equal(section_id, "Fibre3D")) new_fibre3d(new_section, command);
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
