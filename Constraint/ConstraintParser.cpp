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

#include "ConstraintParser.h"
#include <Constraint/Constraint>
#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Recorder/OutputType.h>

void new_bc(unique_ptr<Constraint>& return_obj, istringstream& command, const bool penalty, const bool group) {
    unsigned bc_id;
    if(!get_input(command, bc_id)) {
        suanpan_error("new_bc() needs a valid tag.\n");
        return;
    }

    string dof_id;
    if(!get_input(command, dof_id)) {
        suanpan_error("new_bc() needs a valid DoF identifier.\n");
        return;
    }

    const auto bc_type = suanpan::to_lower(dof_id[0]);

    if(!is_equal(bc_type, 'p') && !is_equal(bc_type, 'e') && !is_equal(bc_type, 'x') && !is_equal(bc_type, 'y') && !is_equal(bc_type, 'z') && !is_equal(bc_type, '1') && !is_equal(bc_type, '2') && !is_equal(bc_type, '3') && !is_equal(bc_type, '4') && !is_equal(bc_type, '5') && !is_equal(bc_type, '6')) {
        suanpan_error("new_bc() needs a valid DoF identifier.\n");
        return;
    }

    uword obj;
    vector<uword> obj_tag;
    while(get_input(command, obj)) obj_tag.push_back(obj);

    penalty ? group ? return_obj = make_unique<GroupPenaltyBC>(bc_id, 0, uvec(obj_tag), bc_type) : return_obj = make_unique<PenaltyBC>(bc_id, 0, uvec(obj_tag), bc_type) : group ? return_obj = make_unique<GroupMultiplierBC>(bc_id, 0, uvec(obj_tag), bc_type) : return_obj = make_unique<MultiplierBC>(bc_id, 0, uvec(obj_tag), bc_type);
}

void new_fixedlength(unique_ptr<Constraint>& return_obj, istringstream& command, const unsigned dof) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_fixedlength() needs a valid tag.\n");
        return;
    }

    uword node_i, node_j;
    if(!get_input(command, node_i) || !get_input(command, node_j)) {
        suanpan_error("new_fixedlength() needs two node tags.\n");
        return;
    }

    return_obj = make_unique<FixedLength>(tag, 0, dof, uvec{node_i, node_j});
}

void new_minimumgap(unique_ptr<Constraint>& return_obj, istringstream& command, const unsigned dof) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_minimumgap() needs a valid tag.\n");
        return;
    }

    uword node_i, node_j;
    if(!get_input(command, node_i) || !get_input(command, node_j)) {
        suanpan_error("new_minimumgap() needs two node tags.\n");
        return;
    }

    double gap;
    if(!get_input(command, gap)) {
        suanpan_error("new_minimumgap() needs a valid minimum gap.\n");
        return;
    }

    return_obj = make_unique<MinimumGap>(tag, 0, dof, gap, uvec{node_i, node_j});
}

void new_maximumgap(unique_ptr<Constraint>& return_obj, istringstream& command, const unsigned dof) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_maximumgap() needs a valid tag.\n");
        return;
    }

    uword node_i, node_j;
    if(!get_input(command, node_i) || !get_input(command, node_j)) {
        suanpan_error("new_maximumgap() needs two node tags.\n");
        return;
    }

    double gap;
    if(!get_input(command, gap)) {
        suanpan_error("new_maximumgap() needs a valid minimum gap.\n");
        return;
    }

    return_obj = make_unique<MaximumGap>(tag, 0, dof, gap, uvec{node_i, node_j});
}

void new_sleeve(unique_ptr<Constraint>& return_obj, istringstream& command, const unsigned dof) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_sleeve() needs a valid tag.\n");
        return;
    }

    uword node_i, node_j;
    if(!get_input(command, node_i) || !get_input(command, node_j)) {
        suanpan_error("new_sleeve() needs two node tags.\n");
        return;
    }

    double min_gap, max_gap;
    if(!get_input(command, min_gap)) {
        suanpan_error("new_sleeve() needs a valid minimum gap.\n");
        return;
    }
    if(!get_input(command, max_gap)) {
        suanpan_error("new_sleeve() needs a valid maximum gap.\n");
        return;
    }

    return_obj = make_unique<Sleeve>(tag, 0, dof, min_gap, max_gap, uvec{node_i, node_j});
}

void new_embed(unique_ptr<Constraint>& return_obj, istringstream& command, const unsigned dof) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_embed() needs a valid tag.\n");
        return;
    }

    unsigned element_tag;
    if(!get_input(command, element_tag)) {
        suanpan_error("new_embed() needs a valid element tag.\n");
        return;
    }

    unsigned node_tag;
    if(!get_input(command, node_tag)) {
        suanpan_error("new_embed() needs a valid node tag.\n");
        return;
    }

    if(2 == dof) return_obj = make_unique<Embed2D>(tag, 0, element_tag, node_tag);
    else return_obj = make_unique<Embed3D>(tag, 0, element_tag, node_tag);
}

void new_mpc(unique_ptr<Constraint>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_mpc() needs a valid tag.\n");
        return;
    }

    unsigned amplitude;
    if(!get_input(command, amplitude)) {
        suanpan_error("new_mpc() needs a valid amplitude tag.\n");
        return;
    }

    double magnitude;
    if(!get_input(command, magnitude)) {
        suanpan_error("new_mpc() needs a valid magnitude.\n");
        return;
    }

    vector<uword> node_tag, dof_tag;
    vector<double> weight_tag;
    while(!command.eof()) {
        double weight;
        uword dof, node;
        if(!get_input(command, node) || !get_input(command, dof) || !get_input(command, weight)) return;
        node_tag.emplace_back(node);
        dof_tag.emplace_back(dof);
        weight_tag.emplace_back(weight);
    }

    return_obj = make_unique<MPC>(tag, 0, amplitude, uvec(node_tag), uvec(dof_tag), vec(weight_tag), magnitude);
}

void new_nodeline(unique_ptr<Constraint>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_nodeline() needs a valid tag.\n");
        return;
    }

    uvec node_tag(3);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_nodeline() needs a valid node tag.\n");
            return;
        }

    return_obj = make_unique<NodeLine>(tag, 0, 0, std::move(node_tag));
}

void new_nodefacet(unique_ptr<Constraint>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_nodefacet() needs a valid tag.\n");
        return;
    }

    uvec node_tag(4);
    for(auto& I : node_tag)
        if(!get_input(command, I)) {
            suanpan_error("new_nodefacet() needs a valid node tag.\n");
            return;
        }

    return_obj = make_unique<NodeFacet>(tag, 0, 0, std::move(node_tag));
}

void new_particlecollision2d(unique_ptr<Constraint>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_particlecollision2d() needs a valid tag.\n");
        return;
    }

    auto space = 1.;
    if(!command.eof() && !get_input(command, space)) {
        suanpan_error("new_particlecollision2d() needs a valid spacing.\n");
        return;
    }

    auto alpha = 1.;
    if(!command.eof() && !get_input(command, alpha)) {
        suanpan_error("new_particlecollision2d() needs a valid multiplier.\n");
        return;
    }

    return_obj = make_unique<ParticleCollision2D>(tag, 0, space, alpha);
}

void new_particlecollision3d(unique_ptr<Constraint>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_particlecollision3d() needs a valid tag.\n");
        return;
    }

    auto space = 1.;
    if(!command.eof() && !get_input(command, space)) {
        suanpan_error("new_particlecollision3d() needs a valid spacing.\n");
        return;
    }

    auto alpha = 1.;
    if(!command.eof() && !get_input(command, alpha)) {
        suanpan_error("new_particlecollision3d() needs a valid multiplier.\n");
        return;
    }

    return_obj = make_unique<ParticleCollision3D>(tag, 0, space, alpha);
}

void new_rigidwall(unique_ptr<Constraint>& return_obj, istringstream& command, const bool finite, const bool penalty) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_rigidwall() needs a valid tag.\n");
        return;
    }

    vec origin(3), norm(3), edge(3);
    for(auto& I : origin) if(!get_input(command, I)) return;
    for(auto& I : norm) if(!get_input(command, I)) return;

    if(finite) for(auto& I : edge) if(!get_input(command, I)) return;

    auto alpha = 1.;
    if(!command.eof() && !get_input(command, alpha)) {
        suanpan_error("new_rigidwall() needs a valid multiplier.\n");
        return;
    }

    return_obj = finite ? penalty ? make_unique<RigidWallPenalty>(tag, 0, 0, std::move(origin), std::move(norm), std::move(edge), alpha) : make_unique<RigidWallMultiplier>(tag, 0, 0, std::move(origin), std::move(norm), std::move(edge), alpha) : penalty ? make_unique<RigidWallPenalty>(tag, 0, 0, std::move(origin), std::move(norm), alpha) : make_unique<RigidWallMultiplier>(tag, 0, 0, std::move(origin), std::move(norm), alpha);
}

int create_new_criterion(const shared_ptr<DomainBase>& domain, istringstream& command) {
    const auto& step_tag = domain->get_current_step_tag();
    if(0 == step_tag) {
        suanpan_error("create_new_criterion() needs a valid step.\n");
        return SUANPAN_SUCCESS;
    }

    string criterion_type;
    if(!get_input(command, criterion_type)) {
        suanpan_error("create_new_criterion() need a criterion type.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("create_new_criterion() requires a tag.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(criterion_type.substr(0, 5), "Logic")) {
        unsigned tag_a, tag_b;
        if(!get_input(command, tag_a) || !get_input(command, tag_b)) {
            suanpan_error("create_new_criterion() requires a valid tag.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(criterion_type, "LogicCriterionAND")) domain->insert(make_shared<LogicCriterionAND>(tag, step_tag, tag_a, tag_b));
        else if(is_equal(criterion_type, "LogicCriterionOR")) domain->insert(make_shared<LogicCriterionOR>(tag, step_tag, tag_a, tag_b));

        return SUANPAN_SUCCESS;
    }

    if(is_equal(criterion_type, "StrainEnergyEvolution")) {
        unsigned incre_level, final_level;
        if(!get_input(command, incre_level)) {
            suanpan_error("create_new_criterion() requires a valid level.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, final_level)) {
            suanpan_error("create_new_criterion() requires a valid level.\n");
            return SUANPAN_SUCCESS;
        }

        auto weight = 1.;
        if(!command.eof() && !get_input(command, weight)) {
            suanpan_error("create_new_criterion() requires a valid weight of central element.\n");
            return SUANPAN_SUCCESS;
        }
        auto iteration = 2;
        if(!command.eof() && !get_input(command, iteration)) {
            suanpan_error("create_new_criterion() requires a valid number of iteration.\n");
            return SUANPAN_SUCCESS;
        }
        auto reactivation = 10;
        if(!command.eof() && !get_input(command, reactivation)) {
            suanpan_error("create_new_criterion() requires a valid number of reactivation ratio.\n");
            return SUANPAN_SUCCESS;
        }
        auto propagation = .5;
        if(!command.eof() && !get_input(command, propagation)) {
            suanpan_error("create_new_criterion() requires a valid propagation factor.\n");
            return SUANPAN_SUCCESS;
        }
        auto tolerance = 1E-5;
        if(!command.eof() && !get_input(command, tolerance)) {
            suanpan_error("create_new_criterion() requires a valid tolerance.\n");
            return SUANPAN_SUCCESS;
        }

        domain->insert(make_shared<StrainEnergyEvolution>(tag, step_tag, incre_level, final_level, weight, iteration, reactivation, propagation, tolerance));

        return SUANPAN_SUCCESS;
    }

    if(is_equal(criterion_type, "MaxHistory")) {
        string type;
        double limit;
        if(!get_input(command, type)) {
            suanpan_error("create_new_criterion() requires a valid type.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, limit)) {
            suanpan_error("create_new_criterion() requires a valid limit.\n");
            return SUANPAN_SUCCESS;
        }

        domain->insert(make_shared<MaxHistory>(tag, step_tag, to_list(type.c_str()), limit));

        return SUANPAN_SUCCESS;
    }

    unsigned node;
    if(!get_input(command, node)) {
        suanpan_error("create_new_criterion() requires a node.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned dof;
    if(!get_input(command, dof)) {
        suanpan_error("create_new_criterion() requires a dof.\n");
        return SUANPAN_SUCCESS;
    }

    double limit;
    if(!get_input(command, limit)) {
        suanpan_error("create_new_criterion() requires a limit.\n");
        return SUANPAN_SUCCESS;
    }

    if(is_equal(criterion_type, "MaxDisplacement")) domain->insert(make_shared<MaxDisplacement>(tag, step_tag, node, dof, limit));
    else if(is_equal(criterion_type, "MinDisplacement")) domain->insert(make_shared<MinDisplacement>(tag, step_tag, node, dof, limit));
    else if(is_equal(criterion_type, "MaxResistance")) domain->insert(make_shared<MaxResistance>(tag, step_tag, node, dof, limit));
    else if(is_equal(criterion_type, "MinResistance")) domain->insert(make_shared<MinResistance>(tag, step_tag, node, dof, limit));

    return SUANPAN_SUCCESS;
}

int create_new_constraint(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string constraint_id;
    if(!get_input(command, constraint_id)) {
        suanpan_error("create_new_constraint() needs constraint type.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Constraint> new_constraint = nullptr;

    if(is_equal(constraint_id, "Embed2D")) new_embed(new_constraint, command, 2);
    else if(is_equal(constraint_id, "Embed3D")) new_embed(new_constraint, command, 3);
    else if(is_equal(constraint_id, "FixedLength2D") || is_equal(constraint_id, "R2D2")) new_fixedlength(new_constraint, command, 2);
    else if(is_equal(constraint_id, "FixedLength3D") || is_equal(constraint_id, "R3D2")) new_fixedlength(new_constraint, command, 3);
    else if(is_equal(constraint_id, "MinimumGap2D") || is_equal(constraint_id, "MinGap2D")) new_minimumgap(new_constraint, command, 2);
    else if(is_equal(constraint_id, "MinimumGap3D") || is_equal(constraint_id, "MinGap3D")) new_minimumgap(new_constraint, command, 3);
    else if(is_equal(constraint_id, "MaximumGap2D") || is_equal(constraint_id, "MaxGap2D")) new_maximumgap(new_constraint, command, 2);
    else if(is_equal(constraint_id, "MaximumGap3D") || is_equal(constraint_id, "MaxGap3D")) new_maximumgap(new_constraint, command, 3);
    else if(is_equal(constraint_id, "Sleeve2D")) new_sleeve(new_constraint, command, 2);
    else if(is_equal(constraint_id, "Sleeve3D")) new_sleeve(new_constraint, command, 3);
    else if(is_equal(constraint_id, "MPC")) new_mpc(new_constraint, command);
    else if(is_equal(constraint_id, "NodeLine")) new_nodeline(new_constraint, command);
    else if(is_equal(constraint_id, "NodeFacet")) new_nodefacet(new_constraint, command);
    else if(is_equal(constraint_id, "ParticleCollision2D")) new_particlecollision2d(new_constraint, command);
    else if(is_equal(constraint_id, "ParticleCollision3D")) new_particlecollision3d(new_constraint, command);
    else if(is_equal(constraint_id, "RigidWallMultiplier")) new_rigidwall(new_constraint, command, false, false);
    else if(is_equal(constraint_id, "RigidWallPenalty")) new_rigidwall(new_constraint, command, false, true);
    else if(is_equal(constraint_id, "FiniteRigidWallMultiplier")) new_rigidwall(new_constraint, command, true, false);
    else if(is_equal(constraint_id, "FiniteRigidWallPenalty")) new_rigidwall(new_constraint, command, true, true);
    else if(is_equal(constraint_id, "Fix") || is_equal(constraint_id, "PenaltyBC")) new_bc(new_constraint, command, true, false);
    else if(is_equal(constraint_id, "GroupPenaltyBC")) new_bc(new_constraint, command, true, true);
    else if(is_equal(constraint_id, "Fix2") || is_equal(constraint_id, "MultiplierBC")) new_bc(new_constraint, command, false, false);
    else if(is_equal(constraint_id, "GroupMultiplierBC")) new_bc(new_constraint, command, false, true);
    else load::object(new_constraint, domain, constraint_id, command);

    if(new_constraint != nullptr) new_constraint->set_start_step(domain->get_current_step_tag());

    if(new_constraint == nullptr || !domain->insert(std::move(new_constraint))) suanpan_error("create_new_constraint() fails to create new constraint.\n");

    return SUANPAN_SUCCESS;
}
