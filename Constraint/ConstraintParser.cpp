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

// ReSharper disable IdentifierTypo
#include "ConstraintParser.h"

#include <Constraint/Constraint>
#include <Domain/ExternalModule.h>
#include <Recorder/OutputType.h>

namespace {
    void new_bc(unique_ptr<Constraint>& return_obj, std::istringstream& command, const bool penalty, const bool group) {
        unsigned bc_id;
        if(!get_input(command, bc_id)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        std::string dof_token;
        if(!get_input(command, dof_token)) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        auto dof_pool = parse_dof(dof_token);
        if(dof_pool.empty()) {
            suanpan_error("A valid dof identifier is required.\n");
            return;
        }

        const auto obj_tag = get_remaining<uword>(command);

        if(penalty) {
            if(group) return_obj = std::make_unique<GroupPenaltyBC>(bc_id, uvec(obj_tag), std::move(dof_pool));
            else return_obj = std::make_unique<PenaltyBC>(bc_id, uvec(obj_tag), std::move(dof_pool));
        }
        else {
            if(group) return_obj = std::make_unique<GroupMultiplierBC>(bc_id, uvec(obj_tag), std::move(dof_pool));
            else return_obj = std::make_unique<MultiplierBC>(bc_id, uvec(obj_tag), std::move(dof_pool));
        }
    }

    template<unsigned dof> void new_fixedlength(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        uword node_i, node_j;
        if(!get_input(command, node_i) || !get_input(command, node_j)) {
            suanpan_error("Two valid nodes are required.\n");
            return;
        }

        return_obj = std::make_unique<FixedLength<dof>>(tag, uvec{node_i, node_j});
    }

    template<unsigned dof> void new_maxforce(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        uword node_i, node_j;
        if(!get_input(command, node_i) || !get_input(command, node_j)) {
            suanpan_error("Two valid nodes are required.\n");
            return;
        }

        double max_force;
        if(!get_input(command, max_force)) {
            suanpan_error("A valid maximum force is required.\n");
            return;
        }

        return_obj = std::make_unique<MaxForce<dof>>(tag, max_force, uvec{node_i, node_j});
    }

    template<unsigned dof> void new_minimumgap(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        uword node_i, node_j;
        if(!get_input(command, node_i) || !get_input(command, node_j)) {
            suanpan_error("Two valid nodes are required.\n");
            return;
        }

        double gap;
        if(!get_input(command, gap)) {
            suanpan_error("A valid minimum gap is required.\n");
            return;
        }

        return_obj = std::make_unique<MinimumGap<dof>>(tag, gap, uvec{node_i, node_j});
    }

    template<unsigned dof> void new_maximumgap(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        uword node_i, node_j;
        if(!get_input(command, node_i) || !get_input(command, node_j)) {
            suanpan_error("Two valid nodes are required.\n");
            return;
        }

        double gap;
        if(!get_input(command, gap)) {
            suanpan_error("A valid minimum gap is required.\n");
            return;
        }

        return_obj = std::make_unique<MaximumGap<dof>>(tag, gap, uvec{node_i, node_j});
    }

    template<unsigned dof> void new_sleeve(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        uword node_i, node_j;
        if(!get_input(command, node_i) || !get_input(command, node_j)) {
            suanpan_error("Two valid nodes are required.\n");
            return;
        }

        double min_gap, max_gap;
        if(!get_input(command, min_gap)) {
            suanpan_error("A valid minimum gap is required.\n");
            return;
        }
        if(!get_input(command, max_gap)) {
            suanpan_error("A valid maximum gap is required.\n");
            return;
        }

        return_obj = std::make_unique<Sleeve<dof>>(tag, min_gap, max_gap, uvec{node_i, node_j});
    }

    template<unsigned dof> void new_embed(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned element_tag;
        if(!get_input(command, element_tag)) {
            suanpan_error("A valid element tag is required.\n");
            return;
        }

        unsigned node_tag;
        if(!get_input(command, node_tag)) {
            suanpan_error("A valid node tag is required.\n");
            return;
        }

        return_obj = std::make_unique<Embed<dof>>(tag, element_tag, node_tag);
    }

    void new_mpc(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        unsigned amplitude;
        if(!get_input(command, amplitude)) {
            suanpan_error("A valid amplitude tag is required.\n");
            return;
        }

        double magnitude;
        if(!get_input(command, magnitude)) {
            suanpan_error("A valid magnitude is required.\n");
            return;
        }

        MPC::Pack pack;
        for(const auto& [node, token, weight] : get_remaining_as_tuple<uword, std::string, double>(command)) {
            const auto dof = parse_dof(token);
            if(1 != dof.size()) {
                suanpan_error("A valid dof is required.\n");
                return;
            }
            pack.emplace_back(node, dof.front(), weight);
        }

        return_obj = std::make_unique<MPC>(tag, amplitude, magnitude, std::move(pack));
    }

    void new_nodeline(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        uvec node_tag(3);
        if(!get_input(command, node_tag)) {
            suanpan_error("A valid node tag is required.\n");
            return;
        }

        return_obj = std::make_unique<NodeLine>(tag, 0, std::move(node_tag));
    }

    void new_nodefacet(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        uvec node_tag(4);
        if(!get_input(command, node_tag)) {
            suanpan_error("A valid node tag is required.\n");
            return;
        }

        return_obj = std::make_unique<NodeFacet>(tag, 0, std::move(node_tag));
    }

    template<unsigned DIM> void new_moleculardynamics(unique_ptr<Constraint>& return_obj, std::istringstream& command) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        if constexpr(2 == DIM)
            return_obj = std::make_unique<MolecularDynamics2D>(tag, get_remaining<uword>(command));
        else
            return_obj = std::make_unique<MolecularDynamics3D>(tag, get_remaining<uword>(command));
    }

    void new_particlecollision(unique_ptr<Constraint>& return_obj, std::istringstream& command, const unsigned dim) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto space = 1.;
        if(!command.eof() && !get_input(command, space)) {
            suanpan_error("A valid spacing is required.\n");
            return;
        }

        auto alpha = 1.;
        if(!command.eof() && !get_input(command, alpha)) {
            suanpan_error("A valid multiplier is required.\n");
            return;
        }

        2 == dim ? return_obj = std::make_unique<ParticleCollision2D>(tag, space, alpha) : return_obj = std::make_unique<ParticleCollision3D>(tag, space, alpha);
    }

    void new_ljpotential(unique_ptr<Constraint>& return_obj, std::istringstream& command, const unsigned) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto space = 1.;
        if(!command.eof() && !get_input(command, space)) {
            suanpan_error("A valid spacing is required.\n");
            return;
        }

        auto alpha = 1.;
        if(!command.eof() && !get_input(command, alpha)) {
            suanpan_error("A valid multiplier is required.\n");
            return;
        }

        return_obj = std::make_unique<LJPotential2D>(tag, space, alpha);
    }

    void new_linearspring(unique_ptr<Constraint>& return_obj, std::istringstream& command, const unsigned) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto space = 1.;
        if(!command.eof() && !get_input(command, space)) {
            suanpan_error("A valid spacing is required.\n");
            return;
        }

        auto alpha = 1.;
        if(!command.eof() && !get_input(command, alpha)) {
            suanpan_error("A valid multiplier is required.\n");
            return;
        }

        return_obj = std::make_unique<LinearSpring2D>(tag, space, alpha);
    }

    void new_rigidwall(unique_ptr<Constraint>& return_obj, std::istringstream& command, const bool finite, const bool penalty) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto p = get_remaining<double>(command);
        if(2 == p.size() || 4 == p.size() || 6 == p.size() || 9 == p.size()) p.emplace_back(1E4);

        switch(p.size()) {
        case 3:
            // 1D origin norm multiplier
            if(penalty) return_obj = std::make_unique<RigidWallPenalty1D>(tag, 0, p[0], p[1], p[2]);
            else return_obj = std::make_unique<RigidWallMultiplier1D>(tag, 0, p[0], p[1], p[2]);
            break;
        case 5:
            if(finite) {
                // 2D origin edge multiplier
                if(penalty) return_obj = std::make_unique<RigidWallPenalty2D>(tag, 0, vec2{p[0], p[1]}, vec3{p[2], p[3], 0.}, p[4]);
                else return_obj = std::make_unique<RigidWallMultiplier2D>(tag, 0, vec2{p[0], p[1]}, vec3{p[2], p[3], 0.}, p[4]);
            }
            else {
                // 2D origin norm multiplier
                if(penalty) return_obj = std::make_unique<RigidWallPenalty2D>(tag, 0, vec2{p[0], p[1]}, vec2{p[2], p[3]}, p[4]);
                else return_obj = std::make_unique<RigidWallMultiplier2D>(tag, 0, vec2{p[0], p[1]}, vec2{p[2], p[3]}, p[4]);
            }
            break;
        case 7:
            // 3D origin norm multiplier
            if(penalty) return_obj = std::make_unique<RigidWallPenalty3D>(tag, 0, vec3{p[0], p[1], p[2]}, vec3{p[3], p[4], p[5]}, p[6]);
            else return_obj = std::make_unique<RigidWallMultiplier3D>(tag, 0, vec3{p[0], p[1], p[2]}, vec3{p[3], p[4], p[5]}, p[6]);
            break;
        case 10:
            // 3D origin edge edge multiplier
            if(penalty) return_obj = std::make_unique<RigidWallPenalty3D>(tag, 0, vec3{p[0], p[1], p[2]}, vec3{p[3], p[4], p[5]}, vec3{p[6], p[7], p[8]}, p[9]);
            else return_obj = std::make_unique<RigidWallMultiplier3D>(tag, 0, vec3{p[0], p[1], p[2]}, vec3{p[3], p[4], p[5]}, vec3{p[6], p[7], p[8]}, p[9]);
            break;
        default:
            suanpan_error("A valid number of parameters is required.\n");
        }
    }

    void new_restitutionwall(unique_ptr<Constraint>& return_obj, std::istringstream& command, const bool finite) {
        unsigned tag;
        if(!get_input(command, tag)) {
            suanpan_error("A valid tag is required.\n");
            return;
        }

        auto p = get_remaining<double>(command);
        if(3 == p.size() || 5 == p.size() || 7 == p.size() || 10 == p.size()) p.emplace_back(1E4);

        switch(p.size()) {
        case 4:
            // 1D origin norm restitution multiplier
            return_obj = std::make_unique<RestitutionWallPenalty1D>(tag, 0, p[0], p[1], p[2], p[3]);
            break;
        case 6:
            if(finite)
                // 2D origin edge restitution multiplier
                return_obj = std::make_unique<RestitutionWallPenalty2D>(tag, 0, vec2{p[0], p[1]}, vec3{p[2], p[3], 0.}, p[4], p[5]);
            else
                // 2D origin norm restitution multiplier
                return_obj = std::make_unique<RestitutionWallPenalty2D>(tag, 0, vec2{p[0], p[1]}, vec2{p[2], p[3]}, p[4], p[5]);
            break;
        case 8:
            // 3D origin norm restitution multiplier
            return_obj = std::make_unique<RestitutionWallPenalty3D>(tag, 0, vec3{p[0], p[1], p[2]}, vec3{p[3], p[4], p[5]}, p[6], p[7]);
            break;
        case 11:
            // 3D origin edge edge restitution multiplier
            return_obj = std::make_unique<RestitutionWallPenalty3D>(tag, 0, vec3{p[0], p[1], p[2]}, vec3{p[3], p[4], p[5]}, vec3{p[6], p[7], p[8]}, p[9], p[10]);
            break;
        default:
            suanpan_error("A valid number of parameters is required.\n");
        }
    }
} // namespace

int create_new_criterion(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    const auto step_tag = domain->get_current_step_tag();
    if(0u == step_tag) {
        suanpan_error("A valid step is required.\n");
        return SUANPAN_SUCCESS;
    }

    std::string criterion_type;
    if(!get_input(command, criterion_type)) {
        suanpan_error("A valid criterion type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    auto flag = true;

    if(is_equal(criterion_type.substr(0, 5), "Logic")) {
        unsigned tag_a, tag_b;
        if(!get_input(command, tag_a, tag_b)) {
            suanpan_error("A valid tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(criterion_type, "LogicCriterionAND")) flag = domain->insert(std::make_shared<LogicCriterionAND>(tag, step_tag, tag_a, tag_b));
        else if(is_equal(criterion_type, "LogicCriterionOR")) flag = domain->insert(std::make_shared<LogicCriterionOR>(tag, step_tag, tag_a, tag_b));
    }
    else if(is_equal(criterion_type, "StrainEnergyEvolution")) {
        unsigned incre_level, final_level;
        if(!get_input(command, incre_level, final_level)) {
            suanpan_error("A valid level is required.\n");
            return SUANPAN_SUCCESS;
        }

        auto weight = 1.;
        if(!command.eof() && !get_input(command, weight)) {
            suanpan_error("A valid weight of central element is required.\n");
            return SUANPAN_SUCCESS;
        }
        auto iteration = 2;
        if(!command.eof() && !get_input(command, iteration)) {
            suanpan_error("A valid number of iteration is required.\n");
            return SUANPAN_SUCCESS;
        }
        auto reactivation = 10;
        if(!command.eof() && !get_input(command, reactivation)) {
            suanpan_error("A valid number of reactivation ratio is required.\n");
            return SUANPAN_SUCCESS;
        }
        auto propagation = .5;
        if(!command.eof() && !get_input(command, propagation)) {
            suanpan_error("A valid propagation factor is required.\n");
            return SUANPAN_SUCCESS;
        }
        auto tolerance = 1E-5;
        if(!command.eof() && !get_input(command, tolerance)) {
            suanpan_error("A valid tolerance is required.\n");
            return SUANPAN_SUCCESS;
        }

        flag = domain->insert(std::make_shared<StrainEnergyEvolution>(tag, step_tag, incre_level, final_level, weight, iteration, reactivation, propagation, tolerance));
    }
    else if(is_equal(criterion_type, "MaxHistory")) {
        std::string type;
        double limit;
        if(!get_input(command, type)) {
            suanpan_error("A valid type is required.\n");
            return SUANPAN_SUCCESS;
        }
        if(!get_input(command, limit)) {
            suanpan_error("A valid limit is required.\n");
            return SUANPAN_SUCCESS;
        }

        flag = domain->insert(std::make_shared<MaxHistory>(tag, step_tag, to_token(type), limit));
    }
    else {
        unsigned node;
        if(!get_input(command, node)) {
            suanpan_error("A valid node tag is required.\n");
            return SUANPAN_SUCCESS;
        }

        unsigned dof;
        if(!get_input(command, dof)) {
            suanpan_error("A valid dof identifier is required.\n");
            return SUANPAN_SUCCESS;
        }

        double limit;
        if(!get_input(command, limit)) {
            suanpan_error("A valid limit is required.\n");
            return SUANPAN_SUCCESS;
        }

        if(is_equal(criterion_type, "MaxDisplacement")) flag = domain->insert(std::make_shared<MaxDisplacement>(tag, step_tag, node, dof, limit));
        else if(is_equal(criterion_type, "MinDisplacement")) flag = domain->insert(std::make_shared<MinDisplacement>(tag, step_tag, node, dof, limit));
        else if(is_equal(criterion_type, "MaxResistance")) flag = domain->insert(std::make_shared<MaxResistance>(tag, step_tag, node, dof, limit));
        else if(is_equal(criterion_type, "MinResistance")) flag = domain->insert(std::make_shared<MinResistance>(tag, step_tag, node, dof, limit));
    }

    if(!flag) suanpan_error("Fail to create new criterion via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_interaction(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string criterion_type;
    if(!get_input(command, criterion_type)) {
        suanpan_error("A valid criterion type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return SUANPAN_SUCCESS;
    }

    auto flag = true;

    if(is_equal(criterion_type, "Hertzian")) flag = domain->insert(std::make_shared<Hertzian>(tag));
    else if(is_equal(criterion_type, "FixedParticle")) {
        double multiplier;
        if(!get_input(command, multiplier)) {
            suanpan_error("A valid multiplier is required.\n");
            return SUANPAN_SUCCESS;
        }

        flag = domain->insert(std::make_shared<FixedParticle>(tag, multiplier, get_remaining_as_set<unsigned>(command)));
    }

    if(!flag) suanpan_error("Fail to create new interaction via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}

int create_new_constraint(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    std::string constraint_id;
    if(!get_input(command, constraint_id)) {
        suanpan_error("A valid constraint type is required.\n");
        return SUANPAN_SUCCESS;
    }

    unique_ptr<Constraint> new_constraint = nullptr;

    if(is_equal(constraint_id, "Embed2D")) new_embed<2u>(new_constraint, command);
    else if(is_equal(constraint_id, "Embed3D")) new_embed<3u>(new_constraint, command);
    else if(is_equal_any(constraint_id, "FiniteRestitutionWall", "FiniteRestitutionWallPenalty")) new_restitutionwall(new_constraint, command, true);
    else if(is_equal_any(constraint_id, "FiniteRigidWall", "FiniteRigidWallPenalty")) new_rigidwall(new_constraint, command, true, true);
    else if(is_equal(constraint_id, "FiniteRigidWallMultiplier")) new_rigidwall(new_constraint, command, true, false);
    else if(is_equal_any(constraint_id, "Fix", "PenaltyBC")) new_bc(new_constraint, command, true, false);
    else if(is_equal_any(constraint_id, "Fix2", "MultiplierBC")) new_bc(new_constraint, command, false, false);
    else if(is_equal_any(constraint_id, "FixedLength2D", "R2D2")) new_fixedlength<2u>(new_constraint, command);
    else if(is_equal_any(constraint_id, "FixedLength3D", "R3D2")) new_fixedlength<3u>(new_constraint, command);
    else if(is_equal(constraint_id, "GroupMultiplierBC")) new_bc(new_constraint, command, false, true);
    else if(is_equal(constraint_id, "GroupPenaltyBC")) new_bc(new_constraint, command, true, true);
    else if(is_equal(constraint_id, "LinearSpring2D")) new_linearspring(new_constraint, command, 2);
    else if(is_equal(constraint_id, "LJPotential2D")) new_ljpotential(new_constraint, command, 2);
    else if(is_equal_any(constraint_id, "MaximumForce2D", "MaxForce2D")) new_maxforce<2u>(new_constraint, command);
    else if(is_equal_any(constraint_id, "MaximumForce3D", "MaxForce3D")) new_maxforce<3u>(new_constraint, command);
    else if(is_equal_any(constraint_id, "MaximumGap2D", "MaxGap2D")) new_maximumgap<2u>(new_constraint, command);
    else if(is_equal_any(constraint_id, "MaximumGap3D", "MaxGap3D")) new_maximumgap<3u>(new_constraint, command);
    else if(is_equal_any(constraint_id, "MinimumGap2D", "MinGap2D")) new_minimumgap<2u>(new_constraint, command);
    else if(is_equal_any(constraint_id, "MinimumGap3D", "MinGap3D")) new_minimumgap<3u>(new_constraint, command);
    else if(is_equal(constraint_id, "MPC")) new_mpc(new_constraint, command);
    else if(is_equal(constraint_id, "NodeFacet")) new_nodefacet(new_constraint, command);
    else if(is_equal(constraint_id, "NodeLine")) new_nodeline(new_constraint, command);
    else if(is_equal(constraint_id, "MolecularDynamics2D")) new_moleculardynamics<2u>(new_constraint, command);
    else if(is_equal(constraint_id, "MolecularDynamics3D")) new_moleculardynamics<3u>(new_constraint, command);
    else if(is_equal(constraint_id, "ParticleCollision2D")) new_particlecollision(new_constraint, command, 2);
    else if(is_equal(constraint_id, "ParticleCollision3D")) new_particlecollision(new_constraint, command, 3);
    else if(is_equal_any(constraint_id, "RestitutionWall", "RestitutionWallPenalty")) new_restitutionwall(new_constraint, command, false);
    else if(is_equal_any(constraint_id, "RigidWall", "RigidWallPenalty")) new_rigidwall(new_constraint, command, false, true);
    else if(is_equal(constraint_id, "RigidWallMultiplier")) new_rigidwall(new_constraint, command, false, false);
    else if(is_equal(constraint_id, "Sleeve2D")) new_sleeve<2u>(new_constraint, command);
    else if(is_equal(constraint_id, "Sleeve3D")) new_sleeve<3u>(new_constraint, command);
    else external_module::object(new_constraint, domain, constraint_id, command);

    if(new_constraint != nullptr) new_constraint->set_start_step(domain->get_current_step_tag());

    if(new_constraint == nullptr || !domain->insert(std::move(new_constraint)))
        suanpan_error("Fail to create new constraint via \"{}\".\n", command.str());

    return SUANPAN_SUCCESS;
}
