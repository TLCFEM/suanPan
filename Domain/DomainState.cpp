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

#include "Domain.h"

#include <Constraint/Constraint.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Element/Element.h>
#include <Load/Load.h>
#include <Recorder/Recorder.h>

void Domain::update_current_resistance() const {
    factory->modify_trial_resistance().zeros();
    if(color_map.empty())
        for(const auto& I : element_pond.get()) factory->assemble_resistance(I->get_current_resistance(), I->get_dof_encoding());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_resistance(I->get_current_resistance(), I->get_dof_encoding());
            });
        });

    factory->commit_resistance();
}

void Domain::update_current_damping_force() const {
    factory->modify_trial_damping_force().zeros();
    if(color_map.empty())
        for(const auto& I : element_pond.get()) factory->assemble_damping_force(I->get_current_damping_force(), I->get_dof_encoding());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_damping_force(I->get_current_damping_force(), I->get_dof_encoding());
            });
        });
    factory->commit_damping_force();
}

void Domain::update_current_nonviscous_force() const {
    factory->modify_trial_nonviscous_force().zeros();
    if(color_map.empty())
        for(const auto& I : element_pond.get()) factory->assemble_nonviscous_force(real(sum(I->get_current_nonviscous_force(), 1)), I->get_dof_encoding());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_nonviscous_force(real(sum(I->get_current_nonviscous_force(), 1)), I->get_dof_encoding());
            });
        });
    factory->commit_nonviscous_force();
}

void Domain::update_current_inertial_force() const {
    factory->modify_trial_inertial_force().zeros();
    if(color_map.empty())
        for(const auto& I : element_pond.get()) factory->assemble_inertial_force(I->get_current_inertial_force(), I->get_dof_encoding());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_inertial_force(I->get_current_inertial_force(), I->get_dof_encoding());
            });
        });
    factory->commit_inertial_force();
}

void Domain::assemble_resistance() const {
    auto& trial_resistance = factory->modify_trial_resistance().zeros();
    if(color_map.empty())
        for(const auto& I : element_pond.get()) factory->assemble_resistance(I->get_trial_resistance(), I->get_dof_encoding());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_resistance(I->get_trial_resistance(), I->get_dof_encoding());
            });
        });

    suanpan::for_all(node_pond.get(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_resistance(trial_resistance(t_node->get_reordered_dof())); });

    // update to sync incre_resistance
    factory->update_trial_resistance(trial_resistance);
}

void Domain::assemble_damping_force() const {
    auto& trial_damping_force = factory->modify_trial_damping_force().zeros();
    if(color_map.empty())
        for(const auto& I : element_pond.get()) factory->assemble_damping_force(I->get_trial_damping_force(), I->get_dof_encoding());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_damping_force(I->get_trial_damping_force(), I->get_dof_encoding());
            });
        });

    suanpan::for_all(node_pond.get(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_damping_force(trial_damping_force(t_node->get_reordered_dof())); });

    // update to sync incre_damping_force
    factory->update_trial_damping_force(trial_damping_force);
}

void Domain::assemble_nonviscous_force() const {
    auto& trial_nonviscous_force = factory->modify_trial_nonviscous_force().zeros();
    if(color_map.empty())
        for(const auto& I : element_pond.get()) factory->assemble_nonviscous_force(real(sum(I->get_trial_nonviscous_force(), 1)), I->get_dof_encoding());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_nonviscous_force(real(sum(I->get_trial_nonviscous_force(), 1)), I->get_dof_encoding());
            });
        });

    suanpan::for_all(node_pond.get(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_nonviscous_force(trial_nonviscous_force(t_node->get_reordered_dof())); });

    // update to sync incre_nonviscous_force
    factory->update_trial_nonviscous_force(trial_nonviscous_force);
}

void Domain::assemble_inertial_force() const {
    auto& trial_inertial_force = factory->modify_trial_inertial_force().zeros();
    if(color_map.empty())
        for(const auto& I : element_pond.get()) factory->assemble_inertial_force(I->get_trial_inertial_force(), I->get_dof_encoding());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_inertial_force(I->get_trial_inertial_force(), I->get_dof_encoding());
            });
        });

    suanpan::for_all(node_pond.get(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_inertial_force(trial_inertial_force(t_node->get_reordered_dof())); });

    // update to sync incre_inertial_force
    factory->update_trial_inertial_force(trial_inertial_force);
}

void Domain::assemble_initial_mass() const {
    factory->clear_mass();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_mass(I->get_initial_mass(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_mass(I->get_initial_mass(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_mass()->csc_condense();
}

void Domain::assemble_current_mass() const {
    factory->clear_mass();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_mass(I->get_current_mass(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_mass(I->get_current_mass(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_mass()->csc_condense();
}

void Domain::assemble_trial_mass() const {
    factory->clear_mass();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_mass(I->get_trial_mass(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_mass(I->get_trial_mass(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_mass()->csc_condense();
}

void Domain::assemble_initial_damping() const {
    factory->clear_damping();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_damping(I->get_initial_viscous(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_damping(I->get_initial_viscous(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_damping()->csc_condense();
}

void Domain::assemble_current_damping() const {
    factory->clear_damping();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_damping(I->get_current_viscous(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_damping(I->get_current_viscous(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_damping()->csc_condense();
}

void Domain::assemble_trial_damping() const {
    factory->clear_damping();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_damping(I->get_trial_viscous(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_damping(I->get_trial_viscous(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_damping()->csc_condense();
}

void Domain::assemble_initial_nonviscous() const {
    if(!factory->is_nonviscous()) return;
    factory->clear_nonviscous();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_nonviscous(I->get_initial_nonviscous(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_nonviscous(I->get_initial_nonviscous(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_nonviscous()->csc_condense();
}

void Domain::assemble_current_nonviscous() const {
    if(!factory->is_nonviscous()) return;
    factory->clear_nonviscous();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_nonviscous(I->get_current_nonviscous(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_nonviscous(I->get_current_nonviscous(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_nonviscous()->csc_condense();
}

void Domain::assemble_trial_nonviscous() const {
    if(!factory->is_nonviscous()) return;
    factory->clear_nonviscous();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_nonviscous(I->get_trial_nonviscous(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_nonviscous(I->get_trial_nonviscous(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_nonviscous()->csc_condense();
}

void Domain::assemble_initial_stiffness() const {
    factory->clear_stiffness();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_stiffness(I->get_initial_stiffness(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_stiffness(I->get_initial_stiffness(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_stiffness()->csc_condense();
}

void Domain::assemble_current_stiffness() const {
    factory->clear_stiffness();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_stiffness(I->get_current_stiffness(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_stiffness(I->get_current_stiffness(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_stiffness()->csc_condense();
}

void Domain::assemble_trial_stiffness() const {
    factory->clear_stiffness();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_stiffness(I->get_trial_stiffness(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_stiffness(I->get_trial_stiffness(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_stiffness()->csc_condense();
}

void Domain::assemble_initial_geometry() const {
    if(!factory->is_nlgeom()) return;

    factory->clear_geometry();
    if(color_map.empty() || is_sparse()) {
        for(const auto& I : element_pond.get())
            if(I->is_nlgeom()) factory->assemble_geometry(I->get_initial_geometry(), I->get_dof_encoding(), I->get_dof_mapping());
    }
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_geometry(I->get_initial_geometry(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_geometry()->csc_condense();
}

void Domain::assemble_current_geometry() const {
    if(!factory->is_nlgeom()) return;

    factory->clear_geometry();
    if(color_map.empty() || is_sparse()) {
        for(const auto& I : element_pond.get())
            if(I->is_nlgeom()) factory->assemble_geometry(I->get_current_geometry(), I->get_dof_encoding(), I->get_dof_mapping());
    }
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_geometry(I->get_current_geometry(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_geometry()->csc_condense();
}

void Domain::assemble_trial_geometry() const {
    if(!factory->is_nlgeom()) return;

    factory->clear_geometry();
    if(color_map.empty() || is_sparse()) {
        for(const auto& I : element_pond.get())
            if(I->is_nlgeom()) factory->assemble_geometry(I->get_trial_geometry(), I->get_dof_encoding(), I->get_dof_mapping());
    }
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_geometry(I->get_trial_geometry(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_geometry()->csc_condense();
}

void Domain::assemble_mass_container() const {
    factory->clear_mass();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_mass(I->get_mass_container(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_mass(I->get_mass_container(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_mass()->csc_condense();
}

void Domain::assemble_stiffness_container() const {
    factory->clear_stiffness();
    if(color_map.empty() || is_sparse())
        for(const auto& I : element_pond.get()) factory->assemble_stiffness(I->get_stiffness_container(), I->get_dof_encoding(), I->get_dof_mapping());
    else
        std::ranges::for_each(color_map, [&](const std::vector<unsigned>& color) {
            suanpan::for_all(color, [&](const unsigned tag) {
                const auto& I = get_element(tag);
                factory->assemble_stiffness(I->get_stiffness_container(), I->get_dof_encoding(), I->get_dof_mapping());
            });
        });

    factory->get_stiffness()->csc_condense();
}

int Domain::update_trial_status() const {
    auto& trial_displacement = factory->get_trial_displacement();
    auto& trial_velocity = factory->get_trial_velocity();
    auto& trial_acceleration = factory->get_trial_acceleration();

    if(AnalysisType::DYNAMICS == factory->get_analysis_type()) suanpan::for_all(node_pond.get(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_status(trial_displacement, trial_velocity, trial_acceleration); });
    else suanpan::for_all(node_pond.get(), [&](const shared_ptr<Node>& t_node) { t_node->update_trial_status(trial_displacement); });

    std::atomic_int code = 0;
    suanpan::for_all(element_pond.get(), [&](const shared_ptr<Element>& t_element) {
        if(t_element->is_local) code += t_element->update_status();
    });

    return allreduce(code.load());
}

int Domain::update_incre_status() const {
    auto& incre_displacement = factory->get_incre_displacement();
    auto& incre_velocity = factory->get_incre_velocity();
    auto& incre_acceleration = factory->get_incre_acceleration();

    if(AnalysisType::DYNAMICS == factory->get_analysis_type()) suanpan::for_all(node_pond.get(), [&](const shared_ptr<Node>& t_node) { t_node->update_incre_status(incre_displacement, incre_velocity, incre_acceleration); });
    else suanpan::for_all(node_pond.get(), [&](const shared_ptr<Node>& t_node) { t_node->update_incre_status(incre_displacement); });

    std::atomic_int code = 0;
    suanpan::for_all(element_pond.get(), [&](const shared_ptr<Element>& t_element) {
        if(t_element->is_local) code += t_element->update_status();
    });

    return allreduce(code.load());
}

int Domain::update_current_status() const {
    const auto analysis_type = factory->get_analysis_type();

    if(analysis_type != AnalysisType::STATICS && analysis_type != AnalysisType::DYNAMICS && analysis_type != AnalysisType::BUCKLE) return SUANPAN_SUCCESS;

    // collect initial nodal quantities into global storage
    if(analysis_type == AnalysisType::DYNAMICS) {
        vec c_g_dsp(factory->get_size(), fill::zeros);
        vec c_g_vel(factory->get_size(), fill::zeros);
        vec c_g_acc(factory->get_size(), fill::zeros);

        for(const auto& I : node_pond.get()) {
            auto& t_dof = I->get_reordered_dof();
            c_g_dsp(t_dof) = I->get_current_displacement();
            c_g_vel(t_dof) = I->get_current_velocity();
            c_g_acc(t_dof) = I->get_current_acceleration();
        }
        factory->update_current_displacement(c_g_dsp);
        factory->update_current_velocity(c_g_vel);
        factory->update_current_acceleration(c_g_acc);
    }
    else {
        vec c_g_dsp(factory->get_size(), fill::zeros);

        for(const auto& I : node_pond.get()) {
            auto& t_dof = I->get_reordered_dof();
            c_g_dsp(t_dof) = I->get_current_displacement();
        }
        factory->update_current_displacement(c_g_dsp);
    }

    std::atomic_int code = 0;
    suanpan::for_all(element_pond.get(), [&](const shared_ptr<Element>& t_element) {
        if(t_element->is_local) code += t_element->update_status();
    });

    if(SUANPAN_SUCCESS != allreduce(code.load())) {
        suanpan_error("Initial conditions cause significant non-linearity, use an analysis step instead.\n");
        return SUANPAN_FAIL;
    }

    // commit element status
    suanpan::for_all(element_pond.get(), [](const shared_ptr<Element>& t_element) {
        t_element->Element::commit_status();
        t_element->commit_status();
    });

    update_current_resistance();
    if(analysis_type == AnalysisType::DYNAMICS) {
        update_current_damping_force();
        update_current_inertial_force();
    }

    return SUANPAN_SUCCESS;
}

void Domain::stage_status() {
    suanpan::for_all(constraint_pond.get(), [&](const shared_ptr<Constraint>& t_constraint) { t_constraint->stage(shared_from_this()); });
}

void Domain::commit_status() const {
    factory->commit_status();

    // element comes before node to account for strain energy update
    suanpan::for_all(element_pond.get(), [](const shared_ptr<Element>& t_element) {
        t_element->Element::commit_status();
        t_element->commit_status();
    });
    suanpan::for_all(node_pond.get(), [](const shared_ptr<Node>& t_node) { t_node->commit_status(); });
    suanpan::for_all(load_pond.get(), [](const shared_ptr<Load>& t_load) { t_load->commit_status(); });
    suanpan::for_all(constraint_pond.get(), [](const shared_ptr<Constraint>& t_constraint) { t_constraint->commit_status(); });
}

void Domain::clear_status() {
    // enable_all();

    updated = false;

    // current tags are treated as state variables
    // current_step_tag = 0;
    // current_converger_tag = 0;
    // current_integrator_tag = 0;
    // current_solver_tag = 0;

    factory->clear_status();

    // element comes before node to account for strain energy update
    suanpan::for_all(element_pond.get(), [](const shared_ptr<Element>& t_element) {
        t_element->Element::clear_status();
        t_element->clear_status();
    });
    suanpan::for_all(node_pond.get(), [](const shared_ptr<Node>& t_node) { t_node->clear_status(); });
    suanpan::for_all(load_pond.get(), [](const shared_ptr<Load>& t_load) { t_load->clear_status(); });
    suanpan::for_all(constraint_pond.get(), [](const shared_ptr<Constraint>& t_constraint) { t_constraint->clear_status(); });

    suanpan::for_all(recorder_pond.get(), [](const shared_ptr<Recorder>& t_recorder) { t_recorder->clear_status(); });
}

void Domain::reset_status() const {
    factory->reset_status();

    // element comes before node to account for strain energy update
    suanpan::for_all(element_pond.get(), [](const shared_ptr<Element>& t_element) {
        t_element->Element::reset_status();
        t_element->reset_status();
    });
    suanpan::for_all(node_pond.get(), [](const shared_ptr<Node>& t_node) { t_node->reset_status(); });
    suanpan::for_all(load_pond.get(), [](const shared_ptr<Load>& t_load) { t_load->reset_status(); });
    suanpan::for_all(constraint_pond.get(), [](const shared_ptr<Constraint>& t_constraint) { t_constraint->reset_status(); });
}
