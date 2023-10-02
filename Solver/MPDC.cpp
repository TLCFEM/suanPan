#include "MPDC.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

MPDC::MPDC(const unsigned T)
    : Solver(T) {}

int MPDC::analyze() {
    auto& C = get_converger();
    auto& G = get_integrator();
    const auto D = G->get_domain();
    auto& W = D->get_factory();

    suanpan_highlight(">> Current Analysis Time: {:.5f}.\n", W->get_trial_time());

    const auto max_iteration = C->get_max_iteration();

    vec samurai;

    // get column index for each nonzero dof
    // uvec load_ref_idx = find(load_ref);
    // for(auto I = 0; I < load_ref_idx.n_elem; ++I) load_ref_idx(I) -= I * load_ref.n_rows;

    const auto idx = to_uvec(W->get_reference_dof());

    if(idx.empty()) {
        suanpan_error("Displacement controlled algorithm is activated but no valid displacement load is applied.\n");
        return SUANPAN_FAIL;
    }

    mat disp_a;

    wall_clock t_clock;

    // iteration counter
    unsigned counter = 0;

    while(true) {
        // update for nodes and elements
        t_clock.tic();
        if(SUANPAN_SUCCESS != G->update_trial_status()) return SUANPAN_FAIL;
        // process modifiers
        if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;
        D->update<Statistics::UpdateStatus>(t_clock.toc());
        // assemble resistance
        t_clock.tic();
        G->assemble_resistance();
        D->update<Statistics::AssembleVector>(t_clock.toc());

        if(constant_matrix()) {
            // some loads may have resistance
            t_clock.tic();
            if(SUANPAN_SUCCESS != G->process_load_resistance()) return SUANPAN_FAIL;
            // some constraints may have resistance
            if(SUANPAN_SUCCESS != G->process_constraint_resistance()) return SUANPAN_FAIL;
            D->update<Statistics::ProcessConstraint>(t_clock.toc());
        }
        else {
            // assemble stiffness
            t_clock.tic();
            G->assemble_matrix();
            D->update<Statistics::AssembleMatrix>(t_clock.toc());
            // process loads
            t_clock.tic();
            if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
            // process constraints
            if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;
            D->update<Statistics::ProcessConstraint>(t_clock.toc());
            // indicate the global matrix has been assembled
            G->set_matrix_assembled_switch(true);
        }

        t_clock.tic();

        // solve ninja
        if(SUANPAN_SUCCESS != G->solve(samurai, G->get_displacement_residual())) return SUANPAN_FAIL;
        // solve reference displacement
        if(SUANPAN_SUCCESS != G->solve(disp_a, G->get_reference_load())) return SUANPAN_FAIL;

        if(const auto n_size = W->get_size(); 0 != W->get_mpc()) {
            mat right, kernel;
            auto& border = W->get_auxiliary_stiffness();
            if(SUANPAN_SUCCESS != G->solve(right, border)) return SUANPAN_FAIL;
            auto& aux_lambda = W->modify_auxiliary_lambda();
            if(!solve(aux_lambda, kernel = border.t() * right.head_rows(n_size), border.t() * samurai.head(n_size) - G->get_auxiliary_residual())) return SUANPAN_FAIL;
            samurai -= right * aux_lambda;
            disp_a -= right * solve(kernel, border.t() * disp_a.head_rows(n_size));
        }

        const vec incre_lambda = solve(mat(disp_a.rows(idx)), W->get_trial_settlement()(idx) - G->get_trial_displacement()(idx) - samurai.rows(idx));

        samurai += disp_a * incre_lambda;

        D->update<Statistics::SolveSystem>(t_clock.toc());

        // avoid machine error accumulation
        G->erase_machine_error(samurai);

        // exit if converged
        if(C->is_converged(counter)) return G->sync_status(true);
        // exit if maximum iteration is hit
        if(++counter > max_iteration) return SUANPAN_FAIL;

        // update internal variable
        G->update_internal(samurai);
        // update trial load factor
        G->update_trial_load_factor(incre_lambda);
        // update trial displacement
        G->update_from_ninja();
        // for tracking
        G->update_load();
        // for tracking
        G->update_constraint();

        if(D->get_attribute(ModalAttribute::LinearSystem)) return G->sync_status(false);
    }
}
