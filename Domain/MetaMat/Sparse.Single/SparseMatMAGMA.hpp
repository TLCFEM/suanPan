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
 * @class SparseMatMAGMA
 * @brief A SparseMatMAGMA class that holds matrices.
 *
 * TODO: improve performance by storing factorization and reusing it
 *
 * @author tlc
 * @date 24/02/2023
 * @version 0.1.0
 * @file SparseMatMAGMA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
// ReSharper disable CppClangTidyBugproneBranchClone
// ReSharper disable CppClangTidyClangDiagnosticMissingFieldInitializers
#ifndef SPARSEMATMAGMA_HPP
#define SPARSEMATMAGMA_HPP

#include "SparseMat.hpp"

#include <magma_v2.h>
#include <magmasparse.h>

// ReSharper disable IdentifierTypo
namespace suanpan::detail::magma {
    template<typename> struct pack {};
    template<> struct pack<double> {
        using opt_type = magma_dopts;
        using mat_type = magma_d_matrix;
        static constexpr auto& csrset = magma_dcsrset;
        static constexpr auto& free = magma_dsolverinfo_free;
        static constexpr auto& init = magma_dsolverinfo_init;
        static constexpr auto& mfree = magma_dmfree;
        static constexpr auto& mtransfer = magma_dmtransfer;
        static constexpr auto& precondsetup = magma_d_precondsetup;
        static constexpr auto& solver = magma_d_solver;
        static constexpr auto& vcopy = magma_dvcopy;
        static constexpr auto& vinit = magma_dvinit;
        static constexpr auto& vset = magma_dvset;
    };
    template<> struct pack<float> {
        using opt_type = magma_sopts;
        using mat_type = magma_s_matrix;
        static constexpr auto& csrset = magma_scsrset;
        static constexpr auto& free = magma_ssolverinfo_free;
        static constexpr auto& init = magma_ssolverinfo_init;
        static constexpr auto& mfree = magma_smfree;
        static constexpr auto& mtransfer = magma_smtransfer;
        static constexpr auto& precondsetup = magma_s_precondsetup;
        static constexpr auto& solver = magma_s_solver;
        static constexpr auto& vcopy = magma_svcopy;
        static constexpr auto& vinit = magma_svinit;
        static constexpr auto& vset = magma_svset;
    };

    template<typename T> auto parse_opts(T& opts, std::istringstream&& command) requires(std::is_same_v<T, magma_dopts> || std::is_same_v<T, magma_sopts>) {
        using ET = decltype(opts.solver_par.atol);

        opts.blocksize = 32;
        opts.alignment = 1;
        opts.scaling = Magma_NOSCALE;
        opts.solver_par.solver = Magma_PGMRES;
        opts.solver_par.atol = std::numeric_limits<ET>::epsilon();
        opts.solver_par.rtol = std::numeric_limits<ET>::epsilon() * 100;
        opts.solver_par.maxiter = 1000;
        opts.solver_par.verbose = 0;
        opts.solver_par.restart = 50;
        opts.precond_par.solver = Magma_ILU;
        opts.precond_par.trisolver = Magma_CUSOLVE;
        opts.precond_par.atol = ET(1);
        opts.precond_par.rtol = std::numeric_limits<ET>::epsilon() * 100;
        opts.precond_par.maxiter = 100;
        opts.precond_par.restart = 10;
        opts.precond_par.levels = 0;
        opts.precond_par.sweeps = 5;
        opts.precond_par.pattern = 1;

        int basic{0};

        std::string token;
        // ReSharper disable StringLiteralTypo
        while(get_input(command, token)) {
            if(is_equal(token, "--mscale")) {
                if(std::string scale; get_input(command, scale)) {
                    if(is_equal("NOSCALE", scale)) opts.scaling = Magma_NOSCALE;
                    else if(is_equal("UNITDIAG", scale)) opts.scaling = Magma_UNITDIAG;
                    else if(is_equal("UNITROW", scale)) opts.scaling = Magma_UNITROW;
                    else if(is_equal("UNITCOL", scale)) opts.scaling = Magma_UNITCOL;
                    else if(is_equal("UNITDIAGCOL", scale)) opts.scaling = Magma_UNITDIAGCOL;
                    else if(is_equal("UNITROWCOL", scale)) opts.scaling = Magma_UNITROWCOL;
                }
            }
            else if(is_equal(token, "--solver")) {
                if(std::string solver; get_input(command, solver)) {
                    if(is_equal("CG", solver)) opts.solver_par.solver = Magma_PCGMERGE;
                    else if(is_equal("PCG", solver)) opts.solver_par.solver = Magma_PCGMERGE;
                    else if(is_equal("BICG", solver)) opts.solver_par.solver = Magma_PBICG;
                    else if(is_equal("PBICG", solver)) opts.solver_par.solver = Magma_PBICG;
                    else if(is_equal("BICGSTAB", solver)) opts.solver_par.solver = Magma_PBICGSTABMERGE;
                    else if(is_equal("PBICGSTAB", solver)) opts.solver_par.solver = Magma_PBICGSTABMERGE;
                    else if(is_equal("QMR", solver)) opts.solver_par.solver = Magma_PQMRMERGE;
                    else if(is_equal("PQMR", solver)) opts.solver_par.solver = Magma_PQMRMERGE;
                    else if(is_equal("TFQMR", solver)) opts.solver_par.solver = Magma_PTFQMRMERGE;
                    else if(is_equal("PTFQMR", solver)) opts.solver_par.solver = Magma_PTFQMRMERGE;
                    else if(is_equal("GMRES", solver)) opts.solver_par.solver = Magma_PGMRES;
                    else if(is_equal("PGMRES", solver)) opts.solver_par.solver = Magma_PGMRES;
                    else if(is_equal("LOBPCG", solver)) opts.solver_par.solver = Magma_LOBPCG;
                    else if(is_equal("LSQR", solver)) opts.solver_par.solver = Magma_LSQR;
                    else if(is_equal("JACOBI", solver)) opts.solver_par.solver = Magma_JACOBI;
                    else if(is_equal("BA", solver)) opts.solver_par.solver = Magma_BAITER;
                    else if(is_equal("BAO", solver)) opts.solver_par.solver = Magma_BAITERO;
                    else if(is_equal("IDR", solver)) opts.solver_par.solver = Magma_PIDRMERGE;
                    else if(is_equal("PIDR", solver)) opts.solver_par.solver = Magma_PIDRMERGE;
                    else if(is_equal("CGS", solver)) opts.solver_par.solver = Magma_PCGSMERGE;
                    else if(is_equal("PCGS", solver)) opts.solver_par.solver = Magma_PCGSMERGE;
                    else if(is_equal("BOMBARDMENT", solver)) opts.solver_par.solver = Magma_BOMBARDMERGE;
                    else if(is_equal("ITERREF", solver)) opts.solver_par.solver = Magma_ITERREF;
                    else if(is_equal("PARDISO", solver)) opts.solver_par.solver = Magma_PARDISO;
                }
            }
            else if(is_equal(token, "--precond")) {
                if(std::string solver; get_input(command, solver)) {
                    if(is_equal("CG", solver)) opts.precond_par.solver = Magma_CGMERGE;
                    else if(is_equal("PCG", solver)) opts.precond_par.solver = Magma_PCG;
                    else if(is_equal("BICGSTAB", solver)) opts.precond_par.solver = Magma_BICGSTABMERGE;
                    else if(is_equal("QMR", solver)) opts.precond_par.solver = Magma_QMRMERGE;
                    else if(is_equal("TFQMR", solver)) opts.precond_par.solver = Magma_TFQMRMERGE;
                    else if(is_equal("PTFQMR", solver)) opts.precond_par.solver = Magma_PTFQMRMERGE;
                    else if(is_equal("GMRES", solver)) opts.precond_par.solver = Magma_GMRES;
                    else if(is_equal("PGMRES", solver)) opts.precond_par.solver = Magma_PGMRES;
                    else if(is_equal("LOBPCG", solver)) opts.precond_par.solver = Magma_LOBPCG;
                    else if(is_equal("JACOBI", solver)) opts.precond_par.solver = Magma_JACOBI;
                    else if(is_equal("BA", solver)) opts.precond_par.solver = Magma_BAITER;
                    else if(is_equal("BAO", solver)) opts.precond_par.solver = Magma_BAITERO;
                    else if(is_equal("IDR", solver)) opts.precond_par.solver = Magma_IDRMERGE;
                    else if(is_equal("PIDR", solver)) opts.precond_par.solver = Magma_PIDRMERGE;
                    else if(is_equal("CGS", solver)) opts.precond_par.solver = Magma_CGSMERGE;
                    else if(is_equal("PCGS", solver)) opts.precond_par.solver = Magma_PCGSMERGE;
                    else if(is_equal("BOMBARDMENT", solver)) opts.precond_par.solver = Magma_BOMBARD;
                    else if(is_equal("ITERREF", solver)) opts.precond_par.solver = Magma_ITERREF;
                    else if(is_equal("ILU", solver) || is_equal("IC", solver)) opts.precond_par.solver = Magma_ILU;
                    else if(is_equal("ILUT", solver) || is_equal("ICT", solver)) opts.precond_par.solver = Magma_ILUT;
                    else if(is_equal("PARILU", solver) || is_equal("AIC", solver)) opts.precond_par.solver = Magma_PARILU;
                    else if(is_equal("PARIC", solver)) opts.precond_par.solver = Magma_PARILU;
                    else if(is_equal("PARICT", solver)) opts.precond_par.solver = Magma_PARICT;
                    else if(is_equal("PARILUT", solver)) opts.precond_par.solver = Magma_PARILUT;
                    else if(is_equal("CUSTOMIC", solver)) opts.precond_par.solver = Magma_CUSTOMIC;
                    else if(is_equal("CUSTOMILU", solver)) opts.precond_par.solver = Magma_CUSTOMILU;
                    else if(is_equal("ISAI", solver)) opts.precond_par.solver = Magma_ISAI;
                    else if(is_equal("NONE", solver)) opts.precond_par.solver = Magma_NONE;
                }
            }
            else if(is_equal(token, "--trisolver")) {
                if(std::string solver; get_input(command, solver)) {
                    if(is_equal("CG", solver)) opts.precond_par.trisolver = Magma_CGMERGE;
                    else if(is_equal("BICGSTAB", solver)) opts.precond_par.trisolver = Magma_BICGSTABMERGE;
                    else if(is_equal("QMR", solver)) opts.precond_par.trisolver = Magma_QMRMERGE;
                    else if(is_equal("TFQMR", solver)) opts.precond_par.trisolver = Magma_TFQMRMERGE;
                    else if(is_equal("GMRES", solver)) opts.precond_par.trisolver = Magma_GMRES;
                    else if(is_equal("JACOBI", solver)) opts.precond_par.trisolver = Magma_JACOBI;
                    else if(is_equal("VBJACOBI", solver)) opts.precond_par.trisolver = Magma_VBJACOBI;
                    else if(is_equal("BA", solver)) opts.precond_par.trisolver = Magma_BAITER;
                    else if(is_equal("BAO", solver)) opts.precond_par.trisolver = Magma_BAITERO;
                    else if(is_equal("IDR", solver)) opts.precond_par.trisolver = Magma_IDRMERGE;
                    else if(is_equal("CGS", solver)) opts.precond_par.trisolver = Magma_CGSMERGE;
                    else if(is_equal("CUSOLVE", solver)) opts.precond_par.trisolver = Magma_CUSOLVE;
                    else if(is_equal("SYNCFREESOLVE", solver)) opts.precond_par.trisolver = Magma_SYNCFREESOLVE;
                    else if(is_equal("ISAI", solver)) opts.precond_par.trisolver = Magma_ISAI;
                    else if(is_equal("NONE", solver)) opts.precond_par.trisolver = Magma_NONE;
                }
            }
            else if(is_equal(token, "--basic")) basic = 1;
            else if(is_equal(token, "--blocksize")) get_input(command, opts.blocksize);
            else if(is_equal(token, "--alignment")) get_input(command, opts.alignment);
            else if(is_equal(token, "--restart")) get_input(command, opts.solver_par.restart);
            else if(is_equal(token, "--atol")) get_input(command, opts.solver_par.atol);
            else if(is_equal(token, "--rtol")) get_input(command, opts.solver_par.rtol);
            else if(is_equal(token, "--maxiter")) get_input(command, opts.solver_par.maxiter);
            else if(is_equal(token, "--verbose")) get_input(command, opts.solver_par.verbose);
            else if(is_equal(token, "--prestart")) get_input(command, opts.precond_par.restart);
            else if(is_equal(token, "--patol")) get_input(command, opts.precond_par.atol);
            else if(is_equal(token, "--prtol")) get_input(command, opts.precond_par.rtol);
            else if(is_equal(token, "--piters")) get_input(command, opts.precond_par.maxiter);
            else if(is_equal(token, "--ppattern")) get_input(command, opts.precond_par.pattern);
            else if(is_equal(token, "--psweeps")) get_input(command, opts.precond_par.sweeps);
            else if(is_equal(token, "--plevels")) get_input(command, opts.precond_par.levels);
        }
        // ReSharper restore StringLiteralTypo

        if(basic == 1) {
            if(opts.solver_par.solver == Magma_PCGMERGE) opts.solver_par.solver = Magma_PCG;
            else if(opts.solver_par.solver == Magma_PBICGSTABMERGE) opts.solver_par.solver = Magma_PBICGSTAB;
            else if(opts.solver_par.solver == Magma_PBICGMERGE) opts.solver_par.solver = Magma_PBICG;
            else if(opts.solver_par.solver == Magma_PTFQMRMERGE) opts.solver_par.solver = Magma_PTFQMR;
            else if(opts.solver_par.solver == Magma_PCGSMERGE) opts.solver_par.solver = Magma_PCGS;
            else if(opts.solver_par.solver == Magma_PQMRMERGE) opts.solver_par.solver = Magma_PQMR;
            else if(opts.solver_par.solver == Magma_PIDRMERGE) opts.solver_par.solver = Magma_PIDR;
            else if(opts.solver_par.solver == Magma_BOMBARDMERGE) opts.solver_par.solver = Magma_BOMBARD;
        }

        if(opts.precond_par.solver == Magma_NONE || opts.precond_par.solver == 0) {
            switch(opts.solver_par.solver) /* NOLINT(clang-diagnostic-switch-enum) */ {
            case Magma_PBICG:
                opts.solver_par.solver = Magma_BICG;
                break;
            case Magma_PBICGMERGE:
                opts.solver_par.solver = Magma_BICGMERGE;
                break;
            case Magma_PBICGSTAB:
                opts.solver_par.solver = Magma_BICGSTAB;
                break;
            case Magma_PBICGSTABMERGE:
                opts.solver_par.solver = Magma_BICGSTABMERGE;
                break;
            case Magma_PCG:
                opts.solver_par.solver = Magma_CG;
                break;
            case Magma_PCGMERGE:
                opts.solver_par.solver = Magma_CGMERGE;
                break;
            case Magma_PCGS:
                opts.solver_par.solver = Magma_CGS;
                break;
            case Magma_PCGSMERGE:
                opts.solver_par.solver = Magma_CGSMERGE;
                break;
            case Magma_PQMR:
                opts.solver_par.solver = Magma_QMR;
                break;
            case Magma_PQMRMERGE:
                opts.solver_par.solver = Magma_QMRMERGE;
                break;
            case Magma_PTFQMR:
                opts.solver_par.solver = Magma_TFQMR;
                break;
            case Magma_PTFQMRMERGE:
                opts.solver_par.solver = Magma_TFQMRMERGE;
                break;
            case Magma_PIDR:
                opts.solver_par.solver = Magma_IDR;
                break;
            case Magma_PIDRMERGE:
                opts.solver_par.solver = Magma_IDRMERGE;
                break;
            case Magma_PGMRES:
                opts.solver_par.solver = Magma_GMRES;
                break;
            default:
                break;
            }
        }

        if(opts.solver_par.solver == Magma_PCG || opts.solver_par.solver == Magma_PCGMERGE) {
            if(opts.precond_par.solver == Magma_ILU) opts.precond_par.solver = Magma_ICC;
            if(opts.precond_par.solver == Magma_PARILU) opts.precond_par.solver = Magma_PARIC;
        }
        if(opts.output_format == Magma_CSR5) {
            if(opts.solver_par.solver == Magma_CGMERGE) opts.solver_par.solver = Magma_CG;
            if(opts.solver_par.solver == Magma_PCGMERGE) opts.solver_par.solver = Magma_PCG;
        }
    }
} // namespace suanpan::detail::magma
// ReSharper restore IdentifierTypo

template<sp_d T> class SparseMatMAGMA final : public SparseMat<T> {
    using magma_s = suanpan::detail::magma::pack<T>;
    using opt_t = typename magma_s::opt_type;
    using mat_t = typename magma_s::mat_type;

    csr_form<T, magma_index_t> csr_mat;

    mat_t A_host{Magma_CSR, Magma_CPU};
    mat_t A_device{Magma_CSR, Magma_DEV};
    mat_t b_host{Magma_DENSE, Magma_CPU};
    mat_t b_device{Magma_DENSE, Magma_DEV};
    mat_t x_host{Magma_DENSE, Magma_CPU};
    mat_t x_device{Magma_DENSE, Magma_DEV};

    opt_t opts{};

    magma_queue_t queue{};

protected:
    int direct_solve(Mat<T>& out_mat, const Mat<T>& in_mat) override { return this->direct_solve(out_mat, Mat<T>(in_mat)); }

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    SparseMatMAGMA(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) {
        magma_init();
        if(SUANPAN_VERBOSE) magma_print_environment();
        magma_queue_create(0, &queue);
        magma_s::init(&opts.solver_par, &opts.precond_par, queue);
    }

    SparseMatMAGMA(const SparseMatMAGMA&);
    SparseMatMAGMA(SparseMatMAGMA&&) = delete;
    SparseMatMAGMA& operator=(const SparseMatMAGMA&) = delete;
    SparseMatMAGMA& operator=(SparseMatMAGMA&&) = delete;

    ~SparseMatMAGMA() override {
        magma_s::mfree(&x_device, queue);
        magma_s::mfree(&b_device, queue);
        magma_s::mfree(&A_device, queue);
        magma_s::mfree(&x_host, queue);
        magma_s::mfree(&b_host, queue);
        magma_s::mfree(&A_host, queue);
        magma_s::free(&opts.solver_par, &opts.precond_par, queue);
        magma_queue_destroy(queue);
        magma_finalize();
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatMAGMA>(*this); }
};

template<sp_d T> SparseMatMAGMA<T>::SparseMatMAGMA(const SparseMatMAGMA& other)
    : SparseMat<T>(other)
    , csr_mat(other.csr_mat) {
    magma_queue_create(0, &queue);
    magma_s::init(&opts.solver_par, &opts.precond_par, queue);
    if(this->factored) {
        magma_s::csrset(csr_mat.n_rows, csr_mat.n_cols, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem(), &A_host, queue);
        magma_s::mtransfer(A_host, &A_device, Magma_CPU, Magma_DEV, queue);
    }
}

template<sp_d T> int SparseMatMAGMA<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    auto b_rows = static_cast<magma_index_t>(B.n_rows), b_cols = static_cast<magma_index_t>(B.n_cols);

    if(!this->factored) {
        csr_mat = csr_form<T, magma_index_t>(this->triplet_mat, SparseBase::ZERO);
        this->factored = true;
        magma_s::csrset(csr_mat.n_rows, csr_mat.n_cols, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem(), &A_host, queue);
        magma_s::mtransfer(A_host, &A_device, Magma_CPU, Magma_DEV, queue);
    }

    suanpan::detail::magma::parse_opts(opts, std::istringstream(this->setting.option));

    magma_s::vset(b_rows, b_cols, (T*)B.memptr(), &b_host, queue);
    magma_s::mtransfer(b_host, &b_device, Magma_CPU, Magma_DEV, queue);

    magma_s::precondsetup(A_device, b_device, &opts.solver_par, &opts.precond_par, queue);

    magma_s::vinit(&x_device, Magma_DEV, b_rows, b_cols, T(0), queue);

    magma_s::solver(A_device, b_device, &x_device, &opts, queue);

    magma_s::mtransfer(x_device, &x_host, Magma_DEV, Magma_CPU, queue);

    X = std::move(B);
    magma_s::vcopy(x_host, &b_rows, &b_cols, X.memptr(), queue);

    return SUANPAN_SUCCESS;
}

#endif

//! @}
