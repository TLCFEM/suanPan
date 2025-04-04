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

#include <argparse/argparse.hpp>
#include <suanPan.h>

void test_mode() {
    argparse::ArgumentParser program("set system_solver mumps", "", argparse::default_arguments::none);

    program.add_argument("--output-error-message").help("[1] output stream for error messages.").default_value(6).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--output-diagnostic-statistics-warning").help("[2] output stream for diagnostic printing and statistics local to each MPI process.").default_value(0).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--output-global-information").help("[3] output stream for global information, collected on the host.").default_value(6).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--printing-level").help("[4] level of printing for error, warning, and diagnostic messages.").default_value(2).choices(0, 1, 2, 3, 4).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--permutation-and-scaling").help("[6] permutes the matrix to a zero-free diagonal and/or scale the matrix.").default_value(7).choices(0, 1, 2, 3, 4, 5, 6, 7).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--symmetric-permutation").help("[7] computes a symmetric permutation (ordering) to determine the pivot order to be used for the factorization in case of sequential analysis.").default_value(7).choices(0, 1, 2, 3, 4, 5, 6, 7).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--scaling-strategy").help("[8] describes the scaling strategy.").default_value(77).choices(-2, -1, 0, 1, 3, 4, 7, 8, 77).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--transpose-matrix").help("[9] computes the solution using A or transpose of A.").default_value(1).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--iterative-refinement").help("[10] applies the iterative refinement to the computed solution.").default_value(0).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--error-analysis").help("[11] computes statistics related to an error analysis of the linear system solved.").default_value(0).choices(0, 1, 2).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--ordering-strategy").help("[12] defines an ordering strategy for symmetric matrices and is used, in conjunction with ICNTL(6), to add constraints to the ordering algorithm.").default_value(0).choices(0, 1, 2, 3).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--root-parallelism").help("[13] controls the parallelism of the root node (enabling or not the use of ScaLAPACK) and also its splitting.").default_value(0).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--working-space-percentage-increase").help("[14] controls the percentage increase in the estimated working space.").default_value(30).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--compression-block-format").help("[15] exploits compression of the input matrix resulting from a block format.").default_value(0).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--openmp-threads").help("[16] controls the setting of the number of OpenMP threads by MUMPS when the setting of multithreading is not possible outside MUMPS.").default_value(0).choices(0, 1, 2, 3).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--distribution-strategy-input").help("[18] defines the strategy for the distributed input matrix.").default_value(0).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--schur-complement").help("[19] computes the Schur complement matrix.").default_value(0).choices(0, 1, 2, 3).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--distribution-strategy-solution").help("[21] determines the distribution (centralized or distributed) of the solution vectors.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--out-of-core").help("[22] controls the in-core/out-of-core (OOC) factorization and solve.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--maximum-working-memory").help("[23] corresponds to the maximum size of the working memory in MB that MUMPS can allocate per working process.").default_value(0).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--null-pivot-row-detection").help("[24] controls the detection of \"null pivot rows\".").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--deficient-and-null-space-basis").help("[25] allows the computation of a solution of a deficient matrix and also of a null space basis.").default_value(0).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--schur-complement-solution").help("[26] drives the solution phase if a Schur complement matrix has been computed.").default_value(0).choices(0, 1, 2).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--rhs-block-size").help("[27] controls the blocking size for multiple right-hand sides.").default_value(-32).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--ordering-computation").help("[28] determines whether a sequential or parallel computation of the ordering is performed.").default_value(0).choices(0, 1, 2).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--inverse-computation").help("[30] computes a user-specified set of entries in the inverse of the original matrix.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--forward-elimination").help("[32] performs the forward elimination of the right-hand sides during the factorization.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--determinant-computation").help("[33] computes the determinant of the input matrix.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--out-of-core-file").help("[34] controls the conservation of the OOC files.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--blr").help("[35] controls the activation of the BLR feature.").default_value(0).choices(0, 1, 2, 3).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--blr-variant").help("[36] controls the choice of BLR factorization variant.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--blr-compression").help("[37] controls the BLR compression of the contribution blocks.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--lu-compression-rate").help("[38] estimates compression rate of LU factors.").default_value(600).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--block-compression-rate").help("[39] estimates compression rate of contribution blocks.").default_value(500).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--tree-parallelism").help("[48] controls multithreading with tree parallelism.").default_value(1).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--compact-working-space").help("[49] compact workarray at the end of factorization phase.").default_value(0).choices(0, 1, 2).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--rank-revealing-factorization").help("[56] detects pseudo-singularities during factorization and factorizes the root node with a rank-revealing method.").default_value(0).choices(0, 1).scan<'i', int>().metavar("INT").nargs(1);
    program.add_argument("--symbolic-factorization").help("[58] defines options for symbolic factorization.").default_value(2).choices(1, 2).scan<'i', int>().metavar("INT").nargs(1);

    std::string command = "--output-error-message 7 --transpose-matrix 1";
    std::istringstream iss(command);
    std::vector<std::string> args;
    args.emplace_back("set system_solver mumps");
    args.insert(args.end(), std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>());

    program.parse_args(args);

    suanpan_info("{}\n", program.get<int>("--output-error-message"));
    suanpan_info("{}\n", program.get<int>("--transpose-matrix"));
    suanpan_info("{}\n", program.help().str());
}
