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
#ifdef SUANPAN_MKL

#include <mkl_cluster_sparse_solver.h>
#include <mpl/mpl.hpp>

int main(int argc, char** argv) {
    const auto& comm_world{mpl::environment::comm_world()};
    const auto& parent = mpl::inter_communicator::parent();
    const auto all = mpl::communicator(parent, mpl::communicator::order_high);

    int config[7]{};
    int iparm[64]{};

    all.bcast(0, config);
    all.bcast(0, iparm);

    const auto mtype = &config[0];
    const auto nrhs = &config[1];
    const auto maxfct = &config[2];
    const auto mnum = &config[3];
    const auto msglvl = &config[4];
    const auto n = &config[5];
    const auto nnz = &config[6];

    const auto float_type = iparm[27];

    const auto np = *n + 1;
    const auto nb = *n * *nrhs;

    std::vector<int> ia(np), ja(*nnz);
    const auto a = std::make_unique<double[]>(*nnz);
    const auto b = std::make_unique<double[]>(nb);
    const auto x = std::make_unique<double[]>(nb);

    if(0 == comm_world.rank()) {
        mpl::irequest_pool requests;

        requests.push(parent.irecv(ia, 0, mpl::tag_t{0}));
        requests.push(parent.irecv(ja, 0, mpl::tag_t{1}));
        if(0 == float_type) {
            requests.push(parent.irecv(a.get(), a.get() + *nnz, 0, mpl::tag_t{2}));
            requests.push(parent.irecv(b.get(), b.get() + nb, 0, mpl::tag_t{3}));
        }
        else {
            requests.push(parent.irecv(reinterpret_cast<float*>(a.get()), reinterpret_cast<float*>(a.get()) + *nnz, 0, mpl::tag_t{2}));
            requests.push(parent.irecv(reinterpret_cast<float*>(b.get()), reinterpret_cast<float*>(b.get()) + nb, 0, mpl::tag_t{3}));
        }

        requests.waitall();
    }

    int64_t pt[64]{};

    const int comm = MPI_Comm_c2f(comm_world.native_handle());

    int error = 0;

    auto finalise = [&] {
        if(0 != comm_world.rank()) return error;

        parent.send(error, 0);
        if(0 != error) return error;

        if(0 == float_type) parent.send(x.get(), x.get() + nb, 0);
        else parent.send(reinterpret_cast<float*>(x.get()), reinterpret_cast<float*>(x.get()) + nb, 0);

        return error;
    };

    int phase = 13;
    cluster_sparse_solver(pt, maxfct, mnum, mtype, &phase, n, a.get(), ia.data(), ja.data(), nullptr, nrhs, iparm, msglvl, b.get(), x.get(), &comm, &error);

    if(0 != error) return finalise();

    phase = -1;
    cluster_sparse_solver(pt, maxfct, mnum, mtype, &phase, n, nullptr, ia.data(), ja.data(), nullptr, nrhs, iparm, msglvl, nullptr, nullptr, &comm, &error);

    return finalise();
}

#else

int main(int, char**) { return 0; }

#endif
