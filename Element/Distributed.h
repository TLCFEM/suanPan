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

#ifndef DISTRIBUTED_H
#define DISTRIBUTED_H

#include <suanPan.h>

class Distributed {
    static constexpr int root_rank{0};

    static auto assign_process(const int obj_tag) { return obj_tag % comm_size; }

    const int tag, process_rank;

public:
    const bool is_local;

    explicit Distributed(const int obj_tag)
        : tag(obj_tag)
        , process_rank(assign_process(tag))
        , is_local(comm_rank == process_rank) {}

#ifdef SUANPAN_DISTRIBUTED
    /**
     * @brief Performs a gather operation on a distributed matrix object.
     *
     * This function initiates a non-blocking gather operation. If the calling process is the root process,
     * it receives data from the specified process. If the calling process is not the root process, it sends
     * its data to the root process.
     *
     * @tparam DT The data type of the matrix elements.
     * @param object The matrix object to be gathered.
     * @return An optional non-blocking request handle for the gather operation.
     */
    template<mpl_data_t DT> std::optional<mpl::irequest> gather(const Mat<DT>& object) {
        if(root_rank == comm_rank) {
            if(!is_local) return comm_world.irecv(const_cast<DT*>(object.memptr()), mpl::contiguous_layout<DT>{object.n_elem}, process_rank, mpl::tag_t{tag});
        }
        else if(is_local) return comm_world.isend(object.memptr(), mpl::contiguous_layout<DT>{object.n_elem}, root_rank, mpl::tag_t{tag});

        return {};
    }
#else
    template<typename T> static auto gather(T&&) {}
#endif
};

#endif // DISTRIBUTED_H
