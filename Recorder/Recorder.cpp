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

#include "Recorder.h"

extern fs::path SUANPAN_OUTPUT;

#ifdef SUANPAN_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

auto Recorder::normalise_size(std::vector<std::vector<vec>>& container) {
    auto cell_size = 0llu, inner_size = 0llu;
    for(const auto& inner : container) {
        if(inner.size() > inner_size) inner_size = inner.size();
        for(const auto& item : inner)
            if(item.n_elem > cell_size) cell_size = item.n_elem;
    }

    for(auto&& inner : container) {
        inner.resize(inner_size);
        for(auto&& item : inner) item.resize(cell_size);
    }

    return std::make_tuple(cell_size, inner_size);
}

const uvec& Recorder::update_tag(const shared_ptr<DomainBase>&) { return object_tag = reference_tag; }

/**
 * \brief ctor
 * \param T unique tag
 * \param B object tags
 * \param L variable type
 * \param I record interval
 * \param H if to use hdf5 format
 */
Recorder::Recorder(const unsigned T, uvec&& B, const OutputType L, const unsigned I, const bool H)
    : UniqueTag(T)
    , use_hdf5(H)
    , original_type(L)
    , variable_type(to_token(to_category(L)))
    , component(to_index(L))
    , interval(I)
    , reference_tag(std::move(B)) {}

void Recorder::initialize(const shared_ptr<DomainBase>& D) { data_pool.resize(update_tag(D).n_elem); }

void Recorder::insert(const double T) { time_pool.emplace_back(T); }

void Recorder::insert(std::vector<vec>&& data, const unsigned index) {
    if(component >= 0) {
        const auto select = static_cast<uword>(component);
        for(auto& item : data) item = vec{item.n_elem > select ? item(select) : 0.};
    }

    data_pool[index].emplace_back(std::move(data));
}
void Recorder::record(const shared_ptr<DomainBase>& D) {
    if(1 == interval || 0 == counter++ % interval) record_impl(D);
}

void Recorder::clear_status() {
    counter = 0u;
    object_tag.reset();
    time_pool.clear();
    data_pool.clear();
}

void Recorder::save() {
    std::ostringstream file_name;
    // ReSharper disable once CppIfCanBeReplacedByConstexprIf
    // ReSharper disable once CppDFAUnreachableCode
    if(comm_size > 1) file_name << 'P' << comm_rank << '-';
    file_name << 'R' << get_tag() << '-' << to_name(original_type);
    const auto origin_name = file_name.str();

#ifdef SUANPAN_HDF5
    if(use_hdf5) {
        file_name << ".h5";

        const auto file_id = H5Fcreate((SUANPAN_OUTPUT / file_name.str()).generic_string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        std::string group_name = "/";
        group_name += origin_name;

        const auto group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for(auto I = 0u; I < data_pool.size(); ++I) {
            auto& object_data = data_pool[I];
            if(object_data.empty()) continue;

            const auto [cell_size, cell_num] = normalise_size(object_data);

            mat data_to_write;
            data_to_write.zeros(cell_size * cell_num + 1, object_data.size());
            data_to_write.row(0) = rowvec{time_pool};

            for(auto J = 0llu; J < data_to_write.n_cols; ++J) {
                auto row = time_pool.empty() ? 0llu : 1llu;
                for(auto& block : object_data[J]) {
                    data_to_write(span(row, row + cell_size - 1), J) = block;
                    row += cell_size;
                }
            }

            hsize_t dimension[2]{data_to_write.n_cols, data_to_write.n_rows};

            std::ostringstream dataset_name;
            dataset_name << origin_name.c_str();
            if(object_tag.size() == data_pool.size()) dataset_name << object_tag(I);
            else dataset_name << "-SUM";

            H5LTmake_dataset(group_id, dataset_name.str().c_str(), 2, dimension, H5T_NATIVE_DOUBLE, data_to_write.mem);
        }

        H5Gclose(group_id);
        H5Fclose(file_id);
    }
    else
#endif
    {
        for(auto I = 0u; I < data_pool.size(); ++I) {
            auto& object_data = data_pool[I];
            if(object_data.empty()) continue;

            const auto [cell_size, cell_num] = normalise_size(object_data);

            mat data_to_write;
            data_to_write.zeros(cell_size * cell_num + 1, object_data.size());
            data_to_write.row(0) = rowvec{time_pool};

            for(auto J = 0llu; J < data_to_write.n_cols; ++J) {
                auto row = time_pool.empty() ? 0llu : 1llu;
                for(auto& block : object_data[J]) {
                    data_to_write(span(row, row + cell_size - 1), J) = block;
                    row += cell_size;
                }
            }

            std::ostringstream dataset_name;
            dataset_name << (SUANPAN_OUTPUT / origin_name).generic_string();
            if(object_tag.size() == data_pool.size()) dataset_name << object_tag(I);
            else dataset_name << "-SUM";

            data_to_write.save(dataset_name.str() + ".txt", raw_ascii);
        }
    }
}

void Recorder::print() {}
