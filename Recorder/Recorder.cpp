/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include <Domain/DomainBase.h>
#include <Toolbox/HDF5Wrapper.h>


extern fs::path SUANPAN_OUTPUT;

auto Recorder::normalise_size(std::vector<std::vector<vec>>& container) {
    uword cell_size{0}, inner_size{0};
    for(const auto& inner : container) {
        if(inner.size() > inner_size) inner_size = static_cast<uword>(inner.size());
        for(const auto& item : inner)
            if(item.n_elem > cell_size) cell_size = item.n_elem;
    }

    for(auto&& inner : container) {
        inner.resize(inner_size);
        for(auto&& item : inner) item.resize(cell_size);
    }

    return std::make_tuple(cell_size, inner_size);
}

std::vector<vec> Recorder::normalise_size(std::vector<vec>&& container) {
    uword max_size{0};
    for(const auto& item : container)
        if(item.n_elem > max_size) max_size = item.n_elem;

    for(auto&& item : container) item.resize(max_size);

    return container;
}

const uvec& Recorder::update_tag(const shared_ptr<DomainBase>&) { return object_tag = reference_tag; }

void Recorder::combine(const std::string_view recorder_name) const {
    if(!is_distributed() || comm_size == 1) return;

    barrier();

    // ReSharper disable once CppIfCanBeReplacedByConstexprIf
    if(0 != comm_rank) return;

#ifdef SUANPAN_HDF5
    std::vector<fs::path> file_paths;
    file_paths.reserve(comm_size);

    for(auto I{0}; I < comm_size; ++I) {
        std::ostringstream distributed_name;
        distributed_name << 'P' << I << '-' << recorder_name;
        file_paths.emplace_back(SUANPAN_OUTPUT / distributed_name.str());
    }

    std::erase_if(file_paths, [](const fs::path& file_path) { return !fs::exists(file_path); });
    if(file_paths.empty()) return;

    const auto target_path = SUANPAN_OUTPUT / recorder_name;

    hdf5::File target(target_path, hdf5::File::Mode::Truncate);
    if(!target) {
        suanpan_warning("Failed to create combined HDF5 file {}.\n", recorder_name);
        return;
    }

    bool success = true;

    for(const auto& file_path : file_paths)
        if(!target.merge(hdf5::File(file_path, hdf5::File::Mode::ReadOnly))) {
            suanpan_warning("Failed to open HDF5 file {} for combination.\n", file_path.generic_string());
            success = false;
            continue;
        }

    target.close();

    if(!success) {
        std::error_code error;
        fs::remove(target_path, error);
        return;
    }

    for(const auto& file_path : file_paths) {
        std::error_code error;
        fs::remove(file_path, error);
        if(error) suanpan_warning("Failed to remove distributed HDF5 file {}: {}.\n", file_path.generic_string(), error.message());
    }
#endif
}

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

void Recorder::initialize(const shared_ptr<DomainBase>& D) { update_tag(D); }

void Recorder::insert(const double T) { time_pool.emplace_back(T); }

void Recorder::insert(std::vector<vec>&& data, const uword tag) {
    if(component >= 0) {
        const auto select = static_cast<uword>(component);
        for(auto& item : data) item = vec{item.n_elem > select ? item(select) : 0.};
    }

    data_pool[tag].emplace_back(std::move(data));
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

void Recorder::save(const shared_ptr<DomainBase>& D) {
    std::ostringstream file_name, recorder_name;
    // ReSharper disable once CppIfCanBeReplacedByConstexprIf
    // ReSharper disable once CppDFAUnreachableCode
    if(comm_size > 1) file_name << 'P' << comm_rank << '-';
    recorder_name << 'R' << get_tag() << '-' << to_name(original_type);
    file_name << recorder_name.str();
    const auto origin_name = file_name.str();

#ifdef SUANPAN_HDF5
    if(use_hdf5) {
        file_name << ".h5";

        hdf5::File file(SUANPAN_OUTPUT / file_name.str(), hdf5::File::Mode::Truncate);

        auto group = file.create_group('/' + origin_name);

        for(auto& [tag, object_data] : data_pool) {
            if(object_data.empty()) continue;

            const auto [cell_size, cell_num] = normalise_size(object_data);

            mat data_to_write(cell_size * cell_num + 1, static_cast<uword>(object_data.size()), fill::zeros);
            data_to_write.row(0) = rowvec{time_pool};

            for(uword J{0}; J < data_to_write.n_cols; ++J) {
                uword row{1};
                for(auto& block : object_data[J]) {
                    data_to_write(span(row, row + cell_size - 1), J) = block;
                    row += cell_size;
                }
            }

            std::ostringstream dataset_name;
            dataset_name << origin_name.c_str() << tag;

            group.create_dataset(dataset_name.str(), data_to_write.mem, {static_cast<hsize_t>(data_to_write.n_cols), static_cast<hsize_t>(data_to_write.n_rows)});
        }

        group.close();
        file.close();

        recorder_name << ".h5";
        D->insert(std::async(std::launch::async, &Recorder::combine, this, recorder_name.str()));
    }
    else
#endif
    {
        for(auto& [tag, object_data] : data_pool) {
            if(object_data.empty()) continue;

            const auto [cell_size, cell_num] = normalise_size(object_data);

            mat data_to_write(static_cast<uword>(object_data.size()), cell_size * cell_num + 1, fill::zeros);
            data_to_write.col(0) = vec{time_pool};

            for(uword J{0}; J < data_to_write.n_rows; ++J) {
                uword col{1};
                for(auto& block : object_data[J]) {
                    data_to_write(J, span(col, col + cell_size - 1)) = block.t();
                    col += cell_size;
                }
            }

            std::ostringstream dataset_name;
            dataset_name << (SUANPAN_OUTPUT / origin_name).generic_string() << tag;

            data_to_write.save(dataset_name.str() + ".txt", raw_ascii);
        }
    }
}

void Recorder::print() {}
