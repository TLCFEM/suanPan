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

extern fs::path SUANPAN_OUTPUT;

#ifdef SUANPAN_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>

namespace {
    std::string join_hdf5_path(const std::string& parent, const std::string& name) {
        if(parent.empty() || "/" == parent) return "/" + name;
        return parent + "/" + name;
    }

    std::vector<std::string> list_hdf5_links(const hid_t file_id, const std::string& group_path) {
        std::vector<std::string> names;

        const auto group_id = H5Gopen2(file_id, group_path.c_str(), H5P_DEFAULT);
        if(group_id < 0) return names;

        H5G_info_t info;
        if(H5Gget_info(group_id, &info) >= 0)
            for(hsize_t I{0}; I < info.nlinks; ++I) {
                const auto length = H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, I, nullptr, 0, H5P_DEFAULT);
                if(length < 0) continue;

                std::vector<char> buffer(static_cast<size_t>(length) + 1, '\0');
                if(H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, I, buffer.data(), buffer.size(), H5P_DEFAULT) >= 0)
                    names.emplace_back(buffer.data(), static_cast<size_t>(length));
            }

        H5Gclose(group_id);

        return names;
    }

    void merge_hdf5_group(const hid_t source_id, const std::string& source_path, const hid_t target_id, const std::string& target_path) {
        for(const auto& name : list_hdf5_links(source_id, source_path)) {
            const auto source_child = join_hdf5_path(source_path, name);
            const auto target_child = join_hdf5_path(target_path, name);

            if(H5Lexists(target_id, target_child.c_str(), H5P_DEFAULT) <= 0) {
                H5Ocopy(source_id, source_child.c_str(), target_id, target_child.c_str(), H5P_DEFAULT, H5P_DEFAULT);
                continue;
            }

            // both objects exist, merge groups recursively and keep the existing dataset otherwise
            const auto target_group = H5Gopen2(target_id, target_child.c_str(), H5P_DEFAULT);
            if(target_group < 0) continue;

            H5Gclose(target_group);
            merge_hdf5_group(source_id, source_child, target_id, target_child);
        }
    }
} // namespace
#endif

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

void Recorder::combine(std::string&& recorder_name) const {
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

    const auto target_id = H5Fcreate(target_path.generic_string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(target_id < 0) {
        suanpan_warning("Failed to create combined HDF5 file {}.\n", recorder_name);
        return;
    }

    bool success = true;

    for(const auto& file_path : file_paths) {
        const auto source_id = H5Fopen(file_path.generic_string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if(source_id < 0) {
            suanpan_warning("Failed to open HDF5 file {} for combination.\n", file_path.generic_string());
            success = false;
            continue;
        }

        merge_hdf5_group(source_id, "/", target_id, "/");

        H5Fclose(source_id);
    }

    H5Fclose(target_id);

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

        const auto file_id = H5Fcreate((SUANPAN_OUTPUT / file_name.str()).generic_string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        std::string group_name = "/";
        group_name += origin_name;

        const auto group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

            hsize_t dimension[2]{data_to_write.n_cols, data_to_write.n_rows};

            std::ostringstream dataset_name;
            dataset_name << origin_name.c_str() << tag;

            H5LTmake_dataset(group_id, dataset_name.str().c_str(), 2, dimension, H5T_NATIVE_DOUBLE, data_to_write.mem);
        }

        H5Gclose(group_id);
        H5Fclose(file_id);

        D->insert(std::async(std::launch::async, &Recorder::combine, (recorder_name << ".h5").str()));
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
