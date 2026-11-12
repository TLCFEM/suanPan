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

#include "HDF5Wrapper.h"

#ifdef SUANPAN_HDF5

#include <utility>

namespace {
    constexpr hid_t invalid_id = -1;

    std::string join_path(const std::string& parent, const std::string& name) {
        if(parent.empty() || "/" == parent || "." == parent) return name;
        return parent + '/' + name;
    }

    hid_t create_group_id(const hid_t location, const std::string& name) {
        if(location < 0 || name.empty()) return invalid_id;
        return H5Gcreate2(location, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t open_group_id(const hid_t location, const std::string& name) {
        if(location < 0 || name.empty()) return invalid_id;
        return H5Gopen2(location, name.c_str(), H5P_DEFAULT);
    }

    hid_t create_dataset_id(const hid_t location, const std::string& name, const double* data, const std::vector<hsize_t>& shape) {
        if(location < 0 || name.empty() || nullptr == data || shape.empty()) return invalid_id;

        const auto space = H5Screate_simple(static_cast<int>(shape.size()), shape.data(), nullptr);
        if(space < 0) return invalid_id;

        const auto dataset = H5Dcreate2(location, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);

        if(dataset < 0) return invalid_id;

        if(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
            H5Dclose(dataset);
            return invalid_id;
        }

        return dataset;
    }

    hid_t open_dataset_id(const hid_t location, const std::string& name) {
        if(location < 0 || name.empty()) return invalid_id;
        return H5Dopen2(location, name.c_str(), H5P_DEFAULT);
    }

    bool exists_impl(const hid_t location, const std::string& name) {
        return location >= 0 && !name.empty() && H5Lexists(location, name.c_str(), H5P_DEFAULT) > 0;
    }

    std::vector<std::string> links_impl(const hid_t location, const std::string& path) {
        std::vector<std::string> names;

        const auto group = open_group_id(location, path);
        if(group < 0) return names;

        H5G_info_t info;
        if(H5Gget_info(group, &info) >= 0)
            for(hsize_t I{0}; I < info.nlinks; ++I) {
                const auto length = H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, I, nullptr, 0, H5P_DEFAULT);
                if(length < 0) continue;

                std::vector<char> buffer(static_cast<std::size_t>(length) + 1, '\0');
                if(H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, I, buffer.data(), buffer.size(), H5P_DEFAULT) >= 0)
                    names.emplace_back(buffer.data(), static_cast<std::size_t>(length));
            }

        H5Gclose(group);

        return names;
    }

    bool is_group_impl(const hid_t location, const std::string& name) {
        const auto group = open_group_id(location, name);
        if(group < 0) return false;
        H5Gclose(group);
        return true;
    }

    bool copy_impl(const hid_t source, const std::string& source_path, const hid_t target, const std::string& target_path) {
        if(source < 0 || target < 0 || source_path.empty() || target_path.empty()) return false;
        return H5Ocopy(source, source_path.c_str(), target, target_path.c_str(), H5P_DEFAULT, H5P_DEFAULT) >= 0;
    }

    bool merge_impl(const hid_t source, const std::string& source_path, const hid_t target, const std::string& target_path) {
        bool success = true;

        for(const auto& name : links_impl(source, source_path)) {
            const auto source_child = join_path(source_path, name), target_child = join_path(target_path, name);

            if(!exists_impl(target, target_child)) {
                if(!copy_impl(source, source_child, target, target_child)) success = false;
            }
            else if(is_group_impl(target, target_child))
                success = merge_impl(source, source_child, target, target_child) && success;
        }

        return success;
    }
} // namespace

namespace hdf5 {
    Dataset::~Dataset() { close(); }

    Dataset::Dataset(Dataset&& other) noexcept
        : id_(std::exchange(other.id_, invalid_id)) {}

    Dataset& Dataset::operator=(Dataset&& other) noexcept {
        if(this != &other) {
            close();
            id_ = std::exchange(other.id_, invalid_id);
        }
        return *this;
    }

    void Dataset::close() noexcept {
        if(id_ >= 0) H5Dclose(id_);
        id_ = invalid_id;
    }

    std::vector<hsize_t> Dataset::shape() const {
        std::vector<hsize_t> dimensions;

        if(id_ < 0) return dimensions;

        const auto space = H5Dget_space(id_);
        if(space < 0) return dimensions;

        if(const auto rank = H5Sget_simple_extent_ndims(space); rank > 0) {
            dimensions.resize(static_cast<std::size_t>(rank));
            if(H5Sget_simple_extent_dims(space, dimensions.data(), nullptr) < 0) dimensions.clear();
        }

        H5Sclose(space);

        return dimensions;
    }

    std::size_t Dataset::size() const {
        if(id_ < 0) return 0;

        std::size_t total = 1;
        for(const auto dimension : shape()) total *= static_cast<std::size_t>(dimension);
        return total;
    }

    bool Dataset::read(double* buffer) const { return id_ >= 0 && nullptr != buffer && H5Dread(id_, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) >= 0; }

    std::vector<double> Dataset::read() const {
        std::vector<double> data;

        if(id_ < 0) return data;

        data.resize(size());
        if(!data.empty() && !read(data.data())) data.clear();

        return data;
    }

    Group::~Group() { close(); }

    Group::Group(Group&& other) noexcept
        : id_(std::exchange(other.id_, invalid_id)) {}

    Group& Group::operator=(Group&& other) noexcept {
        if(this != &other) {
            close();
            id_ = std::exchange(other.id_, invalid_id);
        }
        return *this;
    }

    void Group::close() noexcept {
        if(id_ >= 0) H5Gclose(id_);
        id_ = invalid_id;
    }

    Group Group::create(const File& parent, const std::string& name) { return Group{create_group_id(parent.id(), name)}; }

    Group Group::create(const Group& parent, const std::string& name) { return Group{create_group_id(parent.id(), name)}; }

    Group Group::open(const File& parent, const std::string& name) { return Group{open_group_id(parent.id(), name)}; }

    Group Group::open(const Group& parent, const std::string& name) { return Group{open_group_id(parent.id(), name)}; }

    Group Group::create_group(const std::string& name) const { return Group{create_group_id(id_, name)}; }

    Group Group::open_group(const std::string& name) const { return Group{open_group_id(id_, name)}; }

    Dataset Group::create_dataset(const std::string& name, const double* data, const std::initializer_list<hsize_t> shape) const { return Dataset{create_dataset_id(id_, name, data, std::vector<hsize_t>{shape})}; }

    Dataset Group::create_dataset(const std::string& name, const double* data, const std::vector<hsize_t>& shape) const { return Dataset{create_dataset_id(id_, name, data, shape)}; }

    Dataset Group::create_dataset(const std::string& name, const std::vector<double>& data) const { return Dataset{create_dataset_id(id_, name, data.data(), {data.size()})}; }

    Dataset Group::open_dataset(const std::string& name) const { return Dataset{open_dataset_id(id_, name)}; }

    bool Group::exists(const std::string& name) const { return exists_impl(id_, name); }

    std::vector<std::string> Group::links() const { return links_impl(id_, "."); }

    std::string Group::name() const {
        if(id_ < 0) return {};

        const auto length = H5Iget_name(id_, nullptr, 0);
        if(length <= 0) return {};

        std::string result(static_cast<std::size_t>(length) + 1, '\0');
        if(H5Iget_name(id_, result.data(), result.size()) < 0) return {};

        result.resize(static_cast<std::size_t>(length));
        return result;
    }

    bool Group::merge(const Group& source) const { return merge_impl(source.id(), ".", id_, "."); }

    File::File(const std::filesystem::path& path, const Mode mode) {
        const auto filename = path.string();

        switch(mode) {
        case Mode::ReadOnly:
            id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            break;
        case Mode::ReadWrite:
            id_ = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            break;
        case Mode::CreateNew:
            id_ = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
            break;
        case Mode::Truncate:
        default:
            id_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            break;
        }
    }

    File::File(const std::string& path, const Mode mode)
        : File(std::filesystem::path{path}, mode) {}

    File::~File() { close(); }

    File::File(File&& other) noexcept
        : id_(std::exchange(other.id_, invalid_id)) {}

    File& File::operator=(File&& other) noexcept {
        if(this != &other) {
            close();
            id_ = std::exchange(other.id_, invalid_id);
        }
        return *this;
    }

    void File::close() noexcept {
        if(id_ >= 0) H5Fclose(id_);
        id_ = invalid_id;
    }

    bool File::flush() const { return id_ >= 0 && H5Fflush(id_, H5F_SCOPE_GLOBAL) >= 0; }

    Group File::root() const { return Group{open_group_id(id_, "/")}; }

    Group File::create_group(const std::string& name) const { return Group{create_group_id(id_, name)}; }

    Group File::open_group(const std::string& name) const { return Group{open_group_id(id_, name)}; }

    Dataset File::create_dataset(const std::string& name, const double* data, const std::initializer_list<hsize_t> shape) const { return Dataset{create_dataset_id(id_, name, data, std::vector<hsize_t>{shape})}; }

    Dataset File::create_dataset(const std::string& name, const double* data, const std::vector<hsize_t>& shape) const { return Dataset{create_dataset_id(id_, name, data, shape)}; }

    Dataset File::create_dataset(const std::string& name, const std::vector<double>& data) const { return Dataset{create_dataset_id(id_, name, data.data(), {data.size()})}; }

    Dataset File::open_dataset(const std::string& name) const { return Dataset{open_dataset_id(id_, name)}; }

    bool File::exists(const std::string& name) const { return exists_impl(id_, name); }

    std::vector<std::string> File::links(const std::string& group_name) const { return links_impl(id_, group_name); }

    bool File::copy(const std::string& source, const File& target, const std::string& destination) const { return copy_impl(id_, source, target.id(), destination); }

    bool File::merge(const File& source) const { return merge_impl(source.id(), "/", id_, "/"); }
} // namespace hdf5

#endif
