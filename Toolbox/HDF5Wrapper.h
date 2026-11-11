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
/**
 * @file HDF5Wrapper.h
 * @brief A small RAII wrapper around the HDF5 C API.
 *
 * The wrapper intentionally exposes only the subset of functionality used by
 * the project: creating and opening files, creating and opening groups,
 * creating and reading double precision datasets, listing links, checking
 * object existence and merging HDF5 files.
 *
 * All handles are closed automatically. Objects are move-only, which prevents
 * accidental double closing of HDF5 identifiers.
 *
 * Example:
 * @code
 * hdf5::File file("output.h5", hdf5::File::Mode::Truncate);
 * auto group = file.create_group("R1-Displacement");
 * group.create_dataset("N1", data.mem, {data.n_cols, data.n_rows});
 * @endcode
 *
 * @author tlc
 * @date 11/11/2026
 * @version 0.1.0
 * @addtogroup Utility
 * @{
 */

#ifndef HDF5WRAPPER_H
#define HDF5WRAPPER_H

#ifdef SUANPAN_HDF5

#include <cstddef>
#include <filesystem>
#include <hdf5.h>
#include <initializer_list>
#include <string>
#include <vector>

namespace hdf5 {
    class File;
    class Group;
    class Dataset;

    /**
     * @class Dataset
     * @brief An HDF5 dataset handle.
     */
    class Dataset {
        friend class File;
        friend class Group;

        hid_t id_ = -1;

        explicit Dataset(const hid_t id)
            : id_(id) {}

    public:
        Dataset() = default;
        ~Dataset();

        Dataset(Dataset&&) noexcept;
        Dataset& operator=(Dataset&&) noexcept;
        Dataset(const Dataset&) = delete;
        Dataset& operator=(const Dataset&) = delete;

        [[nodiscard]] bool valid() const noexcept { return id_ >= 0; }
        explicit operator bool() const noexcept { return valid(); }
        [[nodiscard]] hid_t id() const noexcept { return id_; }

        void close() noexcept;

        [[nodiscard]] std::vector<hsize_t> shape() const;
        [[nodiscard]] std::size_t size() const;

        [[nodiscard]] bool read(double*) const;
        [[nodiscard]] std::vector<double> read() const;
    };

    /**
     * @class Group
     * @brief An HDF5 group handle.
     */
    class Group {
        friend class File;

        hid_t id_ = -1;

        explicit Group(const hid_t id)
            : id_(id) {}

    public:
        Group() = default;
        ~Group();

        Group(Group&&) noexcept;
        Group& operator=(Group&&) noexcept;
        Group(const Group&) = delete;
        Group& operator=(const Group&) = delete;

        [[nodiscard]] static Group create(const File&, const std::string&);
        [[nodiscard]] static Group create(const Group&, const std::string&);
        [[nodiscard]] static Group open(const File&, const std::string&);
        [[nodiscard]] static Group open(const Group&, const std::string&);

        [[nodiscard]] bool valid() const noexcept { return id_ >= 0; }
        explicit operator bool() const noexcept { return valid(); }
        [[nodiscard]] hid_t id() const noexcept { return id_; }

        void close() noexcept;

        [[nodiscard]] Group create_group(const std::string&) const;
        [[nodiscard]] Group open_group(const std::string&) const;

        Dataset create_dataset(const std::string&, const double*, std::initializer_list<hsize_t>) const;
        Dataset create_dataset(const std::string&, const double*, const std::vector<hsize_t>&) const;
        Dataset create_dataset(const std::string&, const std::vector<double>&) const;
        [[nodiscard]] Dataset open_dataset(const std::string&) const;

        [[nodiscard]] bool exists(const std::string&) const;
        [[nodiscard]] std::vector<std::string> links() const;
        [[nodiscard]] std::string name() const;

        bool merge(const Group&) const;
    };

    /**
     * @class File
     * @brief An HDF5 file handle.
     */
    class File {
    public:
        enum class Mode {
            ReadOnly,
            ReadWrite,
            Truncate,
            CreateNew
        };

    private:
        hid_t id_ = -1;

        explicit File(const hid_t id)
            : id_(id) {}

    public:
        File() = default;
        explicit File(const std::filesystem::path&, Mode = Mode::Truncate);
        explicit File(const std::string&, Mode = Mode::Truncate);
        ~File();

        File(File&&) noexcept;
        File& operator=(File&&) noexcept;
        File(const File&) = delete;
        File& operator=(const File&) = delete;

        [[nodiscard]] bool valid() const noexcept { return id_ >= 0; }
        explicit operator bool() const noexcept { return valid(); }
        [[nodiscard]] hid_t id() const noexcept { return id_; }

        void close() noexcept;
        [[nodiscard]] bool flush() const;

        [[nodiscard]] Group root() const;
        [[nodiscard]] Group create_group(const std::string&) const;
        [[nodiscard]] Group open_group(const std::string&) const;

        Dataset create_dataset(const std::string&, const double*, std::initializer_list<hsize_t>) const;
        Dataset create_dataset(const std::string&, const double*, const std::vector<hsize_t>&) const;
        Dataset create_dataset(const std::string&, const std::vector<double>&) const;
        [[nodiscard]] Dataset open_dataset(const std::string&) const;

        [[nodiscard]] bool exists(const std::string&) const;
        [[nodiscard]] std::vector<std::string> links(const std::string& = "/") const;

        /**
         * @brief Copy one object (dataset or group) from this file into another file.
         *
         * @param source source path relative to this file, may be absolute.
         * @param target destination file.
         * @param destination destination path relative to the target file.
         */
        bool copy(const std::string&, File&, const std::string&) const;

        /**
         * @brief Recursively merge all objects of a source file into this file.
         *
         * Existing datasets are kept, existing groups are merged recursively.
         */
        bool merge(const File&) const;
    };
} // namespace hdf5

#endif

#endif

//! @}
