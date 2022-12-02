/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
#include <hdf5/hdf5.h>
#include <hdf5/hdf5_hl.h>
#endif

/**
 * \brief ctor
 * \param T unique tag
 * \param B object tags
 * \param L variable type
 * \param I record interval
 * \param R if to record time stamp
 * \param H if to use hdf5 format
 */
Recorder::Recorder(const unsigned T, uvec&& B, const OutputType L, const unsigned I, const bool R, const bool H)
    : Tag(T)
    , object_tag(std::forward<uvec>(B))
    , variable_type(L)
    , data_pool(object_tag.n_elem)
    , record_time(R)
    , use_hdf5(H)
    , interval(I) { suanpan_debug("Recorder %u ctor() called.\n", T); }

Recorder::~Recorder() { suanpan_debug("Recorder %u dtor() called.\n", get_tag()); }

void Recorder::initialize(const shared_ptr<DomainBase>&) {}

void Recorder::set_object_tag(const uvec& T) { object_tag = T; }

const uvec& Recorder::get_object_tag() const { return object_tag; }

void Recorder::set_variable_type(const OutputType T) { variable_type = T; }

const OutputType& Recorder::get_variable_type() const { return variable_type; }

bool Recorder::if_hdf5() const { return use_hdf5; }

bool Recorder::if_record_time() const { return record_time; }

bool Recorder::if_perform_record() { return 1 == interval || 0 == std::remainder(counter++, interval); }

void Recorder::insert(const double T) { time_pool.emplace_back(T); }

void Recorder::insert(const std::vector<vec>& D, const unsigned I) { data_pool[I].emplace_back(D); }

const std::vector<std::vector<std::vector<vec>>>& Recorder::get_data_pool() const { return data_pool; }

const std::vector<double>& Recorder::get_time_pool() const { return time_pool; }

void Recorder::save() {
    if(time_pool.empty() || data_pool.empty() || data_pool.cbegin()->empty() || data_pool.cbegin()->cbegin()->empty() || data_pool.cbegin()->cbegin()->cbegin()->is_empty()) return;

    ostringstream file_name;
    file_name << 'R' << get_tag() << '-' << to_char(variable_type);
    const auto origin_name = file_name.str();

    unsigned idx = 0;

#ifdef SUANPAN_HDF5
    if(use_hdf5) {
        file_name << ".h5";

        const auto file_id = H5Fcreate((SUANPAN_OUTPUT / file_name.str()).generic_string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        string group_name = "/";
        group_name += origin_name;

        const auto group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for(const auto& s_data_pool : data_pool) {
            auto max_size = 0llu;
            for(const auto& I : s_data_pool[0]) if(I.n_elem > max_size) max_size = I.n_elem;

            mat data_to_write(s_data_pool.cbegin()->size() * max_size + 1, time_pool.size(), fill::zeros);

            for(size_t I = 0; I < time_pool.size(); ++I) {
                data_to_write(0, I) = time_pool[I];
                unsigned L = 1;
                for(const auto& J : s_data_pool[I]) for(unsigned K = 0; K < J.n_elem; ++K) data_to_write(L++, I) = J[K];
            }

            hsize_t dimension[2] = {data_to_write.n_cols, data_to_write.n_rows};

            ostringstream dataset_name;
            dataset_name << origin_name.c_str();
            dataset_name << object_tag(idx++);

            H5LTmake_dataset(group_id, dataset_name.str().c_str(), 2, dimension, H5T_NATIVE_DOUBLE, data_to_write.mem);
        }

        H5Gclose(group_id);
        H5Fclose(file_id);
    }
    else {
        for(const auto& s_data_pool : data_pool) {
            auto max_size = 0llu;
            for(const auto& I : s_data_pool[0]) if(I.n_elem > max_size) max_size = I.n_elem;

            mat data_to_write(s_data_pool.cbegin()->size() * max_size + 1, time_pool.size(), fill::zeros);

            for(size_t I = 0; I < time_pool.size(); ++I) {
                data_to_write(0, I) = time_pool[I];
                unsigned L = 1;
                for(const auto& J : s_data_pool[I]) for(unsigned K = 0; K < J.n_elem; ++K) data_to_write(L++, I) = J[K];
            }

            ostringstream dataset_name;
            dataset_name << (SUANPAN_OUTPUT / origin_name).generic_string();
            dataset_name << object_tag(idx++);

            mat(data_to_write.t()).save(dataset_name.str() + ".txt", raw_ascii);
        }
    }
#else
    for(const auto& s_data_pool : data_pool) {
        auto max_size = 0llu;
        for(const auto& I : s_data_pool[0]) if(I.n_elem > max_size) max_size = I.n_elem;

        mat data_to_write(s_data_pool.cbegin()->size() * max_size + 1, time_pool.size(), fill::zeros);

        for(size_t I = 0; I < time_pool.size(); ++I) {
            data_to_write(0, I) = time_pool[I];
            auto L = 1;
            for(const auto& J : s_data_pool[I]) for(unsigned K = 0; K < J.n_elem; ++K) data_to_write(L++, I) = J[K];
        }

        ostringstream dataset_name;
        dataset_name << (SUANPAN_OUTPUT / origin_name).generic_string();
        dataset_name << object_tag(idx++);

        mat(data_to_write.t()).save(dataset_name.str() + ".txt", raw_ascii);
    }
#endif
}

void Recorder::print() { suanpan_info("print() needs to be overwritten.\n"); }
