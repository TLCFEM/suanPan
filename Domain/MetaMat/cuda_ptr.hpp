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
 * @class cuda_ptr
 * @author tlc
 * @date 17/05/2025
 * @version 0.1.0
 * @file cuda_ptr.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef CUDA_PTR_HPP
#define CUDA_PTR_HPP

#include <cusolverSp.h>
#include <cusparse.h>

class cuda_ptr {
    void* ptr{};

    size_t total_size() const { return unit * size; }

public:
    size_t unit{}, size{};

    explicit cuda_ptr(const size_t in_unit = 0, const size_t in_size = 0)
        : unit(in_unit)
        , size(in_size) {
        if(total_size() > 0) cudaMalloc(&ptr, total_size());
    }
    cuda_ptr(const cuda_ptr& other)
        : cuda_ptr(other.unit, other.size) {}
    cuda_ptr(cuda_ptr&&) = delete;
    cuda_ptr& operator=(const cuda_ptr&) = delete;
    cuda_ptr& operator=(cuda_ptr&& other) noexcept {
        if(this != &other) {
            cudaFree(ptr);
            ptr = other.ptr;
            unit = other.unit;
            size = other.size;
            other.ptr = nullptr;
            other.unit = 0;
            other.size = 0;
        }
        return *this;
    }
    ~cuda_ptr() { cudaFree(ptr); }

    auto operator&() { return ptr; }

    auto copy_from(const void* src, cudaStream_t s) { cudaMemcpyAsync(ptr, src, total_size(), cudaMemcpyHostToDevice, s); }

    auto copy_to(void* dest, cudaStream_t s) {
        cudaMemcpyAsync(dest, ptr, total_size(), cudaMemcpyDeviceToHost, s);
        cudaStreamSynchronize(s);
    }
};

#endif

//! @}
