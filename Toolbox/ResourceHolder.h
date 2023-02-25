/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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
 * @class ResourceHolder
 * @author tlc
 * @date 15/02/2023
 * @file ResourceHolder.h
 * @addtogroup Utility
 * @{
 */

#ifndef RESOURCEHOLDER_H
#define RESOURCEHOLDER_H

#include <memory>

template<typename T> class ResourceHolder final {
    std::unique_ptr<T> object = nullptr;

public:
    ResourceHolder() = default;

    ResourceHolder& operator=(const std::shared_ptr<T>& original_object) {
        object = original_object->make_copy();
        return *this;
    }

    ResourceHolder(const ResourceHolder& old_holder)
        : object(old_holder.object ? old_holder.object->make_copy() : nullptr) {}

    ResourceHolder(ResourceHolder&&) noexcept = delete;
    ResourceHolder& operator=(const ResourceHolder&) = delete;
    ResourceHolder& operator=(ResourceHolder&&) noexcept = delete;
    ~ResourceHolder() = default;

    T* operator->() const { return object.get(); }

    explicit operator bool() const { return object != nullptr; }
};

#endif

//! @}
