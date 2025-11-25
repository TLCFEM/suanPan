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

template<typename T> requires requires(T* copyable) { copyable->get_copy(); }
class ResourceHolder final {
    std::unique_ptr<T> object = nullptr;

public:
    ResourceHolder() = default;

    ResourceHolder(const ResourceHolder& old)
        : object(old.object ? old.object->get_copy() : nullptr) {}

    ResourceHolder(ResourceHolder&& old) noexcept { object = std::move(old.object); }

    ResourceHolder& operator=(const ResourceHolder&) = delete;
    ResourceHolder& operator=(ResourceHolder&&) = delete;
    ~ResourceHolder() = default;

    explicit ResourceHolder(std::unique_ptr<T>&& old)
        : object(std::move(old)) {}

    template<typename U> requires std::is_base_of_v<T, U> ResourceHolder& operator=(U&& other) {
        object = std::make_unique<U>(std::forward<U>(other));
        return *this;
    }

    ResourceHolder& operator=(const std::shared_ptr<T>& old) {
        if(old) object = old->get_copy();
        return *this;
    }

    T* operator->() const { return object.get(); }

    explicit operator bool() const { return object != nullptr; }

    bool operator==(std::nullptr_t null) const { return object == null; }
};

#endif

//! @}
