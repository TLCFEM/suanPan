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
 * @brief A registry that maps command tokens to values or factory functions.
 *
 * Tokens are compared case-insensitively against the whole input, matching the
 * behaviour of is_equal().
 *
 * @file TokenRegistry.h
 * @addtogroup Utility
 * @{
 */

#ifndef TOKENREGISTRY_H
#define TOKENREGISTRY_H

#include <Toolbox/utility.h>

/**
 * @brief An immutable registry of values keyed by command token.
 *
 * The registry is normally instantiated as a function local static so that the
 * construction is thread safe and no static initialisation order issue exists.
 *
 * @tparam V value type, which may also be a factory function pointer
 */
template<typename V> class TokenRegistry {
    std::unordered_map<std::string, V> pool;

public:
    TokenRegistry(std::initializer_list<std::pair<std::string_view, V>> entries) {
        pool.reserve(entries.size());
        for(const auto& [token, value] : entries) pool.emplace(suanpan::to_lower(token), value);
    }

    [[nodiscard]] V fetch(const std::string_view token) const {
        if(const auto iter = pool.find(suanpan::to_lower(token)); iter != pool.cend()) return iter->second;
        return V{};
    }

    [[nodiscard]] auto size() const { return pool.size(); }
};

#endif

//! @}
