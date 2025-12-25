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
 * @author tlc
 * @date 01/04/2018
 * @version 0.2.0
 * @file Tag.h
 */

#ifndef TAG_H
#define TAG_H

// ReSharper disable once CppUnusedIncludeDirective
#include <suanPan.h> // for derived classes

/**
 * @class Tag
 * @brief A base Tag class.
 *
 * The `Tag` class is a base class which stores the object's `unique_tag` and
 * `class_tag`. Additionally, the Tag class defines status of an object, which
 * is stored in variable `alive`. By testing its value, the object can be
 * removed or added to the global system.
 *
 */
class Tag {
    bool alive = true; // status flag
    bool guarded = false;

    const unsigned unique_tag; // unique tag of the object
public:
    explicit Tag(unsigned = 0);
    Tag(const Tag&) = default;
    Tag(Tag&&) noexcept = default;
    Tag& operator=(const Tag&) = delete;
    Tag& operator=(Tag&&) = delete;
    virtual ~Tag() = default;

    void set_tag(unsigned) const;
    [[nodiscard]] unsigned get_tag() const;

    void enable();
    void disable();

    void guard();
    void unguard();

    [[nodiscard]] bool is_active() const;
    [[nodiscard]] bool is_guarded() const;

    virtual void print();
};

/**
 * @class CopyableTag
 * @brief Label objects that can be copied.
 *
 */
class CopyableTag : public Tag {
public:
    using Tag::Tag;

    CopyableTag(const CopyableTag&) = default;
    CopyableTag(CopyableTag&&) = default;
    CopyableTag& operator=(const CopyableTag&) = delete;
    CopyableTag& operator=(CopyableTag&&) = delete;
    ~CopyableTag() override = default;
};

/**
 * @class UniqueTag
 * @brief Label objects that cannot be copied.
 */
class UniqueTag : public Tag {
public:
    using Tag::Tag;

    UniqueTag(const UniqueTag&) = delete;
    UniqueTag(UniqueTag&&) = delete;
    UniqueTag& operator=(const UniqueTag&) = delete;
    UniqueTag& operator=(UniqueTag&&) = delete;
    ~UniqueTag() override = default;
};

#endif
