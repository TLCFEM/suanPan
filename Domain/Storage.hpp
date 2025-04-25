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
 * @class Storage
 * @brief A candidate Storage class that stores FEM objects.
 * @author T
 * @date 12/20/2018
 * @version 0.3.0
 * @file Storage.hpp
 * @{
 */

#ifndef STORAGE_H
#define STORAGE_H

#include <Toolbox/container.h>

class Amplitude;
class Expression;
class Constraint;
class Converger;
class Criterion;
class Database;
class Domain;
class DomainBase;
class Element;
class Group;
class Integrator;
class Load;
class Material;
class Modifier;
class Node;
class Orientation;
class Recorder;
class Section;
class Solver;

template<typename T> const char* StorageType() { return "Unknown"; }

template<> inline const char* StorageType<Amplitude>() { return "Amplitude"; }

template<> inline const char* StorageType<Expression>() { return "Expression"; }

template<> inline const char* StorageType<Constraint>() { return "Constraint"; }

template<> inline const char* StorageType<Converger>() { return "Converger"; }

template<> inline const char* StorageType<Criterion>() { return "Criterion"; }

template<> inline const char* StorageType<Database>() { return "Database"; }

template<> inline const char* StorageType<Domain>() { return "Domain"; }

template<> inline const char* StorageType<DomainBase>() { return "Domain"; }

template<> inline const char* StorageType<Element>() { return "Element"; }

template<> inline const char* StorageType<Group>() { return "Group"; }

template<> inline const char* StorageType<Integrator>() { return "Integrator"; }

template<> inline const char* StorageType<Load>() { return "Load"; }

template<> inline const char* StorageType<Material>() { return "Material"; }

template<> inline const char* StorageType<Modifier>() { return "Modifier"; }

template<> inline const char* StorageType<Node>() { return "Node"; }

template<> inline const char* StorageType<Orientation>() { return "Orientation"; }

template<> inline const char* StorageType<Recorder>() { return "Recorder"; }

template<> inline const char* StorageType<Section>() { return "Section"; }

template<> inline const char* StorageType<Solver>() { return "Solver"; }

template<typename T> class Storage : public std::enable_shared_from_this<Storage<T>> {
    const char* type = StorageType<object_type>();

    std::vector<shared_ptr<T>> fish; /**< data storage */

    const shared_ptr<T> empty = nullptr;

    using const_iterator = typename suanpan::unordered_map<unsigned, shared_ptr<T>>::const_iterator;
    using iterator = typename suanpan::unordered_map<unsigned, shared_ptr<T>>::iterator;

    suanpan::unordered_set<unsigned> bait;                /**< data storage */
    suanpan::unordered_map<unsigned, shared_ptr<T>> pond; /**< data storage */

public:
    using object_type = T;

    Storage() = default;
    Storage(const Storage&) = delete;
    Storage(Storage&&) noexcept = delete;
    Storage& operator=(const Storage&) = delete;
    Storage& operator=(Storage&&) noexcept = delete;
    ~Storage() = default;

    const_iterator cbegin() const;
    const_iterator cend() const;
    iterator begin();
    iterator end();

    bool insert(const shared_ptr<T>&);
    shared_ptr<T>& operator[](unsigned);
    const shared_ptr<T>& at(unsigned) const;

    const std::vector<shared_ptr<T>>& get() const;

    [[nodiscard]] bool find(unsigned) const;
    bool erase(unsigned);
    void enable(unsigned);
    void disable(unsigned);

    void update();
    void enable();
    void reset();
    void clear();

    [[nodiscard]] size_t size() const;
};

template<typename T> typename Storage<T>::const_iterator Storage<T>::cbegin() const { return pond.cbegin(); }

template<typename T> typename Storage<T>::const_iterator Storage<T>::cend() const { return pond.cend(); }

template<typename T> typename Storage<T>::iterator Storage<T>::begin() { return pond.begin(); }

template<typename T> typename Storage<T>::iterator Storage<T>::end() { return pond.end(); }

template<typename T> bool Storage<T>::insert(const shared_ptr<T>& I) {
    auto flag = pond.insert({I->get_tag(), I}).second;
    if(!flag)
        suanpan_warning("Fail to insert {} {}.\n", type, I->get_tag());
    return flag;
}

template<typename T> shared_ptr<T>& Storage<T>::operator[](const unsigned L) { return pond[L]; }

template<typename T> const shared_ptr<T>& Storage<T>::at(const unsigned L) const { return pond.contains(L) ? pond.at(L) : empty; }

template<typename T> const std::vector<shared_ptr<T>>& Storage<T>::get() const { return fish; }

template<typename T> bool Storage<T>::find(const unsigned L) const { return pond.contains(L); }

template<typename T> bool Storage<T>::erase(const unsigned L) {
#ifdef SUANPAN_MT
    return pond.unsafe_erase(L) == 1;
#else
    return pond.erase(L) == 1;
#endif
}

template<typename T> void Storage<T>::enable(const unsigned L) {
    if(find(L)) pond[L]->enable();
}

template<typename T> void Storage<T>::disable(const unsigned L) {
    if(find(L)) pond[L]->disable();
}

template<typename T> void Storage<T>::update() {
    reset();
    fish.reserve(size());
    for(const auto& [tag, obj] : pond)
        if(obj->is_active()) fish.push_back(obj);
        else bait.insert(tag);
}

template<typename T> void Storage<T>::enable() {
    for(const auto& I : pond) I.second->enable();
}

template<typename T> void Storage<T>::reset() {
    fish.clear();
    bait.clear();
}

template<typename T> void Storage<T>::clear() {
    pond.clear();
    reset();
}

template<typename T> size_t Storage<T>::size() const { return pond.size(); }

template<typename T> typename Storage<T>::const_iterator cbegin(const Storage<T>& S) { return S.cbegin(); }

template<typename T> typename Storage<T>::const_iterator cend(const Storage<T>& S) { return S.cend(); }

template<typename T> typename Storage<T>::iterator begin(Storage<T>& S) { return S.begin(); }

template<typename T> typename Storage<T>::iterator end(Storage<T>& S) { return S.end(); }

template<typename T> using dual = std::pair<unsigned, shared_ptr<T>>;

using AmplitudeStorage = Storage<Amplitude>;
using ExpressionStorage = Storage<Expression>;
using ConstraintStorage = Storage<Constraint>;
using ConvergerStorage = Storage<Converger>;
using CriterionStorage = Storage<Criterion>;
using DatabaseStorage = Storage<Database>;
using DomainStorage = Storage<Domain>;
using DomainBaseStorage = Storage<DomainBase>;
using ElementStorage = Storage<Element>;
using GroupStorage = Storage<Group>;
using IntegratorStorage = Storage<Integrator>;
using LoadStorage = Storage<Load>;
using MaterialStorage = Storage<Material>;
using ModifierStorage = Storage<Modifier>;
using NodeStorage = Storage<Node>;
using OrientationStorage = Storage<Orientation>;
using RecorderStorage = Storage<Recorder>;
using SectionStorage = Storage<Section>;
using SolverStorage = Storage<Solver>;

#endif

//! @}
