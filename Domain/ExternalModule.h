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
 * @class ExternalModule
 * @brief A ExternalModule class handles communication between the main program
 * and external library.
 *
 * @author tlc
 * @date 28/09/2017
 * @version 0.1.1
 * @file ExternalModule.h
 * @addtogroup Utility
 * @{
 */

#ifndef EXTERNALMODULE_H
#define EXTERNALMODULE_H

#include <Domain/DomainBase.h>
#include <Toolbox/utility.h>

class Element;
class Load;
class Material;
class Section;
class Solver;
class Amplitude;
class Modifier;
class Constraint;

class ExternalModule {
    void* ext_library = nullptr;
    void* ext_creator = nullptr;

    bool locate_module(std::string);

public:
    const std::string library_name;

    explicit ExternalModule(std::string);
    ExternalModule(const ExternalModule&) = delete;
    ExternalModule(ExternalModule&&) = delete;
    ExternalModule& operator=(const ExternalModule&) = delete;
    ExternalModule& operator=(ExternalModule&&) = delete;
    ~ExternalModule();

    bool locate_c_module(const std::string&);
    bool locate_cpp_module(const std::string&);

    void new_object(unique_ptr<Element>&, std::istringstream&) const;
    void new_object(unique_ptr<Load>&, std::istringstream&) const;
    void new_object(unique_ptr<Material>&, std::istringstream&) const;
    void new_object(unique_ptr<Section>&, std::istringstream&) const;
    void new_object(unique_ptr<Solver>&, std::istringstream&) const;
    void new_object(unique_ptr<Amplitude>&, std::istringstream&) const;
    void new_object(unique_ptr<Modifier>&, std::istringstream&) const;
    void new_object(unique_ptr<Constraint>&, std::istringstream&) const;

    void new_adapter(unique_ptr<Element>&, std::istringstream&) const;
    void new_adapter(unique_ptr<Load>&, std::istringstream&) const;
    void new_adapter(unique_ptr<Material>&, std::istringstream&) const;
    void new_adapter(unique_ptr<Section>&, std::istringstream&) const;
    void new_adapter(unique_ptr<Solver>&, std::istringstream&) const;
    void new_adapter(unique_ptr<Amplitude>&, std::istringstream&) const;
    void new_adapter(unique_ptr<Modifier>&, std::istringstream&) const;
    void new_adapter(unique_ptr<Constraint>&, std::istringstream&) const;
};

class load {
public:
    template<typename T> static void object(unique_ptr<T>&, const shared_ptr<DomainBase>&, const std::string&, std::istringstream&);
};

template<typename T> void load::object(unique_ptr<T>& new_object, const shared_ptr<DomainBase>& domain, const std::string& id, std::istringstream& command) {
    // check if the library is already loaded
    auto loaded = false;
    for(const auto& I : domain->get_external_module_pool())
        if(is_equal(I->library_name, id) || I->locate_cpp_module(id) || I->locate_c_module(id)) {
            loaded = true;
            break;
        }

    // not loaded then try load it
    // if loaded find corresponding function
    if(loaded || domain->insert(std::make_shared<ExternalModule>(id)))
        for(const auto& I : domain->get_external_module_pool()) {
            if(I->locate_cpp_module(id)) I->new_object(new_object, command);
            if(new_object != nullptr) break;
            if(I->locate_c_module(id)) I->new_adapter(new_object, command);
            if(new_object != nullptr) break;
        }
}

#endif
