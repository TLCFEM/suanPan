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

// ReSharper disable CppFunctionalStyleCast
#include "ExternalModule.h"

#include <Material/ExternalMaterial.h>
#include <algorithm>

#if defined(SUANPAN_WIN)
#include <Windows.h>
#elif defined(SUANPAN_UNIX)
#include <dlfcn.h>
#endif

using element_creator = void (*)(unique_ptr<Element>&, std::istringstream&);
using load_creator = void (*)(unique_ptr<Load>&, std::istringstream&);
using material_creator = void (*)(unique_ptr<Material>&, std::istringstream&);
using section_creator = void (*)(unique_ptr<Section>&, std::istringstream&);
using solver_creator = void (*)(unique_ptr<Solver>&, std::istringstream&);
using amplitude_creator = void (*)(unique_ptr<Amplitude>&, std::istringstream&);
using modifier_creator = void (*)(unique_ptr<Modifier>&, std::istringstream&);
using constraint_creator = void (*)(unique_ptr<Constraint>&, std::istringstream&);

using external_handler = void (*)(ExternalMaterialData*, int*);

bool ExternalModule::locate_module(std::string module_name) {
    if(ext_library == nullptr) return false;

    suanpan::to_lower(module_name);

#ifdef SUANPAN_WIN
    ext_creator = reinterpret_cast<void*>(GetProcAddress(HINSTANCE(ext_library), LPCSTR(module_name.c_str())));
#elif defined(SUANPAN_UNIX)
    ext_creator = dlsym(ext_library, module_name.c_str());
#endif

    return ext_creator != nullptr;
}

ExternalModule::ExternalModule(std::string L)
    : library_name(std::move(L)) {
#ifdef SUANPAN_WIN
    auto file_name = library_name + ".dll";
    auto gnu_name = "lib" + file_name;

    // no prefix
    ext_library = LoadLibraryA(file_name.c_str());
    if(ext_library == nullptr) {
        suanpan::to_lower(file_name);
        ext_library = LoadLibraryA(file_name.c_str());
    }
    if(ext_library == nullptr) {
        suanpan::to_upper(file_name);
        ext_library = LoadLibraryA(file_name.c_str());
    }

    // linux prefix style
    if(ext_library == nullptr) ext_library = LoadLibraryA(gnu_name.c_str());
    if(ext_library == nullptr) {
        suanpan::to_lower(gnu_name);
        ext_library = LoadLibraryA(gnu_name.c_str());
    }
    if(ext_library == nullptr) {
        suanpan::to_upper(gnu_name);
        ext_library = LoadLibraryA(gnu_name.c_str());
    }
    if(ext_library == nullptr)
        suanpan_error("Cannot load the library with the name \"{}\".\n", file_name);
#elif defined(SUANPAN_UNIX)
    auto file_name = "./lib" + library_name + ".so";
    ext_library = dlopen(file_name.c_str(), RTLD_NOW);
    if(ext_library == nullptr) {
        file_name = "./" + library_name + ".so";
        ext_library = dlopen(file_name.c_str(), RTLD_NOW);
    }
    if(ext_library == nullptr) suanpan_error("Cannot load the library with the name \"{}\".\n", file_name);
#endif
}

ExternalModule::~ExternalModule() {
#ifdef SUANPAN_WIN
    if(ext_library != nullptr) FreeLibrary(HINSTANCE(ext_library));
#elif defined(SUANPAN_UNIX)
    if(ext_library != nullptr) dlclose(ext_library);
#endif
}

bool ExternalModule::locate_c_module(const std::string& module_name) { return locate_module(module_name + "_handler"); }

bool ExternalModule::locate_cpp_module(const std::string& module_name) { return locate_module("new_" + module_name); }

void ExternalModule::new_object(unique_ptr<Element>& return_obj, std::istringstream& command) const { (element_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Load>& return_obj, std::istringstream& command) const { (load_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Material>& return_obj, std::istringstream& command) const { (material_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Section>& return_obj, std::istringstream& command) const { (section_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Solver>& return_obj, std::istringstream& command) const { (solver_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Amplitude>& return_obj, std::istringstream& command) const { (amplitude_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Modifier>& return_obj, std::istringstream& command) const { (modifier_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Constraint>& return_obj, std::istringstream& command) const { (constraint_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_adapter(unique_ptr<Element>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Load>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Material>& return_obj, std::istringstream& command) const {
    unsigned tag;

    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    std::vector<double> pool;

    double para;

    while(!command.eof() && get_input(command, para)) pool.emplace_back(para);

    if(auto ext_obj = std::make_unique<ExternalMaterial>(tag, std::move(pool), ext_creator); ext_obj->validate()) return_obj = std::move(ext_obj);
}

void ExternalModule::new_adapter(unique_ptr<Section>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Solver>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Amplitude>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Modifier>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Constraint>&, std::istringstream&) const {}
