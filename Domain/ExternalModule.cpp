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

#ifdef SUANPAN_WIN
#include <Windows.h>
#else
#include <dlfcn.h>
#endif

bool ExternalModule::locate_module(std::string module_name) {
    if(ext_library == nullptr) return false;

    suanpan::to_lower(module_name);

#ifdef SUANPAN_WIN
    ext_creator = reinterpret_cast<void*>(GetProcAddress(HINSTANCE(ext_library), LPCSTR(module_name.c_str())));
#else
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
#else
#ifdef SUANPAN_UNIX
    static constexpr auto suffix = ".so";
#else
    static constexpr auto suffix = ".dylib";
#endif
    auto file_name = "lib" + library_name + suffix;
    ext_library = dlopen(file_name.c_str(), RTLD_NOW);
    if(ext_library == nullptr) {
        file_name = library_name + suffix;
        ext_library = dlopen(file_name.c_str(), RTLD_NOW);
    }
    if(ext_library == nullptr) suanpan_error("Cannot load the library with the name \"{}\".\n", file_name);
#endif
}

ExternalModule::~ExternalModule() {
#ifdef SUANPAN_WIN
    if(ext_library != nullptr) FreeLibrary(HINSTANCE(ext_library));
#else
    if(ext_library != nullptr) dlclose(ext_library);
#endif
}

bool ExternalModule::locate_c_module(const std::string& module_name) { return locate_module(module_name + "_handler"); }

bool ExternalModule::locate_cpp_module(const std::string& module_name) { return locate_module("new_" + module_name); }

void ExternalModule::new_adapter(unique_ptr<Element>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Load>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Material>& return_obj, std::istringstream& command) const {
    unsigned tag;

    if(!get_input(command, tag)) {
        suanpan_error("A valid tag is required.\n");
        return;
    }

    if(auto ext_obj = std::make_unique<ExternalMaterial>(tag, get_remaining<double>(command), ext_creator); ext_obj->validate()) return_obj = std::move(ext_obj);
}

void ExternalModule::new_adapter(unique_ptr<Section>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Solver>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Amplitude>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Modifier>&, std::istringstream&) const {}

void ExternalModule::new_adapter(unique_ptr<Constraint>&, std::istringstream&) const {}
