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
 * @fn command
 * @brief Parsing commands.
 *
 * @author tlc
 * @date 05/01/2022
 * @version 0.1.0
 * @file command.h
 * @addtogroup Utility
 * @{
 */

#ifndef COMMAND_H
#define COMMAND_H

#include <suanPan.h>

class Bead;

int process_command(const std::shared_ptr<Bead>&, std::istringstream&);

bool normalise_command(std::string&, std::string&);

int process_file(const std::shared_ptr<Bead>&, const char*);

int execute_command(std::istringstream&);

fs::path get_history_path();

#endif

//! @}
