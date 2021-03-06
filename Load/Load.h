/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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
 * @class Load
 * @brief A Load class.
 *
 * The Load class is in charge of returning load level according to given time
 * increment.
 *
 * @author tlc
 * @date 01/10/2017
 * @version 0.2.0
 * @file Load.h
 * @addtogroup Load
 * @{
 */

#ifndef LOAD_H
#define LOAD_H

#include <Domain/Tag.h>

class Amplitude;
class DomainBase;

class Load : public Tag {
protected:
	static const double multiplier;

	const bool initialized = false;

	unsigned start_step, amplitude_tag;

	unsigned end_step = static_cast<unsigned>(-1);

	const uvec nodes, dofs;

	const double pattern;

	shared_ptr<Amplitude> magnitude;

	friend void set_load_multiplier(double);
public:
	explicit Load(unsigned = 0, unsigned = 0, unsigned = 0, uvec&& = {}, uvec&& = {}, double = 0.);
	~Load() override;

	virtual int initialize(const shared_ptr<DomainBase>&);

	virtual int process(const shared_ptr<DomainBase>&) = 0;

	void set_initialized(bool) const;

	void set_start_step(unsigned);
	[[nodiscard]] unsigned get_start_step() const;

	void set_end_step(unsigned);
	[[nodiscard]] unsigned get_end_step() const;

	[[nodiscard]] bool validate_step(const shared_ptr<DomainBase>&) const;

	// some loads may manage state
	virtual void commit_status();
	virtual void clear_status();
	virtual void reset_status();
};

void set_load_multiplier(double);

#endif

//! @}
