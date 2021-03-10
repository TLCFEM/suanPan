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
 * @class ConditionalModifier
 * @brief A ConditionalModifier class.
 *
 * The ConditionalModifier class.
 *
 * @author tlc
 * @date 07/03/2021
 * @version 0.1.0
 * @file ConditionalModifier.h
 * @addtogroup ConditionalModifier
 * @{
 */

#ifndef CONDITIONALMODIFIER_H
#define CONDITIONALMODIFIER_H

#include <Domain/Tag.h>

class DomainBase;
class Amplitude;

class ConditionalModifier : public Tag {
protected:
	const bool initialized = false;

	const bool connected = false;

	unsigned start_step, end_step = static_cast<unsigned>(-1);

	const unsigned amplitude_tag;

	uvec node_encoding; /**< node encoding */

	shared_ptr<Amplitude> magnitude;
public:
	ConditionalModifier(unsigned, unsigned, unsigned, uvec&&);

	virtual int initialize(const shared_ptr<DomainBase>&);

	virtual int process(const shared_ptr<DomainBase>&) = 0;
	virtual int process_resistance(const shared_ptr<DomainBase>&);

	[[nodiscard]] const uvec& get_node_encoding() const;

	void set_initialized(bool) const;
	[[nodiscard]] bool is_initialized() const;

	void set_start_step(unsigned);
	[[nodiscard]] unsigned get_start_step() const;

	void set_end_step(unsigned);
	[[nodiscard]] unsigned get_end_step() const;

	void set_connected(bool) const;
	[[nodiscard]] bool is_connected() const;

	[[nodiscard]] bool validate_step(const shared_ptr<DomainBase>&) const;

	// some may manage state
	virtual void update_status(const vec&);
	virtual void commit_status();
	virtual void clear_status();
	virtual void reset_status();
};

#endif

//! @}
