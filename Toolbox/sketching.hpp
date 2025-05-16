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

#ifndef SKETCHING_HPP
#define SKETCHING_HPP

#include <suanPan.h>

template<typename T> concept row_sketching = requires(T t, const uword col) {
    *t;
    t->n_rows;
    t->n_cols;
    t->head_rows(col);
    t->tail_rows(col);
    t->rows(col, col);
};

template<row_sketching T> mat frequent_row_directions(T&& target, unsigned l) {
    l = std::max(1u, std::min(l, static_cast<unsigned>(target->n_rows / 2)));

    auto current = 2 * l; // how many rows processed

    mat sketch = target->head_rows(current), V;

    const auto process = [&] {
        mat U;
        vec s;
        svd_econ(U, s, V, sketch, "right");
        for(auto I = 0u; I < l; ++I) V.col(I) *= std::sqrt((s(I) + s(l)) * (s(I) - s(l)));
    };

    while(current < target->n_rows) {
        process();

        const auto next = current + l;
        if(next > target->n_rows) {
            const auto remaining = target->n_rows - current;
            sketch.set_size(l + remaining, target->n_cols);
            sketch.tail_rows(remaining) = target->tail_rows(remaining);
        }
        else sketch.tail_rows(l) = target->rows(current, next - 1);
        current = next;

        sketch.head_rows(l) = V.head_cols(l).t();
    }

    process();

    return V.head_cols(l).t();
}

template<typename T> concept col_sketching = requires(T t, const uword col) {
    *t;
    t->n_rows;
    t->n_cols;
    t->head_cols(col);
    t->tail_cols(col);
    t->cols(col, col);
};

template<col_sketching T> mat frequent_col_directions(T&& target, unsigned l) {
    l = std::max(1u, std::min(l, static_cast<unsigned>(target->n_cols / 2)));

    auto current = 2 * l; // how many cols processed

    mat sketch = target->head_cols(current), U;

    const auto process = [&] {
        mat V;
        vec s;
        svd_econ(U, s, V, sketch, "left");
        for(auto I = 0u; I < l; ++I) U.col(I) *= std::sqrt((s(I) + s(l)) * (s(I) - s(l)));
    };

    while(current < target->n_cols) {
        process();

        const auto next = current + l;
        if(next > target->n_cols) {
            const auto remaining = target->n_cols - current;
            sketch.set_size(target->n_rows, l + remaining);
            sketch.tail_cols(remaining) = target->tail_cols(remaining);
        }
        else sketch.tail_cols(l) = target->cols(current, next - 1);
        current = next;

        sketch.head_cols(l) = U.head_cols(l);
    }

    process();

    return U.head_cols(l);
}

template<typename T> concept row_sketching_obj = row_sketching<std::remove_cvref_t<T>*>;
template<typename T> concept col_sketching_obj = col_sketching<std::remove_cvref_t<T>*>;

template<row_sketching_obj T> mat frequent_row_directions(T&& target, unsigned l) { return frequent_row_directions(&target, l); }
template<col_sketching_obj T> mat frequent_col_directions(T&& target, unsigned l) { return frequent_col_directions(&target, l); }

#endif

//! @}
