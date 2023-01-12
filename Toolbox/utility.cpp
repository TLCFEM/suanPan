/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "utility.h"
#include <suanPan.h>

unsigned long long suanpan::binomial(unsigned long long n, unsigned long long k) {
    auto c = 1llu;

    if(k > n - k) k = n - k;

    for(auto i = 1llu; i <= k; i++, n--) {
        if(c / i > std::numeric_limits<unsigned long long>::max() / n) return 0llu;

        c = c / i * n + c % i * n / i;
    }

    return c;
}

char suanpan::to_upper(const char U) { return static_cast<char>(std::toupper(static_cast<int>(U))); }

char suanpan::to_lower(const char L) { return static_cast<char>(std::tolower(static_cast<int>(L))); }

void suanpan::to_upper(string& U) { std::ranges::for_each(U, [](char& C) { C = static_cast<char>(std::toupper(static_cast<int>(C))); }); }

void suanpan::to_lower(string& U) { std::ranges::for_each(U, [](char& C) { C = static_cast<char>(std::tolower(static_cast<int>(C))); }); }

string suanpan::to_upper(const string& U) {
    string C(U);

    to_upper(C);

    return C;
}

string suanpan::to_lower(const string& U) {
    string C(U);

    to_lower(C);

    return C;
}

string suanpan::to_upper(string&& U) {
    to_upper(U);
    return std::move(U);
}

string suanpan::to_lower(string&& U) {
    to_lower(U);
    return std::move(U);
}

std::vector<std::string> suanpan::expression::split(const std::string& variable_string) {
    std::vector<std::string> variable_list;
    auto I = variable_string.cbegin(), J = variable_string.cbegin();
    while(I != variable_string.cend()) {
        if(',' == *I || '"' == *I) {
            if(I != J) variable_list.emplace_back(J, I);
            J = ++I;
        }
        else ++I;
    }

    if(I != J) variable_list.emplace_back(J, I);

    return variable_list;
}

void ignore_whitespace(istringstream& I) {
    while(true)
        if(const auto peek_value = I.peek(); is_equal(peek_value, '\t') || is_equal(peek_value, ' ')) I.ignore();
        else break;
}

bool is_equal(const char* A, const char* B) { return _strcmpi(A, B) == 0; }

bool is_equal(const char A, const char B) { return tolower(static_cast<int>(A)) == tolower(static_cast<int>(B)); }

bool is_equal(const int A, const char B) { return tolower(A) == tolower(static_cast<int>(B)); }

bool is_equal(const string& A, const char* B) { return is_equal(A.c_str(), B); }

bool is_equal(const string& A, const string& B) { return is_equal(A.c_str(), B.c_str()); }

bool if_contain(const string& A, const char* B) { return A.find(B) != string::npos; }

bool if_contain(const string& A, const string& B) { return A.find(B) != string::npos; }

bool if_contain(string&& A, string&& B) { return if_contain(A, B); }

bool is_true(const char* S) { return is_equal(S, "On") || is_equal(S, "True") || is_equal(S, "T") || is_equal(S, "1") || is_equal(S, "Yes") || is_equal(S, "Y"); }

bool is_false(const char* S) { return is_equal(S, "Off") || is_equal(S, "False") || is_equal(S, "F") || is_equal(S, "0") || is_equal(S, "No") || is_equal(S, "N"); }

bool is_true(const string& S) { return is_true(S.c_str()); }

bool is_false(const string& S) { return is_false(S.c_str()); }
