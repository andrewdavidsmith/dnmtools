/* Copyright (C) 2023 Andrew D. Smith
 *
 * Authors: Andrew Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "dnmt_logger.hpp"

using std::chrono;
using std::string;
using std::ostream;

using std::endl;

chrono::steady_clock
log_event(string message) {
  out << prefix << message << endl;
  return chrono::steady_clock::now();
}

chrono::steady_clock
log_event(string message, chrono::steady_clock then) {
  const auto now = steady_clock::now();
  out << prefix << message
    // oss.setf(std::ios::fixed);
    // oss.precision(2);
      << "[delta time: "
      << duration<double>(then - now).count() << "s]"
      << endl;
  return chrono::steady_clock::now();
}

chrono::steady_clock
log_data(string key, string value) {
  out << prefix << key << "=" << value
      << endl;
  return chrono::steady_clock::now();
}

chrono::steady_clock
log_data(string key, string value, chrono::steady_clock sc) {
  out << prefix << key << "=" << value
      << endl;
  return chrono::steady_clock::now();
}
