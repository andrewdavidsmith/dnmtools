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

#include <iomanip>
using std::ostream;
using std::string;

using std::chrono::duration;
using std::chrono::steady_clock;
using std::chrono::time_point;

using std::endl;

time_point<steady_clock>
dnmt_logger::log_event(string message) {
  log_stream << prefix << message << endl;
  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_event(string message, time_point<steady_clock> then) {
  const auto now = steady_clock::now();
  log_stream << prefix << message << " "
             << "[delta time: " << std::fixed << std::setprecision(2)
             << duration<double>(now - then).count() << "s]" << endl;
  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_data(string key, string value) {
  log_stream << prefix << key << "=" << value << endl;
  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_data(string key, string value, time_point<steady_clock> then) {
  const auto now = steady_clock::now();
  log_stream << prefix << key << "=" << value << " "
             << "[delta time: " << std::fixed << std::setprecision(2)
             << duration<double>(now - then).count() << "s]" << endl;
  return steady_clock::now();
}
