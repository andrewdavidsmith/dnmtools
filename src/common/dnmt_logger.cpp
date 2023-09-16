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
#include <vector>
#include <string>
#include <ctime>

using std::ostream;
using std::string;
using std::vector;

using std::chrono::duration;
using std::chrono::steady_clock;
using std::chrono::time_point;
using std::gmtime;
using std::timespec;
using std::timespec_get;
using std::strftime;

static inline string
format_current_time() {
  static constexpr auto time_template = "[%D %T]";
  timespec ts;
  timespec_get(&ts, TIME_UTC);
  vector<char> buf(128);
  strftime(buf.data(), std::size(buf), time_template, std::gmtime(&ts.tv_sec));
  return buf.data();
}

time_point<steady_clock>
dnmt_logger::log_event(string message) {
  static constexpr auto event_tag = "[event]";

  log_stream << prefix << ' '
             << format_current_time() << ' '
             << event_tag << ' '
             << message << '\n';

  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_event(string message, time_point<steady_clock> then) {
  static constexpr auto event_tag = "[event]";
  static constexpr auto delta_time_template = "[delta time: %.3fs]";

  const auto now = steady_clock::now();
  const auto delta_time = duration<double>(now - then).count();

  string formatted_delta_time(128, '\0');
  std::snprintf(formatted_delta_time.data(), std::size(formatted_delta_time),
                delta_time_template, delta_time);
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << event_tag << ' '
             << message << ' '
             << formatted_delta_time << '\n';
  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_data(string key, string value) {
  static constexpr auto data_tag = "[data]";
  // clang-format off
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << data_tag << ' '
             << key << '=' << value
             << '\n';
  // clang-format on
  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_data(string key, string value, time_point<steady_clock> then) {
  static constexpr auto data_tag = "[data]";
  static constexpr auto delta_time_template = "[delta time: %.3fs]";

  const auto now = steady_clock::now();
  const auto delta_time = duration<double>(now - then).count();

  string formatted_delta_time(128, '\0');
  std::snprintf(formatted_delta_time.data(), std::size(formatted_delta_time),
                delta_time_template, delta_time);
  // clang-format off
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << data_tag << ' '
             << key << '=' << value << ' '
             << formatted_delta_time << '\n';
  // clang-format on
  return steady_clock::now();
}
