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

using std::endl;

static inline string
format_current_time() {
  static constexpr auto time_format = "%D %T";
  timespec ts;
  timespec_get(&ts, TIME_UTC);
  vector<char> buf(128);
  strftime(buf.data(), std::size(buf), time_format, std::gmtime(&ts.tv_sec));
  return buf.data();
}

time_point<steady_clock>
dnmt_logger::log_event(string message) {
  static constexpr auto event_tag = "[event] ";
  log_stream << prefix
             << event_tag
             << message
             << '\n';
  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_event(string message, time_point<steady_clock> then) {
  static constexpr auto event_tag = "[event] ";
  static constexpr auto time_tag_left = " [delta time: ";
  static constexpr auto time_tag_right = "s]";
  const auto now = steady_clock::now();
  const auto delta_time = duration<double>(now - then).count();

  string formatted_data_and_time(128);
  snprintf(formatted_data_and_time.data(),
           128ul, time_format_template);

  log_stream << prefix << event_tag << message << time_tag_left << std::fixed
             << std::setprecision(2) << delta_time << time_tag_right
             << " " << format_current_time() << '\n';
  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_data(string key, string value) {
  static constexpr auto data_tag = "[data] ";
  log_stream << prefix << data_tag << key << "=" << value << '\n';
  return steady_clock::now();
}

time_point<steady_clock>
dnmt_logger::log_data(string key, string value, time_point<steady_clock> then) {
  static constexpr auto data_tag = "[data] ";
  static constexpr auto time_tag_left = " [delta time: ";
  static constexpr auto time_tag_right = "s]";
  const auto now = steady_clock::now();
  const auto delta_time = duration<double>(now - then).count();
  log_stream << prefix << data_tag << key << '=' << value << time_tag_left
             << std::fixed << std::setprecision(2) << delta_time
             << time_tag_right << '\n';
  return steady_clock::now();
}


// static void
// print_with_time(const string &s) {
//   std::timespec ts;
//   std::timespec_get(&ts, TIME_UTC);
//   vector<char> buf(128, '\0');
//   std::strftime(buf, std::size(buf), "%D %T", std::gmtime(&ts.tv_sec));
//   std::cout << "Current time: " << buf << '.' << ts.tv_nsec << " UTC\n";

//   auto tmp = system_clock::to_time_t(system_clock::now());
//   string time_fmt(std::ctime(&tmp));
//   time_fmt.pop_back();
//   cerr << "[" << time_fmt << "] " << s << endl;
// }
