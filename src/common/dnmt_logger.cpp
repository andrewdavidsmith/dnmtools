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

static inline auto
format_current_time() -> string {
  static constexpr auto time_template = "[%D %T]";
  timespec ts;
  timespec_get(&ts, TIME_UTC);
  vector<char> buf(128);
  strftime(buf.data(), std::size(buf), time_template, std::gmtime(&ts.tv_sec));
  return buf.data();
}

static inline auto
format_elapsed_time(time_point<steady_clock> then) -> string {
  static constexpr auto delta_time_template = "[delta time: %.3fs]";
  static constexpr auto max_delta_time_size = 128;
  const auto now = steady_clock::now();
  const auto delta_time = duration<double>(now - then).count();
  string formatted_delta_time(max_delta_time_size, '\0');
  const auto delta_time_size =
    std::snprintf(formatted_delta_time.data(), std::size(formatted_delta_time),
                  delta_time_template, delta_time);
  formatted_delta_time.resize(delta_time_size);
  return formatted_delta_time;
}

auto
dnmt_log::screen::event(string message) -> time_point<steady_clock> {
  static constexpr auto event_tag = "[event]";
  // clang-format off
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << event_tag << ' '
             << message
             << '\n';
  // clang-format on
  return steady_clock::now();
}

auto
dnmt_log::screen::event(string message, time_point<steady_clock> then) -> time_point<steady_clock> {
  static constexpr auto event_tag = "[event]";
  // clang-format off
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << event_tag << ' '
             << message << ' '
             << format_elapsed_time(then)
             << '\n';
  // clang-format on
  return steady_clock::now();
}

auto
dnmt_log::screen::data(string key, string value) -> time_point<steady_clock> {
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

auto
dnmt_log::screen::data(string key, string value,
                       time_point<steady_clock> then) -> time_point<steady_clock> {
  static constexpr auto data_tag = "[data]";
  // clang-format off
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << data_tag << ' '
             << key << '=' << value << ' '
             << format_elapsed_time(then) << '\n';
  // clang-format on
  return steady_clock::now();
}

auto
dnmt_log::logfile::event(string message) -> time_point<steady_clock> {
  static constexpr auto event_tag = "[event]";
  // clang-format off
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << event_tag << ' '
             << message
             << '\n';
  // clang-format on
  return steady_clock::now();
}

auto
dnmt_log::logfile::event(string message, time_point<steady_clock> then) -> time_point<steady_clock> {
  static constexpr auto event_tag = "[event]";
  // clang-format off
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << event_tag << ' '
             << message << ' '
             << format_elapsed_time(then)
             << '\n';
  // clang-format on
  return steady_clock::now();
}

auto
dnmt_log::logfile::data(string key, string value) -> time_point<steady_clock> {
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

auto
dnmt_log::logfile::data(string key, string value, time_point<steady_clock> then) -> time_point<steady_clock> {
  static constexpr auto data_tag = "[data]";
  // clang-format off
  log_stream << prefix << ' '
             << format_current_time() << ' '
             << data_tag << ' '
             << key << '=' << value << ' '
             << format_elapsed_time(then)
             << '\n';
  // clang-format on
  return steady_clock::now();
}

auto
dnmt_log::event(std::string message)
  -> std::chrono::time_point<std::chrono::steady_clock> {
  screen::get().event(message);
  return logfile::get().event(message);
}

auto
dnmt_log::event(std::string message,
                std::chrono::time_point<std::chrono::steady_clock> scl)
  -> std::chrono::time_point<std::chrono::steady_clock> {
  screen::get().event(message, scl);
  return logfile::get().event(message, scl);
}

auto
dnmt_log::data(std::string key, std::string value)
  -> std::chrono::time_point<std::chrono::steady_clock> {
  screen::get().data(key, value);
  return logfile::get().data(key, value);
}

auto
dnmt_log::data(std::string key, std::string value,
     std::chrono::time_point<std::chrono::steady_clock> scl)
  -> std::chrono::time_point<std::chrono::steady_clock> {
  screen::get().data(key, value, scl);
  return logfile::get().data(key, value, scl);
}

auto
dnmt_log::init(std::string p, std::string fn, std::ostream &ls) -> void {
  dnmt_log::screen::get(p, ls);
  dnmt_log::logfile::get(p, fn);
}

auto
dnmt_log::init(std::string p, std::ostream &ls) -> void {
  dnmt_log::screen::get(p, ls);
  dnmt_log::logfile::get(p, string());
}
