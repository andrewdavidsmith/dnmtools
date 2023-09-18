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

#ifndef DNMT_LOGGER_HPP
#define DNMT_LOGGER_HPP

#include <chrono>
#include <iostream>
#include <ostream>
#include <string>
#include <fstream>

namespace dnmt_log {

struct screen {
  static auto get(std::string p = "", std::ostream &ls = std::clog)
    -> screen & {
    static screen instance(p, ls);
    return instance;
  }

  auto event(std::string message)
    -> std::chrono::time_point<std::chrono::steady_clock>;

  auto event(std::string message,
             std::chrono::time_point<std::chrono::steady_clock> scl)
    -> std::chrono::time_point<std::chrono::steady_clock>;

  auto data(std::string key, std::string value)
    -> std::chrono::time_point<std::chrono::steady_clock>;

  auto data(std::string key, std::string value,
            std::chrono::time_point<std::chrono::steady_clock> scl)
    -> std::chrono::time_point<std::chrono::steady_clock>;
private:
  std::string prefix;
  std::ostream &log_stream;

  screen(std::string p, std::ostream &ls): prefix{p}, log_stream{ls} {}
};

typedef screen &screen_ref;

struct logfile {
  static auto get(std::string p = "", std::string fn = "") -> logfile & {
    static logfile instance(p, fn);
    return instance;
  }

  auto event(std::string message)
    -> std::chrono::time_point<std::chrono::steady_clock>;

  auto event(std::string message,
             std::chrono::time_point<std::chrono::steady_clock> scl)
    -> std::chrono::time_point<std::chrono::steady_clock>;

  auto data(std::string key, std::string value)
    -> std::chrono::time_point<std::chrono::steady_clock>;

  auto data(std::string key, std::string value,
            std::chrono::time_point<std::chrono::steady_clock> scl)
    -> std::chrono::time_point<std::chrono::steady_clock>;
private:
  std::string prefix;
  std::ofstream log_stream;

  logfile(std::string p, std::string fn)
    : prefix{p}, log_stream{std::ofstream(fn, std::ios_base::app)} {}
};

typedef logfile &logfile_ref;

auto
init(std::string p, std::string fn = "", std::ostream &ls = std::clog) -> void;

auto
init(std::string p, std::ostream &ls = std::clog) -> void;

auto
event(std::string message)
  -> std::chrono::time_point<std::chrono::steady_clock>;

auto
event(std::string message,
      std::chrono::time_point<std::chrono::steady_clock> scl)
  -> std::chrono::time_point<std::chrono::steady_clock>;

auto
data(std::string key, std::string value)
  -> std::chrono::time_point<std::chrono::steady_clock>;

auto
data(std::string key, std::string value,
     std::chrono::time_point<std::chrono::steady_clock> scl)
  -> std::chrono::time_point<std::chrono::steady_clock>;

};  // namespace dnmt_log
#endif
