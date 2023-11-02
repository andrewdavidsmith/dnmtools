/* guessprotocol: a program for guessing whether a wgbs protocol is
 * wgbs, pbat or random pbat; also checks if the protocol is RRBS
 *
 * Copyright (C) 2019-2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
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

#include <bamxx.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>

#include "OptionParser.hpp"
#include "numerical_utils.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::min;
using std::runtime_error;
using std::string;
using std::vector;

using bamxx::bgzf_file;

static inline auto
entropy(const vector<double> &v) -> double {
  return
    transform_reduce(cbegin(v), cend(v), 0.0, std::plus<double>(),
                     [](const double x) -> double {
                       return x > 0.0 ? x*log(x)/log(2.0) : x;
                     });
}

static inline auto
is_cpg_pos(const string &s, const uint32_t idx) -> bool {
  return std::toupper(s[idx]) == 'C' &&
    (idx < size(s) - 1 && std::toupper(s[idx + 1]) == 'G');
}

static inline auto
is_cpg_neg(const string &s, const uint32_t idx) -> bool {
  return std::toupper(s[idx]) == 'G' &&
    (idx > 0 && std::toupper(s[idx - 1]) == 'C');
}

static inline uint8_t
complement_index(const uint8_t c) {
  return 3u - c;
}

constexpr int nuc_to_idx[] = {
 /*  0*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 16*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 32*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 48*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 64*/  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 80*/  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 96*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /*112*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /*128*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};


static inline auto
get_five_prime_kmer(const string &s, const uint32_t k) -> uint32_t {
  uint32_t r = 0;
  for (auto i = 0u; i < k; ++i)
    r = (r << 2) | nuc_to_idx[static_cast<uint8_t>(s[i])];
  return r;
}


struct nucleotide_model {
  vector<double> pr{};
  vector<double> lpr{};
  double bisulfite_conversion_rate{};
  bool is_t_rich{};

  nucleotide_model(const vector<double> &bc, const double conv_rate,
                   const bool is_t_rich)
      : pr{bc}, bisulfite_conversion_rate{conv_rate}, is_t_rich{is_t_rich} {
    auto nuc_from = is_t_rich ? 1 : 2;
    auto nuc_to = is_t_rich ? 3 : 0;
    pr[nuc_to] += bisulfite_conversion_rate * pr[nuc_from];
    pr[nuc_from] *= (1.0 - bisulfite_conversion_rate);
    assert(reduce(cbegin(pr), cend(pr), 0.0) == 1.0);
    lpr.resize(std::size(pr));
    transform(cbegin(pr), cend(pr), begin(lpr),
              [](const double x) { return std::log(x); });
  }

  double operator()(const string &s) const {
    return accumulate(cbegin(s), cend(s), 0.0,
                      [&](const double x, const char c) {
                        const auto i = nuc_to_idx[static_cast<uint8_t>(c)];
                        return i == 4 ? x : x + lpr[i];
                      });
  };

  string tostring() const {
    std::ostringstream oss;
    oss << "pr:\n";
    for (const auto i : pr) oss << i << '\n';
    oss << "log pr:\n";
    for (const auto i : lpr) oss << i << '\n';
    oss << bisulfite_conversion_rate << '\n' << is_t_rich;
    return oss.str();
  }
};

static inline auto
apply_nucleotide_model(const uint32_t k,
                       const vector<uint32_t> &c,
                       const nucleotide_model &m) -> double {
  double tot = 0.0;
  for (auto i = 0u; i < size(c); ++i) {
    double tmp = 0.0;
    auto x = i;
    for (auto j = 0u; j < k; ++j) {
      tmp += m.lpr[x & 3];
      x /= 4;
    }
    tot += c[i]*tmp;
  }
  return tot;
}

template<typename T> static inline auto
apply_nucleotide_model(const uint32_t k, const nucleotide_model &m)
  -> vector<double> {
  const auto dim = (1 << 2 * k);
  auto u = vector<double>(dim, 0.0);
  for (auto i = 0u; i < dim; ++i) {
    double tot = 0.0;
    auto x = i;
    for (auto j = 0u; j < k; ++j) {
      tot += m.lpr[x & 3];
      x /= 4;
    }
    u[i] = tot;
  }
  return u;
}

static inline auto
add_pseudocount(const double pseudocount, vector<double> &p) -> void {
  transform(cbegin(p), cend(p), begin(p),
            [=](const double x) { return x + pseudocount; });
  const auto tot = reduce(cbegin(p), cend(p), 0.0);
  transform(cbegin(p), cend(p), begin(p),
            [=](const double x) { return x / tot; });
}


struct rrbs_model {
  string name{};
  vector<string> consensus{};
  vector<vector<double>> pr{};
  vector<vector<double>> lpr{};

  rrbs_model(const string &name, const vector<string> &consensus,
             const uint32_t alphabet_size, const double pseudocount)
    : name{name}, consensus{consensus} {
    // alternate model
    pr.resize(size(consensus.front()), vector<double>(alphabet_size, 0.0));
    for (auto i = 0u; i < size(consensus); ++i) {
      for (auto j = 0u; j < size(consensus[i]); ++j) {
        const uint8_t b = consensus[i][j];
        pr[j][nuc_to_idx[b]] += 1.0;
      }
    }

    for (auto &p : pr) add_pseudocount(pseudocount, p);

    lpr = pr;
    for (auto &p : lpr)
      transform(cbegin(p), cend(p), begin(p),
                [](const double x) { return std::log(x); });
  }

  // auto operator()(const string &s) const -> double {
  //   auto lpr_itr = begin(lpr);
  //   return accumulate(cbegin(s), cbegin(s) + size(consensus), 0.0,
  //                     [&](const double x, const char c) {
  //                       const auto i = nuc_to_idx[static_cast<uint8_t>(c)];
  //                       return i == 4 ? x : x + (*lpr_itr++)[i];
  //                     });

  auto operator()(const vector<uint32_t> &c) const -> double {
    const auto k = size(consensus.front());
    double tot = 0.0;
    for (auto i = 0u; i < size(c); ++i) {
      double kmer_tot = 0.0;
      auto position = i;
      for (auto j = 0u; j < k; ++j) {
        kmer_tot += lpr[k - j - 1][position & 3];
        position >>= 2;
      }
      tot += c[i]*kmer_tot;
    }
    return tot;
  }

  auto tostring() const -> string {
    constexpr auto precision_val = 4u;
    std::ostringstream oss;
    oss.precision(precision_val);
    oss.setf(std::ios_base::fixed, std::ios_base::floatfield);
    oss << "name: " << name << '\n';
    oss << "consensus:\n";
    for (const auto &c : consensus)
      oss << c << ",\n";
    oss << "pr:\n";
    for (const auto &i : pr) {
      for (const auto j : i) oss << j << ' ';
      oss << '\n';
    }
    oss << "lpr:\n";
    for (const auto &i : lpr) {
      for (const auto j : i) oss << j << ' ';
      oss << '\n';
    }
    return oss.str();
  }
};



struct guessprotocol_summary {

  static constexpr auto wgbs_cutoff_confident = 0.99;
  static constexpr auto wgbs_cutoff_unconfident = 0.9;
  static constexpr auto rpbat_cutoff_confident_high = 0.8;
  static constexpr auto rpbat_cutoff_confident_low = 0.2;
  static constexpr auto pbat_cutoff_unconfident = 0.1;
  static constexpr auto pbat_cutoff_confident = 0.01;

  // protocol is the guessed protocol (wgbs, pbat, rpbat, or inconclusive)
  // based on the content of the reads.
  string protocol;
  // confidence indicates the level of confidence in the guess for the
  // protocol.
  string confidence;
  // layout indicates whether the reads are paired or single-ended.
  string layout;
  // n_reads_wgbs is the average number of reads (for single-ended reads) or
  // read pairs (for paired reads) where read1 is T-rich.
  double n_reads_wgbs{};
  // n_reads is the number of evaluated reads or read pairs.
  uint64_t n_reads{};
  // wgbs_fraction is the probability that a read (for single-ended reads) or
  // the read1 of a read pair (for paired reads) is T-rich.
  double wgbs_fraction{};
  // rrbs_fraction is the sum over all reads of the probability that
  // the read is from RRBS as indicated by the prefix of the read.
  double rrbs_fraction{};


  void evaluate() {

    const auto frac = n_reads_wgbs / n_reads;
    protocol = "inconclusive";

    // assigning wgbs (near one)
    if (frac > wgbs_cutoff_confident) {
      protocol = "wgbs";
      confidence = "high";
    }
    else if (frac > wgbs_cutoff_unconfident) {
      protocol = "wgbs";
      confidence = "low";
    }
    // assigning pbat (near zero)
    else if (frac < pbat_cutoff_confident) {
      protocol = "pbat";
      confidence = "high";
    }
    else if (frac < pbat_cutoff_unconfident) {
      protocol = "pbat";
      confidence = "low";
    }
    // assigning rpbat (towards middle)
    else if (frac > rpbat_cutoff_confident_low &&
             frac < rpbat_cutoff_confident_high) {
      protocol = "rpbat";
      confidence = "high";
    }
    else {
      protocol = "rpbat";
      confidence = "low";
    }

    wgbs_fraction = frac;
    rrbs_fraction = rrbs_fraction / n_reads;
  }

  string tostring() const {
    std::ostringstream oss;
    oss << "protocol: " << protocol << '\n'
        << "confidence: " << confidence << '\n'
        << "wgbs_fraction: " << wgbs_fraction  << '\n'
        << "n_reads_wgbs: " << n_reads_wgbs << '\n'
        << "n_reads: " << n_reads << '\n'
        << "rrbs_fraction: " << rrbs_fraction;
    return oss.str();
  }
};

// store each read from one end
struct FASTQRecord {
  string name;
  string seq;
};

// see if two reads from two ends match to each other (they should
// have the same name)
static bool
mates(const size_t to_ignore_at_end, // in case names have #0/1 name ends
      const FASTQRecord &a, const FASTQRecord &b) {
  assert(to_ignore_at_end < std::size(a.name));
  return equal(cbegin(a.name), cend(a.name) - to_ignore_at_end,
               cbegin(b.name));
}

// Read 4 lines one time from fastq and fill in the FASTQRecord structure
static bgzf_file &
operator>>(bgzf_file &s, FASTQRecord &r) {
  constexpr auto n_error_codes = 5u;

  enum err_code { none, bad_name, bad_seq, bad_plus, bad_qual };

  static const array<runtime_error, n_error_codes> error_msg = {
    runtime_error(""), runtime_error("failed to parse fastq name line"),
    runtime_error("failed to parse fastq sequence line"),
    runtime_error("failed to parse fastq plus line"),
    runtime_error("failed to parse fastq qual line")
  };

  err_code ec = err_code::none;

  if (!getline(s, r.name)) return s;

  if (r.name.empty() || r.name[0] != '@') ec = err_code::bad_name;

  const auto nm_end = r.name.find_first_of(" \t");
  const auto nm_sz = (nm_end == string::npos ? r.name.size() : nm_end) - 1;
  r.name.erase(copy_n(cbegin(r.name) + 1, nm_sz, begin(r.name)), cend(r.name));

  if (!getline(s, r.seq)) ec = err_code::bad_seq;

  string tmp;
  if (!getline(s, tmp)) ec = err_code::bad_plus;

  if (!getline(s, tmp)) ec = err_code::bad_qual;

  if (ec != err_code::none) throw error_msg[ec];

  return s;
}

int
main_guessprotocol(int argc, const char **argv) {

  try {

    constexpr auto alphabet_size = 4u;
    static const vector<double> human_base_comp = {0.295, 0.205, 0.205, 0.295};
    static const vector<double> flat_base_comp = {0.25, 0.25, 0.25, 0.25};

    constexpr auto description = "guess bisulfite protocol for a library";

    bool verbose = false;
    bool use_human;
    string outfile;
    size_t reads_to_check = 1000000;
    size_t name_suffix_len = 0;
    double bisulfite_conversion_rate = 0.98;
    double pseudocount = 0.1;
    uint32_t kmer = 3;

    namespace fs = std::filesystem;
    const string cmd_name = std::filesystem::path(argv[0]).filename();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(cmd_name, description,
                           "<end1-fastq> [<end2-fastq>]");
    opt_parse.add_opt("nreads", 'n', "number of reads in initial check",
                      false, reads_to_check);
    opt_parse.add_opt("ignore", 'i', "length of read name suffix "
                      "to ignore when matching", false, name_suffix_len);
    opt_parse.add_opt("bisulfite", 'b', "bisulfite conversion rate",
                      false, bisulfite_conversion_rate);
    opt_parse.add_opt("human", 'H', "assume human genome",
                      false, use_human);
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("verbose", 'v',
                      "report available information during the run",
                      false, verbose);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested() || leftover_args.size() > 2) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const vector<string> reads_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    auto base_comp = flat_base_comp;
    if (use_human) base_comp = human_base_comp;

    nucleotide_model t_rich_model(base_comp, bisulfite_conversion_rate, true);
    nucleotide_model a_rich_model(base_comp, bisulfite_conversion_rate, false);

    guessprotocol_summary summary;
    summary.layout = reads_files.size() == 2 ? "paired" : "single";

    if (verbose) {
      if (reads_files.size() == 2)
        cerr << "data layout: "
             << "paired" << '\n'
             << "read1 file: " << reads_files.front() << '\n'
             << "read2 file: " << reads_files.back() << '\n';
      else
        cerr << "data layout: "
             << "single" << '\n'
             << "read file: " << reads_files.front() << '\n';
      cerr << "reads to check: " << reads_to_check << '\n'
           << "read name suffix length: " << name_suffix_len << '\n'
           << "bisulfite conversion: " << bisulfite_conversion_rate << '\n';
    }

    const auto n_kmers = (1u << 2*kmer);
    vector<uint32_t> kmer_count(n_kmers, 0u);

    if (reads_files.size() == 2) {

      // input: paired-end reads with end1 and end2
      bgzf_file in1(reads_files.front(), "r");
      if (!in1)
        throw runtime_error("cannot open file: " + reads_files.front());

      bgzf_file in2(reads_files.back(), "r");
      if (!in2)
        throw runtime_error("cannot open file: " + reads_files.back());

      FASTQRecord r1, r2;
      while (in1 >> r1 && in2 >> r2 && summary.n_reads < reads_to_check) {
        summary.n_reads++;

        if (!mates(name_suffix_len, r1, r2))
          throw runtime_error("expected mates: " + r1.name + ", " + r2.name);

        const double ta = t_rich_model(r1.seq) + a_rich_model(r2.seq);
        const double at = a_rich_model(r1.seq) + t_rich_model(r2.seq);

        const auto prob_read1_t_rich = exp(ta - log_sum_log(ta, at));
        summary.n_reads_wgbs += prob_read1_t_rich;

        if (r1.seq.find_first_not_of("ACGT") >= kmer)
          ++kmer_count[get_five_prime_kmer(r1.seq, kmer)];
      }
    }
    else {

      // input: single-end reads
      bgzf_file in(reads_files.front(), "r");
      if (!in)
        throw runtime_error("cannot open file: " + reads_files.front());

      FASTQRecord r;
      while (in >> r && summary.n_reads < reads_to_check) {
        summary.n_reads++;

        const double t = t_rich_model(r.seq);
        const double a = a_rich_model(r.seq);

        const auto prob_t_rich = exp(t - log_sum_log(t, a));
        summary.n_reads_wgbs += prob_t_rich;

        if (r.seq.find_first_not_of("ACGT") >= kmer)
          ++kmer_count[get_five_prime_kmer(r.seq, kmer)];
      }
    }
    summary.evaluate();

    auto MspI = rrbs_model("MspI", {"CGG", "TGG"}, alphabet_size, pseudocount);

    cerr << MspI.tostring() << endl;

    const double mspi_evidence = MspI(kmer_count);
    const double null_evidence = summary.wgbs_fraction > 0.5 ?
      apply_nucleotide_model(kmer, kmer_count, t_rich_model) :
      apply_nucleotide_model(kmer, kmer_count, a_rich_model);
    cerr << mspi_evidence << endl;
    cerr << null_evidence << endl;

    const uint64_t tot = reduce(cbegin(kmer_count), cend(kmer_count));

    vector<double> kmer_freq;
    for (const auto p : kmer_count)
      kmer_freq.push_back(p/static_cast<double>(tot));
    cerr << entropy(kmer_freq) << endl;

    const double prob_rrbs = exp(mspi_evidence - log_sum_log(mspi_evidence, null_evidence));
    cout << prob_rrbs << endl;
    if (!outfile.empty()) {
      std::ofstream out(outfile);
      if (!out) throw runtime_error("failed to open: " + outfile);
      out << summary.tostring() << endl;
    }
    else cout << summary.tostring() << endl;
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
