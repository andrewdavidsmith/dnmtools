/* methstates: a program for converting read sequences in SAM format
 * files into methylation states at CpGs covered by those reads
 *
 * Copyright (C) 2011-2022 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith and Masaru Nakajima
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

#include <algorithm>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <charconv>
#include <vector>

#include "OptionParser.hpp"
#include "bam_record_utils.hpp"
#include "cigar_utils.hpp"
#include "dnmt_error.hpp"
#include "htslib_wrapper.hpp"
#include "sam_record.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::lower_bound;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::distance;
using std::min;
using std::swap;

using bamxx::bam_rec;

static const char b2c[] = "TNGNNNCNNNNNNNNNNNNA";

template<class BidirIt, class OutputIt>
// constexpr // since C++20
OutputIt
revcomp_copy(BidirIt first, BidirIt last, OutputIt d_first) {
  for (; first != last; ++d_first) *d_first = b2c[*(--last) - 'A'];
  return d_first;
}

static inline bool
is_cpg(const string &s, const uint64_t idx) {
  return s[idx] == 'C' && s[idx + 1] == 'G';
}

static void
collect_cpgs(const string &s, vector<uint64_t> &cpgs) {
  cpgs.clear();
  const uint64_t lim = std::size(s) - 1;
  for (auto i = 0u; i < lim; ++i)
    if (is_cpg(s, i)) cpgs.push_back(i);
}

static bool
convert_meth_states_pos(const vector<uint64_t> &cpgs,
                        const bamxx::bam_header &hdr, const bam_rec &aln,
                        uint64_t &first_cpg_index, string &states) {
  states.clear();

  const uint64_t seq_start = get_pos(aln);
  const uint64_t width = rlen_from_cigar(aln);
  const uint64_t seq_end = seq_start + width;

  string seq_str;
  get_seq_str(aln, seq_str);
  apply_cigar(aln, seq_str, 'N');

  if (std::size(seq_str) != width)
    throw dnmt_error("bad sam record format: " + to_string(hdr, aln));

  // get the first cpg site equal to or large than seq_start
  auto cpg_itr = lower_bound(begin(cpgs), end(cpgs), seq_start);
  auto first_cpg_itr = end(cpgs);

  if (cpg_itr == end(cpgs)) return false;

  for (; cpg_itr != end(cpgs) && *cpg_itr < seq_end; cpg_itr++) {
    const char x = seq_str[*cpg_itr - seq_start];
    states += (x == 'T') ? 'T' : ((x == 'C') ? 'C' : 'N');
    if (first_cpg_itr == end(cpgs)) first_cpg_itr = cpg_itr;
  }

  if (first_cpg_itr != end(cpgs))
    first_cpg_index = distance(begin(cpgs), first_cpg_itr);

  return states.find_first_of("CT") != string::npos;
}

static bool
convert_meth_states_neg(const vector<uint64_t> &cpgs,
                        const bamxx::bam_header &hdr, const bam_rec &aln,
                        uint64_t &first_cpg_index, string &states) {
  /* ADS: the "revcomp" on the read sequence is needed for the cigar
     to be applied, since the cigar is relative to the genome
     coordinates and not the read's sequence. But the read sequence
     may is assumed to have been T-rich to begin with, so it becomes
     A-rich. And the position of the C in the CpG becomes the G
     position.
   */

  states.clear();

  const uint64_t seq_start = get_pos(aln);
  const uint64_t width = rlen_from_cigar(aln);
  const uint64_t seq_end = seq_start + width;

  string orig_seq;
  get_seq_str(aln, orig_seq);

  string seq_str;
  seq_str.resize(orig_seq.size());
  revcomp_copy(begin(orig_seq), end(orig_seq), begin(seq_str));
  apply_cigar(aln, seq_str, 'N');

  if (seq_str.size() != width)
    throw dnmt_error("bad sam record format: " + to_string(hdr, aln));

  // get the first cpg site equal to or large than seq_start - 1
  // the -1 is because we look for G in the read corresponding to a
  // CpG in chromosome, which are indexed in cpgs based on the position of C
  auto cpg_itr =
    lower_bound(begin(cpgs), end(cpgs), seq_start > 0 ? seq_start - 1 : 0);
  auto first_cpg_itr = end(cpgs);

  if (cpg_itr == end(cpgs)) { return false; }
  else {
    for (; cpg_itr != end(cpgs) && *cpg_itr < seq_end - 1; cpg_itr++) {
      const char x = seq_str[*cpg_itr - seq_start + 1];
      states += (x == 'G') ? 'C' : ((x == 'A') ? 'T' : 'N');
      if (first_cpg_itr == end(cpgs)) first_cpg_itr = cpg_itr;
    }
  }

  if (first_cpg_itr != end(cpgs)) {
    first_cpg_index = distance(begin(cpgs), first_cpg_itr);
  }

  return states.find_first_of("CT") != string::npos;
}

static void
get_chrom(const string &chrom_name, const vector<string> &all_chroms,
          const unordered_map<string, uint64_t> &chrom_lookup, string &chrom) {
  auto the_chrom = chrom_lookup.find(chrom_name);
  if (the_chrom == end(chrom_lookup))
    throw dnmt_error("could not find chrom: " + chrom_name);

  chrom = all_chroms[the_chrom->second];
  if (chrom.empty()) throw dnmt_error("problem with chrom: " + chrom_name);
}

struct state_set {
  int32_t tid{};
  uint64_t pos{};
  string seq{};
  state_set(int32_t tid, uint64_t pos, string &&seq)
    : tid{tid}, pos{pos}, seq{move(seq)} {}
  uint32_t serialize(const bamxx::bam_header &hdr, char *buf, const uint32_t buf_size) const {
    const auto b_end = buf + buf_size;
    auto b_itr = buf;

    auto n_itr = sam_hdr_tid2name(hdr.h, tid);
    while (b_itr != b_end && *n_itr != '\0') *b_itr++ = *n_itr++;
    *b_itr++ = '\t';

    auto [ptr, ec] = std::to_chars(b_itr, b_end, pos);
    b_itr = ptr;
    *b_itr++ = '\t';

    const auto n = size(seq);
    copy_n(cbegin(seq), n, b_itr);
    b_itr += n;
    *b_itr++ = '\n';

    return distance(buf, b_itr);
  }
};

int
main_methstates(int argc, const char **argv) {
  try {

    constexpr uint32_t buf_size = 1024;

    const string description =
      "Convert mapped reads in SAM format into a format that indicates binary \
      sequences of methylation states in each read, indexed by the identity   \
      of the CpG they cover, along with the chromosome. Only reads that       \
      cover a CpG site are included in the output. All output is relative to  \
      the positive reference strand. This format is used as input to other    \
      tools, and is not intended to be human-interpretable. All chromosome    \
      sequences are loaded at once.";

    bool VERBOSE = false;
    bool compress_output = false;

    string chrom_file;
    string outfile("-");
    int n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], description, "<sam-file>");
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("chrom", 'c', "fasta format reference genome file", true,
                      chrom_file);
    opt_parse.add_opt("threads", 't', "threads to use for reading input", false,
                      n_threads);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (n_threads < 0) throw dnmt_error("thread count cannot be negative");

    /* first load in all the chromosome sequences and names, and make
       a map from chromosome name to the location of the chromosome
       itself */
    vector<string> all_chroms, chrom_names;
    read_fasta_file_short_names(chrom_file, chrom_names, all_chroms);
    for (auto &&i : all_chroms)
      transform(begin(i), end(i), begin(i),
                [](const char c) { return std::toupper(c); });

    unordered_map<string, uint64_t> chrom_lookup;
    for (uint64_t i = 0; i < chrom_names.size(); ++i)
      chrom_lookup[chrom_names[i]] = i;

    if (VERBOSE) cerr << "n_chroms: " << all_chroms.size() << endl;

    bamxx::bam_tpool tp(n_threads);  // declared first; destroyed last

    bamxx::bam_in in(mapped_reads_file);
    if (!in) throw dnmt_error("cannot open input file " + mapped_reads_file);
    bamxx::bam_header hdr(in);
    if (!hdr) throw dnmt_error("cannot read heade" + mapped_reads_file);

    /* set the threads for the input file decompression */
    if (n_threads > 1) tp.set_io(in);

    // open the output file
    const string output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out) throw dnmt_error("error opening output file: " + outfile);

    vector<uint64_t> cpgs;

    unordered_set<string> chroms_seen;
    string chrom_name;
    string chrom;

    vector<char> buf(buf_size, '\0');
    vector<state_set> out1;
    vector<state_set> out2;

    // If the next position precedes both pos1 and pos2, we have moved
    // to a new tid, and we output everything in out1 and out2. If the
    // next position matches one of pos1 or pos2, we enqueue the
    // record for the appropriate among out1 or out2. If the next
    // position follows pos1 and pos2, then we will never again see
    // the earlier of pos1 and pos2, so all records matching the
    // earlier of pos1 or pos2 can be written.
    uint64_t pos1{}, pos2{};

    // iterate over records/reads in the SAM file, sequentially
    // processing each before considering the next
    bam_rec aln;
    while (in.read(hdr, aln)) {
      // get the correct chrom if it has changed
      if (string(sam_hdr_tid2name(hdr, aln)) != chrom_name) {
        chrom_name = sam_hdr_tid2name(hdr, aln);

        // make sure all reads from same chrom are contiguous in the file
        if (chroms_seen.find(chrom_name) != end(chroms_seen))
          throw dnmt_error("chroms out of order (check SAM file sorted)");

        if (VERBOSE) cerr << "processing " << chrom_name << endl;

        get_chrom(chrom_name, all_chroms, chrom_lookup, chrom);
        collect_cpgs(chrom, cpgs);
      }

      uint64_t first_cpg_index = std::numeric_limits<uint64_t>::max();
      string seq;

      const bool has_cpgs =
        bam_is_rev(aln)
          ? convert_meth_states_neg(cpgs, hdr, aln, first_cpg_index, seq)
          : convert_meth_states_pos(cpgs, hdr, aln, first_cpg_index, seq);

      if (has_cpgs) {
        if (first_cpg_index == pos2)
          out2.emplace_back(get_tid(aln), first_cpg_index, move(seq));
        else if (first_cpg_index == pos1)
          out1.emplace_back(get_tid(aln), first_cpg_index, move(seq));
        else {
          for (const auto &x : out1) {
            const auto n_bytes = x.serialize(hdr, buf.data(), buf_size);
            if (!out.write(buf.data(), n_bytes)) {
              cerr << "failure writing output" << endl;
              return EXIT_FAILURE;
            }
          }
          out1.clear();
          pos1 = 0; // if this is equal to pos2, as happens on advance
                    // of tid, out2 will be used by order of
                    // conditions above, and we might flush an empty
                    // queue once per tid; otherwise this value will
                    // never be used

          if (first_cpg_index < min(pos1, pos2)) {
            // advance tid: just flushed out1, but we need to also
            // flush out2 and then out2 will be ready
            for (const auto &x : out2) {
              const auto n_bytes = x.serialize(hdr, buf.data(), buf_size);
              if (!out.write(buf.data(), n_bytes)) {
                cerr << "failure writing output" << endl;
                return EXIT_FAILURE;
              }
            }
            out2.clear();
          }
          else {
            // advancing pos: just flushed out1, so swap and out2 is
            // ready for current aln
            swap(pos1, pos2);
            swap(out1, out2);
          }

          pos2 = first_cpg_index;
          out2.emplace_back(get_tid(aln), first_cpg_index, move(seq));
        }
      }
    }
    for (const auto &x : out1) {
      const auto n_bytes = x.serialize(hdr, buf.data(), buf_size);
      if (!out.write(buf.data(), n_bytes)) {
        cerr << "failure writing output" << endl;
        return EXIT_FAILURE;
      }
    }
    for (const auto &x : out2) {
      const auto n_bytes = x.serialize(hdr, buf.data(), buf_size);
      if (!out.write(buf.data(), n_bytes)) {
        cerr << "failure writing output" << endl;
        return EXIT_FAILURE;
      }
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
