#pragma once

#include <string> // std::string
#include <vector> // std::vector<T>


namespace paw
{


struct FastaRecord
{
  std::string id;
  std::string seq;
};


class Fasta
{
public:
  using Records = std::vector<FastaRecord>;

  std::vector<std::string> ids;
  std::vector<std::string> seqs;

  Fasta(){}

  void add_record(std::string id, std::string seq);
  void load(std::string const & fn);
  void store(std::string const & fn);

  std::string get_sequence(std::size_t const index) const;

};


} // namespace paw


#ifdef IMPLEMENT_PAW


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


namespace io = boost::iostreams;


namespace paw
{

void
Fasta::add_record(std::string id, std::string seq)
{
  if (id.size() == 0 || id[0] != '>')
    id = ">" + id;

  ids.push_back(id);
  seqs.push_back(seq);
}


std::string
Fasta::get_sequence(std::size_t const index) const
{
  assert(index < seqs.size());
  return seqs[index];
}


void
Fasta::load(std::string const & fn)
{
  std::ifstream file(fn, std::ios_base::in | std::ios_base::binary);
  io::filtering_istream in;

  // If filename ends with ".gz", assume we should decompress with gzip
  if (fn.size() >= 3 && fn.substr(fn.size() - 3, 3) == ".gz")
    in.push(io::gzip_decompressor()); // gzip file

  // Read file
  in.push(file);

  std::stringstream ss; // sequence stream
  //FastaRecord rec = {"", ""}; // Empty FASTA record
  std::string id;
  std::string seq;
  std::getline(in, id); // Read first header

  for (std::string line; std::getline(in, line); )
  {
    // Check if the next line is a header line
    if (line[0] == '>')
    {
      seq = ss.str();

      // Force uppercase
      std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

      // Clear stringstream
      ss.str(std::string());

      // Add completed record
      add_record(std::move(id), std::move(seq));

      // Create a new record, where line is the ID
      id = std::move(line);
      seq.clear();
    }

    // Add line to sequence
    ss << line;
  }

  // Add final sequence
  seq = ss.str();

  // Force uppercase to final sequence
  std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

  // Add final record
  ids.push_back(id);
  seqs.push_back(seq);
}


void
Fasta::store(std::string const & fn)
{
  std::stringstream ss;

  for (std::size_t i = 0; i < ids.size(); ++i)
  {
    auto const & seq = seqs[i];

    // Print ID
    ss << ids[i] << "\n";

    // Print sequence
    std::size_t constexpr MAX_WIDTH = 90;

    for (std::size_t j = 0; j < seq.size(); j += MAX_WIDTH)
      ss << seq.substr(j, MAX_WIDTH) << "\n";
  }

  if (fn == "-")
  {
    std::cout << ss.rdbuf();
  }
  else
  {
    std::ofstream ofile(fn, std::ios_base::out | std::ios_base::binary);
    io::filtering_streambuf<boost::iostreams::input> out;

    if (boost::algorithm::ends_with(fn, ".gz"))
      out.push(boost::iostreams::gzip_compressor());

    out.push(ss);
    boost::iostreams::copy(out, ofile);
  }
}


} // namespace paw


#endif // IMPLEMENT_PAW
