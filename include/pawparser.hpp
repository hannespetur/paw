#pragma once
#include <exception> // std::exception
#include <sstream> // std::stringstream
#include <unordered_map> // std::unordered_map
#include <unordered_set> // std::unordered_set
#include <vector> // std::vector


namespace paw
{

/** TYPE DEFINITIONS */
using TParserFlags = std::unordered_map<std::string, std::string>;
using TPositional = std::vector<std::string>;

/** READ ONLY VARIABLES */
char const NO_SHORT_OPTION = ' ';
char const LIST_DELIMITER = ',';

/** EXCEPTIONS */
class MissingValueException :
  public std::exception
{
private:
  std::string msg;

public:
  MissingValueException(char const shrt, std::string const & lng)
  {
    std::stringstream ss;

    if (shrt != NO_SHORT_OPTION)
      ss << "[pawparser::MissingValue] Flag '--" << lng << "=value' (or '-" << shrt << "value') requires a value, but was not passed any.";
    else
      ss << "[pawparser::MissingValue] Flag '--" << lng << "=value' was not passed any value.";

    msg = ss.str();
  }


  virtual const char *
  what() const throw()
  {
    return msg.c_str();
  }


};

class MissingPositionalException :
  public std::exception
{
private:
  std::string msg;

public:
  MissingPositionalException(size_t const index, std::string const & meta_string)
  {
    std::stringstream ss;
    ss << "[pawparser::MissingPositional] Positional argument at index " << index << " '" << meta_string << "' was expected but is missing.";
    msg = ss.str();
  }


  virtual const char *
  what() const throw()
  {
    return msg.c_str();
  }


};


class ExtraOptionException :
  public std::exception
{
private:
  std::string msg;

public:
  ExtraOptionException(std::string const & unkown_option)
  {
    std::stringstream ss;
    ss << "[pawparser::ExtraOption] No option '" << unkown_option << "' is available.";
    msg = ss.str();
  }


  virtual const char *
  what() const throw()
  {
    return msg.c_str();
  }


};


struct Argument
{
  char shrt;
  std::string lng;
  std::string description;
  std::string meta_string;
  std::string default_value;
};


class Parser
{
private:
  /** CONTAINERS */
  unsigned major_version = 0;
  unsigned minor_version = 0;

  // Keys that maps options (or 'flags') to their given values as given by the user
  TParserFlags flags;

  // Vector of all values of positional arguments, in the same order as they were inserted
  TPositional positional;
  TPositional::iterator next_positional; // Iterator pointing to the next positional element to parse

  // Vector of all arguments
  std::vector<Argument> args;

  /** PRIVATE METHODS */
  TParserFlags::iterator inline check_flag(char const shrt, std::string const & lng); // Checks if option was set


public:
  Parser(Parser const & p) = delete; // Disallow copy contructor
  Parser(int argc, char ** argv); // Constructor uses the argc and argv variables passed to the 'main' function
  Parser(std::vector<std::string> const & args);

  void inline
  check_for_extra_options();

  template <typename T>
  void inline
  parse_option(T & val, char const shrt, std::string const & lng, std::string const & description);

  template <typename T>
  void inline
  parse_option_list(T & list, char const shrt, std::string const & lng, std::string const & description, char const delimiter = LIST_DELIMITER);

  template <typename T>
  void inline
  parse_positional_argument(T & val, std::string const & meta_string, std::string const & description);

  template <typename T>
  void inline
  parse_remaining_positional_arguments(T & list, std::string const & meta_string, std::string const & description);

  void inline
  set_version(unsigned MAJOR, unsigned MINOR);

};



/** PRIVATE METHODS */
TParserFlags::iterator inline
Parser::check_flag(char const shrt, std::string const & lng)
{
  // First check long options
  auto flag_it = flags.find(lng);

  if (flag_it != flags.end())
    return flag_it; // Found the flag

  // Then check short option
  return flags.find(std::string(1, shrt));
}


/** PUBLIC METHODS */
Parser::Parser(int argc, char ** argv)
  : Parser(std::vector<std::string>(argv + 1, argv + argc))
{}


Parser::Parser(std::vector<std::string> const & arguments)
{
  for (auto arg_it = arguments.cbegin(); arg_it != arguments.cend(); ++arg_it)
  {
    if (arg_it->size() <= 1 || (*arg_it)[0] != '-')
    {
      positional.push_back(std::move(*arg_it));
    }
    else if (arg_it->size() == 2 && (*arg_it)[1] == '-')
    {
      // Force all remaining arguments to be positional arguments
      std::move(arg_it + 1, arguments.end(), std::back_inserter(positional));
      break;
    }
    else if ((*arg_it)[1] == '-')
    {
      // long option flag was used
      auto it = std::find(arg_it->begin(), arg_it->end(), '=');
      std::string flag(arg_it->begin() + 2, it);

      if (it != arg_it->end())
        flags[std::move(flag)] = std::string(it + 1, arg_it->end()); // Equal sign used
      else
        flags[std::move(flag)] = ""; // No equal sign
    }
    else
    {
      // short option flag was used, e.g. '-d' or with value '-dtab'
      flags[std::string(arg_it->begin() + 1, arg_it->begin() + 2)] = std::string(arg_it->begin() + 2, arg_it->end());
    }
  }

  next_positional = positional.begin();
}


void inline
Parser::check_for_extra_options()
{
  // Find all known options
  std::unordered_set<std::string> known_options;

  for (auto const & arg : args)
  {
    if (arg.shrt != NO_SHORT_OPTION)
      known_options.insert(std::string(1, arg.shrt));

    known_options.insert(arg.lng);
  }

  // Check if any parsed arguments are not defined
  for (auto it = flags.begin(); it != flags.end(); ++it)
  {
    if (known_options.find(it->first) == known_options.end())
      throw ExtraOptionException(it->first);
  }
}


void inline
Parser::set_version(unsigned const major, unsigned const minor)
{
  major_version = major;
  minor_version = minor;
}


template <typename T>
void inline
Parser::parse_option(T & val, char const shrt, std::string const & lng, std::string const & description)
{
  {
    std::stringstream ss;
    ss << val;
    Argument arg = {shrt, lng, description, "value", ss.str()};
    args.push_back(std::move(arg));
  }

  auto flag_it = check_flag(shrt, lng);

  if (flag_it != flags.end())
  {
    if (flag_it->second.size() == 0)
      throw MissingValueException(shrt, lng);

    std::stringstream ss;
    ss << flag_it->second;
    ss >> val;
  }
}


/** bool specialization */
template <>
void inline
Parser::parse_option(bool & val, char const shrt, std::string const & lng, std::string const & description)
{
  args.push_back({shrt, lng, description, "", ""});
  auto flag_it = check_flag(shrt, lng);

  if (flag_it != flags.end())
    val ^= true; // Flips values
}


template <typename T>
void inline
Parser::parse_option_list(T & list, char const shrt, std::string const & lng, std::string const & description, char const delimiter)
{
  using TVal = typename T::value_type;
  args.push_back({shrt, lng, description, "list...", ""});
  auto flag_it = check_flag(shrt, lng);

  if (flag_it != flags.end())
  {
    if (flag_it->second.size() == 0)
      throw MissingValueException(shrt, lng);

    auto it = flag_it->second.begin();
    auto find_it = std::find(it, flag_it->second.end(), delimiter);

    while (it != flag_it->second.end())
    {
      // insert value
      {
        TVal val;
        std::stringstream ss;
        ss << std::string(it, find_it);
        ss >> val;
        list.insert(list.end(), val);
      }

      if (find_it == flag_it->second.end())
        break;

      it = find_it + 1;
      find_it = std::find(it, flag_it->second.end(), delimiter);
    }
  }
}


template <typename T>
void inline
Parser::parse_positional_argument(T & val, std::string const & meta_string, std::string const & description)
{
  args.push_back({NO_SHORT_OPTION, "", description, meta_string, ""});

  if (next_positional == positional.end())
    throw MissingPositionalException(std::distance(positional.begin(), next_positional), meta_string);

  std::stringstream ss;
  ss << *next_positional;
  ss >> val;
  ++next_positional;
}


template <typename T>
void inline
Parser::parse_remaining_positional_arguments(T & list, std::string const & meta_string, std::string const & description)
{
  using TVal = typename T::value_type;

  while (next_positional != positional.end())
  {
    TVal val;
    std::stringstream ss;
    ss << *next_positional;
    ss >> val;
    list.insert(list.end(), val);
    ++next_positional;
  }
}


} // namespace paw
