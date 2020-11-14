#pragma once

#include <algorithm> // std::find
#include <cassert>
#include <iostream> //
#include <map> // std::multimap
#include <unordered_set> // std::unordered_set
#include <sstream> // std::istringstream, std::ostringstream
#include <string> // std::string
#include <utility> // std::pair
#include <vector> // std::vector<T>

#include <paw/parser/exception.hpp>
#include <paw/parser/internal.hpp>


/** \file paw/parser.hpp
 * The main header file of paw::Parser.
 */

/** Top-level namespace of the paw library.*/
namespace paw
{

/** Main Data Structure for paw parser.*/
class Parser
{
public:
  /** \brief Data structure for command line arguments.
   *  \details Arguments in this data structure are ready for display in the 'help' page of the
   *          program.
   */
  struct Arg
  {
public:
    /** Short option for the 'help' page.
     * No option shown if it is paw::internal::NO_SHORT_OPTION.
     */
    char shrt;
    std::string lng;           /**< Long option for the 'help' page.*/
    std::string description;   /**< Description of the argument to display in the 'help' page.*/
    std::string meta_string;   /**< String to represent the value with in the 'help' page.*/
    std::string default_value; /**< The default value to display in the 'help' page.*/
  };

  /* TYPE DEFINITIONS */
  /** Type definition for the hash table container to use for parsing arguments.*/
  using FlagMap = std::multimap<std::string, std::string>;

  /** Type definition of the container used to store the arguments.*/
  using Args = std::vector<Arg>;

public:
  /** \brief Construction of a paw parser instance.
   * \details This constructor uses the argc and argv variables that are passed to the 'main'
   *          function of the program.
   * \param[in] argc length of the argv array.
   * \param[in] argv arguments passed, including the name of the executable.
   * \returns paw parser object.
   * \exception None.
   */
  Parser(int argc, char ** argv);

  /** \brief Construction of a paw Parser instance.
   * \details This constructor parses a vector of strings that should include the programs name.
   * \param[in] args Dynamic vector with all arguments.
   * \returns paw Parser object.
   * \exception None.
   */
  Parser(std::vector<std::string> const & args);

  /** \brief Adds an subcommand to the parse.
   * \details Add subcommand and use the paw::Parser::parse_subcommand() to parse them.
   * \param[in] subcommand_name Name of the subcommand.
   * \param[in] description Description of the subcommand.
   * \exception None.
   */
  void add_subcommand(std::string const & subcommand_name, std::string const & description);

  /** Checks if the user passed invalid options.
   * \exception paw::exception::invalid_option thrown if user passed an invalid option.
   */
  void check_for_invalid_options();

  /** \brief Finalizes paw Parser.
   * \details Parses help flag and checks for invalid options.
   * \param[in] allow_no_arguments Set if the parser should not throw help when no arguments
   *                               are passed.
   * \exception paw::exception::invalid_option thrown if user passed an invalid option.
   */
  void finalize(bool const allow_no_arguments = true);


  /** Get a reference to the multimap containing all the flags the user passed.
   * \returns Multimap with flags.
   * \exception None.
   */
  FlagMap const &
  get_flag_map_reference() const;

  /**
   * Get a read-only reference to the list of subcommands.
   * \returns std::vector of Subcommands
   * \exception None.
   */
  std::vector<std::pair<std::string, std::string> > const &
  get_subcommands_reference() const;

  /**
   * Get a read-only reference to the list of positional arguments
   * \returns std::vector<std::string> containing the arguments in order
   * \exception None.                                            \
   */
  std::vector<std::string> const &
  get_positional_arguments_reference() const;

  /** Parses an option that was passed by the user.
   * If 'T' is not a boolean, we expect that the argument has a single value.
   * \param[in,out] val Reference to the value to change if the option was parsed.
   * \param[in] shrt Short option.
   * \param[in] lng Long option.
   * \param[in] description Short description of the option.
   * \param[in] meta_string Meta string for the option's value.
   * \exception paw::exception::missing_value thrown when option was parsed
   *            but without a value.
   */
  template <typename T>
  inline void
  parse_option(T & val,
               char const shrt,
               std::string const & lng,
               std::string const & description,
               std::string const & meta_string = "value");

  /** Parses an advanced/hidden option that was passed by the user.
   * If 'T' is not a boolean, we expect that the argument has a single value.
   * \param[in,out] val Reference to the value to change if the option was parsed.
   * \param[in] shrt Short option.
   * \param[in] lng Long option.
   * \param[in] description Short description of the option.
   * \param[in] meta_string Meta string for the option's value.
   * \exception paw::exception::missing_value thrown when option was parsed
   *            but without a value.
   */
  template <typename T>
  inline void
  parse_advanced_option(T & val,
                        char const shrt,
                        std::string const & lng,
                        std::string const & description,
                        std::string const & meta_string = "value");

  /** Parses an argument with a list of values.
   * \param[in,out] list Reference to the list to change if the option was parsed.
   * \param[in] shrt Short option.
   * \param[in] lng Long option.
   * \param[in] description Short description of the option.
   * \param[in] delimiter Delimiter between items in the list.
   * \param[in] meta_string Meta string for the option's value.
   * \exception paw::exception::missing_value thrown when option list was parsed,
   *            but without a value.
   */
  template <typename T>
  inline void
  parse_option_list(T & list,
                    char const shrt,
                    std::string const & lng,
                    std::string const & description,
                    char const delimiter = paw::internal::DEFAULT_LIST_DELIMITER,
                    std::string const & meta_string = "value1,...");

  /** Parses the next positional argument.
   * \param[in,out] val Reference to the value to change if the option was parsed.
   * \param[in] meta_string Meta string for the option's value.
   * \param[in] description Short description of the option.
   * \exception paw::exception::missing_positional_argument thrown when there was no
   *            positional argument passed.
   */
  template <typename T>
  inline void
  parse_positional_argument(T & val,
                            std::string const & meta_string,
                            std::string const & description);

  /** Parse all remaining positional arguments.
   * \param[in,out] list Reference to the list to change if the option was parsed.
   * \param[in] meta_string Meta string for the option's values.
   * \param[in] description Short description of the option.
   * \exception None.
   */
  template <typename T>
  inline void
  parse_remaining_positional_arguments(T & list,
                                       std::string const & meta_string,
                                       std::string const & description);

  /** Parses the subcommand for the program.
   * \param[in,out] subcommand Subcommand to parse.
   * \exception paw::exception::invalid_subcommand when the users does not pass a valid
   *            subcommand.
   */
  void
  parse_subcommand(std::string & subcommand);

  /** Changes the name of the program.
   * \param[in] name Name to set.
   * \exception None.
   */
  void
  set_name(std::string const & name);

  /** Sets the version of the program using two unsigned integers.
   * \param[in] major Major version number.
   * \param[in] minor Minor version number.
   * \param[in] patch Patch version number.
   * \exception None.
   */
  void
  set_version(unsigned const major,
              unsigned const minor,
              unsigned const patch);

  /** Sets the version of the program using a std::string.
   * \param[in] version Version of the program.
   * \exception None.
   */
  void set_version(std::string const & version);

  /** Generates a help message to display to the user.
   * \exception None.
   */
  std::string generate_help_message() const;

  /** Throw exceptions with help message */
  void throw_help() const;


private:
  /** Vector of all optional arguments.*/
  Args opt_args;

  /** Vector of all advanced/hidden arguments.*/
  Args opt_adv_args;

  /** Vector of all positional arguments. */
  Args pos_args;

  /** The name of the program. */
  std::string program_name{};

  /** The version of the program. */
  std::string version{};

  /** Keys that maps options (or 'flags') to their given values as given by the user.*/
  FlagMap flag_map;

  /** Indicator whether there was a missing positional argument. */
  bool missing_positional_argument{false};

  /** Vector of all values of positional arguments, in the same order as they were inserted.*/
  std::vector<std::string> positional{};

  /** Iterator that points to the next positional argument.*/
  long next_positional{0};

  /** Copy of the arguments parsed.*/
  std::vector<std::string> raw_args;

  /** Vector of all subcommands of the program. The subcommands are added using the
   * paw::Parser::add_subcommands() method.
   */
  std::vector<std::pair<std::string, std::string> > subcommands;

  /** Copy of the subcommand used */
  std::string subcommand;

  /** \brief Finds flag in parsed argument options.
   * \details Returns iterator pointing to the end if not found.
   * \returns Iterator to the flag
   * \exception None.
   */
  FlagMap::iterator find_flag(char const shrt, std::string const & lng);

};


/**********************
 * TEMPLATE FUNCTIONS *
 **********************/


template <typename T>
inline void
Parser::parse_option(T & val,
                     char const shrt,
                     std::string const & lng,
                     std::string const & description,
                     std::string const & meta_string)
{
  {
    std::ostringstream ss;
    ss << val;
    paw::Parser::Arg arg = {shrt, lng, description, meta_string, ss.str()};
    opt_args.push_back(std::move(arg));
  }

  auto flag_it = find_flag(shrt, lng);

  if (flag_it != flag_map.end())
  {
    if (flag_it->second.size() == 0)
      throw paw::exception::missing_value(shrt, lng, this->generate_help_message());

    std::string value;

    // values that could have been positional start with 2 (integer)
    if (flag_it->second[0] == static_cast<char>(2))
    {
      std::string index_str(std::next(flag_it->second.begin()), flag_it->second.end());
      long index{0};

      std::istringstream ss_i {
        index_str
      };
      ss_i >> index;

      assert(index < static_cast<long>(positional.size()));

      value = positional[index];
      positional[index].clear();
    }
    else
    {
      value = flag_it->second;
    }

    assert(value.size() > 0);

    std::istringstream ss {
      value
    };
    ss >> val;

    // Check if there were any logical errors
    if (ss.fail() || !ss.eof())
    {
      throw paw::exception::invalid_option_value(shrt,
                                                 lng,
                                                 ss.str(),
                                                 this->generate_help_message());
    }
  }
}


/** Explicit specialization of paw::Parser::parse_option() with T as type bool.
 * Here we only check if the argument was passed or not and except no value.
 */
template <>
inline void
Parser::parse_option(bool & val,
                     char const shrt,
                     std::string const & lng,
                     std::string const & description,
                     std::string const & /*mega_string makes no sense for booleans*/)
{
  Arg arg; // = {shrt, lng, description, meta_string, std::string("")};
  arg.shrt = shrt;
  arg.lng = lng;
  arg.description = description;
  arg.meta_string = paw::internal::OPTION_HAS_NO_VALUE;
  arg.default_value = "";
  opt_args.push_back(std::move(arg));
  auto flag_it = find_flag(shrt, lng);

  if (flag_it != flag_map.end())
    val ^= true; // Flips value
}


template <typename T>
inline void
Parser::parse_advanced_option(T & val,
                              char const shrt,
                              std::string const & lng,
                              std::string const & description,
                              std::string const & meta_string)
{
  {
    std::ostringstream ss;
    ss << val;
    paw::Parser::Arg arg = {shrt, lng, description, meta_string, ss.str()};
    opt_adv_args.push_back(std::move(arg));
  }

  auto flag_it = find_flag(shrt, lng);

  if (flag_it != flag_map.end())
  {
    if (flag_it->second.size() == 0)
      throw paw::exception::missing_value(shrt, lng, this->generate_help_message());

    std::string value;

    // values that could have been positional start with 2 (integer)
    if (flag_it->second[0] == static_cast<char>(2))
    {
      std::string index_str(std::next(flag_it->second.begin()), flag_it->second.end());
      long index{0};

      std::istringstream ss_i {
        index_str
      };
      ss_i >> index;

      assert(index < static_cast<long>(positional.size()));

      value = positional[index];
      positional[index].clear();
    }
    else
    {
      value = flag_it->second;
    }

    assert(value.size() > 0);

    std::istringstream ss {
      value
    };
    ss >> val;

    // Check if there were any logical errors
    if (ss.fail() || !ss.eof())
    {
      throw paw::exception::invalid_option_value(shrt,
                                                 lng,
                                                 ss.str(),
                                                 this->generate_help_message());
    }
  }
}


template <>
inline void
Parser::parse_advanced_option(bool & val,
                              char const shrt,
                              std::string const & lng,
                              std::string const & description,
                              std::string const & /*mega_string makes no sense for booleans*/)
{
  Arg arg; // = {shrt, lng, description, meta_string, std::string("")};
  arg.shrt = shrt;
  arg.lng = lng;
  arg.description = description;
  arg.meta_string = paw::internal::OPTION_HAS_NO_VALUE;
  arg.default_value = "";
  opt_adv_args.push_back(std::move(arg));
  auto flag_it = find_flag(shrt, lng);

  if (flag_it != flag_map.end())
    val ^= true; // Flips value
}


template <typename T>
inline void
Parser::parse_option_list(T & list,
                          char const shrt,
                          std::string const & lng,
                          std::string const & description,
                          char const delimiter,
                          std::string const & meta_string)
{
  /* Lambda help function that takes in a range and adds all items in the range to a container
   * with an 'insert' method. */
  auto parse_option_range_to_list_lambda =
    [&](std::pair<FlagMap::iterator, FlagMap::iterator> it_range)
    {
      while (it_range.first != it_range.second)
      {
        if (it_range.first->second.size() == 0)
          throw paw::exception::missing_value(shrt, lng, this->generate_help_message());

        std::string::iterator begin = it_range.first->second.begin();
        std::string::iterator end = it_range.first->second.end();
        std::string::iterator delim = std::find(begin, end, delimiter);

        while (begin != end)
        {
          {
            typename T::value_type val;
            std::istringstream ss {
              std::string(begin, delim)
            };
            ss >> val;

            // Check if there were any logical errors
            if (ss.fail() || !ss.eof())
            {
              throw paw::exception::invalid_option_value(shrt,
                                                         lng,
                                                         ss.str(),
                                                         this->generate_help_message());
            }

            list.insert(list.end(), val); // insert new value back to list
          }

          if (delim == end)
            break;

          begin = delim + 1;
          delim = std::find(begin, end, delimiter);
        }

        ++it_range.first;
      }
    };

  opt_args.push_back({shrt, lng, description, meta_string, ""});

  // Parse long option flag
  parse_option_range_to_list_lambda(flag_map.equal_range(lng));

  // Parse short option flag
  parse_option_range_to_list_lambda(flag_map.equal_range(std::string(1, shrt)));
}


template <typename T>
inline void
Parser::parse_positional_argument(T & val,
                                  std::string const & meta_string,
                                  std::string const & description)
{
  pos_args.push_back({paw::internal::NO_SHORT_OPTION, "", description, meta_string, ""});

  while (next_positional < static_cast<long>(positional.size()) &&
         positional[next_positional].size() == 0)
  {
    ++next_positional;
  }

  if (next_positional >= static_cast<long>(positional.size()))
  {
    missing_positional_argument = true;
    return;
  }

  std::istringstream ss {
    positional[next_positional]
  };
  ss >> val;

  // Check if there were any logical errors
  if (ss.fail() || !ss.eof())
    throw paw::exception::invalid_positional(ss.str(), this->generate_help_message());

  ++next_positional;
}


template <typename T>
inline void
Parser::parse_remaining_positional_arguments(T & list,
                                             std::string const & meta_string,
                                             std::string const & description)
{
  pos_args.push_back({paw::internal::NO_SHORT_OPTION,
                      "",
                      description,
                      meta_string,
                      ""});

  while (next_positional < positional.size())
  {
    if (positional[next_positional].size() == 0)
    {
      ++next_positional;
      continue;
    }

    typename T::value_type val;
    std::istringstream ss {
      positional[next_positional]
    };
    ss >> val;
    list.insert(list.end(), val);

    // Check if there were any logical errors
    if (ss.fail() || !ss.eof())
      throw paw::exception::invalid_positional(ss.str(), this->generate_help_message());

    ++next_positional;
  }
}


} // namespace paw


#if defined(IMPLEMENT_PAW)


namespace paw
{

/******************
 * IMPLEMENTATION *
 ******************/

/* PUBLIC METHODS */
Parser::Parser(int argc, char ** argv) :
  Parser(std::vector<std::string>(argv, argv + argc))
{}


Parser::Parser(std::vector<std::string> const & arguments)
  : raw_args(arguments)
{
  // raw_args is assumed to be never empty
  if (raw_args.size() == 0)
    raw_args.push_back("<program>");

  for (auto arg_it = raw_args.cbegin() + 1; arg_it != raw_args.cend(); ++arg_it)
  {
    if (arg_it->size() <= 1 || (*arg_it)[0] != '-')
    {
      positional.push_back(std::move(*arg_it));
    }
    else if (arg_it->size() == 2 && (*arg_it)[1] == '-')
    {
      // Force all remaining arguments to be positional arguments
      std::move(std::next(arg_it), raw_args.cend(), std::back_inserter(positional));
      break;
    }
    else if ((*arg_it)[1] == '-')
    {
      // long option flag was used
      auto it = std::find(arg_it->begin(), arg_it->end(), '=');
      std::string flag(arg_it->begin() + 2, it);

      if (it != arg_it->end())
      {
        // Equal sign was used
        flag_map.emplace(std::move(flag), std::string(it + 1, arg_it->end()));
      }
      else
      {
        // No equal sign used, check if the next argument starts with '-'
        auto next_it = std::next(arg_it);

        if (next_it == raw_args.cend() ||
            next_it->size() == 0 ||
            (next_it->size() > 1 && (*next_it)[0] == '-'))
        {
          // next argument starts with '-', there is no value
          flag_map.emplace(std::move(flag), std::string(""));
        }
        else
        {
          // Add the next argument as a positional, but prepare it to be a value as well
          std::string prefix(1, static_cast<char>(2));
          assert(prefix.size() == 1);
          assert(prefix[0] == 2);
          prefix.append(std::to_string(positional.size()));
          flag_map.emplace(flag, prefix);
          positional.push_back(*next_it);
          arg_it = next_it;
        }
      }
    }
    else
    {
      std::string flag(std::string(arg_it->begin() + 1, arg_it->begin() + 2));

      // short option flag was used, e.g. '-d' or with value '-dtab'
      if (arg_it->size() > 2)
      {
        flag_map.emplace(std::move(flag), std::string(arg_it->begin() + 2, arg_it->end()));
      }
      else
      {
        // No equal sign used, check if the next argument starts with '-'
        auto next_it = std::next(arg_it);

        if (next_it == raw_args.cend() ||
            next_it->size() == 0 ||
            (next_it->size() > 1 && (*next_it)[0] == '-'))
        {
          flag_map.emplace(std::move(flag), std::string(""));
        }
        else
        {
          std::string prefix(1, static_cast<char>(2));
          assert(prefix.size() == 1);
          assert(prefix[0] == 2);
          prefix.append(std::to_string(positional.size()));
          flag_map.emplace(flag, prefix);

          positional.push_back(*next_it);
          arg_it = next_it;
        }
      }
    }
  }
}


void
Parser::add_subcommand(std::string const & subcommand_name, std::string const & description)
{
  subcommands.push_back({subcommand_name, description});
}


void
Parser::check_for_invalid_options()
{
  // Find all known options
  std::unordered_set<std::string> known_options;

  for (auto const & arg : opt_args)
  {
    if (arg.shrt != paw::internal::NO_SHORT_OPTION)
      known_options.insert(std::string(1, arg.shrt));

    known_options.insert(arg.lng);
  }

  // Check if any parsed arguments are not defined
  for (auto it = flag_map.begin(); it != flag_map.end(); ++it)
  {
    if (known_options.find(it->first) == known_options.end())
    {
      throw paw::exception::invalid_option(it->first, this->generate_help_message());
    }
  }
}


Parser::FlagMap const &
Parser::get_flag_map_reference() const
{
  return flag_map;
}


std::vector<std::pair<std::string, std::string> > const &
Parser::get_subcommands_reference() const
{
  return subcommands;
}


std::vector<std::string> const &
Parser::get_positional_arguments_reference() const
{
  return positional;
}


void
Parser::parse_subcommand(std::string & subcommand)
{
  parse_positional_argument(subcommand,
                            "subcommand",
                            "Subcommand to execute. "
                            "List of available subcommands are shown in the following section.");

  if (subcommands.size() == 0)
    throw paw::exception::invalid_subcommand(subcommand, this->generate_help_message());

  this->subcommand = subcommand;
}


void
Parser::set_name(std::string const & name)
{
  program_name = name;
}


void
Parser::set_version(unsigned const major, unsigned const minor, unsigned const patch)
{
  std::ostringstream ss;
  ss << major << '.' << minor << '.' << patch;
  version = ss.str();
}


void
Parser::set_version(std::string const & version)
{
  this->version = version;
}


void
Parser::finalize(bool const allow_no_arguments)
{
  // Parse help argument
  bool help_flag = false;
  this->parse_option(help_flag, 'h', "help", "Show this help.");

  if (help_flag || (!allow_no_arguments && raw_args.size() <= 1))
    throw_help();

  // Check if there are missing subcommands
  if (subcommands.size() > 0 && subcommand.size() > 0)
  {
    auto find_subcommand_it =
      std::find_if(subcommands.cbegin(),
                   subcommands.cend(),
                   [this](std::pair<std::string, std::string> const & cmd)
      {
        return cmd.first == this->subcommand;
      });

    if (find_subcommand_it == subcommands.cend())
    {
      throw paw::exception::invalid_subcommand(subcommand,
                                               this->generate_help_message());
    }
  }

  if (missing_positional_argument)
  {
    assert(next_positional < static_cast<long>(pos_args.size()));
    throw paw::exception::missing_positional_argument(pos_args[next_positional].meta_string,
                                                      this->generate_help_message());
  }

  this->check_for_invalid_options();
}


std::string
Parser::generate_help_message() const
{

  using paw::internal::print_string;

  std::ostringstream ss;
  std::string const & binary_name = this->raw_args[0];
  std::size_t constexpr INDENT = 3;
  std::size_t constexpr MAX_WIDTH = 80;

  // Name section
  if (program_name.size() > 0)
  {
    ss << "\nNAME";
    print_string(ss, program_name, INDENT, MAX_WIDTH);
    ss << "\n\n";
  }

  // Usage section
  ss << "USAGE";

  {
    std::ostringstream usage_ss;
    usage_ss << binary_name;

    if (subcommands.size() > 0)
    {
      if (subcommand.size() > 0)
        usage_ss << ' ' << subcommand;
      else
        usage_ss << " <subcommand>";
    }

    if (pos_args.size() > 0)
    {
      auto pos_arg_it = pos_args.cbegin();

      if (subcommands.size() > 0)
        ++pos_arg_it;

      while (pos_arg_it != pos_args.cend())
      {
        usage_ss << " <" << pos_arg_it->meta_string << ">";
        ++pos_arg_it;
      }
    }

    usage_ss << " [OPTIONS]";
    print_string(ss, usage_ss.str(), INDENT, MAX_WIDTH);
  }

  ss << "\n";

  // Positional arguments
  if (pos_args.size() > 0)
  {
    auto pos_arg_it = pos_args.cbegin();

    // If a subcommand was used, skip the first positional argument
    if (subcommand.size() > 0)
      ++pos_arg_it;

    while (pos_arg_it != pos_args.end())
    {
      if (pos_arg_it->description.size() == 0)
        continue;

      ss << '\n' << std::string(INDENT, ' ') << '<' << pos_arg_it->meta_string << '>';
      print_string(ss, pos_arg_it->description, INDENT * 2, MAX_WIDTH);
      ss << '\n';
      ++pos_arg_it;
    }
  }

  // Subcommands section
  if (subcommands.size() > 0)
  {
    ss << "\nSUBCOMMANDS";

    long max_subcmd_size = 0;

    for (auto const & subcmd : subcommands)
    {
      if (static_cast<long>(subcmd.first.size()) > max_subcmd_size)
        max_subcmd_size = subcmd.first.size();
    }

    for (auto const & subcmd : subcommands)
    {
      long const num_spaces = max_subcmd_size + 1 - subcmd.first.size();
      print_string(ss, subcmd.first + std::string(num_spaces, ' ') + subcmd.second, INDENT, MAX_WIDTH);
    }

    ss << '\n';
  }

  // Options section
  if (opt_args.size() > 0)
  {
    ss << "\nOPTIONS";

    for (auto const & arg : opt_args)
    {
      std::ostringstream opt_ss;
      opt_ss << "--" << arg.lng;

      if (arg.meta_string.size() > 0 && arg.meta_string != paw::internal::OPTION_HAS_NO_VALUE)
        opt_ss << "=" << arg.meta_string;

      if (arg.shrt != paw::internal::NO_SHORT_OPTION)
      {
        opt_ss << " or -" << arg.shrt;

        if (arg.meta_string.size() > 0 && arg.meta_string != paw::internal::OPTION_HAS_NO_VALUE)
          opt_ss << arg.meta_string;
      }

      if (arg.default_value.size() > 0)
        opt_ss << " [default: " << arg.default_value << "]";

      print_string(ss, opt_ss.str(), INDENT, MAX_WIDTH);
      print_string(ss, arg.description, INDENT * 2, MAX_WIDTH);
      ss << '\n';
    }
  }

  // Version number
  if (version.size() > 0)
    ss << "\nVERSION\n" << std::string(INDENT, ' ') << version << "\n";

  return ss.str();
}


void
Parser::throw_help() const
{
  throw paw::exception::help(this->generate_help_message());
}


/* PRIVATE METHODS */
Parser::FlagMap::iterator
Parser::find_flag(char const shrt, std::string const & lng)
{
  // First check long options
  auto flag_it = flag_map.find(lng);

  if (flag_it != flag_map.end())
    return flag_it; // Found the flag

  // Then check short option
  return flag_map.find(std::string(1, shrt));
}


} // namespace paw


#endif // #define IMPLEMENT_PAW
