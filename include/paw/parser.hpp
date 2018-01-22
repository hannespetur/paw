#pragma once

#include <algorithm> // std::find
#include <exception> // std::exception
#include <iomanip> // std::setw
#include <ios> // std::ios::failbit
#include <map> // std::multimap
#include <unordered_set> // std::unordered_set
#include <sstream> // std::istringstream, std::ostringstream
#include <string> // std::string
#include <vector> // std::vector

/** \mainpage
 * paw parser is a small, simple, and explicit command-line argument parser.
 * \author Hannes P Eggertsson
 * \copyright GNU GPLv3
 */

/** \file paw/parser.hpp
 * A single header file with the entire paw parser.
 */

/** Top-level namespace of the paw library.*/
namespace paw
{

/** Namespace for internal variables and function used by the paw library. */
namespace internal
{

void inline
print_string(std::ostringstream & ss,
             std::string const & str,
             std::size_t const base_indent,
             std::size_t const MAX_WIDTH
             );

} // internal


/** Main Data Structure for paw parser.*/
class parser
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
     * No option shown if it is paw::parser::NO_SHORT_OPTION.
     */
    char shrt;
    std::string lng;     /**< Long option for the 'help' page.*/
    std::string description;     /**< Description of the argument to display in the 'help' page.*/
    std::string meta_string;     /**< String to represent the value with in the 'help' page.*/
    std::string default_value;     /**< The default value to display in the 'help' page.*/
  };


  /** Type definition of the container used to store the arguments.*/
  using Args = std::vector<Arg>;
  using Subcommands = std::vector<std::pair<std::string, std::string> >;

  /* READ ONLY VARIABLES */
  /** Specifies how to represent arguments that have no short option.*/
  char static const NO_SHORT_OPTION = ' ';

  /** Specifies how the default seperator that splits argument lists.*/
  char static const DEFAULT_LIST_DELIMITER = ',';

  /** Specifies how to indicate the option has no value. */
  std::string const OPTION_HAS_NO_VALUE = "__OPTIONS_HAS_NO_VALUE__";

  /* EXCEPTIONS */
  /** Exception indicating that the user passed the help argument.*/
  class help_exception : public std::exception
  {
private:
    std::string help_message;     /**< Message to display when exception is thrown.*/

    /** Constructs the help exception message.
      * \param[in] program_name Name of the program.
      * \param[in] binary_name Name of the binary.
      * \param[in] opt_args Options to display in help page.
      * \param[in] pos_args Positional arguments to display in help page.
      * \param[in] version Version of the program.
      * \param[in] subcommands Available subcommands.
      * \param[in] subcommand Subcommand used, if any.
      */
public:
    help_exception(std::string const & help_message);

    /** Gets the help message.
     * \returns Help message.
     */
    virtual const char * what() const throw();
  };

  /** Exception indicating that the user passed an argument that isn't available.
   * E.g. user parser a "-t" option but no "-t" argument is defiend by the program.
   */
  class invalid_option_exception : public std::exception
  {
private:
    std::string error_message;     /**< Message to display when exception is thrown.*/

public:
    /** Constructs an invalid option exception.
     * \param[in] invalid_option the parsed invalid option.
     * \param[in] help_message the help message to display with the exception.
     */
    invalid_option_exception(std::string const & invalid_option,
                             std::string const & help_message
                             );

    /** Gets the error message of this exception.
     * \returns Error message.
     */
    virtual const char * what() const throw();
  };

  /** Exception indicating that the user passed an argument with value of the incorrect type.
   * E.g. user parser a "-ta" option but t takes in an integer.
   */
  class invalid_option_value_exception : public std::exception
  {
private:
    std::string error_message;     /**< Message to display when exception is thrown. */

public:
    /** Constructs an invalid option value exception.
     * \param[in] shrt Short option.
     * \param[in] lng Long option.
     * \param[in] invalid_option_value the parsed invalid option value.
     * \param[in] help_message the help message to display with the exception.
     */
    invalid_option_value_exception(char const shrt,
                                   std::string const & lng,
                                   std::string const & invalid_option_value,
                                   std::string const & help_message
                                   );
    /** Gets the error message of this exception.
     * \returns Error message.
     */
    virtual const char * what() const throw();
  };

  /** Exception indicating that the user passed a positional argument of the incorrect type.*/
  class invalid_positional_exception : public std::exception
  {
private:
    std::string error_message;     /**< Message to display when exception is thrown. */

public:
    /** Constructs an invalid positional exception.
     * \param[in] invalid_positional the parsed invalid positional value.
     * \param[in] help_message the help message to display with the exception.
     */
    invalid_positional_exception(std::string const & invalid_positional,
                                 std::string const & help_message
                                 );

    /** Gets the error message of this exception.
     * \returns Error message.
     */
    virtual const char * what() const throw();
  };

  /** Exception indicating an invalid subcommand
   * Thrown when the user passes a subcommand which had not been defined.
   */
  class invalid_subcommand_exception : public std::exception
  {
    std::string error_message;   /**< Message to display when exception is thrown. */

public:
    /** Constructs an invalid subcommand exception
     * \param[in] subcommand Subcommand the user passed.
     * \param[in] help_message the help message to display with the exception.
     */
    invalid_subcommand_exception(std::string const & subcommand,
                                 std::string const & help_message
                                 );

    /** Gets the error message of this exception.
     * \returns Error message.
     */
    virtual const char * what() const throw();
  };

  /** Exception indicating a missing value.
   * Thrown when the user passes an argument without a value, but the argument requires it.
   */
  class missing_value_exception : public std::exception
  {
    std::string error_message;   /**< Message to display when exception is thrown. */

public:
    /** Constructs a missing value exception.
     * \param[in] shrt short option.
     * \param[in] lng long option.
     * \param[in] help_message the help message to display with the exception.
     */
    missing_value_exception(char const shrt, std::string const & lng,
                            std::string const & help_message
                            );

    /** Gets the error message of this exception.
     * \returns Error message.
     */
    virtual const char * what() const throw();
  };


  /** Exception indicating that the parser was expecting a postional argument,
   * but was not passed any.
   */
  class missing_positional_argument_exception : public std::exception
  {
private:
    std::string error_message;     /**< Message to display when exception is thrown. */

public:
    /** Constructs a missing positional argument exception
     * \param[in] meta_string Meta string of the positional argument.
     * \param[in] help_message the help message to display with the exception.
     */
    missing_positional_argument_exception(std::string const & meta_string,
                                          std::string const & help_message
                                          );

    /** Gets the error message of this exception.
     * \returns Error message.
     */
    virtual const char * what() const throw();

  };

  /* TYPE DEFINITIONS */
  /** Type definition for the hash table container to use for parsing arguments.*/
  using FlagMap = std::multimap<std::string, std::string>;

public:
  /** \brief Construction of a paw parser instance.
   * \details This constructor uses the argc and argv variables that are passed to the 'main'
   *          function of the program.
   * \param[in] argc length of the argv array.
   * \param[in] argv arguments passed, including the name of the executable.
   * \returns paw parser object.
   * \exception None.
   */
  parser(int argc, char ** argv);

  /** \brief Construction of a paw parser instance.
   * \details This constructor parses a vector of strings that should include the programs name.
   * \param[in] args Dynamic vector with all arguments.
   * \returns paw parser object.
   * \exception None.
   */
  parser(std::vector<std::string> const & args);

  /** \brief Adds an subcommand to the parse.
   * \details Add subcommand and use the paw::parser::parse_subcommand() to parse them.
   * \param[in] subcommand_name Name of the subcommand.
   * \param[in] description Description of the subcommand.
   * \exception None.
   */
  void add_subcommand(std::string const & subcommand_name, std::string const & description);

  /** Checks if the user passed invalid options.
   * \exception paw::parser::invalid_option_exception thrown if user passed an invalid option.
   */
  void check_for_invalid_options();

  /** \brief Finalizes the parser.
   * \details Parses help flag and checks for invalid options.
   * \param[in] allow_no_arguments Set if the parser should not throw help when no arguments
   *                               are passed.
   * \exception paw::parser::invalid_option_exception thrown if user passed an invalid option.
   */
  void finalize(bool const allow_no_arguments = true);


  /** Get a reference to the multimap containing all the flags the user passed.
   * \returns Multimap with flags.
   * \exception None.
   */
  FlagMap const &
  get_flag_map_reference();

  /** Parses an option that was passed by the user.
   * If 'T' is not a boolean, we expect that the argument has a single value.
   * \param[in,out] val Reference to the value to change if the option was parsed.
   * \param[in] shrt Short option.
   * \param[in] lng Long option.
   * \param[in] description Short description of the option.
   * \param[in] meta_string Meta string for the option's value.
   * \exception paw::parser::missing_value_exception thrown when option was parsed
   *            but without a value.
   */
  template <typename T>
  inline void
  parse_option(T & val,
               char const shrt,
               std::string const & lng,
               std::string const & description,
               std::string const & meta_string = "value"
               );

  /** Parses an argument with a list of values.
   * \param[in,out] list Reference to the list to change if the option was parsed.
   * \param[in] shrt Short option.
   * \param[in] lng Long option.
   * \param[in] description Short description of the option.
   * \param[in] delimiter Delimiter between items in the list.
   * \param[in] meta_string Meta string for the option's value.
   * \exception paw::parser::missing_value_exception thrown when option list was parsed,
   *            but without a value.
   */
  template <typename T>
  inline void
  parse_option_list(T & list,
                    char const shrt,
                    std::string const & lng,
                    std::string const & description,
                    char const delimiter = DEFAULT_LIST_DELIMITER,
                    std::string const & meta_string = "value1,..."
                    );

  /** Parses the next positional argument.
   * \param[in,out] val Reference to the value to change if the option was parsed.
   * \param[in] meta_string Meta string for the option's value.
   * \param[in] description Short description of the option.
   * \exception paw::parser::missing_positional_argument_exception thrown when there was no
   *            positional argument passed.
   */
  template <typename T>
  inline void
  parse_positional_argument(T & val,
                            std::string const & meta_string,
                            std::string const & description
                            );

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
                                       std::string const & description
                                       );

  /** Parses the subcommand for the program.
   * \param[in,out] subcommand Subcommand to parse.
   * \exception paw::parser::invalid_subcommand_exception when the users does not pass a valid
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
   * \param[in] MAJOR major version.
   * \param[in] MINOR minor version.
   * \exception None.
   */
  void
  set_version(unsigned const MAJOR,
              unsigned const MINOR
              );

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
  /** Type definition for the container to use for positional arguments.*/
  using Positional = std::vector<std::string>;

  /** Vector of all optional arguments.*/
  Args opt_args;

  /** Vector of all positional arguments. */
  Args pos_args;

  /** The name of the program. */
  std::string program_name;

  /** The version of the program. */
  std::string version;

  /** Keys that maps options (or 'flags') to their given values as given by the user.*/
  FlagMap flag_map;

  /** Indicator whether there was a missing positional argument. */
  bool missing_positional_argument = false;

  /** Indicator whether the program has subcommands. */
  bool has_subcommands = false;

  /** Vector of all values of positional arguments, in the same order as they were inserted.*/
  Positional positional;

  /** Iterator that points to the next positional argument.*/
  Positional::iterator next_positional;

  /** Copy of the arguments parsed.*/
  std::vector<std::string> raw_args;

  /** Vector of all subcommands of the program. The subcommands are added using the
   * paw::parser::add_subcommands() method.
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
inline
void
parser::parse_option(T & val,
                     char const shrt,
                     std::string const & lng,
                     std::string const & description,
                     std::string const & meta_string
                     )
{
  {
    std::ostringstream ss;
    ss << val;
    paw::parser::Arg arg = {shrt, lng, description, meta_string, ss.str()};
    opt_args.push_back(std::move(arg));
  }

  auto flag_it = find_flag(shrt, lng);

  if (flag_it != flag_map.end())
  {
    if (flag_it->second.size() == 0)
      throw paw::parser::missing_value_exception(shrt, lng, this->generate_help_message());

    std::istringstream ss {
      flag_it->second
    };
    ss >> val;

    // Check if there were any logical errors
    if (ss.fail() || !ss.eof())
    {
      throw paw::parser::invalid_option_value_exception(shrt,
                                                        lng,
                                                        ss.str(),
                                                        this->generate_help_message()
                                                        );
    }
  }
}


/** Explicit specialization of paw::parser::parse_option() with T as type bool.
 * Here we only check if the argument was passed or not and except no value.
 */
template <>
inline void
parser::parse_option(bool & val,
                     char const shrt,
                     std::string const & lng,
                     std::string const & description,
                     std::string const & /*mega_string makes no sense for booleans*/
                     )
{
  Arg arg; // = {shrt, lng, description, meta_string, std::string("")};
  arg.shrt = shrt;
  arg.lng = lng;
  arg.description = description;
  arg.meta_string = paw::parser::OPTION_HAS_NO_VALUE;
  arg.default_value = "";
  opt_args.push_back(std::move(arg));
  auto flag_it = find_flag(shrt, lng);

  if (flag_it != flag_map.end())
    val ^= true; // Flips value
}


template <typename T>
inline void
parser::parse_option_list(T & list,
                          char const shrt,
                          std::string const & lng,
                          std::string const & description,
                          char const delimiter,
                          std::string const & meta_string
                          )
{
  /* Lambda help function that takes in a range and adds all items in the range to a container
   * with an 'insert' method. */
  auto parse_option_range_to_list_lambda =
    [&](std::pair<FlagMap::iterator, FlagMap::iterator> it_range)
    {
      while (it_range.first != it_range.second)
      {
        if (it_range.first->second.size() == 0)
          throw paw::parser::missing_value_exception(shrt, lng, this->generate_help_message());

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
              throw paw::parser::invalid_option_value_exception(shrt,
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
parser::parse_positional_argument(T & val,
                                  std::string const & meta_string,
                                  std::string const & description
                                  )
{
  pos_args.push_back({paw::parser::NO_SHORT_OPTION, "", description, meta_string, ""});

  if (next_positional == positional.end())
  {
    missing_positional_argument = true;
    return;
  }

  std::istringstream ss {
    *next_positional
  };
  ss >> val;

  // Check if there were any logical errors
  if (ss.fail() || !ss.eof())
    throw paw::parser::invalid_positional_exception(ss.str(), this->generate_help_message());

  ++next_positional;
}


template <typename T>
inline void
parser::parse_remaining_positional_arguments(T & list,
                                             std::string const & meta_string,
                                             std::string const & description
                                             )
{
  pos_args.push_back({paw::parser::NO_SHORT_OPTION,
                      "",
                      description,
                      meta_string,
                      ""}
                     );

  while (next_positional != positional.end())
  {
    typename T::value_type val;
    std::istringstream ss {
      *next_positional
    };
    ss >> val;
    list.insert(list.end(), val);

    // Check if there were any logical errors
    if (ss.fail() || !ss.eof())
      throw paw::parser::invalid_positional_exception(ss.str(), this->generate_help_message());

    ++next_positional;
  }
}


#ifdef IMPLEMENT

/******************
 * IMPLEMENTATION *
 ******************/

namespace internal
{

void
print_string(std::ostringstream & ss,
             std::string const & str,
             std::size_t const base_indent,
             std::size_t const MAX_WIDTH
             )
{
  // Write indentation
  ss << "\n" << std::string(base_indent, ' ');
  std::size_t LINE_SIZE = base_indent;
  auto str_it = str.begin();

  while (str_it != str.end())
  {
    auto whitespace_it = std::find(str_it, str.end(), ' ');
    auto newline_it = std::find(str_it, str.end(), '\n');

    // If newline comes first, print until the newline
    if (std::distance(str_it, newline_it) <
        std::distance(str_it, whitespace_it)
        )
    {
      whitespace_it = newline_it;
    }

    if (std::distance(str_it, whitespace_it) + LINE_SIZE < MAX_WIDTH)
    {
      ss << std::string(str_it, whitespace_it);
      LINE_SIZE += std::distance(str_it, whitespace_it);
    }
    else
    {
      ss << '\n' << std::string(base_indent, ' ')
         << std::string(str_it, whitespace_it);
      LINE_SIZE = base_indent + std::distance(str_it, whitespace_it);
    }

    str_it = whitespace_it;

    if (str_it != str.end())
    {
      if (*str_it == '\n')
      {
        ss << '\n' << std::string(base_indent, ' ');
        LINE_SIZE = base_indent;
      }
      else
      {
        ss << ' ';
        ++LINE_SIZE;
      }

      ++str_it;
    }
  }
}


} // paw::internal


/* EXCEPTIONS */
parser::help_exception::help_exception(std::string const & _help_message)
{
  help_message = _help_message;
}


const char *
parser::help_exception::what() const throw()
{
  return help_message.c_str();
}


/* Invalid option exception */
parser::invalid_option_exception::invalid_option_exception(std::string const & invalid_option,
                                                           std::string const & help_message
                                                           )
{
  std::ostringstream ss;
  ss << "[paw::parser::InvalidOption] ERROR: Unknown option '"
     << invalid_option
     << "' was passed.\n"
     << help_message;
  error_message = ss.str();
}


const char *
parser::invalid_option_exception::what() const throw()
{
  return error_message.c_str();
}


/* Invalid option value exception */
parser::invalid_option_value_exception::invalid_option_value_exception(
  char const shrt,
  std::string const & lng,
  std::string const & invalid_option_value,
  std::string const & help_message
  )
{
  std::ostringstream ss;

  if (shrt != paw::parser::NO_SHORT_OPTION)
  {
    ss << "[paw::parser::InvalidOptionValue] Value '"
       << invalid_option_value
       << "' for option '"
       << lng
       << "' (or '-"
       << shrt
       << "') is invalid.";
  }
  else
  {
    ss << "[paw::parser::InvalidOptionValue] Value '"
       << invalid_option_value
       << "' for option '"
       << lng
       << "' is invalid.";
  }

  ss << "\n" << help_message;
  error_message = ss.str();
}


const char *
parser::invalid_option_value_exception::what() const throw()
{
  return error_message.c_str();
}


/* Invalid positional exception */
parser::invalid_positional_exception::invalid_positional_exception(
  std::string const & invalid_positional,
  std::string const & help_message
  )
{
  std::ostringstream ss;
  ss << "[paw::parser::InvalidPositional] Positional argument <"
     << invalid_positional
     << "> is of invalid type.\n"
     << help_message;

  error_message = ss.str();
}


const char *
parser::invalid_positional_exception::what() const throw()
{
  return error_message.c_str();
}


/* Invalid subcommand exception */
parser::invalid_subcommand_exception::invalid_subcommand_exception(
  std::string const & subcommand,
  std::string const & help_message
  )
{
  std::ostringstream ss;
  ss << "[paw::parser::InvalidSubcommand] Subcommand '"
     << subcommand
     << "' is invalid.\n"
     << help_message;

  error_message = ss.str();
}


const char *
parser::invalid_subcommand_exception::what() const throw()
{
  return error_message.c_str();
}


/* Missing value exception */
parser::missing_value_exception::missing_value_exception(char const shrt,
                                                         std::string const & lng,
                                                         std::string const & help_message
                                                         )
{
  std::ostringstream ss;

  if (shrt != NO_SHORT_OPTION)
  {
    ss << "[paw::parser::MissingValue] Option '--"
       << lng
       << "' (or '-"
       << shrt
       << "') requires a value, but was not passed any.";
  }
  else
  {
    ss << "[paw::parser::MissingValue] Flag '--"
       << lng
       << "' was not passed any value.";
  }

  ss << "\n" << help_message;
  error_message = ss.str();
}


const char *
parser::missing_value_exception::what() const throw()
{
  return error_message.c_str();
}


/* Missing positional argument exception */
parser::missing_positional_argument_exception::missing_positional_argument_exception(
  std::string const & meta_string,
  std::string const & help_message
  )
{
  std::ostringstream ss;
  ss << "[paw::parser::MissingPositionalArgument] ERROR: "
     << "Positional argument <"
     << meta_string
     << "> is missing.\n"
     << help_message;

  error_message = ss.str();
}


const char *
parser::missing_positional_argument_exception::what() const throw()
{
  return error_message.c_str();
}


/* PUBLIC METHODS */
parser::parser(int argc, char ** argv) :
  parser(std::vector<std::string>(argv, argv + argc))
{}


parser::parser(std::vector<std::string> const & arguments)
{
  raw_args = std::vector<std::string>(arguments);

  // raw_args is assumed to be never empty
  if (raw_args.size() == 0)
    raw_args.push_back("<program>");

  for (auto arg_it = arguments.cbegin() + 1; arg_it != arguments.cend(); ++arg_it)
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
      {
        flag_map.emplace(std::move(flag), std::string(it + 1, arg_it->end())); // Equal sign used
      }
      else
      {
        // No equal sign used, check if the next argument starts with '-'
        auto next_it = std::next(arg_it);

        if (next_it == arguments.cend() ||
            next_it->size() == 0 ||
            (next_it->size() > 1 && (*next_it)[0] == '-')
            )
        {
          flag_map.emplace(std::move(flag), std::string(""));
        }
        else
        {
          flag_map.emplace(std::move(flag), std::string(*next_it));
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

        if (next_it == arguments.cend() ||
            next_it->size() == 0 ||
            (next_it->size() > 1 && (*next_it)[0] == '-')
            )
        {
          flag_map.emplace(std::move(flag), std::string(""));
        }
        else
        {
          flag_map.emplace(std::move(flag), std::string(*next_it));
          arg_it = next_it;
        }
      }
    }
  }

  next_positional = positional.begin();
}


void
parser::add_subcommand(std::string const & subcommand_name, std::string const & description)
{
  subcommands.push_back({subcommand_name, description});
}


void
parser::check_for_invalid_options()
{
  // Find all known options
  std::unordered_set<std::string> known_options;

  for (auto const & arg : opt_args)
  {
    if (arg.shrt != NO_SHORT_OPTION)
      known_options.insert(std::string(1, arg.shrt));

    known_options.insert(arg.lng);
  }

  // Check if any parsed arguments are not defined
  for (auto it = flag_map.begin(); it != flag_map.end(); ++it)
  {
    if (known_options.find(it->first) == known_options.end())
      throw paw::parser::invalid_option_exception(it->first, this->generate_help_message());
  }
}


parser::FlagMap const &
parser::get_flag_map_reference()
{
  return flag_map;
}


void
parser::parse_subcommand(std::string & subcommand)
{
  parse_positional_argument(subcommand,
                            "subcommand",
                            "Subcommand to execute. "
                            "List of available subcommands are shown in the following section."
                            );

  if (subcommands.size() == 0)
    throw paw::parser::invalid_subcommand_exception(subcommand, this->generate_help_message());

  this->subcommand = subcommand;
}


void
parser::set_name(std::string const & name)
{
  program_name = name;
}


void
parser::set_version(unsigned const major, unsigned const minor)
{
  std::ostringstream ss;
  ss << major << '.' << minor;
  version = ss.str();
}


void
parser::set_version(std::string const & version)
{
  this->version = version;
}


void
parser::finalize(bool const allow_no_arguments)
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
      throw paw::parser::invalid_subcommand_exception(subcommand, this->generate_help_message());
  }

  std::size_t const N = std::distance(positional.begin(), next_positional);

  if (missing_positional_argument)
  {
    throw paw::parser::missing_positional_argument_exception(pos_args[N].meta_string,
                                                             this->generate_help_message()
                                                             );
  }

  this->check_for_invalid_options();
}


std::string
parser::generate_help_message() const
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
  if (subcommand.size() > 0)
  {
    ss << "\nSUBCOMMANDS";

    for (auto const & subcmd : subcommands)
    {
      print_string(ss, subcmd.first, INDENT, MAX_WIDTH);
      print_string(ss, subcmd.second, INDENT * 2, MAX_WIDTH);
      ss << '\n';
    }
  }

  // Options section
  if (opt_args.size() > 0)
  {
    ss << "\nOPTIONS";

    for (auto const & arg : opt_args)
    {
      std::ostringstream opt_ss;
      opt_ss << "--" << arg.lng;

      if (arg.meta_string.size() > 0 && arg.meta_string != paw::parser::OPTION_HAS_NO_VALUE)
        opt_ss << "=" << arg.meta_string;

      if (arg.shrt != paw::parser::NO_SHORT_OPTION)
      {
        opt_ss << " or -" << arg.shrt;

        if (arg.meta_string.size() > 0 && arg.meta_string != paw::parser::OPTION_HAS_NO_VALUE)
          opt_ss << arg.meta_string;
      }

      if (arg.default_value.size() > 0)
        opt_ss << " [default: " << arg.default_value << "]";

      print_string(ss, opt_ss.str(), INDENT, MAX_WIDTH);
      print_string(ss, arg.description, INDENT * 2, MAX_WIDTH);
      ss << "\n";
    }
  }

  // Version number
  if (version.size() > 0)
    ss << "\nVERSION\n" << std::string(INDENT, ' ') << version << "\n";

  return ss.str();
}


void
parser::throw_help() const
{
  throw paw::parser::help_exception(this->generate_help_message());
}


/* PRIVATE METHODS */
parser::FlagMap::iterator
parser::find_flag(char const shrt, std::string const & lng)
{
  // First check long options
  auto flag_it = flag_map.find(lng);

  if (flag_it != flag_map.end())
    return flag_it; // Found the flag

  // Then check short option
  return flag_map.find(std::string(1, shrt));
}


#endif // #define IMPLEMENT

} // namespace paw
