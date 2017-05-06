#pragma once
#include <algorithm> // std::find
#include <exception> // std::exception
#include <iomanip> // std::setw
#include <ios> // std::ios::failbit
#include <map> // std::multimap
#include <unordered_set> // std::unordered_set
#include <sstream> // std::stringstream, std::ostringstream
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

/** Main Data Structure for paw parser.*/
class parser
{
  public:
    /** \brief Data structure for command line arguments.
     * \details Arguments in this data structure are ready for display in the 'help' page of the
     *          program.
     */
    struct Arg
    {
      public:
        char shrt; /**< Short option for the 'help' page.
                   *< No option shown if it is paw::parser::NO_SHORT_OPTION.
                   */
        std::string lng; /**< Long option for the 'help' page.*/
        std::string description; /**< Description of the argument to display in the 'help' page.*/
        std::string meta_string; /**< String to represent the value with in the 'help' page.*/
        std::string default_value; /**< The default value to display in the 'help' page.*/

        /** Default constructor is disabled.*/
        Arg() = delete;
    };


    /** Type definition of the container used to store the arguments.*/
    using Args = std::vector<Arg>;

    /* READ ONLY VARIABLES */
    /** Specifies how to represent arguments that have no short option.*/
    static const char NO_SHORT_OPTION = ' ';
    /** Specifies how the default seperator that splits argument lists.*/
    static const char DEFAULT_LIST_DELIMITER = ',';

    /* EXCEPTIONS */
    /** Exception indicating that the user passed the help argument.*/
    class help_exception : public std::exception
    {
      private:
        std::string help_message; /**< Message to display when exception is thrown.*/

        /** Constructs the help exception message.
          * \param program_name Name of the program.
          * \param binary_name Name of the binary.
          * \param args Arguments to display in help page.
          */
      public:
        help_exception(const std::string& program_name,
                       const std::string& binary_name,
                       const Args& args
                       );

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
        std::string error_message; /**< Message to display when exception is thrown.*/

      public:
        /** Constructs an invalid option exception.
         * \param[in] invalid_option the parsed invalid option.
         */
        invalid_option_exception(const std::string& invalid_option);
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
        std::string error_message; /**< Message to display when exception is thrown. */

      public:
        /** Constructs an invalid option value exception.
         * \param[in] invalid_option_value the parsed invalid option value.
         */
        invalid_option_value_exception(const char shrt,
                                       const std::string& lng,
                                       const std::string& invalid_option_value
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
        std::string error_message; /**< Message to display when exception is thrown. */

      public:
        /** Constructs an invalid positional exception.
         * \param[in] invalid_positional the parsed invalid positional value.
         */
        invalid_positional_exception(const std::string& invalid_positional);
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
      std::string error_message; /**< Message to display when exception is thrown. */

      public:
        /** Constructs a missing value exception.
         * \param[in] shrt short option.
         * \param[in] lng long option.
         */
        missing_value_exception(const char shrt, const std::string& lng);
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
        std::string error_message; /**< Message to display when exception is thrown. */

      public:
        /** Constructs a missing positional argument exception
         * \param[in] index index of the positional argument.
         * \param[in] meta_string Meta string of the positional argument.
         */
        missing_positional_argument_exception(const size_t index,
                                              const std::string& meta_string
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
    /** Constructor uses the argc and argv variables passed to the 'main' function of the program.
     * \param[in] argc length of the argv array.
     * \param[in] argv arguments passed, including the name of the executable.
     * \returns paw parser object.
     * \exception None.
     */
    parser(int argc, char ** argv);
    /** Constructor that parses a vector of strings.
     * \param[in] args Dynamic vector with all arguments.
     * \returns paw parser object.
     * \exception None.
     */
    parser(const std::vector<std::string>& args);

    /** Checks if the user passed invalid options.
     * \exception paw::parser::invalid_option_exception thrown if user passed an invalid option.
     */
    inline void
    check_for_invalid_options();

    /** \brief Finalizes the parser.
     * \details Parses help flag and checks for invalid options.
     * \exception paw::parser::invalid_option_exception thrown if user passed an invalid option.
     */
    inline void
    finalize();

    /** Get a reference to the multimap containing all the flags the user passed.
     * \returns Multimap with flags.
     * \exception None.
     */
    inline const FlagMap &
    get_flag_map_reference();

    /** Parses an option that was passed by the user.
     * If 'T' is not a boolean, we expect that the argument has a single value.
     * \param[in,out] val reference to the value to change if the option was parsed.
     * \param[in] shrt short option.
     * \param[in] lng long option.
     * \param[in] description Short description of the option.
     * \exception paw::parser::missing_value_exception thrown when option was parsed
     *            but without a value.
     */
    template <typename T>
    inline void
    parse_option(T& val,
                 const char shrt,
                 const std::string& lng,
                 const std::string& description
                 );

    /** Parses an argument with a list of values.
     * \exception paw::parser::missing_value_exception thrown when option list was parsed,
     *            but without a value.
     */
    template <typename T>
    inline void
    parse_option_list(T& list,
                      const char shrt,
                      const std::string& lng,
                      const std::string& description,
                      const char delimiter = DEFAULT_LIST_DELIMITER
                      );

    /** Parses the next positional argument.
     * \exception paw::parser::missing_positional_argument_exception thrown when there was no
     *            positional argument passed.
     */
    template <typename T>
    inline void
    parse_positional_argument(T& val,
                              const std::string& meta_string,
                              const std::string& description
                              );

    /** Parse all remaining positional arguments.
     * \exception None.
     */
    template <typename T>
    inline void
    parse_remaining_positional_arguments(T& list,
                                         const std::string& meta_string,
                                         const std::string& description
                                         );

    /** Sets the name of the program
     * \param name Name to set.
     * \exception None.
     */
    inline void
    set_name(const std::string& name);

    /** Sets the version of the program using two unsigned integers.
     * \param MAJOR major version.
     * \param MINOR minor version.
     * \exception None.
     */
    inline void
    set_version(const unsigned MAJOR,
                const unsigned MINOR
                );

    /** Sets the version of the program using a std::string.
     * \param version version of the program.
     * \exception None.
     */
    inline void
    set_version(const std::string& version);


  private:
    /** Type definition for the container to use for positional arguments.*/
    using Positional = std::vector<std::string>;

    /** Vector of all arguments.*/
    Args args;

    /** The name of the program. */
    std::string program_name;

    /** The version of the program. */
    std::string version;

    /** Keys that maps options (or 'flags') to their given values as given by the user.*/
    FlagMap flag_map;

    /** Vector of all values of positional arguments, in the same order as they were inserted.*/
    Positional positional;

    /** Iterator that points to the next positional argument.*/
    Positional::iterator next_positional;

    /** Copy of the arguments parsed.*/
    std::vector<std::string> raw_args;

    /** \brief Finds flag in parsed argument options.
     * \details Returns iterator pointing to the end if not found.
     * \returns Iterator to the flag
     * \exception None.
     */
    inline FlagMap::iterator
    find_flag(const char shrt, const std::string& lng);

};


/* EXCEPTIONS */
parser::help_exception::help_exception(const std::string& program_name,
                                       const std::string& binary_name,
                                       const Args& arguments)
{
  std::ostringstream ss;
  size_t constexpr INDENT = 4;
  size_t constexpr MAX_WIDTH = 100;

  auto print_string =
    [&](const std::string& str, const size_t base_indent)
    {
      // Write indentation
      ss << "\n" << std::string(base_indent, ' ');
      size_t LINE_SIZE = base_indent;
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
          ss << "\n" << std::string(base_indent, ' ')
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
    };

  // Name section
  if (program_name.size() > 0)
  {
    ss << "NAME";
    print_string(program_name, INDENT);
    ss << "\n\n";
  }

  // Usage section
  ss << "USAGE";

  {
    std::ostringstream usage_ss;
    const std::size_t slash = binary_name.rfind('/');

    if (slash == std::string::npos)
      usage_ss << binary_name;
    else
      usage_ss << std::string(binary_name.begin() + slash + 1, binary_name.end());

    usage_ss << " [OPTIONS]";
    print_string(usage_ss.str(), INDENT);
    ss << '\n';
  }

  // Options section
  ss << "\nOPTIONS";

  for (const auto& arg : arguments)
  {
    ss << "\n" << std::string(INDENT, ' ');

    if (arg.shrt != paw::parser::NO_SHORT_OPTION)
    {
      ss << '-' << arg.shrt;

      if (arg.meta_string.size() > 0)
        ss << arg.meta_string;

      ss << ", ";
    }

    ss << "--" << arg.lng;

    if (arg.meta_string.size() > 0)
      ss << "=" << arg.meta_string;

    print_string(arg.description, INDENT * 2);
    ss << "\n";
  }

  help_message = ss.str();
}

const char *
parser::help_exception::what() const throw()
{
  return help_message.c_str();
}

/* Invalid option exception */
parser::invalid_option_exception::invalid_option_exception(std::string const& invalid_option)
{
  std::ostringstream ss;
  ss << "[paw::parser::InvalidOption] No option '" << invalid_option << "' is invalid.";
  error_message = ss.str();
}

const char *
parser::invalid_option_exception::what() const throw()
{
  return error_message.c_str();
}

/* Invalid option value exception */
parser::invalid_option_value_exception::invalid_option_value_exception(
  const char shrt,
  const std::string& lng,
  std::string const& invalid_option_value
  )
{
  std::ostringstream ss;

  if (shrt != paw::parser::NO_SHORT_OPTION)
  {
    ss << "[paw::parser::InvalidOptionValue] Value '" << invalid_option_value
       << "' for option '" << lng << "=value' (or '-" << shrt << "value') is invalid.";
  }
  else
  {
    ss << "[paw::parser::InvalidOptionValue] Value '" << invalid_option_value
       << "' for option '" << lng << "=value' is invalid.";
  }

  error_message = ss.str();
}

const char *
parser::invalid_option_value_exception::what() const throw()
{
  return error_message.c_str();
}

/* Invalid positional exception */
parser::invalid_positional_exception::invalid_positional_exception(
  std::string const& invalid_positional
  )
{
  std::ostringstream ss;
  ss << "[paw::parser::InvalidPositional] No positional option '" << invalid_positional
     << "' is of invalid type.";
  error_message = ss.str();
}

const char *
parser::invalid_positional_exception::what() const throw()
{
  return error_message.c_str();
}

/* Missing value exception */
parser::missing_value_exception::missing_value_exception(const char shrt,
                                                         const std::string& lng
                                                         )
{
  std::ostringstream ss;

  if (shrt != NO_SHORT_OPTION)
  {
    ss << "[paw::parser::MissingValue] Option '--"
       << lng << "=value' (or '-"
       << shrt << "value') requires a value, but was not passed any.";
  }
  else
  {
    ss << "[paw::parser::MissingValue] Flag '--"
       << lng << "=value' was not passed any value.";
  }

  error_message = ss.str();
}

const char *
parser::missing_value_exception::what() const throw()
{
  return error_message.c_str();
}

/* Missing positional argument exception */
parser::missing_positional_argument_exception::missing_positional_argument_exception
  (const size_t index,
  const std::string& meta_string
  )
{
  std::ostringstream ss;
  ss << "[paw::parser::MissingPositionalArgument] Positional argument at index "
     << index << " '" << meta_string << "' was expected but is missing.";
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


parser::parser(const std::vector<std::string>& arguments)
{
  raw_args = std::vector<std::string>(arguments);

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
        flag_map.emplace(std::move(flag), std::string(it + 1, arg_it->end())); // Equal sign used
      else
        flag_map.emplace(std::move(flag), std::string("")); // No equal sign used
    }
    else
    {
      // short option flag was used, e.g. '-d' or with value '-dtab'
      flag_map.emplace(std::string(arg_it->begin() + 1, arg_it->begin() + 2),
                       std::string(arg_it->begin() + 2, arg_it->end())
                       );
    }
  }

  next_positional = positional.begin();
}

inline void
parser::check_for_invalid_options()
{
  // Find all known options
  std::unordered_set<std::string> known_options;

  for (const auto& arg : args)
  {
    if (arg.shrt != NO_SHORT_OPTION)
      known_options.insert(std::string(1, arg.shrt));

    known_options.insert(arg.lng);
  }

  // Check if any parsed arguments are not defined
  for (auto it = flag_map.begin(); it != flag_map.end(); ++it)
  {
    if (known_options.find(it->first) == known_options.end())
      throw paw::parser::invalid_option_exception(it->first);
  }
}

inline const parser::FlagMap &
parser::get_flag_map_reference()
{
  return flag_map;
}

template <typename T>
inline void
parser::parse_option(T& val,
                     const char shrt,
                     const std::string& lng,
                     const std::string& description
                     )
{
  {
    std::stringstream ss;
    ss << val;
    paw::parser::Arg arg = {shrt, lng, description, "value", ss.str()};
    args.push_back(std::move(arg));
  }

  auto flag_it = find_flag(shrt, lng);

  if (flag_it != flag_map.end())
  {
    if (flag_it->second.size() == 0)
      throw paw::parser::missing_value_exception(shrt, lng);

    std::stringstream ss;
    ss << flag_it->second;
    ss >> val;

    // Check if there were any logical errors
    if ((ss.rdstate() & std::ios::failbit) != 0)
      throw paw::parser::invalid_option_value_exception(shrt, lng, ss.str());
  }
}

/** Explicit specialization of paw::parser::parse_option() with T as type bool.
 * Here we only check if the argument was passed or not and except no value.
 */
template <>
inline void
parser::parse_option(bool& val,
                     const char shrt,
                     const std::string& lng,
                     const std::string& description
                     )
{
  args.push_back({shrt, lng, description, "", ""});
  auto flag_it = find_flag(shrt, lng);

  if (flag_it != flag_map.end())
    val ^= true; // Flips value
}

template <typename T>
inline void
parser::parse_option_list(T& list,
                          char const shrt,
                          const std::string& lng,
                          const std::string& description,
                          const char delimiter
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
          throw paw::parser::missing_value_exception(shrt, lng);

        std::string::iterator begin = it_range.first->second.begin();
        std::string::iterator end = it_range.first->second.end();
        std::string::iterator delim = std::find(begin, end, delimiter);

        while (begin != end)
        {
          {
            typename T::value_type val;
            std::stringstream ss;
            ss << std::string(begin, delim);
            ss >> val;

            // Check if there were any logical errors
            if ((ss.rdstate() & std::ios::failbit) != 0)
              throw paw::parser::invalid_option_value_exception(shrt, lng, ss.str());

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

  args.push_back({shrt, lng, description, "list...", ""});

  // Parse long option flag
  parse_option_range_to_list_lambda(flag_map.equal_range(lng));

  // Parse short option flag
  parse_option_range_to_list_lambda(flag_map.equal_range(std::string(1, shrt)));
}

template <typename T>
inline void
parser::parse_positional_argument(T& val,
                                  const std::string& meta_string,
                                  const std::string& description
                                  )
{
  args.push_back({paw::parser::NO_SHORT_OPTION, "", description, meta_string, ""});
  const std::size_t N = std::distance(positional.begin(), next_positional);

  if (next_positional == positional.end())
  {
    throw paw::parser::missing_positional_argument_exception(N, meta_string);
  }

  std::stringstream ss;
  ss << *next_positional;
  ss >> val;

  // Check if there were any logical errors
  if ((ss.rdstate() & std::ios::failbit) != 0)
  {
    throw paw::parser::invalid_option_value_exception(paw::parser::NO_SHORT_OPTION,
                                                      std::to_string(N),
                                                      ss.str()
                                                      );
  }

  ++next_positional;
}

template <typename T>
inline void
parser::parse_remaining_positional_arguments(T& list,
                                             const std::string& meta_string,
                                             const std::string& description
                                             )
{
  while (next_positional != positional.end())
  {
    typename T::value_type val;
    std::stringstream ss;
    ss << *next_positional;
    ss >> val;
    list.insert(list.end(), val);
    ++next_positional;
  }
}

inline void
parser::set_name(const std::string& name)
{
  program_name = name;
}

inline void
parser::set_version(const unsigned major, const unsigned minor)
{
  std::ostringstream ss;
  ss << major << "." << minor;
  version = ss.str();
}

inline void
parser::set_version(const std::string& version)
{
  this->version = version;
}

inline void
parser::finalize()
{
  // Parse help argument
  bool help_flag = false;
  this->parse_option(help_flag, 'h', "help", "Show this help.");

  if (help_flag)
    throw paw::parser::help_exception(this->program_name, this->raw_args[0], this->args);

  this->check_for_invalid_options();
}

/* PRIVATE METHODS */
inline parser::FlagMap::iterator
parser::find_flag(const char shrt, const std::string& lng)
{
  // First check long options
  auto flag_it = flag_map.find(lng);

  if (flag_it != flag_map.end())
    return flag_it; // Found the flag

  // Then check short option
  return flag_map.find(std::string(1, shrt));
}

} // namespace paw
