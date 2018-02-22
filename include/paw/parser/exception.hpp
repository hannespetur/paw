#pragma once

#include <exception> // std::exception
#include <string> // std::string

#include <paw/parser/internal.hpp>


namespace paw
{
namespace exception
{

/* EXCEPTIONS */
/** Exception indicating that the user passed the help argument.*/
class help : public std::exception
{
private:
  std::string help_message;       /**< Message to display when exception is thrown.*/

  /** Constructs the help exception message.
    * \param[in] help_message The help message to display with the exception.
    */
public:
  help(std::string const & help_message);

  /** Gets the help message.
   * \returns Help message.
   */
  virtual const char * what() const throw();
};

/** Exception indicating that the user passed an argument that isn't available.
 * E.g. user parser a "-t" option but no "-t" argument is defiend by the program.
 */
class invalid_option : public std::exception
{
private:
  std::string error_message;       /**< Message to display when exception is thrown.*/

public:
  /** Constructs an invalid option exception.
   * \param[in] invalid_option the parsed invalid option.
   * \param[in] help_message the help message to display with the exception.
   */
  invalid_option(std::string const & invalid_option,
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
class invalid_option_value : public std::exception
{
private:
  std::string error_message;       /**< Message to display when exception is thrown. */

public:
  /** Constructs an invalid option value exception.
   * \param[in] shrt Short option.
   * \param[in] lng Long option.
   * \param[in] invalid_option_value the parsed invalid option value.
   * \param[in] help_message the help message to display with the exception.
   */
  invalid_option_value(char const shrt,
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
class invalid_positional : public std::exception
{
private:
  std::string error_message;       /**< Message to display when exception is thrown. */

public:
  /** Constructs an invalid positional exception.
   * \param[in] invalid_positional the parsed invalid positional value.
   * \param[in] help_message the help message to display with the exception.
   */
  invalid_positional(std::string const & invalid_positional,
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
class invalid_subcommand : public std::exception
{
  std::string error_message;     /**< Message to display when exception is thrown. */

public:
  /** Constructs an invalid subcommand exception
   * \param[in] subcommand Subcommand the user passed.
   * \param[in] help_message the help message to display with the exception.
   */
  invalid_subcommand(std::string const & subcommand,
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
class missing_value : public std::exception
{
  std::string error_message;     /**< Message to display when exception is thrown. */

public:
  /** Constructs a missing value exception.
   * \param[in] shrt short option.
   * \param[in] lng long option.
   * \param[in] help_message the help message to display with the exception.
   */
  missing_value(char const shrt, std::string const & lng,
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
class missing_positional_argument : public std::exception
{
private:
  std::string error_message;       /**< Message to display when exception is thrown. */

public:
  /** Constructs a missing positional argument exception
   * \param[in] meta_string Meta string of the positional argument.
   * \param[in] help_message the help message to display with the exception.
   */
  missing_positional_argument(std::string const & meta_string,
                              std::string const & help_message
                              );

  /** Gets the error message of this exception.
   * \returns Error message.
   */
  virtual const char * what() const throw();

};

} // namespace exception
} // namespace paw


#if defined(IMPLEMENT_PAW) || defined(__JETBRAINS_IDE__)


namespace paw
{
namespace exception
{
/* EXCEPTIONS */
help::help(std::string const & _help_message)
  : help_message(_help_message)
{}


const char *
help::what() const throw()
{
  return help_message.c_str();
}


/* Invalid option exception */
invalid_option::invalid_option(std::string const & invalid_option,
                               std::string const & help_message
                               )
  : error_message("")
{
  std::ostringstream ss;
  ss << "[paw::Parser::InvalidOption] ERROR: Unknown option '"
     << invalid_option
     << "' was passed.\n"
     << help_message;
  error_message = ss.str();
}


const char *
invalid_option::what() const throw()
{
  return error_message.c_str();
}


/* Invalid option value exception */
invalid_option_value::invalid_option_value(
  char const shrt,
  std::string const & lng,
  std::string const & invalid_option_value,
  std::string const & help_message
  )
{
  std::ostringstream ss;

  if (shrt != paw::internal::NO_SHORT_OPTION)
  {
    ss << "[paw::Parser::InvalidOptionValue] Value '"
       << invalid_option_value
       << "' for option '"
       << lng
       << "' (or '-"
       << shrt
       << "') is invalid.";
  }
  else
  {
    ss << "[paw::Parser::InvalidOptionValue] Value '"
       << invalid_option_value
       << "' for option '"
       << lng
       << "' is invalid.";
  }

  ss << "\n" << help_message;
  error_message = ss.str();
}


const char *
invalid_option_value::what() const throw()
{
  return error_message.c_str();
}


/* Invalid positional exception */
invalid_positional::invalid_positional(
  std::string const & invalid_positional,
  std::string const & help_message
  )
{
  std::ostringstream ss;
  ss << "[paw::Parser::InvalidPositional] Positional argument <"
     << invalid_positional
     << "> is of invalid type.\n"
     << help_message;

  error_message = ss.str();
}


const char *
invalid_positional::what() const throw()
{
  return error_message.c_str();
}


/* Invalid subcommand exception */
invalid_subcommand::invalid_subcommand(
  std::string const & subcommand,
  std::string const & help_message
  )
{
  std::ostringstream ss;
  ss << "[paw::Parser::InvalidSubcommand] Subcommand '"
     << subcommand
     << "' is invalid.\n"
     << help_message;

  error_message = ss.str();
}


const char *
invalid_subcommand::what() const throw()
{
  return error_message.c_str();
}


/* Missing value exception */
missing_value::missing_value(char const shrt,
                             std::string const & lng,
                             std::string const & help_message
                             )
{
  std::ostringstream ss;

  if (shrt != paw::internal::NO_SHORT_OPTION)
  {
    ss << "[paw::Parser::MissingValue] Option '--"
       << lng
       << "' (or '-"
       << shrt
       << "') requires a value, but was not passed any.";
  }
  else
  {
    ss << "[paw::Parser::MissingValue] Flag '--"
       << lng
       << "' was not passed any value.";
  }

  ss << "\n" << help_message;
  error_message = ss.str();
}


const char *
missing_value::what() const throw()
{
  return error_message.c_str();
}


/* Missing positional argument exception */
missing_positional_argument::missing_positional_argument(
  std::string const & meta_string,
  std::string const & help_message
  )
{
  std::ostringstream ss;
  ss << "[paw::Parser::MissingPositionalArgument] ERROR: "
     << "Positional argument <"
     << meta_string
     << "> is missing.\n"
     << help_message;

  error_message = ss.str();
}


const char *
missing_positional_argument::what() const throw()
{
  return error_message.c_str();
}


} // namespace exception
} // namespace paw

#endif //IMPLEMENT_PAW
