#pragma once

#include <iterator> // std::distance
#include <sstream> // std::ostringstream
#include <string> // std::string


namespace paw
{
/** Namespace for internal variables and function used by the paw library. */
namespace internal
{

/* READ ONLY VARIABLES */
/** Specifies how to represent arguments that have no short option.*/
char static const NO_SHORT_OPTION = ' ';

/** Specifies how the default seperator that splits argument lists.*/
char static const DEFAULT_LIST_DELIMITER = ',';

/** Specifies how to indicate the option has no value. */
std::string const OPTION_HAS_NO_VALUE = "__OPTIONS_HAS_NO_VALUE__";


/** Prints a string to stream with a fixed maximum width. */
void inline
print_string(std::ostringstream & ss,
             std::string const & str,
             std::size_t const base_indent,
             std::size_t const MAX_WIDTH
             );

} // namespace internal
} // namespace paw


#ifdef IMPLEMENT_PAW

/******************
 * IMPLEMENTATION *
 ******************/

#include <algorithm> // std::find
#include <iomanip> // std::setw


namespace paw
{
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


} // namespace internal
} // namespace paw

#endif // IMPLEMENT_PAW
