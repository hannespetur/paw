#include <algorithm> // std::sort
#include <iostream> // std::cout, std::cerr
#include <iterator> // std::ostream_iterator
#include <numeric> // std::accumulate
#include <string> // std::string
#include <vector> // std::vector

#include <paw/parser.hpp>  // paw::parser


struct Options
{
  std::string subcmd;
  bool use_double = false;
  std::vector<long> my_ints;
  std::vector<double> my_doubles;
};


template <typename T>
void
print_sum(T first, T last)
{
  using TVal = typename T::value_type;
  TVal sum = std::accumulate(first, last, static_cast<TVal>(0));
  std::cout << sum << "\n";
}


template <typename T>
void
print_sorted_vector(T first, T last)
{
  using TVal = typename T::value_type;
  std::sort(first, last);
  std::copy(first, last, std::ostream_iterator<TVal>(std::cout, " "));
  std::cout << "\n";
}


int
main(int argc, char ** argv)
{
  Options options;
  paw::Parser parser(argc, argv);

  try
  {
    parser.set_name("Example 3 - Program with subcommands.");
    parser.set_version("1.0");
    parser.add_subcommand("sum", "Subcommand that finds the sum of stuff.");
    parser.add_subcommand("sort", "Subcommand that sorts stuff.");
    parser.parse_subcommand(options.subcmd);
    parser.parse_option(options.use_double, 'd', "double", "Parse numbers as doubles.");

    // Parse as ints if no '--double' option was passed
    if (options.use_double)
    {
      parser.parse_remaining_positional_arguments(options.my_doubles,
                                                  "args...",
                                                  "Other remaining stuff."
                                                  );
    }
    else
    {
      parser.parse_remaining_positional_arguments(options.my_ints,
                                                  "args...",
                                                  "Other remaining stuff."
                                                  );
    }

    parser.finalize();
  }
  catch (const std::exception & e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  if (options.subcmd == "sum")
  {
    if (options.use_double)
      print_sum(options.my_doubles.begin(), options.my_doubles.end());
    else
      print_sum(options.my_ints.begin(), options.my_ints.end());
  }
  else if (options.subcmd == "sort")
  {
    if (options.use_double)
      print_sorted_vector(options.my_doubles.begin(), options.my_doubles.end());
    else
      print_sorted_vector(options.my_ints.begin(), options.my_ints.end());
  }

  return EXIT_SUCCESS;
}
