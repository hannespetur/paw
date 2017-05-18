#include <iostream> // std::cout, std::cerr
#include <numeric> // std::accumulate
#include <string> // std::string
#include <vector> // std::vector

#include <paw/parser.hpp>  // paw::parser


struct Options
{
  std::string subcmd;
  bool use_double = false;
  bool is_right_to_left = false;
  std::vector<long> my_ints;
  std::vector<double> my_doubles;
};


int
main(int argc, char ** argv)
{
  Options options;
  paw::parser parser(argc, argv);

  try
  {
    parser.set_name("Example 3 - Program with subcommands.");
    parser.set_version("1.0");
    parser.parse_option(options.use_double, 'd', "double", "Parse numbers as doubles.");
    parser.add_subcommand("add", "Subcommand that adds stuff.");
    parser.add_subcommand("subtract", "Subcommand that subtracts stuff.");
    parser.add_subcommand("divide", "Subcommand that divides stuff.");
    parser.add_subcommand("sort", "Subcommand that sorts stuff.");
    parser.parse_subcommand(options.subcmd);

    // Options exclusive to the subtract subcommand
    if (options.subcmd == "subtract")
      parser.parse_option(options.is_right_to_left, 'r', "--right-to-left",
                          "Parse right to left.");

    // Parse as ints if no '--double' option was passed
    if (options.use_double)
    {
      parser.parse_remaining_positional_arguments(options.my_doubles,
                                                  "remaining args...",
                                                  "Other remaining stuff."
                                                  );
    }
    else
    {
      parser.parse_remaining_positional_arguments(options.my_ints,
                                                  "remaining args...",
                                                  "Other remaining stuff."
                                                  );
    }

    parser.finalize();
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  if (options.subcmd == "add")
  {
    if (options.use_double)
    {
      std::cout << std::accumulate(options.my_doubles.begin(), options.my_doubles.end(), 0.0)
                << "\n";
    }
    else
    {
      std::cout << std::accumulate(options.my_ints.begin(), options.my_ints.end(), 0l)
                << "\n";
    }
  }

  return EXIT_SUCCESS;
}
