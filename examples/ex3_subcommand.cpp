#include <iostream> // std::cout, std::cerr
#include <string> // std::string
#include <vector> // std::vector
#include <paw/parser.hpp>  // paw::parser


struct Options
{
  std::string subcmd;
  bool my_bool = false;
  int my_int = 0;
  unsigned my_uint = 0;
  double my_double = 0.0;
  std::string my_string = "";
  std::vector<int> my_ints;
  std::vector<std::string> my_strings;
  std::string my_first_pos_argument;
  std::string my_second_pos_argument;
  std::vector<std::string> my_remaining_pos_arguments;
};

int
main(int argc, char ** argv)
{
  Options options;

  try
  {
    paw::parser parser(argc, argv);
    parser.set_name("Example 3 - Program with subcommands.");
    parser.set_version("1.0");
    parser.add_subcommand("add", "Subcommand that adds");
    parser.add_subcommand("divide", "Subcommand that divides");
    parser.parse_subcommand(options.subcmd);
    parser.finalize();
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
