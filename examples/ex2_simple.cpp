#include <iostream> // std::cout, std::cerr
#include <string> // std::string
#include <vector> // std::vector
#include <paw/parser.hpp>  // paw::parser


struct Options
{
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
  paw::parser parser(argc, argv);
  parser.set_name("Example 2 - A simple program that uses Paw parser.");
  parser.set_version("3.14.15");

  // Parse arguments
  try
  {
    parser.parse_option(options.my_bool,
                        'b',
                        "bool",
                        "Test boolean value. Lorem ipsum dolor sit amet, consectetur adipiscing "
                        "elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. "
                        "Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi "
                        "ut aliquip ex ea commodo consequat. Duis aute irure dolor in "
                        "reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla "
                        "pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa "
                        "qui officia deserunt mollit anim id est laborum.");
    parser.parse_option(options.my_int, 'i', "int", "Test int value.");
    parser.parse_option(options.my_uint, 'u', "uint", "Test unsigned int value.");
    parser.parse_option(options.my_double, 'd', "double", "Test double value.");
    parser.parse_option(options.my_string, 's', "string", "Test string value.");
    parser.parse_option_list(options.my_strings, 'S', "strings", "Test a list of strings.");
    parser.parse_option_list(options.my_ints, 'I', "ints", "Test a list of ints.");
    parser.parse_positional_argument(options.my_first_pos_argument,
                                     "first",
                                     "First positional argument."
                                     );
    parser.parse_positional_argument(options.my_second_pos_argument,
                                     "second",
                                     "Second positional argument."
                                     );
    parser.parse_remaining_positional_arguments(options.my_remaining_pos_arguments,
                                                "remaining args...",
                                                "Other remaining stuff."
                                                );
    parser.finalize();
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "my_bool = " << options.my_bool << "\n";
  std::cout << "my_int = " << options.my_int << "\n";
  std::cout << "my_uint = " << options.my_uint << "\n";
  std::cout << "my_double = " << options.my_double << "\n";
  std::cout << "my_string = \"" << options.my_string << "\"\n";
  std::cout << "my_ints = [ ";


  for (const auto& my_int : options.my_ints)
    std::cout << my_int << " ";

  std::cout << "]\n";
  std::cout << "my_strings = [ ";

  for (const auto& my_string : options.my_strings)
    std::cout << my_string << " ";
  std::cout << "]\n";

  std::cout << "my_first_pos_argument = " << options.my_first_pos_argument << "\n";
  std::cout << "my_second_pos_argument = " << options.my_second_pos_argument << "\n";
  std::cout << "my_remaining_pos_arguments = [ ";

  for (const auto& my_pos_args : options.my_remaining_pos_arguments)
    std::cout << my_pos_args << " ";
  std::cout << "]\n";
  return EXIT_SUCCESS;
}
