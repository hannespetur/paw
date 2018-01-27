#include <iostream> // std::cout, std::cerr
#include <string> // std::string

#include <paw/parser.hpp> // paw::parser


struct Options
{
  bool my_bool = false;
  unsigned my_uint = 0;
  std::string my_pos_arg;
};


int
main(int argc, char ** argv)
{
  Options options;
  paw::Parser parser(argc, argv);
  parser.set_name("Example 1 - A very basic program that uses Paw parser.");
  parser.set_version("1.0");

  try
  {
    parser.parse_positional_argument(options.my_pos_arg, "positional", "Description of pos arg.");
    parser.parse_option(options.my_bool, 'b', "bool", "Description of the option.");
    parser.parse_option(options.my_uint, 'u', "uint", "Description of the option.", "N");
    parser.finalize();
  }
  catch (const std::exception & e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "my_bool = " << options.my_bool << "\n";
  std::cout << "my_uint = " << options.my_uint << "\n";
  std::cout << "my_pos_arg = " << options.my_pos_arg << "\n";
  return EXIT_SUCCESS;
}
