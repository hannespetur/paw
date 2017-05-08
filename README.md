## Paw
The paw library is a header-only C++11 library.
Its main design goal is to be very simple and convenient to use.

### Parser
Paw parser is a small header-only file for parsing command-line arguments defined in `parser.hpp`.
Its most distinctive feature that it automatically parses argument values to their respective types (see example below).

#### Example
This example is taken from `examples/ex1_basic_usage.cpp`.
It shows how to create a dummy program which has two options, one flag (option without a value) and another that has an unsigned integer.
```cpp
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
  paw::parser parser(argc, argv);
  parser.set_name("Example 1 - A very basic program that uses Paw parser.");
  parser.set_version("1.0");

  try
  {
    parser.parse_option(options.my_bool, 'b', "bool", "Description of the option.");
    parser.parse_option(options.my_uint, 'u', "uint", "Description of the option.", "N");
    parser.parse_positional_argument(options.my_pos_arg, "positional", "Description of pos arg.");
    parser.finalize();
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "my_bool = " << options.my_bool << "\n";
  std::cout << "my_uint = " << options.my_uint << "\n";
  std::cout << "my_pos_arg = " << options.my_pos_arg << "\n";
  return EXIT_SUCCESS;
}
```

### More to come...
Hopefully.

### Author
Hannes P Eggertsson

### License
GNU GPLv3
