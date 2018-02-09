#include <exception> // std::exception
#include <iostream> // std::cerr
#include <string> // std::string
#include <vector> // std::vector

#include <boost/filesystem.hpp>

#include <paw/parser.hpp>  // paw::parser


namespace fs = boost::filesystem; // You could also use the STL filesystem (c++17 feature)


struct Options
{
  fs::path input;
  fs::path output;
};


int
main(int argc, char ** argv)
{
  Options options;
  paw::Parser parser(argc, argv);

  try
  {
    parser.set_name("Example 4 - Program using Boost filesystem.");
    parser.set_version("1.0");
    parser.parse_positional_argument(options.input, "input", "Input filename.");
    parser.parse_positional_argument(options.output, "output", "Output filename.");
    parser.finalize();
  }
  catch (const std::exception & e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  if (!fs::exists(options.input))
    std::cerr << "!!FAIL!! File '" << options.input << "' does not exist.\n";
  else
    std::cerr << "OK. File '" << options.input << "' was found.\n";

  fs::path const & out_dir = options.output.parent_path();
  std::cout << "Output directory is " << out_dir << "\n";

  if (!out_dir.empty() && !fs::exists(out_dir))
  {
    std::cerr << "Creating output directory.\n";

    if (!fs::create_directory(out_dir))
    {
      std::cerr << "Failed to create output directory.\n";
    }
  }


  return EXIT_SUCCESS;
}
