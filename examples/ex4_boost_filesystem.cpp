#include <iostream> // std::cerr
#include <string> // std::string
#include <vector> // std::vector

#include <boost/filesystem.hpp>

#include <paw/parser.hpp>  // paw::parser


struct Options
{
  boost::filesystem::path input;
  boost::filesystem::path output;
};


int
main(int argc, char ** argv)
{
  namespace fs = boost::filesystem;

  Options options;
  paw::parser parser(argc, argv);

  try
  {
    parser.set_name("Example 4 - Program using Boost filesystem.");
    parser.set_version("1.0");
    parser.parse_option(options.input, 'i', "input", "Input filename.");
    parser.parse_option(options.output, 'o', "output", "Output filename.");

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

    parser.finalize();
  }
  catch (const std::exception & e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}
