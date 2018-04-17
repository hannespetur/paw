#include <paw/align.hpp>


int
main(int, char **)
{
  std::string current_archs = paw::get_current_arch();
  std::cerr << "Current archs are: " << current_archs << "\n";

  paw::align("AC", "GT");
}
