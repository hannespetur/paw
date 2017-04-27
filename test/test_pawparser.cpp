#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "pawparser.hpp"

TEST_CASE("Parse argc argv")
{
  std::vector<std::string> args = {"program"};

  SECTION("Short option used")
  {
    args.push_back("-t");
    paw::Parser pawparser(args);
    paw::TParserFlags const & flags = pawparser.get_flags_reference();

    REQUIRE(flags.count("t") == 1); // Make sure option was registered
    REQUIRE(flags.size() == 1); // Make sure it is the only item registered
  }
}
