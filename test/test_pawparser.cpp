#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "pawparser.hpp"


TEST_CASE("No options parsed")
{
  paw::parser pawparser({"program"});
  const paw::parser::FlagMap& flag_map = pawparser.get_flag_map_reference();
  REQUIRE(flag_map.size() == 0);
}


TEST_CASE("Parse argc argv with short options")
{
  SECTION("Short option used once")
  {
    paw::parser pawparser({"program", "-t"});
    const paw::parser::FlagMap& flag_map = pawparser.get_flag_map_reference();

    REQUIRE(flag_map.count("t") == 1);
    REQUIRE(flag_map.size() == 1);
  }

  SECTION("Short option used twice")
  {
    paw::parser pawparser({"program", "-t", "-t"});
    const paw::parser::FlagMap& flag_map = pawparser.get_flag_map_reference();
    REQUIRE(flag_map.count("t") == 2);
    REQUIRE(flag_map.size() == 2);
  }

  SECTION("Short option with value")
  {
    paw::parser pawparser({"program", "-t4"});
    const paw::parser::FlagMap& flag_map = pawparser.get_flag_map_reference();
    REQUIRE(flag_map.count("t") == 1);
    REQUIRE(flag_map.size() == 1);
    auto range_it = flag_map.equal_range("t");
    REQUIRE(std::distance(range_it.first, range_it.second) == 1);
    REQUIRE(range_it.first->first == std::string("t"));
    REQUIRE(range_it.first->second == std::string("4"));
    REQUIRE(std::next(range_it.first) == range_it.second);
  }

  SECTION("Multiple short options with values")
  {
    paw::parser pawparser({"program", "-t4", "-t84", "-t42"});
    const paw::parser::FlagMap& flag_map = pawparser.get_flag_map_reference();
    REQUIRE(flag_map.count("t") == 3);
    REQUIRE(flag_map.size() == 3);
    auto range_it = flag_map.equal_range("t");
    REQUIRE(std::distance(range_it.first, range_it.second) == 3);
    REQUIRE(range_it.first->first == std::string("t"));
    REQUIRE(range_it.first->second == std::string("4"));
    ++range_it.first;
    REQUIRE(range_it.first->first == std::string("t"));
    REQUIRE(range_it.first->second == std::string("84"));
    ++range_it.first;
    REQUIRE(range_it.first->first == std::string("t"));
    REQUIRE(range_it.first->second == std::string("42"));
    REQUIRE(std::next(range_it.first) == range_it.second);
  }
}


TEST_CASE("Parse argc argv with long options")
{
  SECTION("Long option used once")
  {
    paw::parser pawparser({"program", "--test"});
    const paw::parser::FlagMap& flag_map = pawparser.get_flag_map_reference();
    REQUIRE(flag_map.count("test") == 1);
    REQUIRE(flag_map.size() == 1);
  }

  SECTION("Long option used multiple times with value")
  {
    paw::parser pawparser({"program", "--test=1", "--test=3", "--test=2"});
    const paw::parser::FlagMap& flag_map = pawparser.get_flag_map_reference();
    REQUIRE(flag_map.count("test") == 3);
    REQUIRE(flag_map.size() == 3);
    auto range_it = flag_map.equal_range("test");
    REQUIRE(std::distance(range_it.first, range_it.second) == 3);
    REQUIRE(range_it.first->first == std::string("test"));
    REQUIRE(range_it.first->second == std::string("1"));
    ++range_it.first;
    REQUIRE(range_it.first->first == std::string("test"));
    REQUIRE(range_it.first->second == std::string("3"));
    ++range_it.first;
    REQUIRE(range_it.first->first == std::string("test"));
    REQUIRE(range_it.first->second == std::string("2"));
    REQUIRE(std::next(range_it.first) == range_it.second);
  }
}


TEST_CASE("Parse argc argv with short and long options")
{
  SECTION("Mixture of multiple short and long options")
  {
    paw::parser pawparser({"program", "-t42", "--test=3", "-t1", "--test=2", "--test=4"});
    const paw::parser::FlagMap& flag_map = pawparser.get_flag_map_reference();
    REQUIRE(flag_map.count("test") == 3);
    REQUIRE(flag_map.count("t") == 2);
    REQUIRE(flag_map.size() == 5);

    auto range_it = flag_map.equal_range("test");
    REQUIRE(std::distance(range_it.first, range_it.second) == 3);
    REQUIRE(range_it.first->first == std::string("test"));
    REQUIRE(range_it.first->second == std::string("3"));
    ++range_it.first;
    REQUIRE(range_it.first->first == std::string("test"));
    REQUIRE(range_it.first->second == std::string("2"));
    ++range_it.first;
    REQUIRE(range_it.first->first == std::string("test"));
    REQUIRE(range_it.first->second == std::string("4"));
    REQUIRE(std::next(range_it.first) == range_it.second);

    range_it = flag_map.equal_range("t");
    REQUIRE(std::distance(range_it.first, range_it.second) == 2);
    REQUIRE(range_it.first->first == std::string("t"));
    REQUIRE(range_it.first->second == std::string("42"));
    ++range_it.first;
    REQUIRE(range_it.first->first == std::string("t"));
    REQUIRE(range_it.first->second == std::string("1"));
    REQUIRE(std::next(range_it.first) == range_it.second);
  }
}


TEST_CASE("Parse optional values")
{
  SECTION("Valid arguments")
  {
    paw::parser pawparser({"program", "--int16=-1", "--uint16=-1", "--str=StRiNg", "--bool"});

  }


}
