#include <iterator> // std::distance, std::next
#include <string> // std::string
#include <vector> // std::vector<T>

#include "../catch.hpp"

#include <paw/parser.hpp>


TEST_CASE("No options parsed")
{
  paw::Parser parser({"program"});
  paw::Parser::FlagMap const & flag_map = parser.get_flag_map_reference();
  REQUIRE(flag_map.size() == 0);
  parser.set_name("program2");
}


TEST_CASE("Parse argc argv with short options")
{
  SECTION("Short option used once")
  {
    paw::Parser pawparser({"program", "-t"});
    paw::Parser::FlagMap const & flag_map = pawparser.get_flag_map_reference();

    REQUIRE(flag_map.count("t") == 1);
    REQUIRE(flag_map.size() == 1);
  }

  SECTION("Short option used twice")
  {
    paw::Parser pawparser({"program", "-t", "-t"});
    paw::Parser::FlagMap const & flag_map = pawparser.get_flag_map_reference();
    REQUIRE(flag_map.count("t") == 2);
    REQUIRE(flag_map.size() == 2);
  }

  SECTION("Short option with value")
  {
    paw::Parser pawparser({"program", "-t4"});
    paw::Parser::FlagMap const & flag_map = pawparser.get_flag_map_reference();
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
    paw::Parser pawparser({"program", "-t4", "-t84", "-t42"});
    paw::Parser::FlagMap const & flag_map = pawparser.get_flag_map_reference();
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
    paw::Parser pawparser({"program", "--test"});
    paw::Parser::FlagMap const & flag_map = pawparser.get_flag_map_reference();
    REQUIRE(flag_map.count("test") == 1);
    REQUIRE(flag_map.size() == 1);
  }

  SECTION("Long option used multiple times with value")
  {
    paw::Parser pawparser({"program", "--test=1", "--test=3", "--test=2"});
    paw::Parser::FlagMap const & flag_map = pawparser.get_flag_map_reference();
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
    paw::Parser pawparser({"program", "-t42", "--test=3", "-t1", "--test=2", "--test=4"});
    paw::Parser::FlagMap const & flag_map = pawparser.get_flag_map_reference();
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


struct Options
{
  bool my_bool = false;
  int my_int = 0;
  unsigned my_uint = 0;
  double my_double = 0.0;
  std::string my_string = "";
  std::vector<int> my_ints;
  std::vector<std::string> my_strings;
};


TEST_CASE("Parse options")
{
  Options options;

  SECTION("Valid options are passed")
  {
    paw::Parser pawparser(
      {"program", "-d-100.5", "--int=-1", "--uint=-1", "--string=StRiNg", "--bool"}
      );
    pawparser.parse_option(options.my_bool, 'b', "bool", "Test boolean value.");
    pawparser.parse_option(options.my_int, 'i', "int", "Test int value.");
    pawparser.parse_option(options.my_uint, 'u', "uint", "Test unsigned int value.");
    pawparser.parse_option(options.my_double, 'd', "double", "Test double value.");
    pawparser.parse_option(options.my_string, 's', "string", "Test string value.");

    REQUIRE(options.my_bool);
    REQUIRE(options.my_int == -1);
    REQUIRE(options.my_uint > 0);
    REQUIRE(options.my_uint == static_cast<unsigned int>(-1));
    REQUIRE(options.my_double == -100.5);
    REQUIRE(options.my_string == "StRiNg");
  }

  SECTION("Pass invalid numbers")
  {
    paw::Parser pawparser(
      {"program", "--double=0.5a", "--int=-1b", "--uint=1."}
      );

    REQUIRE_THROWS_AS(
      pawparser.parse_option(options.my_double, 'd', "double", "Test double value."),
      paw::exception::invalid_option_value
      );

    REQUIRE_THROWS_AS(
      pawparser.parse_option(options.my_int, 'i', "int", "Test int value."),
      paw::exception::invalid_option_value
      );

    REQUIRE_THROWS_AS(
      pawparser.parse_option(options.my_int, 'u', "uint", "Test uint value."),
      paw::exception::invalid_option_value
      );
  }

  SECTION("Option with missing value")
  {
    paw::Parser pawparser({"program", "-d", "-b"});
    REQUIRE_THROWS_AS(
      pawparser.parse_option(options.my_double, 'd', "double", "Test double value."),
      paw::exception::missing_value
      );

    REQUIRE_NOTHROW(pawparser.parse_option(options.my_bool, 'b', "bool", "Test bool value."));
    REQUIRE(options.my_bool);
  }

  SECTION("Valid option string list")
  {
    paw::Parser pawparser({"program", "--files=f1.txt,f2.txt", "-ff3.txt", "main.txt"});
    pawparser.parse_option_list(options.my_strings, 'f', "files", "Test file list");
    REQUIRE(options.my_strings.size() == 3);
    REQUIRE(options.my_strings[0] == std::string("f1.txt"));
    REQUIRE(options.my_strings[1] == std::string("f2.txt"));
    REQUIRE(options.my_strings[2] == std::string("f3.txt"));
  }

  SECTION("Option list with missing value")
  {
    paw::Parser pawparser({"program", "--files", "f1.txt"});
  }
}
