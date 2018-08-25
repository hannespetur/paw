#include <paw/align.hpp>
#include <paw/parser.hpp>
#include <paw/internal/config.hpp>

#include <unordered_map>


/// A very simple graph structure (only support for adding vertices/edges but not removing them)
class Graph
{
  std::vector<std::string> vertices; // ID of vertex is its index in this vector
  std::unordered_map<long, std::vector<long> > from_to_edges; // from->to edge, where from is key
  std::unordered_map<long, std::vector<long> > to_from_edges; // from->to edge, where to is key

public:
  inline void
  add_vertex(std::string && seq)
  {
    vertices.push_back(std::move(seq));
  }


  inline void
  add_edge(long const from, long const to)
  {
    from_to_edges[from].push_back(to);
    to_from_edges[to].push_back(from);
  }


  inline std::string
  get_sequence(long index) const
  {
    return vertices.at(index);
  }


  inline std::vector<long>
  get_indexes_inbound(long index) const
  {
    if (to_from_edges.count(index) == 0)
      return std::vector<long>();

    return to_from_edges.at(index);
  }


  inline std::vector<long>
  get_indexes_outbound(long index) const
  {
    if (from_to_edges.count(index) == 0)
      return std::vector<long>();

    return from_to_edges.at(index);
  }
};


std::vector<Graph>
create_test_graphs()
{
  std::vector<Graph> graphs;

  // Test graph 1
  {
    // AAAAAAAAAA
    // A|T
    // AAAAAAAAA
    Graph g;
    g.add_vertex("A");
    g.add_vertex("A");
    g.add_vertex("T");
    g.add_vertex("AAAAAAAAA");
    g.add_edge(0, 1);
    g.add_edge(0, 2);
    g.add_edge(1, 3);
    g.add_edge(2, 3);

    assert(g.get_indexes_inbound(0).size() == 0);
    assert(g.get_indexes_outbound(0).size() == 2);
    assert(g.get_indexes_inbound(1).size() == 1);
    assert(g.get_indexes_outbound(1).size() == 1);
    assert(g.get_indexes_inbound(2).size() == 1);
    assert(g.get_indexes_outbound(2).size() == 1);
    assert(g.get_indexes_inbound(3).size() == 2);
    assert(g.get_indexes_outbound(3).size() == 0);
    assert(g.get_sequence(0) == "A");
    assert(g.get_sequence(1) == "A");
    assert(g.get_sequence(2) == "T");
    assert(g.get_sequence(3) == "AAAAAAAAA");
    graphs.push_back(std::move(g));
  }

  return graphs;
}


int
main(int argc, char ** argv)
{
  std::vector<int> tests_to_run;

  try
  {
    paw::Parser parser(argc, argv);
    parser.set_name("Alignment example 5 - Graph tests.");
    parser.parse_remaining_positional_arguments(tests_to_run,
                                                "list of tests to run...",
                                                "List of all tests to run (all if not specified)."
      );
    parser.finalize();
  }
  catch (std::exception const & e)
  {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }

  std::vector<Graph> graphs = create_test_graphs();

  // Graph 1
  auto const & g = graphs[0];
  paw::AlignmentOptions<uint8_t> opts;
  opts.continuous_alignment = true;
  //opts.left_column_free = true;
  opts.right_column_free = true;
  opts.set_match(1).set_mismatch(-2).set_gap_open(-5).set_gap_extend(-1);
  //std::string query = ;
  auto ar = paw::global_alignment(std::string("AA"), g.get_sequence(0), opts);
  auto ar2 = paw::global_alignment(std::string("AA"), g.get_sequence(0), opts);
  std::cout << "Score = " << ar.score << "\n";
  std::cout << "Score 2 = " << ar2.score << "\n";

  std::cout << "Num test graphs = " << graphs.size() << "\n";
  return EXIT_SUCCESS;
}