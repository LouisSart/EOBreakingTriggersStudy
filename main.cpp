#include "DRM-ARM-stats/cubie_cube.hpp"
#include <cassert>
#include <deque>
#include <iomanip>
#include <iostream>

/*
R F :
  - in Eslice : {1E, 2NE, 1fNE}
  - out Eslice : {2fE, 1E, 1fNE}
  - 3+ 3-
  5 644 800 RZP states
R U2 F :
    - in Eslice : {1E, 2NE, 1fNE}
    - out Eslice : {1E, 2fE, 1fNE}
    - 6+
    auto roots = init_roots({{F3, U2, R3},
                        {F3, U2, R},
                        {B3, U2, L},
                        {B3, U2, L3},
                        {F3, D2, L3},
                        {F3, D2, L},
                        {B3, D2, R},
                        {B3, D2, R3}});
R U F :
    - in Eslice : {1E, 1fE, 1NE, 1fNE}
    - out Eslice : {2fE}
    - 3+
    451584 cases
    auto roots = init_roots({{F3, U3, R3},
                        {F3, U3, R},
                        {B3, U3, L},
                        {B3, U3, L3},
                        {F3, D3, L3},
                        {F3, D3, L},
                        {B3, D3, R},
                        {B3, D3, R3}});
F R F :
    - in Eslice : {3E, 1NE}
    - out Eslice : {1E, 2fNE}
    - 3+
    37632 cases
R U' F:
  auto roots = init_roots({{F3, U, R3},
                           {F3, U, R},
                           {B3, U, L},
                           {B3, U, L3},
                           {F3, D, L3},
                           {F3, D, L},
                           {B3, D, R},
                           {B3, D, R3}});
*/

unsigned N_EO = ipow(2, NE - 1);
unsigned N_CO = ipow(3, NC - 1);
unsigned N_ESL = binomial(NE, 4);

bool is_E_edge(const Cubie &e) {
  return e == RF || e == RB || e == LF || e == LB;
}

unsigned dr_coord(CubieCube &cc) {
  static unsigned layout[NE] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (unsigned i = 0; i < NE; i++) {
    layout[i] = (int)is_E_edge(cc.ep[i]);
  }

  auto elc = layout_coord(layout, NE);
  auto eoc = eo_coord(cc.eo, NE - 1);
  auto coc = co_coord(cc.co, NC - 1);
  return (elc * N_EO + eoc) * N_CO + coc;
}

class DRNode {
public:
  CubieCube state;
  unsigned depth;

  DRNode(CubieCube state = CubieCube(), unsigned depth = 0)
      : state(state), depth(depth) {}

  DRNode make_child(Move m) {
    CubieCube new_state = state;
    new_state.apply(elementary_transformations[m]);
    return DRNode(new_state, depth + 1);
  }
};

auto generator(const std::deque<DRNode> &roots) {

  std::set<unsigned> visited;
  std::vector<unsigned> visited_counts;
  std::deque<DRNode> queue = roots;
  unsigned current_depth = 0;
  std::cout << "Searching at depth: 0" << std::endl;

  while (queue.size() > 0) {
    auto node = queue.back();
    unsigned coord = dr_coord(node.state);

    if (node.depth > current_depth) {
      std::cout << "Visited: " << visited.size() << std::endl;
      visited_counts.push_back(visited.size());
      current_depth = node.depth;
      std::cout << "Searching at depth: " << node.depth << std::endl;
    }

    if (!visited.contains(coord)) {
      visited.insert(coord);
      for (Move m : {U, U3, U2, D, D2, D3, R2, L2, F2, B2}) {
        auto child = node.make_child(m);
        queue.push_front(child);
      }
    }
    queue.pop_back();
  }
  std::cout << "Visited: " << visited.size() << std::endl;
  visited_counts.push_back(visited.size());

  return visited_counts;
}

void display_distribution(const std::vector<unsigned> &visited_counts) {
  std::cout << std::setprecision(3);
  std::cout << "[length][cases][% of <= length]" << std::endl;
  auto N = visited_counts.back();
  for (unsigned i = 0; i < visited_counts.size(); i++) {
    std::cout << "Depth " << i << ": " << visited_counts[i] << " "
              << ((double)visited_counts[i] / N) * 100 << "%" << std::endl;
  }
}

auto init_roots(std::vector<Algorithm> triggers) {
  std::deque<DRNode> roots;
  for (auto trigger : triggers) {
    auto cc = CubieCube();
    cc.apply(trigger);
    auto node = DRNode(cc);
    roots.push_back(node);
  }
  return roots;
}

int main() {
  auto roots = init_roots({{F3, U, R3},
                           {F3, U, R},
                           {B3, U, L},
                           {B3, U, L3},
                           {F3, D, L3},
                           {F3, D, L},
                           {B3, D, R},
                           {B3, D, R3}});
  auto distribution = generator(roots);
  display_distribution(distribution);

  return 0;
}