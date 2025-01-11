#include "DRM-ARM-stats/cubie_cube.hpp"
#include <iostream>
#include <cassert>

/* 
R F : 
  - in Eslice : {1E, 2NE, 1fNE}
  - out Eslice : {2fE, 1E, 1fNE}
  - 3+ 3-
R U2 F :
    - in Eslice : {1E, 2NE, 1fNE}
    - out Eslice : {1E, 2fE, 1fNE}
    - 6+
R U F :
    - in Eslice : {1E, 1fE, 1NE, 1fNE}
    - out Eslice : {2fE}
    - 3+
*/

unsigned N_EO = ipow(2, NE);
unsigned N_CO = ipow(3, NC);
unsigned N_ESL = binomial(NE, 4);

bool is_E_edge(const Cubie &e){
    return e == RF || e == RB || e == LF || e == LB;
}

unsigned dr_coord(CubieCube &cc) {
    static unsigned layout[NE] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    for (unsigned i = 0; i < NE; i++) {
        layout[i] = (int)is_E_edge(cc.ep[i]);
    }

    auto elc = layout_coord(layout, NE);
    auto eoc = eo_coord(cc.eo, NE);
    auto coc = co_coord(cc.co, NC);
    return (elc * N_EO + eoc) * N_CO + coc;
}

int main() {
    std::srand(std::time(nullptr));
    std::set<unsigned> visited;

    for (unsigned k = 0; k < 1000; ++k){
        auto cc = CubieCube::random_state();
        auto n = dr_coord(cc);
        visited.insert(n);
    }
    assert(*prev(visited.end()) < 4294967295);
    assert(visited.size() == 1000);

    return 0;
}