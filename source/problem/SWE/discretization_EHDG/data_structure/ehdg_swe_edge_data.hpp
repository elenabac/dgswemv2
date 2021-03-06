#ifndef EHDG_SWE_EDGE_DATA_HPP
#define EHDG_SWE_EDGE_DATA_HPP

#include "ehdg_swe_edge_data_state.hpp"
#include "ehdg_swe_edge_data_internal.hpp"

namespace SWE {
namespace EHDG {
struct EdgeData {
    EdgeState edge_state;
    EdgeInternal edge_internal;

    void initialize() {
        this->edge_state    = EdgeState(this->ndof);
        this->edge_internal = EdgeInternal(this->ngp);
    }

    uint get_ndof() { return this->ndof; }
    uint get_ngp() { return this->ngp; }

    void set_ndof(const uint ndof) { this->ndof = ndof; }
    void set_ngp(const uint ngp) { this->ngp = ngp; }

  private:
    uint ndof;
    uint ngp;
};
}
}

#endif
