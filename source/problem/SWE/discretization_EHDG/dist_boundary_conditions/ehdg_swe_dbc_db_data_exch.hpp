#ifndef EHDG_SWE_DBC_DB_DATA_EXCH_HPP
#define EHDG_SWE_DBC_DB_DATA_EXCH_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"

namespace SWE {
namespace EHDG {
namespace DBC {
struct DBIndex {
    uint q_in;
    uint Fn_in;

    uint q_ex;
    uint Fn_ex;
};

class DBDataExchanger {
  private:
    DBIndex index;

    std::vector<double>& send_preproc_buffer;
    std::vector<double>& receive_preproc_buffer;

    std::vector<double>& send_buffer;
    std::vector<double>& receive_buffer;

    std::vector<double>& send_postproc_buffer;
    std::vector<double>& receive_postproc_buffer;

  public:
    DBDataExchanger() = default;
    DBDataExchanger(const DBIndex& index,
                    std::vector<double>& send_preproc_buffer,
                    std::vector<double>& receive_preproc_buffer,
                    std::vector<double>& send_buffer,
                    std::vector<double>& receive_buffer,
                    std::vector<double>& send_postproc_buffer,
                    std::vector<double>& receive_postproc_buffer);

    void SetEX(const std::vector<Vector<double, SWE::n_variables>>& q_in,
               const std::vector<Vector<double, SWE::n_variables>>& Fn_in);

    void GetEX(const uint gp, Vector<double, SWE::n_variables>& q_ex, Vector<double, SWE::n_variables>& Fn_ex);
};

DBDataExchanger::DBDataExchanger(const DBIndex& index,
                                 std::vector<double>& send_preproc_buffer,
                                 std::vector<double>& receive_preproc_buffer,
                                 std::vector<double>& send_buffer,
                                 std::vector<double>& receive_buffer,
                                 std::vector<double>& send_postproc_buffer,
                                 std::vector<double>& receive_postproc_buffer)
    : index(index),
      send_preproc_buffer(send_preproc_buffer),
      receive_preproc_buffer(receive_preproc_buffer),
      send_buffer(send_buffer),
      receive_buffer(receive_buffer),
      send_postproc_buffer(send_postproc_buffer),
      receive_postproc_buffer(receive_postproc_buffer) {}

void DBDataExchanger::SetEX(const std::vector<Vector<double, SWE::n_variables>>& q_in,
                            const std::vector<Vector<double, SWE::n_variables>>& Fn_in) {
    for (uint gp = 0; gp < q_in.size(); gp++) {
        for (uint var = 0; var < SWE::n_variables; var++) {
            this->send_buffer[this->index.q_in + SWE::n_variables * gp + var]  = q_in[gp][var];
            this->send_buffer[this->index.Fn_in + SWE::n_variables * gp + var] = Fn_in[gp][var];
        }
    }
}

void DBDataExchanger::GetEX(const uint gp,
                            Vector<double, SWE::n_variables>& q_ex,
                            Vector<double, SWE::n_variables>& Fn_ex) {
    for (uint var = 0; var < SWE::n_variables; var++) {
        q_ex[SWE::n_variables - var - 1]  = this->receive_buffer[this->index.q_ex - SWE::n_variables * gp - var];
        Fn_ex[SWE::n_variables - var - 1] = this->receive_buffer[this->index.Fn_ex - SWE::n_variables * gp - var];
    }
}
}
}
}

#endif
