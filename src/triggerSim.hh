#ifndef triggerSim_hh
#define triggerSim_hh

//C, C++
#include <vector>
#include <array>

//root
#include <TROOT.h>

class sipmCameraHist;
class TH1D;

class triggerSim {

public :
  //
  triggerSim(const sipmCameraHist* simphist);
  ~triggerSim();

  //std::vector<std::array<int, 2>> get_trigger(const std::vector<std::vector<int>> &wf);
  std::vector<std::vector<int>> get_trigger_test(const std::vector<std::vector<int>> &wf);
  std::vector<std::vector<int>> get_trigger_test();
  std::vector<std::vector<unsigned int>> get_trigger(const std::vector<std::vector<int>> &wf,
						     TH1D *h1_digital_sum,
						     TH1D *h1_digital_sum_3ns,
						     TH1D *h1_digital_sum_5ns,
						     TH1D *h1_fadc_val);
  static void print_trigger_vec(const std::vector<std::array<int, 2>> &trg_vector);
  static void print_trigger_vec(const std::vector<std::vector<unsigned int>> &trg_vector);

private:

  int get_flower_digital_sum(const unsigned int ch_i, const unsigned int wf_j, const std::vector<std::vector<int>> &wf, Int_t w_l, Int_t w_r, Bool_t norm_yes);
  int get_digital_sum(const unsigned int ch_i, const unsigned int wf_j, const std::vector<std::vector<int>> &wf, Int_t w_l, Int_t w_r, Bool_t norm_yes, Int_t sum_type);

  const sipmCameraHist* _simphist;

};

#endif
