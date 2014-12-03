#include <exception>

#include "r_interface/check_interrupt.h"
#include "r_interface/error.h"

#include "Models/Glm/PoissonRegressionModel.hpp"
#include "Models/Glm/PosteriorSamplers/PoissonRegressionSpikeSlabSampler.hpp"

#include "Models/ChisqModel.hpp"
#include "Models/IndependentMvnModelGivenScalarSigma.hpp"

#include "r_interface/boom_r_tools.hpp"
#include "r_interface/handle_exception.hpp"
#include "r_interface/list_io.hpp"
#include "r_interface/print_R_timestamp.hpp"
#include "r_interface/prior_specification.hpp"
#include "r_interface/seed_rng_from_R.hpp"

#include "cpputil/Ptr.hpp"

namespace {
  using namespace BOOM;  // NOLINT
  using namespace BOOM::RInterface;  // NOLINT
  Ptr<PoissonRegressionModel> SpecifyPoissonRegression(
      SEXP r_design_matrix,
      SEXP r_integer_response_vector,
      SEXP r_exposure_vector,
      SEXP r_spike_slab_prior,
      SEXP r_nthreads,
      RListIoManager *io_manager) {
    Matrix design_matrix(ToBoomMatrix(r_design_matrix));
    std::vector<int> response(ToIntVector(r_integer_response_vector));
    Vector exposure(ToBoomVector(r_exposure_vector));

    NEW(PoissonRegressionModel, model)(design_matrix.ncol());
    int n = response.size();
    for (int i = 0; i < n; ++i) {
      NEW(PoissonRegressionData, data_point)(
          response[i],
          design_matrix.row(i),
          exposure[i]);
      model->add_data(data_point);
    }
    SpikeSlabGlmPrior prior_spec(r_spike_slab_prior);
    NEW(PoissonRegressionSpikeSlabSampler, sampler)(
        model.get(),
        prior_spec.slab(),
        prior_spec.spike(),
        Rf_asInteger(r_nthreads));
    if (prior_spec.max_flips() > 0) {
      sampler->limit_model_selection(prior_spec.max_flips());
    }
    model->set_method(sampler);

    io_manager->add_list_element(new GlmCoefsListElement(
        model->coef_prm(),
        "beta"));
    return model;
  }
}  // namespace

extern "C" {
  using BOOM::RErrorReporter;
  using BOOM::RCheckInterrupt;
  using namespace BOOM;  // NOLINT

  SEXP poisson_regression_spike_slab(
      SEXP r_design_matrix,
      SEXP r_integer_response_vector,
      SEXP r_exposure_vector,
      SEXP r_spike_slab_prior,
      SEXP r_niter,
      SEXP r_ping,
      SEXP r_nthreads,
      SEXP r_initial_beta,
      SEXP r_seed) {
    RErrorReporter error_reporter;
    try {
      seed_rng_from_R(r_seed);
      RListIoManager io_manager;
      Ptr<PoissonRegressionModel> model = SpecifyPoissonRegression(
          r_design_matrix,
          r_integer_response_vector,
          r_exposure_vector,
          r_spike_slab_prior,
          r_nthreads,
          &io_manager);
      int niter = Rf_asInteger(r_niter);
      int ping = Rf_asInteger(r_ping);
      SEXP ans;
      PROTECT(ans = io_manager.prepare_to_write(niter));
      for (int i = 0; i < niter; ++i) {
        print_R_timestamp(i, ping);
        model->sample_posterior();
        io_manager.write();
      }
      UNPROTECT(1);
      return ans;
    } catch (std::exception &e) {
      handle_exception(e);
    } catch (...) {
      handle_unknown_exception();
    }
    return R_NilValue;
  }

}  // extern "C"
