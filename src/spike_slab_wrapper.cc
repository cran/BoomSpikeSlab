// Copyright 2010 Google Inc. All Rights Reserved.
// Author: stevescott@google.com (Steve Scott)

#include <exception>

#include "r_interface/check_interrupt.h"
#include "r_interface/error.h"

#include "Models/Glm/PosteriorSamplers/BregVsSampler.hpp"
#include "Models/Glm/PosteriorSamplers/SpikeSlabDaRegressionSampler.hpp"
#include "Models/Glm/RegressionModel.hpp"

#include "Models/ChisqModel.hpp"
#include "Models/IndependentMvnModelGivenScalarSigma.hpp"
#include "Models/MvnGivenScalarSigma.hpp"

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

  Ptr<RegressionModel> SpecifyRegressionModel(
      SEXP r_design_matrix,
      SEXP r_response_vector,
      SEXP r_spike_slab_prior,
      SEXP r_bma_method,
      SEXP r_oda_options,
      BOOM::RListIoManager *io_manager) {
    Vector y(ToBoomVector(r_response_vector));
    Matrix X(ToBoomMatrix(r_design_matrix));
    Ptr<RegressionModel> model(new RegressionModel(X.ncol()));
    int n = X.nrow();
    for (int i = 0; i < n; ++i) {
      model->add_data(Ptr<RegressionData>(new RegressionData(y[i], X.row(i))));
    }

    model->coef().drop_all();
    model->coef().add(0);  // start with the intercept

    std::string bma_method = ToString(r_bma_method);
    if (bma_method == "SSVS") {
      BOOM::RInterface::RegressionConjugateSpikeSlabPrior prior(
          r_spike_slab_prior, model->Sigsq_prm());
      NEW(BregVsSampler, ssvs_sampler)(
          model.get(),
          prior.slab(),
          prior.siginv_prior(),
          prior.spike());
      ssvs_sampler->set_sigma_upper_limit(prior.sigma_upper_limit());
      if (prior.max_flips() > 0) {
        ssvs_sampler->limit_model_selection(prior.max_flips());
      }
      model->set_method(ssvs_sampler);
    } else if (bma_method == "ODA") {
      BOOM::RInterface::IndependentRegressionSpikeSlabPrior prior(
          r_spike_slab_prior, model->Sigsq_prm());
      double eigenvalue_fudge_factor = .01;
      double fallback_probability = 0.0;
      if (!Rf_isNull(r_oda_options)) {
        eigenvalue_fudge_factor = Rf_asReal(getListElement(
            r_oda_options,
            "eigenvalue.fudge.factor"));
        fallback_probability = Rf_asReal(getListElement(
            r_oda_options,
            "fallback.probability"));
      }
      NEW(SpikeSlabDaRegressionSampler, oda_sampler)(
              model.get(),
              prior.slab(),
              prior.siginv_prior(),
              prior.prior_inclusion_probabilities(),
              eigenvalue_fudge_factor,
              fallback_probability);
      oda_sampler->set_sigma_upper_limit(prior.sigma_upper_limit());
      model->set_method(oda_sampler);
    }
    io_manager->add_list_element(
        new GlmCoefsListElement(model->coef_prm(), "beta"));
    io_manager->add_list_element(
        new StandardDeviationListElement(model->Sigsq_prm(), "sigma"));
    return(model);
  }
}  // namespace

extern "C" {
  using BOOM::RErrorReporter;
  using BOOM::RCheckInterrupt;
  using namespace BOOM;  // NOLINT

  SEXP do_spike_slab(SEXP r_design_matrix,
                                       SEXP r_response_vector,
                                       SEXP r_spike_slab_prior,
                                       SEXP r_niter,
                                       SEXP r_ping,
                                       SEXP r_bma_method,
                                       SEXP r_oda_options,
                                       SEXP r_seed) {
    RErrorReporter error_reporter;
    try {
      seed_rng_from_R(r_seed);
      RListIoManager io_manager;
      Ptr<RegressionModel> model = SpecifyRegressionModel(
          r_design_matrix,
          r_response_vector,
          r_spike_slab_prior,
          r_bma_method,
          r_oda_options,
          &io_manager);

      int niter = Rf_asInteger(r_niter);
      int ping = Rf_asInteger(r_ping);
      SEXP ans;
      PROTECT(ans = io_manager.prepare_to_write(niter));
      for (int i = 0; i < niter; ++i) {
        if (RCheckInterrupt()) {
          error_reporter.SetError("Canceled by user.");
          return R_NilValue;
        }
        print_R_timestamp(i, ping);
        model->sample_posterior();
        io_manager.write();
      }
      UNPROTECT(1);
      return ans;
    } catch(std::exception &e) {
      handle_exception(e);
    } catch (...) {
      handle_unknown_exception();
    }
    return R_NilValue;
  }
}
