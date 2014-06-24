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
#include "r_interface/seed_rng_from_R.hpp"

#include "cpputil/Ptr.hpp"

namespace {
  using namespace BOOM;  // NOLINT
  using namespace BOOM::RInterface;  // NOLINT

  Ptr<IndependentMvnModelGivenScalarSigma>
  ExtractIndependentBetaPrior(SEXP prior, Ptr<UnivParams> sigsq) {
    Vector prior_mean(ToBoomVector(getListElement(prior, "mu")));
    Vector prior_variance_diagonal(
        ToBoomVector(getListElement(prior, "prior.variance.diagonal")));
    return new IndependentMvnModelGivenScalarSigma(
        new VectorParams(prior_mean),
        new VectorParams(prior_variance_diagonal),
        sigsq);
  }

  Ptr<MvnGivenScalarSigma>
  ExtractZellnerPrior(SEXP prior, Ptr<UnivParams> sigsq) {
    Vector prior_mean(ToBoomVector(getListElement(prior, "mu")));
    SpdMatrix prior_precision(ToBoomSpd(getListElement(prior, "siginv")));
    return new MvnGivenScalarSigma(prior_mean, prior_precision, sigsq);
  }

  Ptr<RegressionModel> SpecifyRegressionModel(
      SEXP r_design_matrix,
      SEXP r_response_vector,
      SEXP r_spike_slab_prior,
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

    Vector prior_inclusion_probabilities =
        ToBoomVector(getListElement(
            r_spike_slab_prior, "prior.inclusion.probabilities"));

    double prior_df = Rf_asReal(getListElement(
        r_spike_slab_prior, "prior.df"));
    double sigma_guess = Rf_asReal(getListElement(
        r_spike_slab_prior, "sigma.guess"));
    Ptr<ChisqModel> siginv_prior(new ChisqModel(prior_df, sigma_guess));

    if (Rf_inherits(r_spike_slab_prior, "SpikeSlabPrior")) {
      Ptr<MvnGivenScalarSigma> beta_prior = ExtractZellnerPrior(
          r_spike_slab_prior, model->Sigsq_prm());
      Ptr<BregVsSampler> sampler(new BregVsSampler(
          model.get(),
          beta_prior,
          siginv_prior,
          new VariableSelectionPrior(prior_inclusion_probabilities)));
      int max_flips = Rf_asInteger(getListElement(
          r_spike_slab_prior, "max.flips"));
      if (max_flips > 0) {
        sampler->limit_model_selection(max_flips);
      }
      model->set_method(sampler);
    } else if (Rf_inherits(r_spike_slab_prior, "IndependentSpikeSlabPrior")) {
      Ptr<IndependentMvnModelGivenScalarSigma> beta_prior =
          ExtractIndependentBetaPrior(r_spike_slab_prior, model->Sigsq_prm());
      Ptr<SpikeSlabDaRegressionSampler> sampler(
          new SpikeSlabDaRegressionSampler(
              model.get(),
              beta_prior,
              siginv_prior,
              prior_inclusion_probabilities));
      model->set_method(sampler);
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
                                       SEXP r_seed) {
    RErrorReporter error_reporter;
    try {
      seed_rng_from_R(r_seed);
      RListIoManager io_manager;
      Ptr<RegressionModel> model = SpecifyRegressionModel(
          r_design_matrix,
          r_response_vector,
          r_spike_slab_prior,
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
