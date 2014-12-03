// Copyright 2011 Google Inc. All Rights Reserved.
// Author: stevescott@google.com (Steve Scott)

#include <exception>
#include <string>

#include "cpputil/Ptr.hpp"
#include "r_interface/print_R_timestamp.hpp"
#include "r_interface/handle_exception.hpp"
#include "r_interface/seed_rng_from_R.hpp"
#include "r_interface/prior_specification.hpp"
#include "r_interface/list_io.hpp"
#include "r_interface/boom_r_tools.hpp"

#include "Models/Glm/BinomialLogitModel.hpp"
#include "Models/Glm/PosteriorSamplers/BinomialLogitCompositeSpikeSlabSampler.hpp"
#include "Models/Glm/PosteriorSamplers/BinomialLogitSpikeSlabSampler.hpp"
#include "Models/Glm/VariableSelectionPrior.hpp"
#include "Models/MvnModel.hpp"

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include "Rinternals.h"
#include "R.h"

namespace {
// Set the coefficients equal to their initial values, and determine
// which coefficients are initially excluded (i.e. forced to zero).
// Args:
//   initial_beta:  Vector containing initial coefficients.
//   prior_inclusion_probabilities: Prior probabilities that each
//     coefficient is nonzero.
//   model:  The model that owns the coefficients.
//   sampler:  The sampler that will make posterior draws for the model.
inline void InitializeCoefficients(
    const BOOM::Vector &initial_beta,
    const BOOM::Vector &prior_inclusion_probabilities,
    BOOM::Ptr<BOOM::BinomialLogitModel> model,
    BOOM::Ptr<BOOM::BinomialLogitCompositeSpikeSlabSampler> sampler) {
  model->set_Beta(initial_beta);
  if (min(prior_inclusion_probabilities) >= 1.0) {
    // Ensure all coefficients are included if you're not going to
    // do model averaging.
    sampler->allow_model_selection(false);
    model->coef().add_all();
  } else {
    // Model averaging is desired.  "Small" coefficients start off
    // excluded from the model.  Large ones start off included.
    // Adding or dropping is idempotent, so no need to worry about
    // dropping an already excluded coefficient.
    for (int i = 0; i < initial_beta.size(); ++i) {
      if (fabs(initial_beta[i]) < 1e-8) {
        model->coef().drop(i);
      } else {
        model->coef().add(i);
      }
    }
  }
}

}  // unnamed namespace

extern "C" {
  using BOOM::SpdMatrix;
  using BOOM::Matrix;
  using BOOM::Vector;
  using BOOM::ConstVectorView;
  using BOOM::VectorView;
  using BOOM::MvnModel;
  using BOOM::VariableSelectionPrior;
  using BOOM::Ptr;

  // This function is a wrapper for spike and slab regression.  It
  // takes input from R and formats it for the appropriate BOOM
  // objects that handle the comptuations.
  SEXP logit_spike_slab_wrapper(
      SEXP r_x,              // design matrix
      SEXP r_y,              // vector of success counts
      SEXP r_ny,             // vector of trial counts
      SEXP r_prior,          // SpikeSlabPrior
      SEXP r_niter,          // number of mcmc iterations
      SEXP r_ping,           // frequency of desired status updates
      SEXP r_nthreads,       // number of imputation threads
      SEXP r_beta0,          // initial value in the MCMC simulation
      SEXP r_clt_threshold,  // see comments in ../R/logit.spike.R
      SEXP r_mh_chunk_size,  // see comments in ../R/logit.spike.R
      SEXP r_seed)  {
    try {
      BOOM::RInterface::seed_rng_from_R(r_seed);
      Matrix design_matrix(BOOM::ToBoomMatrix(r_x));
      std::vector<int> successes(BOOM::ToIntVector(r_y));
      std::vector<int> trials(BOOM::ToIntVector(r_ny));
      Ptr<BOOM::BinomialLogitModel> model(new BOOM::BinomialLogitModel(
          design_matrix.ncol()));
      int number_of_observations = successes.size();
      for (int i = 0; i < number_of_observations; ++i) {
        Ptr<BOOM::BinomialRegressionData>
            dp(new BOOM::BinomialRegressionData(
                successes[i],
                trials[i],
                design_matrix.row(i)));
        model->add_data(dp);
      }

      BOOM::RInterface::SpikeSlabGlmPrior prior(r_prior);

      double proposal_degrees_of_freedom = 3;
      int max_tim_chunk_size = Rf_asInteger(r_mh_chunk_size);
      int max_rwm_chunk_size = 1;
      double rwm_variance_scale_factor = .025;
      Ptr<BOOM::BinomialLogitCompositeSpikeSlabSampler> sampler(
          new BOOM::BinomialLogitCompositeSpikeSlabSampler(
              model.get(),
              prior.slab(),
              prior.spike(),
              Rf_asInteger(r_clt_threshold),
              proposal_degrees_of_freedom,
              max_tim_chunk_size,
              max_rwm_chunk_size,
              rwm_variance_scale_factor));
      int nthreads = Rf_asInteger(r_nthreads);
      if (nthreads > 1) {
        sampler->set_number_of_workers(nthreads);
      }

      if (prior.max_flips() > 0) {
        sampler->limit_model_selection(prior.max_flips());
      }

      InitializeCoefficients(BOOM::ToBoomVector(r_beta0),
                             prior.prior_inclusion_probabilities(),
                             model,
                             sampler);

      int niter = Rf_asInteger(r_niter);
      BOOM::RListIoManager io_manager;
      io_manager.add_list_element(
          new BOOM::GlmCoefsListElement(
              model->coef_prm(),
              "beta"));
      SEXP ans;
      PROTECT(ans = io_manager.prepare_to_write(niter));
      int ping = Rf_asInteger(r_ping);

      for (int i = 0; i < niter; ++i) {
        R_CheckUserInterrupt();
        BOOM::print_R_timestamp(i, ping);
        sampler->draw();
        io_manager.write();
      }
      UNPROTECT(1);
      return ans;
    } catch(std::exception &e) {
      BOOM::RInterface::handle_exception(e);
    } catch(...) {
      BOOM::RInterface::handle_unknown_exception();
    }
    return R_NilValue;
  }
}
