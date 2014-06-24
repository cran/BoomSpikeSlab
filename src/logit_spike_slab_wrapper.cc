// Copyright 2011 Google Inc. All Rights Reserved.
// Author: stevescott@google.com (Steve Scott)

#include <exception>
#include <string>

#include "cpputil/Ptr.hpp"
#include "r_interface/print_R_timestamp.hpp"
#include "r_interface/handle_exception.hpp"

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
//   beta0:  C array containing initial coefficients.
//   number_of_coefficients:  Length of beta0.
//   prior_inclusion_probabilities: Prior probabilities that each
//     coefficient is nonzero.
//   model:  The model that owns the coefficients.
//   sampler:  The sampler that will make posterior draws for the model.
inline void InitializeCoefficients(
    const double *beta0, int number_of_coefficients,
    const BOOM::Vec & prior_inclusion_probabilities,
    BOOM::Ptr<BOOM::BinomialLogitModel> model,
    BOOM::Ptr<BOOM::BinomialLogitCompositeSpikeSlabSampler> sampler) {
  BOOM::Vec initial_beta(beta0, beta0 + number_of_coefficients);
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
    for (int i = 0; i < number_of_coefficients; ++i) {
      if (fabs(beta0[i]) < 1e-8) {
        model->coef().drop(i);
      } else {
        model->coef().add(i);
      }
    }
  }
}

}  // unnamed namespace

extern "C" {
  using BOOM::Spd;
  using BOOM::Mat;
  using BOOM::Vec;
  using BOOM::ConstVectorView;
  using BOOM::VectorView;
  using BOOM::MvnModel;
  using BOOM::VariableSelectionPrior;
  using BOOM::BinomialLogitModel;
  using BOOM::BinomialRegressionData;
  using BOOM::BinomialLogitSpikeSlabSampler;
  using BOOM::BinomialLogitCompositeSpikeSlabSampler;
  using BOOM::Ptr;

  // This function is a wrapper for spike and slab regression.  It
  // takes input from R and formats it for the appropriate BOOM
  // objects that handle the comptuations.
  void logit_spike_slab_wrapper(
      const double *x,              // design matrix
      const int *y,                 // vector of success counts
      const int *ny,                // vector of trial counts
      const double *mu,             // prior mean vector
      const double *sigma_inverse,  // prior precision matrix
      const double *pi,             // prior inclusion probabilities
      const int *max_flips,         // max number of indicators to sample
      const int *n,                 // length of y, nrow(x)
      const int *p,                 // ncol(x), ncol(ans)
      const int *niter,             // number of mcmc iterations
      const int *ping,              // frequency of desired status updates
      const int *nthreads,          // number of imputation threads
      const double *beta0,          // initial value in the MCMC simulation
      const int *clt_threshold,     // see comments in ../R/logit.spike.R
      const int *mh_chunk_size,     // see comments in ../R/logit.spike.R
      const int *seed,
      double *beta_draws            // output:  mcmc draws of beta
                                ) {
    try {
      BOOM::GlobalRng::rng.seed(*seed);
      Ptr<BinomialLogitModel>  model(new BinomialLogitModel(*p));
      for (int i = 0; i < *n; ++i) {
        ConstVectorView xi(x+i, *p, *n);
        Ptr<BinomialRegressionData>
            dp(new BinomialRegressionData(y[i], ny[i], xi));
        model->add_data(dp);
      }

      // Copy data to BOOM structures in preparation for building the
      // model
      Vec Mu(mu, mu + (*p));
      Spd Siginv(Mat(sigma_inverse,
                     sigma_inverse + (*p)*(*p),
                     *p,
                     *p));
      Vec prior_inclusion_probabilities(pi, pi + (*p));

      Ptr<MvnModel> beta_prior(new MvnModel(Mu, Siginv, true));
      Ptr<VariableSelectionPrior> inclusion_prior(
          new VariableSelectionPrior(
              prior_inclusion_probabilities));

      double proposal_degrees_of_freedom = 3;
      int max_chunk_size = *mh_chunk_size;
      Ptr<BOOM::BinomialLogitCompositeSpikeSlabSampler> sampler(
          new BinomialLogitCompositeSpikeSlabSampler(
              model.get(),
              beta_prior,
              inclusion_prior,
              *clt_threshold,
              proposal_degrees_of_freedom,
              max_chunk_size));
      if (*nthreads > 1) {
        // TODO(stevescott): Add threading support back into
        // BOOM::BinomialLogitSpikeSlabSampler.
      }
      sampler->limit_model_selection(*max_flips);

      InitializeCoefficients(beta0, *p, prior_inclusion_probabilities,
                             model, sampler);

      for (int i = 0; i < *niter; ++i) {
        // TODO(stevescott): Find a way to allow user interrupts that
        // ensures smart-pointer protected resources are freed.
        R_CheckUserInterrupt();
        BOOM::print_R_timestamp(i, *ping);
        sampler->draw();
        // Copy the draw to row i of the matrix beta_draws.
        // beta_view refers to row i of beta_draws.
        // It has length *p, with stride *niter.
        VectorView beta_view(beta_draws + i, *p, *niter);
        beta_view = model->Beta();
      }
    } catch(std::exception &e) {
      BOOM::RInterface::handle_exception(e);
    } catch(...) {
      BOOM::RInterface::handle_unknown_exception();
    }
  }
}
