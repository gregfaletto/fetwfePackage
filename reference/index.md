# Package index

## Estimators

Fit fused extended two-way fixed effects (FETWFE) and its comparison
estimators for difference-in-differences with staggered adoptions.

- [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  : Fused extended two-way fixed effects
- [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
  : Extended two-way fixed effects
- [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)
  : Bridge-penalized extended two-way fixed effects
- [`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
  : Two-way fixed effects with covariates and a single pooled treatment
  effect per cohort

## Inference

Confidence intervals and debiased overall-ATT inference.

- [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  : Debiased overall-ATT with a uniformly-valid standard error
- [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
  : Parametric simultaneous (1 - alpha) confidence intervals over a
  family of treatment effects

## Treatment-effect summaries

Extract cohort, event-study, and cohort-by-time treatment effects.

- [`getTes()`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md)
  : Compute True Treatment Effects
- [`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md)
  : Per-cohort average treatment effects
- [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
  : Compute pooled event-time treatment-effect estimates
- [`cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortTimeATTs.md)
  : Per-(cohort, time) average treatment effects

## Simulation

Generate coefficients and panel data for simulation studies.

- [`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
  : Generate Coefficient Vector for Data Generation
- [`genCoefsCore()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefsCore.md)
  : Generate Coefficient Vector for Data Generation (core)
- [`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md)
  : Generate Random Panel Data for FETWFE Simulations
- [`simulateDataCore()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateDataCore.md)
  : Generate Random Panel Data for FETWFE Simulations (core)

## Fit on simulated data

Convenience drivers that fit directly on simulateData() output.

- [`fetwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfeWithSimulatedData.md)
  : Run FETWFE on Simulated Data
- [`etwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfeWithSimulatedData.md)
  : Run ETWFE on Simulated Data
- [`betwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfeWithSimulatedData.md)
  : Run BETWFE on Simulated Data
- [`twfeCovsWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovsWithSimulatedData.md)
  : Run twfeCovs on Simulated Data
- [`debiasedATTWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATTWithSimulatedData.md)
  : Run debiasedATT() on simulated data

## Data conversion

Convert data frames from other DiD packages into the FETWFE format.

- [`attgtToFetwfeDf()`](https://gregfaletto.github.io/fetwfePackage/reference/attgtToFetwfeDf.md)
  :

  Convert data formatted for `att_gt()` to a dataframe suitable for
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  /
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)

- [`etwfeToFetwfeDf()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfeToFetwfeDf.md)
  :

  Convert data prepared for `etwfe::etwfe()` to the format required by
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  and
  [`fetwfe::etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)

## Methods

broom (tidy / glance / augment) and base S3 methods (plot, print) for
the fitted objects.

- [`tidy(`*`<FETWFE_tes>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.FETWFE_tes.md)
  :

  Tidy a `FETWFE_tes` simulation truth object

- [`tidy(`*`<betwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.betwfe.md)
  :

  Tidy a `betwfe` fitted object

- [`tidy(`*`<cohortStudy>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortStudy.md)
  :

  Tidy a `cohortStudy` object

- [`tidy(`*`<cohortTimeATTs>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortTimeATTs.md)
  :

  Tidy a `cohortTimeATTs` object

- [`tidy(`*`<debiased_att>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.debiased_att.md)
  : Tidy a debiased overall-ATT estimate

- [`tidy(`*`<etwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.etwfe.md)
  :

  Tidy an `etwfe` fitted object

- [`tidy(`*`<eventStudy>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.eventStudy.md)
  :

  Tidy an `eventStudy` object

- [`tidy(`*`<fetwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.fetwfe.md)
  :

  Tidy an `fetwfe` fitted object

- [`tidy(`*`<simultaneous_cis>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.simultaneous_cis.md)
  : Tidy a simultaneous-confidence-interval object

- [`tidy(`*`<twfeCovs>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.twfeCovs.md)
  :

  Tidy a `twfeCovs` fitted object

- [`glance(`*`<betwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/glance.betwfe.md)
  :

  Glance a `betwfe` fitted object

- [`glance(`*`<etwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/glance.etwfe.md)
  :

  Glance an `etwfe` fitted object

- [`glance(`*`<fetwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/glance.fetwfe.md)
  :

  Glance an `fetwfe` fitted object

- [`glance(`*`<twfeCovs>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/glance.twfeCovs.md)
  :

  Glance a `twfeCovs` fitted object

- [`augment(`*`<betwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/augment.betwfe.md)
  : Augment user-supplied data with fitted values and residuals from a
  betwfe fit

- [`augment(`*`<etwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/augment.etwfe.md)
  : Augment user-supplied data with fitted values and residuals from an
  etwfe fit

- [`augment(`*`<fetwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/augment.fetwfe.md)
  : Augment user-supplied data with fitted values and residuals from a
  fetwfe fit

- [`augment(`*`<twfeCovs>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/augment.twfeCovs.md)
  : Augment is not defined for a twfeCovs fit (documented omission)

- [`plot(`*`<betwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/plot.betwfe.md)
  : Plot CATT or event-study estimates from a fitted BETWFE

- [`plot(`*`<etwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/plot.etwfe.md)
  : Plot CATT or event-study estimates from a fitted ETWFE

- [`plot(`*`<fetwfe>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/plot.fetwfe.md)
  : Plot CATT or event-study estimates from a fitted FETWFE / ETWFE /
  BETWFE

- [`plot(`*`<simultaneous_cis>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/plot.simultaneous_cis.md)
  : Plot a simultaneous-confidence-interval object

- [`plot(`*`<twfeCovs>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/plot.twfeCovs.md)
  : Plot is not defined for a twfeCovs fit (documented omission)

- [`print(`*`<debiased_att>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/print.debiased_att.md)
  : Print a debiased overall-ATT estimate

- [`print(`*`<simultaneous_cis>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/print.simultaneous_cis.md)
  : Print a simultaneous-confidence-interval object

## Classes and low-level accessors

Class definitions and the `catt_df` element accessors.

- [`FETWFE_coefs-class`](https://gregfaletto.github.io/fetwfePackage/reference/FETWFE_coefs-class.md)
  : FETWFE Coefficient-Vector Class

- [`FETWFE_simulated-class`](https://gregfaletto.github.io/fetwfePackage/reference/FETWFE_simulated-class.md)
  : Simulated Panel-Data Class

- [`FETWFE_tes-class`](https://gregfaletto.github.io/fetwfePackage/reference/FETWFE_tes-class.md)
  : Compute True Treatment Effects Output Class

- [`betwfe-class`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe-class.md)
  : Bridge-Penalized Extended Two-Way Fixed Effects Output Class

- [`etwfe-class`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe-class.md)
  : Extended Two-Way Fixed Effects Output Class

- [`fetwfe-class`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe-class.md)
  : Fused Extended Two-Way Fixed Effects Output Class

- [`twfeCovs-class`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs-class.md)
  : TWFE-With-Covariates Output Class

- [`` `$`( ``*`<catt_df>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/cash-.catt_df.md)
  :

  Dollar-sign access on a `catt_df` object

- [`` `$<-`( ``*`<catt_df>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/cash-set-.catt_df.md)
  :

  Dollar-sign assignment on a `catt_df` object

- [`` `[`( ``*`<catt_df>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/sub-.catt_df.md)
  :

  Single-bracket access on a `catt_df` object

- [`` `[[`( ``*`<catt_df>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/sub-sub-.catt_df.md)
  :

  Double-bracket access on a `catt_df` object

- [`` `[[<-`( ``*`<catt_df>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/sub-subset-.catt_df.md)
  :

  Double-bracket assignment on a `catt_df` object

- [`` `[<-`( ``*`<catt_df>`*`)`](https://gregfaletto.github.io/fetwfePackage/reference/subset-.catt_df.md)
  :

  Single-bracket assignment on a `catt_df` object
