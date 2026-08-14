# print.fetwfe output is stable

    Code
      print(fit_fetwfe)
    Output
      Fused Extended Two-Way Fixed Effects Results
      ===========================================
      
      Overall Average Treatment Effect (ATT):
        Estimate:   0.0000
        Std. Error: 0.0000
        P-value:    NA
        Selected:   FALSE
        95% CI:    [0.0000, 0.0000]
      
      Cohort Average Treatment Effects (CATT) [simultaneous 95% CI]:
       cohort estimate se ci_low ci_high p_value selected
            2        0  0      0       0      NA    FALSE
            3        0  0      0       0      NA    FALSE
      
      Event-Study Average Treatment Effects (per event time) [simultaneous 95% CI]:
       event_time n_cohorts estimate se ci_low ci_high p_value
                0         2        0 NA     NA      NA      NA
                1         2        0 NA     NA      NA      NA
                2         2        0 NA     NA      NA      NA
                3         1        0 NA     NA      NA      NA
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (G) : 2
        Covariates (d)      : 2
        Features (p)        : 41
        Selected size       : 0
        Lambda*             : 0.0299

# print.summary.fetwfe output is stable

    Code
      print(summary(fit_fetwfe))
    Output
      Summary of Fused Extended Two-Way Fixed Effects
      ================================================
      
      Overall ATT: 0.0000  (SE = 0.0000, p = NA, 95% CI = [0.0000, 0.0000])
      Selected: FALSE
      
      CATT (preview) [simultaneous 95% CI]:
       cohort estimate se ci_low ci_high p_value selected
            2        0  0      0       0      NA    FALSE
            3        0  0      0       0      NA    FALSE
      
      Event Study (preview) [simultaneous 95% CI]:
       event_time n_cohorts estimate se ci_low ci_high p_value
                0         2        0 NA     NA      NA      NA
                1         2        0 NA     NA      NA      NA
                2         2        0 NA     NA      NA      NA
                3         1        0 NA     NA      NA      NA
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (G) : 2
        Covariates (d)      : 2
        Features (p)        : 41
        Selected size       : 0
        Lambda*             : 0.0299

# print.etwfe output is stable

    Code
      print(fit_etwfe)
    Output
      Extended Two-Way Fixed Effects Results
      =====================================
      
      Overall Average Treatment Effect (ATT):
        Estimate:   -0.0456
        Std. Error: 0.3843
        P-value:    0.9055
        95% CI:    [-0.7988, 0.7075]
      
      Cohort Average Treatment Effects (CATT) [simultaneous 95% CI]:
       cohort     estimate        se    ci_low   ci_high   p_value
            2  0.009060593 0.4879617 -1.079644 1.0977651 0.9997735
            3 -0.127673929 0.4655589 -1.166395 0.9110471 0.9518819
      
      Event-Study Average Treatment Effects (per event time) [simultaneous 95% CI]:
       event_time n_cohorts    estimate        se    ci_low   ci_high   p_value
                0         2 -0.50429204 0.4341539 -1.563527 0.5549433 0.5918871
                1         2 -0.05370444 0.4826019 -1.231142 1.1237328 0.9998858
                2         2  0.15743661 0.4826578 -1.020137 1.3350103 0.9923936
                3         1  0.44849429 0.6429757 -1.120218 2.0172067 0.8916439
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (G) : 2
        Covariates (d)      : 2
        Features (p)        : 41

# print.summary.etwfe output is stable

    Code
      print(summary(fit_etwfe))
    Output
      Summary of Extended Two-Way Fixed Effects
      ========================================
      
      Overall ATT: -0.0456  (SE = 0.3843, p = 0.9055, 95% CI = [-0.7988, 0.7075])
      
      CATT (preview) [simultaneous 95% CI]:
       cohort     estimate        se    ci_low   ci_high   p_value
            2  0.009060593 0.4879617 -1.079644 1.0977651 0.9997735
            3 -0.127673929 0.4655589 -1.166395 0.9110471 0.9518819
      
      Event Study (preview) [simultaneous 95% CI]:
       event_time n_cohorts    estimate        se    ci_low   ci_high   p_value
                0         2 -0.50429204 0.4341539 -1.563527 0.5549433 0.5918871
                1         2 -0.05370444 0.4826019 -1.231142 1.1237328 0.9998858
                2         2  0.15743661 0.4826578 -1.020137 1.3350103 0.9923936
                3         1  0.44849429 0.6429757 -1.120218 2.0172067 0.8916439
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (G) : 2
        Covariates (d)      : 2
        Features (p)        : 41

# print.betwfe output is stable

    Code
      print(fit_betwfe)
    Output
      Bridge-Penalized Extended Two-Way Fixed Effects Results
      =======================================================
      
      Overall Average Treatment Effect (ATT):
        Estimate:   0.0000
        Std. Error: 0.0000
        P-value:    NA
        Selected:   FALSE
        95% CI:    [0.0000, 0.0000]
      
      Cohort Average Treatment Effects (CATT) [simultaneous 95% CI]:
       cohort estimate se ci_low ci_high p_value selected
            2        0  0      0       0      NA    FALSE
            3        0  0      0       0      NA    FALSE
      
      Event-Study Average Treatment Effects (per event time) [simultaneous 95% CI]:
       event_time n_cohorts estimate se ci_low ci_high p_value
                0         2        0 NA     NA      NA      NA
                1         2        0 NA     NA      NA      NA
                2         2        0 NA     NA      NA      NA
                3         1        0 NA     NA      NA      NA
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (G) : 2
        Covariates (d)      : 2
        Features (p)        : 41
        Selected size       : 0
        Lambda*             : 0.0265

# print.summary.betwfe output is stable

    Code
      print(summary(fit_betwfe))
    Output
      Summary of Bridge-Penalized Extended Two-Way Fixed Effects
      ==========================================================
      
      Overall ATT: 0.0000  (SE = 0.0000, p = NA, 95% CI = [0.0000, 0.0000])
      Selected: FALSE
      
      CATT (preview) [simultaneous 95% CI]:
       cohort estimate se ci_low ci_high p_value selected
            2        0  0      0       0      NA    FALSE
            3        0  0      0       0      NA    FALSE
      
      Event Study (preview) [simultaneous 95% CI]:
       event_time n_cohorts estimate se ci_low ci_high p_value
                0         2        0 NA     NA      NA      NA
                1         2        0 NA     NA      NA      NA
                2         2        0 NA     NA      NA      NA
                3         1        0 NA     NA      NA      NA
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (G) : 2
        Covariates (d)      : 2
        Features (p)        : 41
        Selected size       : 0
        Lambda*             : 0.0265

# print.twfeCovs output is stable

    Code
      print(fit_twfecovs)
    Output
      TWFE (with covariates) Results
      ==============================
      
      Overall Average Treatment Effect (ATT):
        Estimate:   -0.3605
        Std. Error: 0.3549
        P-value:    0.3098
        95% CI:    [-1.0561, 0.3351]
      
      Cohort Average Treatment Effects (CATT) [simultaneous 95% CI]:
       cohort   estimate        se    ci_low   ci_high   p_value
            2 -0.3090865 0.4615481 -1.339713 0.7215401 0.7493122
            3 -0.4375667 0.4317387 -1.401630 0.5264962 0.5202684
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (G) : 2
        Covariates (d)      : 2
        Features (p)        : 10

# print.summary.twfeCovs output is stable

    Code
      print(summary(fit_twfecovs))
    Output
      Summary of TWFE (with covariates)
      =================================
      
      Overall ATT: -0.3605  (SE = 0.3549, p = 0.3098, 95% CI = [-1.0561, 0.3351])
      
      CATT (preview) [simultaneous 95% CI]:
       cohort   estimate        se    ci_low   ci_high   p_value
            2 -0.3090865 0.4615481 -1.339713 0.7215401 0.7493122
            3 -0.4375667 0.4317387 -1.401630 0.5264962 0.5202684
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (G) : 2
        Covariates (d)      : 2
        Features (p)        : 10

# print.fetwfe output is stable on a distinct-dimension fit

    Code
      print(fit_distinct)
    Output
      Fused Extended Two-Way Fixed Effects Results
      ===========================================
      
      Overall Average Treatment Effect (ATT):
        Estimate:   -0.2039
        Std. Error: 0.1290
        P-value:    0.1139
        Selected:   TRUE
        95% CI:    [-0.4566, 0.0489]
      
      Cohort Average Treatment Effects (CATT) [simultaneous 95% CI]:
       cohort   estimate        se     ci_low    ci_high   p_value selected
            2 -0.2038581 0.1289615 -0.4566691 0.04895278 0.1139313     TRUE
            3 -0.2038581 0.1289615 -0.4566691 0.04895278 0.1139313     TRUE
            4 -0.2038581 0.1289615 -0.4566691 0.04895278 0.1139313     TRUE
      
      Event-Study Average Treatment Effects (per event time) [simultaneous 95% CI]:
       event_time n_cohorts   estimate        se     ci_low    ci_high   p_value
                0         3 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
                1         3 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
                2         3 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
                3         2 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
                4         1 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
      
      Model Details:
        Units (N)           : 40
        Time periods (T)    : 6
        Treated cohorts (G) : 3
        Covariates (d)      : 1
        Features (p)        : 41
        Selected size       : 1
        Lambda*             : 0.0396

# print.summary.fetwfe output is stable on a distinct-dimension fit

    Code
      print(summary(fit_distinct))
    Output
      Summary of Fused Extended Two-Way Fixed Effects
      ================================================
      
      Overall ATT: -0.2039  (SE = 0.1290, p = 0.1139, 95% CI = [-0.4566, 0.0489])
      Selected: TRUE
      
      CATT (preview) [simultaneous 95% CI]:
       cohort   estimate        se     ci_low    ci_high   p_value selected
            2 -0.2038581 0.1289615 -0.4566691 0.04895278 0.1139313     TRUE
            3 -0.2038581 0.1289615 -0.4566691 0.04895278 0.1139313     TRUE
            4 -0.2038581 0.1289615 -0.4566691 0.04895278 0.1139313     TRUE
      
      Event Study (preview) [simultaneous 95% CI]:
       event_time n_cohorts   estimate        se     ci_low    ci_high   p_value
                0         3 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
                1         3 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
                2         3 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
                3         2 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
                4         1 -0.2038581 0.1289615 -0.4566706 0.04895432 0.1139313
      
      Model Details:
        Units (N)           : 40
        Time periods (T)    : 6
        Treated cohorts (G) : 3
        Covariates (d)      : 1
        Features (p)        : 41
        Selected size       : 1
        Lambda*             : 0.0396

