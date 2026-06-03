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
        Selected size       : 1
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
        Selected size       : 1
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
        Selected size       : 1
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
        Selected size       : 1
        Lambda*             : 0.0265

