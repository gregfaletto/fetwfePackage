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
      
      Cohort Average Treatment Effects (CATT):
       cohort estimate se ci_low ci_high p_value selected
            2        0  0      0       0      NA    FALSE
            3        0  0      0       0      NA    FALSE
      
      Event-Study Average Treatment Effects (per event time):
       event_time n_cohorts estimate se ci_low ci_high p_value
                0         2        0 NA     NA      NA      NA
                1         2        0 NA     NA      NA      NA
                2         2        0 NA     NA      NA      NA
                3         1        0 NA     NA      NA      NA
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (R) : 2
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
      
      CATT (preview):
       cohort estimate se ci_low ci_high p_value selected
            2        0  0      0       0      NA    FALSE
            3        0  0      0       0      NA    FALSE
      
      Event Study (preview):
       event_time n_cohorts estimate se ci_low ci_high p_value
                0         2        0 NA     NA      NA      NA
                1         2        0 NA     NA      NA      NA
                2         2        0 NA     NA      NA      NA
                3         1        0 NA     NA      NA      NA
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (R) : 2
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
      
      Cohort Average Treatment Effects (CATT):
       cohort     estimate        se     ci_low   ci_high   p_value
            2  0.009060593 0.4879617 -0.9473267 0.9654479 0.9851855
            3 -0.127673929 0.4655589 -1.0401525 0.7848047 0.7839017
      
      Event-Study Average Treatment Effects (per event time):
       event_time n_cohorts    estimate        se     ci_low   ci_high   p_value
                0         2 -0.50429204 0.4341539 -1.3552181 0.3466340 0.2454178
                1         2 -0.05370444 0.4826019 -0.9995867 0.8921778 0.9113935
                2         2  0.15743661 0.4826578 -0.7885553 1.1034285 0.7442830
                3         1  0.44849429 0.6429757 -0.8117150 1.7087036 0.4854717
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (R) : 2
        Covariates (d)      : 2
        Features (p)        : 41

# print.summary.etwfe output is stable

    Code
      print(summary(fit_etwfe))
    Output
      Summary of Extended Two-Way Fixed Effects
      ========================================
      
      Overall ATT: -0.0456  (SE = 0.3843, p = 0.9055, 95% CI = [-0.7988, 0.7075])
      
      CATT (preview):
       cohort     estimate        se     ci_low   ci_high   p_value
            2  0.009060593 0.4879617 -0.9473267 0.9654479 0.9851855
            3 -0.127673929 0.4655589 -1.0401525 0.7848047 0.7839017
      
      Event Study (preview):
       event_time n_cohorts    estimate        se     ci_low   ci_high   p_value
                0         2 -0.50429204 0.4341539 -1.3552181 0.3466340 0.2454178
                1         2 -0.05370444 0.4826019 -0.9995867 0.8921778 0.9113935
                2         2  0.15743661 0.4826578 -0.7885553 1.1034285 0.7442830
                3         1  0.44849429 0.6429757 -0.8117150 1.7087036 0.4854717
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (R) : 2
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
      
      Cohort Average Treatment Effects (CATT):
       cohort estimate se ci_low ci_high p_value selected
            2        0  0      0       0      NA    FALSE
            3        0  0      0       0      NA    FALSE
      
      Event-Study Average Treatment Effects (per event time):
       event_time n_cohorts estimate se ci_low ci_high p_value
                0         2        0 NA     NA      NA      NA
                1         2        0 NA     NA      NA      NA
                2         2        0 NA     NA      NA      NA
                3         1        0 NA     NA      NA      NA
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (R) : 2
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
      
      CATT (preview):
       cohort estimate se ci_low ci_high p_value selected
            2        0  0      0       0      NA    FALSE
            3        0  0      0       0      NA    FALSE
      
      Event Study (preview):
       event_time n_cohorts estimate se ci_low ci_high p_value
                0         2        0 NA     NA      NA      NA
                1         2        0 NA     NA      NA      NA
                2         2        0 NA     NA      NA      NA
                3         1        0 NA     NA      NA      NA
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (R) : 2
        Covariates (d)      : 2
        Features (p)        : 41
        Selected size       : 1
        Lambda*             : 0.0265

