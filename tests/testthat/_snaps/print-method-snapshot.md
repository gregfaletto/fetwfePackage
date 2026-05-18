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
       Cohort Estimated TE SE ConfIntLow ConfIntHigh P_value selected
            2            0  0          0           0      NA    FALSE
            3            0  0          0           0      NA    FALSE
      
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
       Cohort Estimated TE SE ConfIntLow ConfIntHigh P_value selected
            2            0  0          0           0      NA    FALSE
            3            0  0          0           0      NA    FALSE
      
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
        Std. Error: 0.4012
        P-value:    0.9094
        95% CI:    [-0.8319, 0.7406]
      
      Cohort Average Treatment Effects (CATT):
       Cohort Estimated TE        SE ConfIntLow ConfIntHigh   P_value
            2  0.009060593 0.4879617 -0.9473267   0.9654479 0.9851855
            3 -0.127673929 0.4655589 -1.0401525   0.7848047 0.7839017
      
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
      
      Overall ATT: -0.0456  (SE = 0.4012, p = 0.9094, 95% CI = [-0.8319, 0.7406])
      
      CATT (preview):
       Cohort Estimated TE        SE ConfIntLow ConfIntHigh   P_value
            2  0.009060593 0.4879617 -0.9473267   0.9654479 0.9851855
            3 -0.127673929 0.4655589 -1.0401525   0.7848047 0.7839017
      
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
       Cohort Estimated TE SE ConfIntLow ConfIntHigh P_value selected
            2            0  0          0           0      NA    FALSE
            3            0  0          0           0      NA    FALSE
      
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
       Cohort Estimated TE SE ConfIntLow ConfIntHigh P_value selected
            2            0  0          0           0      NA    FALSE
            3            0  0          0           0      NA    FALSE
      
      Model Details:
        Units (N)           : 30
        Time periods (T)    : 5
        Treated cohorts (R) : 2
        Covariates (d)      : 2
        Features (p)        : 41
        Selected size       : 1
        Lambda*             : 0.0265

