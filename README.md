# Reassessing the ARRA Employment Multiplier: A Continuous Treatment Approach

This repository contains an original econometric research paper that completes a comprehensive reassessment of the employment effects that the American Recovery and Reinvestment Act (ARRA) spending policy carried across the fifty United States between February 2009 and December 2010.

## The Paper
Reassessing the ARRA Employment Multiplier 2026.pdf

## The Code
The analysis evaluates the evolutionary and compounding effects of the job multiplier using panel data econometric techniques, specifically applying continuous treatment estimators. The code covers the main tables and dynamic event studies presented in the paper.

### Repository Structure:
* `Reassessing the ARRA Employment Multiplier 2026.pdf`: The final output describing the continuous treatment difference-in-difference methodology and the results.
* `Scripts`: Contains the code for each step of the panel data analysis.

### Script Guide:
To replicate specific tables and figures, refer to the corresponding scripts *(Note: update the file names if yours differ)*:
* **Table 1:** `Table 1.R` - *ARRA Employment Multiplier Estimates comparing TWFE models and the CGBS ATT estimator*
* **Table 2:** `Table 2.R` - *Cost per Job Multipliers conversion across different estimators*
* **Figure 1:** `Figure 1.R` - *Pre-Period Parallel Trends tests for Median Split & Quarters*
* **Figure 2:** `Figure 2.R` - *Dynamic Event Study applying the Callaway, Goodman-Bacon, and Sant'Anna (CGBS) estimator*
* **Figure 3:** `Figure 3.R` - *ARRA Spending Pace vs. Marginal Employment Effect overtime*

### Data Gathering
* The primary data used for this paper relies on the harmonized state-level dataset created by Chodorow-Reich (2019), spanning from December 2007 up to January 2011. Granular state-level statistics reported by federal agencies were sourced from the Recovery.gov website.
* Current employment statistics (CES) and quarterly census of employment and wages (QCEW) were provided by the Bureau of Labor Statistics (BLS). Additional macro-indicative data was extracted from the Bureau of Economic Analysis (BEA).

**Tools used:** R, Two-Way Fixed Effects (TWFE), 2SLS, Callaway-Goodman-Bacon-Sant'Anna (CGBS) continuous treatment Difference-in-Differences

---
*This project was conducted as the Final Paper for Panel Data Econometrics at Maastricht University (Academic Year 2025-2026).*
