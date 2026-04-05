# Reassessing the ARRA Employment Multiplier: A Continuous Treatment Approach

This repository contains an original econometric research paper that completes a comprehensive reassessment of the employment effects that the American Recovery and Reinvestment Act (ARRA) spending policy carried across the fifty United States between February 2009 and December 2010.

## The Paper
Reassessing the ARRA Employment Multiplier 2026.pdf

## The Code
The analysis evaluates the evolutionary and compounding effects of the job multiplier using panel data econometric techniques, specifically applying continuous treatment estimators. The code covers the main tables and dynamic event studies presented in the paper.

### Repository Structure:
* `Tommaso Maccaferri-Reassessing the ARRA Employment Multiplier 2026.pdf`: The final output describing the continuous treatment difference-in-difference methodology and the results.
* `Script: Reassessing the ARRA Employment Multiplier.R`: Contains the code for each step of the panel data analysis.

### Script Guide:
The Reassessing the ARRA Employment Multiplier.R script contains:
- *Pre-Period Parallel Trends tests for Median Split & Quarters*
- *Chodorow-Reich (2019) original estimation using IV*
- *Dynamic Event Study applying the Callaway and Sant'Anna (2021) estimator*
- *Two-way Fixed Effects estimation*
- *Callaway, Goodman-Bacon, and Sant'Anna (CGBS) estimation*
- *Dynamic Event Study applying CGBS estimator*
- *Cost per Job Multipliers conversion across different estimators*
- *ARRA Employment Multiplier Estimates comparing baseline IV model, TWFE models and the CGBS ATT estimator*

### Data Gathering
* The primary data used for this paper relies on the harmonized state-level dataset created by Chodorow-Reich (2019), spanning from December 2007 up to January 2011. Granular state-level statistics reported by federal agencies were sourced from the Recovery.gov website.
* Current employment statistics (CES) and quarterly census of employment and wages (QCEW) were provided by the Bureau of Labor Statistics (BLS). Additional macro-indicative data was extracted from the Bureau of Economic Analysis (BEA).

**Tools used:** R, Two-Way Fixed Effects (TWFE), 2SLS, Callaway-Goodman-Bacon-Sant'Anna (CGBS) continuous treatment Difference-in-Differences

---
*This project was conducted as the Final Paper for Panel Data Econometrics at Maastricht University (Academic Year 2025-2026).*
