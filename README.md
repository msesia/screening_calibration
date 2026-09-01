# Distribution-Free Selection of Low-Risk Oncology Patients for Survival Beyond a Time Horizon


This repository contains R code accompanying the following paper:

[Distribution-Free Selection of Low-Risk Oncology Patients for Survival Beyond a Time Horizon](https://arxiv.org/pdf/2512.18118v1)


## Paper abstract

We study how to select a subset of patients who are unlikely to experience an adverse event within a given time horizon, by calibrating a screening rule based on the output of any survival model. We consider two complementary frameworks. The first extends the classical idea of estimating the event rate among selected patients using a hold-out dataset, integrating it with the Learn-Then-Test method. This provides approximate high-probability guarantees that are comparable to those obtainable from simultaneous confidence bands while being often less conservative. The second takes a different perspective by reformulating the problem in terms of multiple hypotheses testing, enabling false discovery rate (FDR) control via the Benjamini-Hochberg procedure applied to selective conformal p-values. This provides approximate guarantees in expectation. We clarify the theoretical relationship between these approaches, explain how both can handle right-censored data and be made doubly robust via augmented inverse probability of censoring weighting, and compare them empirically using simulations and oncology data from the Flatiron Health Research Database. Our results reveal a practical trade-off: FDR-based screening is typically more powerful, while high-probability calibration is more conservative but offers stronger guarantees, especially when few patients are selected. We also provide practical guidance on implementation and tuning.

