# SJM-LongitudinalSurrogate

This repository contains code to implement the methods proposed in: Wang, X., Zhou, J., Parast, L. and Greene, T., 2024. Semiparametric Joint Modeling to Estimate the Treatment Effect on a Longitudinal Surrogate with Application to Chronic Kidney Disease Trials. arXiv preprint arXiv:2412.14124. 

These methods use Semiparametric Joint Modeling (SJM) to examine the treatment effect on a longitudinal surrogate marker. Specifically, these functions estimate the treatment effect on a longitudinal surrogate, specifically the slope of the longitudinal outcome, by jointly modeling the longitudinal outcome and the informative terminal event, where the model for the longitudinal outcome is semiparametric, the relationship between the longitudinal outcome and the terminal event is nonparametric, and the terminal event is modeled via a semiparametric Cox model. 

The main file is mainSJM.R. The funsSJM.R file and the two .rda data files are called from within the mainSJM.R file. Download the repository, then open mainSJM.R and run the code. The est.linear function estimates the treatment effect on a (assumed to be linear) longitudinal surrogate, while est.nonlinear relaxes the linear assumption allows for a nonlinear trajectory. The .rda data files are used to illustrate these two functions and the resulting parameter estimates and standard error estimation. 


