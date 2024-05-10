# Implementation for "Combining Experimental and Historical Data for Policy Evaluation"

This repository contains the implementation for the paper "Combining Experimental and Historical Data for Policy Evaluation" in R.

## Summary of the paper

This paper studies policy evaluation with multiple data sources, especially in scenarios that involve one experimental dataset with two arms, complemented by a historical dataset generated under a single control arm. We propose novel data integration methods that linearly integrate base policy value estimators constructed based on the experimental and historical data, with weights optimized to minimize the mean square error (MSE) of the resulting combined estimator. We further apply the pessimistic principle to obtain more robust estimators, and extend these developments to sequential decision making. Theoretically, we establish non-asymptotic error bounds for the MSEs of our proposed estimators, and derive their oracle, efficiency and robustness properties across a broad spectrum of reward shift scenarios. Numerical experiments and real-data-based analyses from a ridesharing company demonstrate the superior performance of the proposed estimators.

<img src=https://github.com/tingstat/Data_Combination/blob/main/single_stage.png alt=Syntheetic width="600">

** Figure 1**.  In this example, we consider a contextual bandit setting where the sample size of the experimental data is $|\mathcal{D}_e|=48$, and the sample size of the historical data is set to be $|\mathcal{D}_h|=m|\mathcal{D}_e|$ with {$m \in \{1,2,3\}$.} A deterministic switchback design is adopted to generate $\mathcal{D}_e$.
We vary the mean reward shift $b_h$ (the average reward under the control between the historical data and the experimental data) within the range from 0 to 1.2, incrementing by 0.1 at each step. We also vary the conditional variance of the reward and use $d$ to characterize this difference. 
Comparison is made among the following ATE estimators.
- NonPessi: The proposed non-pessimistic estimator.
- Pessi: The proposed pessimistic estimator.
- EDO: The doubly robust estimator $\widehat \tau_e$ constructed based on the experimental data only (see \eqref{eqn:taue}).
- Lasso: A weighted estimator $\widehat \tau_{Lasso}=w\widehat \tau_e + (1-w) \widehat \tau_h$ that linearly combines the ATE estimator $\widehat \tau_e$ based on experimental data and $\widehat \tau_h$ based on historical data, where the weight $w$ is chosen to minimize the estimated variance of the final ATE estimator with the Lasso penalty \cite{cheng2021adaptive},
- SPE: The Semiparametrically Efficient Estimator proposed by \cite{li2023improving} is developed under the assumption of no reward shift between the experimental and historical data, i.e., $r_e(0,s)=r_h(s)$ for any $s$. 


## File Overview

- Est_CombiningData.R: the estimation functions and the data generating function.
- simu_CombiningData.R: demo codes to run the experiment with respect to Example 6.1.
