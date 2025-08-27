# Advanced Statistical Methods Explained

**Note:** This document provides a detailed explanation of some of the key multivariate statistical methods used in this project. It is compiled and expanded from two internal reports (`PRC.pdf` and `Power estimation for the multivariate microbiome data.pdf`), which are now stored in the `./archive` directory for historical reference.

---

## 1. Principal Response Curves (PRC)

### Theoretical Background

Principal Response Curves (PRC) are a specialized form of multivariate analysis designed specifically for longitudinal studies where the goal is to assess the effects of a treatment over time, relative to an untreated control group. It is a direct extension of Redundancy Analysis (RDA) and directly answers the question: **"How does the microbial community in the treatment group deviate from the control group at each time point?"**

The method works by modeling the community composition data (e.g., ASV abundances) with respect to two main constraints: **time** and **treatment**. The resulting plot is powerful because it distills this complex, high-dimensional data into a simple and intuitive visualization:

*   The **x-axis** represents the different time points of the study.
*   The **y-axis** (the first canonical axis) represents the principal treatment effect. By definition, the **control group is represented by the horizontal line at y=0**. The plotted line shows the trajectory of the treatment group's community composition *relative to the control*. A significant deviation from the zero line indicates a treatment effect.
*   The **right-hand axis** shows the "species scores" (in our case, ASV scores). These scores indicate how much each individual ASV contributes to the overall treatment effect shown in the main plot.

### How to Interpret the PRC Plot

Interpreting a PRC plot requires understanding the relationship between the main trajectory line and the species scores.

1.  **Main Effect Line:** If the line for the test group is, for example, in the negative region of the y-axis, it means the community has changed relative to the control. The farther the line is from zero, the larger the difference.
2.  **ASV/Species Scores:** The interpretation of these scores depends on the sign of the main effect line.
    *   If the main effect line is **negative**, an ASV with a **negative** score is predicted to be *more abundant* in the test group compared to the control. Conversely, an ASV with a **positive** score is predicted to be *less abundant*.
    *   If the main effect line were positive, the interpretation would be reversed.

In our study, the PRC plot (Figure 4 in the paper) shows a clear effect of the OMT after Week 2. For instance, ASV40 (Chlorobi G-1) has a large negative score (-1.5), and since the main PRC line is also negative, we infer that ASV40's abundance increased significantly in the test group after OMT.

### Implementation in this Project

The PRC analysis was performed to visualize the temporal impact of the Oral Microbiota Transplant (OMT) on the canine oral microbiome. The analysis used log-transformed abundance data to focus on the differences between the test and control groups over the 12-week study period.

*   **Script for Execution:** The code for this analysis is located in `R/04_statistical_and_clinical_analysis.R`.

---

## 2. Retrospective Power Analysis for PERMANOVA

### Theoretical Background

Analyzing differences in entire microbial communities is complex. We cannot simply use a t-test. Instead, we first calculate the "distance" or "dissimilarity" between every pair of samples (e.g., using Bray-Curtis dissimilarity). To test if there is a significant difference between experimental groups (e.g., "Control" vs. "Test"), we use a **Permutational Multivariate Analysis of Variance (PERMANOVA)**. PERMANOVA is a powerful non-parametric method that can determine if the distance *between* groups is significantly larger than the distance *within* groups.

After finding a significant result, a crucial follow-up question is: **"Did our study have enough samples to reliably detect this effect?"** This is a question of statistical power. A low-power study might find a non-significant result simply because it didn't have enough samples, not because there was no real effect (a Type II error).

Power estimation for PERMANOVA is not straightforward. We used the method developed by Kelly et al. (2015), implemented in the R package `micropower`. This framework works via simulation and bootstrapping:

1.  It takes our actual distance matrix (e.g., Bray-Curtis).
2.  It repeatedly creates new, smaller distance matrices by randomly sampling (with replacement) from our data.
3.  For each of these new "bootstrap" matrices, it runs a PERMANOVA test.
4.  **Power is calculated as the proportion of these simulated tests that yield a significant p-value (e.g., p < 0.05).**

A power of 0.88, for example, means that if we were to repeat the experiment under the same conditions, we would have an 88% chance of detecting a true effect.

### Effect Size: R² vs. Omega-Squared (ω²)

When reporting the results of an ANOVA-type test, it's important to report the **effect size**—a measure of how much of the total variation is explained by our grouping factor. The coefficient of determination (R²) is a common measure, but it is known to be biased upwards. A less biased and more conservative measure is **Omega-squared (ω²)**, which accounts for mean-squared error. We calculated ω² to provide a more accurate estimate of the effect size in our PERMANOVA tests.

### Implementation in this Project

A retrospective power analysis was performed on our Bray-Curtis dissimilarity matrix to evaluate whether the study design had sufficient statistical power. As shown in the results table in the original report, the power to detect the overall difference between the Control and Test groups at Week 4 was very high (0.99), validating the robustness of our findings.

*   **Script for Execution:** The primary PERMANOVA tests are performed in `R/03_diversity_and_ordination_analysis.R`. The specific retrospective power analysis calculations are performed in `R/04_statistical_and_clinical_analysis.R`.

---

### References

1.  Kelly, B.J., Gross, R., Bittinger, K., et al. (2015). Power and sample-size estimation for microbiome studies using pairwise distances and PERMANOVA. *Bioinformatics*, 31(15), 2461-2468.
2.  Van den Brink P.J., & ter Braak C.J.F. (1999). Principal response curves: analysis of time-dependent multivariate responses of a biological community to stress. *Environmental Toxicology and Chemistry*, 18(2), 138-148.