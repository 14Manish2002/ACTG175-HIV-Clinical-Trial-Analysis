# Clinical Trial Survival Analysis and Treatment Efficacy Evaluation  
## Using ACTG 175 HIV Clinical Trial Data

**Student:** Manish Pandey  
**Course:** Biostatistics / Clinical Data Analysis  
**Dataset:** AIDS Clinical Trials Group Study 175  
**Tools Used:** R (survival, survminer, tableone, ggplot2, pROC), Python (pandas, scikit-learn, matplotlib)

---

## Table of Contents
1. [Introduction](#1-introduction)
2. [Dataset Description](#2-dataset-description)
3. [Data Preprocessing](#3-data-preprocessing)
4. [Exploratory Data Analysis](#4-exploratory-data-analysis)
5. [Survival Analysis](#5-survival-analysis)
6. [Cox Proportional Hazards Model](#6-cox-proportional-hazards-model)
7. [Machine Learning — Logistic Regression](#7-machine-learning--logistic-regression)
8. [Key Findings & Conclusions](#8-key-findings--conclusions)
9. [Limitations & Future Work](#9-limitations--future-work)
10. [References](#10-references)

---

## 1. Introduction

HIV/AIDS remains one of the most significant public health challenges globally. Identifying which treatment regimens slow disease progression and which clinical factors predict poor outcomes are central questions in HIV clinical research.

The **ACTG 175 study** was a landmark randomised controlled trial comparing four antiretroviral treatment strategies in HIV-positive adults with CD4 counts of 200–500 cells/mm³. This analysis re-examines the ACTG 175 data to:

- Estimate and compare **event-free survival** across treatment groups using Kaplan-Meier curves
- Quantify the independent effect of treatment and clinical variables on disease progression via a **Cox Proportional Hazards model**
- Build a **logistic regression classifier** to predict individual patient outcomes

The primary endpoint (`label`) is a composite of AIDS-defining event or death.

---

## 2. Dataset Description

| Attribute | Detail |
|---|---|
| Total observations | 2,139 patients |
| Variables | 25 |
| Missing values | **None** |
| Duplicate records | **None** |
| Events (label = 1) | 521 (24.4%) |
| Censored (label = 0) | 1,618 (75.6%) |
| Follow-up range | 14 — 1,231 days |

### Variable Dictionary

| Variable | Type | Description |
|---|---|---|
| `time` | Numeric | Days to event or censoring |
| `label` | Binary | 1 = AIDS event/death, 0 = censored |
| `trt` | Integer (0–3) | Treatment arm (4 groups) |
| `treat` | Binary | 0 = ZDV only, 1 = combination therapy |
| `age` | Numeric | Age in years |
| `wtkg` | Numeric | Body weight in kg |
| `cd40` | Numeric | Baseline CD4 count (cells/mm³) |
| `cd420` | Numeric | CD4 count at 20 weeks |
| `cd80` | Numeric | Baseline CD8 count |
| `cd820` | Numeric | CD8 count at 20 weeks |
| `karnof` | Integer | Karnofsky performance score (70–100) |
| `symptom` | Binary | 1 = symptomatic at baseline |
| `gender` | Binary | 0 = female, 1 = male |
| `hemo` | Binary | Hemophilia status |
| `homo` | Binary | Homosexual activity |
| `drugs` | Binary | IV drug use history |
| `race` | Binary | Race indicator |
| `preanti` | Numeric | Days of prior antiretroviral therapy |
| `str2` | Binary | Antiretroviral history |
| `strat` | Integer | CD4 stratum at screening |

---

## 3. Data Preprocessing

```r
# Convert categorical variables to factors
actg$gender  <- as.factor(actg$gender)
actg$homo    <- as.factor(actg$homo)
actg$drugs   <- as.factor(actg$drugs)
actg$hemo    <- as.factor(actg$hemo)
actg$symptom <- as.factor(actg$symptom)
actg$treat   <- as.factor(actg$treat)
```

- No imputation was required — dataset had **zero missing values**
- No records were removed — dataset had **zero duplicates**
- The binary `treat` variable was used for survival analysis (rather than the 4-level `trt`) to answer the core clinical question: does combination therapy outperform monotherapy?

---

## 4. Exploratory Data Analysis

### Figure 1 — EDA Overview

![EDA Overview](plots/01_EDA_overview.png)

**Observations:**
- **Age** ranged from 12 to 70 years; mean = 35.25, approximately right-skewed
- **Weight** had a roughly normal distribution centred around 75 kg
- **Baseline CD4** was right-skewed; median = 340 cells/mm³, consistent with inclusion criteria (200–500)
- **Treatment groups:** 532 (24.9%) patients received ZDV monotherapy; 1,607 (75.1%) received combination therapy — reflecting the trial's predominant enrolment into combination arms
- **Gender:** 82.8% male, 17.2% female — typical of early HIV clinical trials
- **Outcome:** 521 events (24.4%) vs 1,618 censored (75.6%)

### Figure 2 — Baseline Characteristics by Treatment Group

![Baseline by Treatment](plots/02_baseline_by_treatment.png)

**Table 1 — Baseline Characteristics (Stratified by Treatment)**

| Variable | ZDV Only (n=532) | Combination (n=1,607) | p-value |
|---|---|---|---|
| Age (mean ± SD) | 35.23 ± 8.85 | 35.26 ± 8.66 | 0.945 |
| Weight kg (mean ± SD) | 76.06 ± 13.22 | 74.82 ± 13.27 | 0.060 |
| Baseline CD4 (mean ± SD) | 353.20 ± 114.11 | 349.61 ± 120.04 | 0.544 |
| Baseline CD8 (mean ± SD) | 987.25 ± 475.22 | 986.42 ± 481.98 | 0.972 |
| Karnofsky (mean ± SD) | 95.43 ± 5.98 | 95.45 ± 5.88 | 0.949 |

All baseline differences were **non-significant (p > 0.05)**, confirming successful randomisation. This means any observed survival differences can be attributed to treatment rather than baseline confounding.

---

## 5. Survival Analysis

### 5.1 Kaplan-Meier Curves

### Figure 3 — Kaplan-Meier Survival Curves

![Kaplan-Meier](plots/03_kaplan_meier_curves.png)

**Interpretation:**
- Patients on **combination therapy** (blue) maintained a consistently higher event-free survival probability throughout the entire follow-up period
- The separation between curves begins early and widens over time
- Neither group reached median survival (median = NA for both groups), meaning more than 50% of patients in both groups remained event-free at study end
- The **number at risk** decreases over time as expected, with more patients at risk in the larger combination therapy group

### 5.2 Log-Rank Test

```
Log-rank Chi-squared = 47.7,  df = 1,  p = 5 × 10⁻¹²
```

**Interpretation:** The log-rank test is strongly significant (p < 0.001). The null hypothesis (no difference in survival between groups) is decisively rejected. Patients receiving combination antiretroviral therapy had significantly better survival compared to ZDV monotherapy.

### 5.3 CD4 and CD8 Changes

### Figure 5 — CD4/CD8 Count Change from Baseline to Week 20

![CD4 CD8 Change](plots/05_cd4_cd8_change.png)

**Observations:**
- The combination therapy group showed a **larger increase in CD4 count** from baseline to week 20, reflecting better immune recovery
- CD8 counts remained relatively stable in both groups
- This immunological improvement is consistent with the observed survival benefit

### 5.4 Event Rates Across Subgroups

### Figure 6 — Disease Progression by Subgroup

![Subgroups](plots/06_event_rates_subgroups.png)

**Key observations:**
- Symptomatic patients had a notably higher event rate (~42%) vs asymptomatic patients (~19%)
- IV drug users had slightly higher event rates
- The treatment group difference is the most pronounced subgroup effect

---

## 6. Cox Proportional Hazards Model

### 6.1 Model Specification

```r
cox_model <- coxph(
  Surv(time, label) ~ age + wtkg + cd40 + cd80 + karnof + treat,
  data = actg
)
```

### 6.2 Results

### Figure 4 — Forest Plot (Hazard Ratios)

![Forest Plot](plots/04_cox_forest_plot.png)

**Full Model Results:**

| Predictor | Coefficient | HR | 95% CI | p-value | Significance |
|---|---|---|---|---|---|
| Age | 0.00789 | 1.008 | [0.998, 1.018] | 0.121 | NS |
| Weight (kg) | 0.00067 | 1.001 | [0.994, 1.007] | 0.844 | NS |
| **Baseline CD4** | **−0.00445** | **0.996** | **[0.995, 0.996]** | **< 2×10⁻¹⁶** | *** |
| **Baseline CD8** | **0.00046** | **1.000** | **[1.000, 1.001]** | **2.3×10⁻⁸** | *** |
| **Karnofsky Score** | **−0.02729** | **0.973** | **[0.960, 0.986]** | **6.4×10⁻⁵** | *** |
| **Combination Therapy** | **−0.68059** | **0.506** | **[0.422, 0.607]** | **1.9×10⁻¹³** | *** |

**Model Fit:**
- Concordance Index = **0.682** (SE = 0.012) — acceptable discriminative ability
- Likelihood Ratio Test: 192.2 on 6 df, p < 2.2×10⁻¹⁶

### 6.3 Interpretation of Hazard Ratios

**Combination Therapy (HR = 0.506):**
This is the most clinically important finding. Patients receiving combination antiretroviral therapy have a **49.4% lower hazard** of disease progression or death compared to ZDV monotherapy, after adjusting for age, weight, CD4, CD8, and Karnofsky score. This is a large, highly significant protective effect.

**Baseline CD4 Count (HR = 0.996 per cell/mm³):**
Higher CD4 at baseline is protective. Each additional 100 cells/mm³ in CD4 count reduces the hazard by approximately 33%. This makes biological sense — a stronger immune system at study entry is associated with better outcomes.

**Baseline CD8 Count (HR = 1.0005 per cell/mm³):**
Slightly elevated hazard with higher CD8 count. This is somewhat counterintuitive, but higher CD8 may reflect a state of chronic immune activation in more advanced HIV disease. The per-unit effect is small but statistically significant due to the large sample.

**Karnofsky Performance Score (HR = 0.973 per unit):**
Better functional status is protective. Each 10-point improvement in Karnofsky score corresponds to approximately a 24% reduction in hazard. Patients who are more functionally capable at baseline have better survival outcomes.

**Age and Weight:** Neither reached statistical significance (p > 0.05), suggesting their prognostic value is captured by the other clinical markers after adjustment.

### 6.4 Proportional Hazards Assumption Check

```
Schoenfeld Residuals Test:

Variable     Chi-sq   df    p-value    Verdict
Age           1.131    1    0.288      ✅ Satisfied
Weight        0.458    1    0.499      ✅ Satisfied  
Baseline CD4 13.143    1    0.0003     ❌ Violated
Baseline CD8  1.606    1    0.205      ✅ Satisfied
Karnofsky     0.566    1    0.452      ✅ Satisfied
Treatment     7.374    1    0.007      ❌ Violated
GLOBAL       28.335    6    8.1×10⁻⁵   ❌ Violated
```

**Interpretation:** The proportional hazards assumption is violated for `cd40` and `treat`. This means their effect on survival is **not constant over time** — the treatment benefit may be strongest in early follow-up and attenuate later. A time-varying Cox model would address this limitation.

---

## 7. Machine Learning — Logistic Regression

### 7.1 Setup

```r
# 80/20 train-test split, seed = 123
# Features: age, wtkg, cd40, cd80, karnof, symptom, treat
# Outcome: event (binary: 0 = no event, 1 = event)
```

### 7.2 Model Coefficients

```
(Intercept)  : 2.332   (p = 0.024)
age          : 0.012   (p = 0.082)
wtkg         : 0.001   (p = 0.841)  ← not significant
cd40         : -0.004  (p < 0.001) ***
cd80         : 0.0005  (p < 0.001) ***
karnof       : -0.027  (p = 0.005) **
symptom1     : 0.471   (p = 0.001) **
treat1       : -0.722  (p < 0.001) ***
```

### 7.3 Model Performance

### Figure 7 — ROC Curve and Confusion Matrix

![ROC and CM](plots/07_roc_confusion_matrix.png)

**Confusion Matrix (Test Set, n = 428):**

```
               Predicted 0   Predicted 1
Actual 0           321            8
Actual 1            84           15
```

**Performance Metrics:**

| Metric | Value | Interpretation |
|---|---|---|
| **Accuracy** | **78.5%** | Correctly classifies 78.5% of patients |
| **AUC-ROC** | **0.731** | Acceptable discrimination ability |
| Sensitivity | 15.2% | Model misses many true events |
| Specificity | 97.6% | Rarely misclassifies non-events |
| PPV | 65.2% | 65% of predicted events are correct |
| NPV | 79.3% | 79% of predicted non-events are correct |

**Note on class imbalance:** Only 24.4% of patients experienced the event. The model is biased toward predicting the majority class (no event), resulting in high specificity but low sensitivity. Techniques like SMOTE, class weighting, or threshold tuning could improve sensitivity at the cost of specificity.

### Figure 8 — Feature Importance

![Feature Importance](plots/08_feature_importance.png)

Consistent with the Cox model, **treatment group** and **baseline CD4** are the strongest predictors in the logistic model. Symptomatic status, Karnofsky score, and CD8 count also contribute significantly.

---

## 8. Key Findings & Conclusions

### Figure 9 — Correlation Heatmap

![Correlation Heatmap](plots/09_correlation_heatmap.png)

Notable correlations:
- `cd40` and `cd420` are moderately correlated (r ≈ 0.50) — baseline predicts week-20 CD4
- `cd80` and `cd820` similarly correlated
- `label` (event) negatively correlates with `cd40` (r ≈ −0.18) — higher CD4 → fewer events
- `treat` shows modest negative correlation with `label` — combination therapy reduces events

### Figure 10 — Summary Dashboard

![Summary Dashboard](plots/10_summary_dashboard.png)

### Summary of All Key Results

| Analysis | Key Finding |
|---|---|
| **Randomisation check** | Groups well-balanced at baseline (all p > 0.05) ✅ |
| **Kaplan-Meier** | Combination therapy has higher event-free survival throughout follow-up |
| **Log-rank test** | p = 5×10⁻¹² — highly significant survival difference |
| **Cox HR (treatment)** | HR = 0.506 — combination therapy reduces hazard by **49.4%** |
| **Cox HR (CD4 baseline)** | HR = 0.996 — higher CD4 is strongly protective (p < 2×10⁻¹⁶) |
| **Cox HR (Karnofsky)** | HR = 0.973 — better function is protective (p < 0.001) |
| **PH assumption** | Violated for CD4 and treatment — time-varying model recommended |
| **Logistic regression AUC** | 0.731 — acceptable model discrimination |
| **Overall conclusion** | **Combination antiretroviral therapy significantly reduces HIV disease progression** |

---

## 9. Limitations & Future Work

### Limitations

1. **Proportional Hazards Violation** — The Cox model assumes constant hazard ratios over time, but this was violated for `cd40` and `treat`. The reported HRs represent average effects and may not capture time-varying dynamics.

2. **Class Imbalance** — Only 24.4% of patients experienced the event. This leads to low sensitivity in the logistic regression model, which may miss true positives in clinical practice.

3. **Binary Treatment Grouping** — Collapsing four distinct treatment arms into binary `treat` simplifies interpretation but loses granularity. Sub-analyses comparing all four arms (ddI alone, ddC+ZDV, ddI+ZDV vs ZDV monotherapy) would add insight.

4. **Single Institution / Era** — The trial was conducted in the mid-1990s; modern antiretroviral regimens are significantly different. Findings reflect the historical treatment landscape.

5. **No External Validation** — The logistic regression model was evaluated on a held-out split of the same dataset, not an independent external cohort.

### Suggested Future Work

- [ ] **Time-varying Cox model** to address PH violations for CD4 and treatment
- [ ] **Four-arm survival comparison** using `trt` (0, 1, 2, 3) instead of binary `treat`
- [ ] **Random Forest / XGBoost** for improved predictive accuracy
- [ ] **SMOTE or class-weighted logistic regression** to address class imbalance
- [ ] **Competing risks analysis** — separate cause-specific hazards (AIDS event vs death)
- [ ] **Stratified Cox model** stratified by `strat` (CD4 stratum at screening)

---

## 10. References

1. Hammer SM, Katzenstein DA, Hughes MD, et al. (1996). A trial comparing nucleoside monotherapy with combination therapy in HIV-infected adults with CD4 cell counts from 200 to 500 per cubic millimeter. *New England Journal of Medicine*, 335(15), 1081–1090.

2. Kaplan EL, Meier P. (1958). Nonparametric estimation from incomplete observations. *Journal of the American Statistical Association*, 53(282), 457–481.

3. Cox DR. (1972). Regression models and life-tables. *Journal of the Royal Statistical Society: Series B*, 34(2), 187–202.

4. Therneau TM, Grambsch PM. (2000). *Modeling Survival Data: Extending the Cox Model*. Springer, New York.

5. UCI Machine Learning Repository — ACTG 175 Dataset. https://archive.ics.uci.edu/ml/datasets/aids+clinical+trials+group+study+175

---

*Analysis performed using R 4.4.x and Python 3.x*  
*Packages: survival, survminer, tableone, ggplot2, pROC (R); pandas, scikit-learn, matplotlib, seaborn (Python)*
