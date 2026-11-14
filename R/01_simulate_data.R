# =============================================================================
# 01_simulate_data.R
# Synthetic overall-survival cohort (n = 160) modelled on a typical
# haematological-oncology registry: age group, cytogenetic risk class,
# smoking history and TP53 mutation status. Survival times come from a
# Weibull proportional-hazards model with administrative censoring.
# =============================================================================

set.seed(2024)
RNGkind(sample.kind = "Rejection")
n <- 160

age <- round(pmin(pmax(rnorm(n, 64, 11), 25), 89))
age_group <- cut(age, c(-Inf, 59, 70, Inf), labels = c("< 60", "60-70", "> 70"))
cytogenetic_risk <- factor(sample(c("Good", "Intermediate", "Poor"), n, TRUE, c(0.15, 0.50, 0.35)),
                           levels = c("Good", "Intermediate", "Poor"))
smoking_history <- factor(rbinom(n, 1, 0.40), labels = c("No", "Yes"))
tp53_mutated <- factor(rbinom(n, 1, ifelse(cytogenetic_risk == "Poor", 0.30, 0.08)), labels = c("No", "Yes"))
sex <- factor(sample(c("Female", "Male"), n, TRUE))

# True log hazard ratios
lp <- 0.00 * (age_group == "60-70") + 0.55 * (age_group == "> 70") +
      0.60 * (cytogenetic_risk == "Intermediate") + 1.20 * (cytogenetic_risk == "Poor") +
      0.10 * (smoking_history == "Yes") + 0.80 * (tp53_mutated == "Yes")

shape <- 0.9; scale <- 60                                   # Weibull baseline (months)
t_event <- scale * (-log(runif(n)) / exp(lp))^(1 / shape)
t_cens  <- runif(n, 6, 72)                                  # staggered entry, 6-72 months follow-up
time_months <- round(pmin(t_event, t_cens), 2)
death <- as.integer(t_event <= t_cens)

cohort <- data.frame(patient_id = sprintf("PT%03d", 1:n), age, age_group, sex,
                     cytogenetic_risk, smoking_history, tp53_mutated,
                     survival_time_months = time_months,
                     survival_time_days = round(time_months * 30.44),
                     death)
dir.create("data", showWarnings = FALSE)
write.csv(cohort, "data/survival_cohort.csv", row.names = FALSE)
cat(sprintf("Saved data/survival_cohort.csv: %d patients, %d deaths (%.0f%%)\n",
            n, sum(death), 100 * mean(death)))
