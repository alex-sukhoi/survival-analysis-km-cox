# =============================================================================
# 02_survival_analysis.R
#   1. Descriptive table: alive vs died, with chi-square / Fisher tests
#   2. Kaplan-Meier curves (overall, by age group, by cytogenetic risk)
#      with numbers at risk, median survival with 95% CI, log-rank tests
#   3. Cox proportional hazards: univariable and multivariable HRs,
#      forest plot, proportional-hazards check (Schoenfeld residuals),
#      concordance
# =============================================================================

suppressPackageStartupMessages({ library(survival); library(ggplot2); library(gridExtra) })
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

d <- read.csv("data/survival_cohort.csv", stringsAsFactors = TRUE)
d$age_group <- factor(d$age_group, levels = c("< 60", "60-70", "> 70"))
d$cytogenetic_risk <- factor(d$cytogenetic_risk, levels = c("Good", "Intermediate", "Poor"))
S <- with(d, Surv(survival_time_months, death))
pal <- c("#2b6cb0", "#dd6b20", "#c53030", "#319795")

# --- 1. Descriptives -----------------------------------------------------------
vars <- c("age_group", "cytogenetic_risk", "smoking_history", "tp53_mutated")
desc <- do.call(rbind, lapply(vars, function(v) {
  tab <- table(d[[v]], factor(d$death, labels = c("Alive", "Died")))
  test <- if (any(chisq.test(tab)$expected < 5)) fisher.test(tab) else chisq.test(tab)
  data.frame(variable = v, level = rownames(tab),
             total = sprintf("%d (%.1f%%)", rowSums(tab), 100 * rowSums(tab) / sum(tab)),
             alive = sprintf("%d (%.1f%%)", tab[, 1], 100 * tab[, 1] / sum(tab[, 1])),
             died  = sprintf("%d (%.1f%%)", tab[, 2], 100 * tab[, 2] / sum(tab[, 2])),
             test = c(sprintf("%s p = %.3f", ifelse(inherits(test, "htest") && grepl("Fisher", test$method),
                                                    "Fisher", "Chi-square"), test$p.value),
                      rep("", nrow(tab) - 1)))
}))
write.csv(desc, "results/01_descriptives_by_outcome.csv", row.names = FALSE)

# --- 2. Kaplan-Meier -----------------------------------------------------------
km_table <- function(fit, label) {
  s <- summary(fit)$table
  if (is.null(dim(s))) s <- t(as.matrix(s))
  data.frame(group = if (is.null(rownames(s))) "All patients" else sub(".*=", "", rownames(s)),
             n = s[, "records"], deaths = s[, "events"],
             median_months = round(s[, "median"], 1),
             ci_95 = sprintf("%s-%s", round(s[, "0.95LCL"], 1), round(s[, "0.95UCL"], 1)),
             analysis = label, row.names = NULL)
}
fit_all  <- survfit(S ~ 1, data = d)
fit_age  <- survfit(S ~ age_group, data = d)
fit_risk <- survfit(S ~ cytogenetic_risk, data = d)
lr_age  <- survdiff(S ~ age_group, data = d)
lr_risk <- survdiff(S ~ cytogenetic_risk, data = d)
p_lr <- function(x) pchisq(x$chisq, length(x$n) - 1, lower.tail = FALSE)
medians <- rbind(km_table(fit_all, "Overall"), km_table(fit_age, "Age group"),
                 km_table(fit_risk, "Cytogenetic risk"))
medians$log_rank_p <- c(NA, rep(signif(p_lr(lr_age), 3), 3), rep(signif(p_lr(lr_risk), 3), 3))
write.csv(medians, "results/02_km_median_survival.csv", row.names = FALSE)

times <- seq(0, 60, 12)
surv_at <- summary(fit_all, times = c(12, 24, 36))
write.csv(data.frame(month = surv_at$time, survival = round(surv_at$surv, 3),
                     ci_low = round(surv_at$lower, 3), ci_high = round(surv_at$upper, 3)),
          "results/03_overall_survival_at_1_2_3_years.csv", row.names = FALSE)

km_plot <- function(fit, title, subtitle, file) {
  s <- summary(fit, censored = TRUE)
  strata <- if (is.null(fit$strata)) rep("All patients", length(fit$time)) else
    rep(sub(".*=", "", names(fit$strata)), fit$strata)
  df <- data.frame(time = fit$time, surv = fit$surv, lower = fit$lower, upper = fit$upper,
                   n.censor = fit$n.censor, strata = strata)
  start <- data.frame(time = 0, surv = 1, lower = 1, upper = 1, n.censor = 0, strata = unique(strata))
  df <- rbind(start, df); df$strata <- factor(df$strata, levels = unique(strata))
  g <- ggplot(df, aes(time, surv, colour = strata, fill = strata))
  if (nlevels(df$strata) == 1)   # confidence band only for the single-curve plot
    g <- g + geom_ribbon(aes(ymin = lower, ymax = upper), stat = "identity", alpha = 0.12,
                         colour = NA, outline.type = "full")
  g <- g + geom_step(linewidth = 0.9) +
    geom_point(data = subset(df, n.censor > 0), shape = 3, size = 2) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50") +
    scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
    scale_x_continuous(breaks = times, limits = c(0, 60)) + scale_y_continuous(labels = scales::percent) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Overall survival", colour = NULL, fill = NULL) +
    theme_minimal(base_size = 12) + theme(legend.position = "top")
  risk <- summary(fit, times = times, extend = TRUE)
  rdf <- data.frame(time = risk$time, n = risk$n.risk,
                    strata = if (is.null(risk$strata)) "All patients" else sub(".*=", "", risk$strata))
  rdf$strata <- factor(rdf$strata, levels = rev(levels(df$strata)))
  r <- ggplot(rdf, aes(time, strata, label = n, colour = strata)) + geom_text(size = 3.6) +
    scale_colour_manual(values = rev(pal[seq_along(levels(df$strata))]), guide = "none") +
    scale_x_continuous(breaks = times, limits = c(0, 60)) +
    labs(x = "Months since diagnosis", y = NULL, subtitle = "Number at risk") +
    theme_minimal(base_size = 11) + theme(panel.grid = element_blank())
  ggsave(file, arrangeGrob(g, r, heights = c(3.2, 0.6 + 0.22 * nlevels(df$strata))),
         width = 7.5, height = 5.4, dpi = 110)
}
km_plot(fit_all, "Overall survival, all patients",
        sprintf("Median %.1f months (95%% CI %.1f-%.1f)", summary(fit_all)$table["median"],
                summary(fit_all)$table["0.95LCL"], summary(fit_all)$table["0.95UCL"]),
        "figures/km_overall.png")
km_plot(fit_age, "Overall survival by age group", sprintf("Log-rank p = %.3f", p_lr(lr_age)),
        "figures/km_age_group.png")
km_plot(fit_risk, "Overall survival by cytogenetic risk class", ifelse(p_lr(lr_risk) < 1e-4, "Log-rank p < 0.0001", sprintf("Log-rank p = %.4f", p_lr(lr_risk))),
        "figures/km_cytogenetic_risk.png")

# --- 3. Cox proportional hazards ---------------------------------------------
hr_table <- function(fit, label) {
  s <- summary(fit)
  data.frame(model = label, term = rownames(s$coefficients),
             beta = round(s$coefficients[, "coef"], 3), se = round(s$coefficients[, "se(coef)"], 3),
             HR = round(s$conf.int[, "exp(coef)"], 2), ci_low = round(s$conf.int[, "lower .95"], 2),
             ci_high = round(s$conf.int[, "upper .95"], 2), p = signif(s$coefficients[, "Pr(>|z|)"], 3),
             row.names = NULL)
}
uni <- do.call(rbind, lapply(c("age_group", "cytogenetic_risk", "smoking_history", "tp53_mutated"), function(v)
  hr_table(coxph(as.formula(paste("S ~", v)), data = d), "Univariable")))
cox <- coxph(S ~ age_group + cytogenetic_risk + smoking_history + tp53_mutated, data = d)
multi <- hr_table(cox, "Multivariable")
write.csv(rbind(uni, multi), "results/04_cox_hazard_ratios.csv", row.names = FALSE)

zph <- cox.zph(cox)
ph <- data.frame(term = rownames(zph$table), chisq = round(zph$table[, "chisq"], 2),
                 df = zph$table[, "df"], p = round(zph$table[, "p"], 3), row.names = NULL)
write.csv(ph, "results/05_proportional_hazards_test.csv", row.names = FALSE)
model_stats <- data.frame(statistic = c("Concordance (C-index)", "Likelihood-ratio test p", "Events", "N"),
                          value = c(sprintf("%.3f", cox$concordance["concordance"]),
                                    format.pval(summary(cox)$logtest["pvalue"], digits = 3),
                                    cox$nevent, cox$n))
write.csv(model_stats, "results/06_cox_model_fit.csv", row.names = FALSE)

labels <- c("age_group60-70" = "Age 60-70 (vs < 60)", "age_group> 70" = "Age > 70 (vs < 60)",
            "cytogenetic_riskIntermediate" = "Intermediate risk (vs good)",
            "cytogenetic_riskPoor" = "Poor risk (vs good)",
            "smoking_historyYes" = "Smoking history", "tp53_mutatedYes" = "TP53 mutated")
fp <- multi; fp$label <- factor(labels[fp$term], levels = rev(labels))
gf <- ggplot(fp, aes(HR, label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, colour = "#2b6cb0") +
  geom_point(size = 3, colour = "#2b6cb0") +
  geom_text(aes(x = max(ci_high) * 1.6, label = sprintf("%.2f (%.2f-%.2f)  %s", HR, ci_low, ci_high, ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))),
            hjust = 0, size = 3.5) +
  scale_x_log10(limits = c(min(fp$ci_low) * 0.8, max(fp$ci_high) * 9)) +
  labs(title = "Multivariable Cox model: hazard ratios",
       subtitle = sprintf("C-index = %.2f; HR > 1 means higher risk of death", cox$concordance["concordance"]),
       x = "Hazard ratio (log scale, 95% CI)", y = NULL) +
  theme_minimal(base_size = 12)
ggsave("figures/cox_forest_plot.png", gf, width = 8.5, height = 4, dpi = 110)

cat("\nMedian survival:\n"); print(medians, row.names = FALSE)
cat("\nSurvival at 1/2/3 years:\n"); print(read.csv("results/03_overall_survival_at_1_2_3_years.csv"))
cat("\nCox (multivariable):\n"); print(multi, row.names = FALSE)
cat("\nPH test:\n"); print(ph, row.names = FALSE)
cat("\nModel:\n"); print(model_stats, row.names = FALSE)
