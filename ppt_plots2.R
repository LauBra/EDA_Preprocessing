library(tidyverse)
library(bestNormalize)   # for BoxCox & YeoJohnson transforms

set.seed(123)

# 1. Simulate CRP (very skewed, realistic clinical distribution)
crp <- rlnorm(500, meanlog = log(5), sdlog = 1)

# 2. Simulate an inflammation outcome that depends on log(CRP)
infl_score <- 50 + 12 * log(crp) + rnorm(500, sd = 5)

# 3. Compute transformations
log_crp  <- log(crp)
sqrt_crp <- sqrt(crp)

# Box–Cox (requires > 0)
bc_obj <- boxcox(crp)
boxcox_crp <- predict(bc_obj)

# Yeo–Johnson (allows zeros/negatives)
yj_obj <- yeojohnson(crp)
yeojohnson_crp <- predict(yj_obj)

# 4. Create long dataset for faceted plotting
df <- tibble(
  Raw = crp,
  `Log(CRP)` = log_crp,
  `Sqrt(CRP)` = sqrt_crp,
  `BoxCox(CRP)` = boxcox_crp,
  `YeoJohnson(CRP)` = yeojohnson_crp,
  inflammation = infl_score
) %>%
  pivot_longer(
    cols = -inflammation,
    names_to = "Transform",
    values_to = "CRP_trans"
  )

# 5. Plot relationships with smoothing curve
aa <- ggplot(df, aes(x = CRP_trans, y = inflammation)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, colour = "red", linewidth = 1) +
  facet_wrap(~ Transform, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "CRP (transformed versions)",
    y = "Inflammation score",
    title = "Effect of CRP transformations on the relationship with outcome"
  )


# 5. Plot relationships with smoothing curve
bb <- ggplot(df, aes(x = CRP_trans)) +
  geom_histogram(bins = 30) + #fill="grey80", color="white") +
  facet_wrap(~ Transform, scales = "free") +
  theme_minimal() +
  labs(
    x = "CRP (transformed versions)",
    y = "Inflammation score",
    title = "Effect of CRP transformations on the relationship with outcome"
  )


pdf("One.pdf", 8, 4)
print(aa)
dev.off()

pdf("Two.pdf", 8, 4)
print(bb)
dev.off()


crp_z      <- as.numeric(scale(crp))                                # standardisation
crp_minmax <- (crp - min(crp)) / (max(crp) - min(crp))              # min–max
crp_robust <- (crp - median(crp)) / IQR(crp)                        # robust (median/IQR)

df_scale <- tibble(
  `Raw CRP`          = crp,
  `Z-score CRP`      = crp_z,
  `Min–max CRP`      = crp_minmax,
  `Robust-scaled CRP`= crp_robust,
  inflammation       = infl_score
) %>%
  pivot_longer(
    cols = -inflammation,
    names_to = "Scale",
    values_to = "CRP_scaled"
  )

# Plot: scalings vs outcome
cc <- ggplot(df_scale, aes(x = CRP_scaled, y = inflammation)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, colour = "red", linewidth = 1) +
  facet_wrap(~ Scale, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "CRP (scaled versions)",
    y = "Inflammation score",
    title = "CRP scalings and relationship with outcome"
  )




dd <- ggplot(df_scale, aes(x = CRP_scaled)) +
  geom_histogram(bins = 30) + #fill="grey80", color="white") +
  facet_wrap(~ Scale, scales = "free") +
  theme_minimal() +
  labs(
    x = "CRP (transformed versions)",
    y = "Inflammation score",
    title = "Effect of CRP transformations on the relationship with outcome"
  )


pdf("Third.pdf", 8, 4)
print(cc)
dev.off()

pdf("Fourth.pdf", 8, 4)
print(dd)
dev.off()
