library(ggplot2)
library(mgcv)
library(DHARMa) 
library(car)     
library(mgcViz)
library(gratia)

# Load data
dataset <- read.delim("clipboard")

# Shapiro-Wilk test for normality
shapiro.test(dataset$canopy_max)
hist(dataset$canopy_max)

gam_model <- gam(canopy_max ~ s(CH, k=4), 
                 family = tw(link = "log"), 
                 data = dataset, method = "REML")
summary(gam_model)

# # # # # DHARMA Diagnostics # # # # # #
png("DHARMa_main_diagnostics_max.png",
    width = 2400, height = 2400, res = 300)

sim_res <- simulateResiduals(fittedModel = gam_model, n = 1000)
plot(sim_res)

dev.off()


# # # # # # DHARMA uniformity # # # # # #
png("DHARMa_uniformity_test_max.png",
    width = 2400, height = 1800, res = 300)

testUniformity(sim_res)

dev.off()


# # # # # # DHARMA dispersion # # # # # #
png("DHARMa_dispersion_test_max.png",
    width = 2400, height = 1800, res = 300)

testDispersion(sim_res)

dev.off()


# Create a new data frame for prediction
newdata <- data.frame(CH = seq(min(dataset$CH), max(dataset$CH), length.out = 200))

# Predict fitted values and standard errors
pred <- predict(gam_model, newdata, type = "link", se.fit = TRUE)

# Add predictions and confidence intervals to the data-frame
newdata$fit <- exp(pred$fit)
newdata$upper <- exp(pred$fit + 1.96 * pred$se.fit)
newdata$lower <- exp(pred$fit - 1.96 * pred$se.fit)


# Plot
p_gam <- ggplot() +
  geom_point(data = dataset, aes(x = CH, y = canopy_max),
             color = "red", size = 2, alpha = 0.7) +
  geom_line(data = newdata, aes(x = CH, y = fit),
            color = "black", linewidth = 1.2) +
  geom_ribbon(data = newdata, aes(x = CH, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.3) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    x = "Cumulative Herbivory Index",
    y = "Canopy cover (%)"
  )

p_gam




######### Derivatives ##########
# Compute derivatives of the smooth using the 'select' argument
deriv <- derivatives(
  gam_model,
  select = "s(CH)",   # instead of term
  type = "central",
  n = 200
)

draw(deriv, add_change=T, change_type = "sizer")

# Check column names
colnames(deriv)

# Create a column indicating where the slope is significantly different from 0
deriv$significant <- (deriv$.lower_ci > 0) | (deriv$.upper_ci < 0)

# Extract CHI values where slope is significant
change_points <- deriv$CH[deriv$significant]

# Add a column for shading the significant slope regions
deriv$shade <- ifelse(deriv$significant, deriv$CH, NA)

# Find the point where derivative is most negative
max_decline_index <- which.min(deriv$.derivative)  # index of most negative slope
max_decline_CHi <- deriv$CH[max_decline_index]     # corresponding CHi
max_decline_slope <- deriv$.derivative[max_decline_index]

cat("Most negative slope occurs at CHi =", max_decline_CHi, 
    "with slope =", max_decline_slope, "\n")

# Plot
p_gam_2 = ggplot() +
  # Highlight significant slope regions
  geom_ribbon(data = deriv, aes(x = CH, ymin = -Inf, ymax = Inf, fill = significant), alpha = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#de2d26", "FALSE" = NA), guide = "none") +
  
  # GAM fit and CI
  geom_ribbon(data = newdata, aes(x = CH, ymin = lower, ymax = upper), 
              fill = "grey70", alpha = 0.3) +
  geom_line(data = newdata, aes(x = CH, y = fit), color = "black", linewidth = 1.2) +
  
  # Actual data points
  geom_point(data = dataset, aes(x = CH, y = canopy_max), color = "red", size = 2, alpha = 0.7) +
  
  coord_cartesian(ylim = c(0, 60)) +
  
  # Labels and theme
  labs(x = "Cumulative Herbivory Index (CH)", 
       y = "Maximum canopy cover (%)",
       title = "") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

p_gam_2

ggsave(
  filename = "GAM_canopy_vs_CH.png",
  plot = p_gam_2,
  width = 6,
  height = 5,
  dpi = 300
)


