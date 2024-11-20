#####################################
# Exercise II
#####################################

# Load libraries
library(ggplot2)
library(dplyr)
library(lme4)

# Read the data
data <- read_csv("allvar.csv")

#------------------------------------------------------------------------------
# i) 
#------------------------------------------------------------------------------

#We transform the variable CD4PCT
data$sqrt_CD4PCT <- sqrt(data$CD4PCT)

#Create a new variable time, defined as difference between time of visit and baseage.
data$time <- data$visage-data$baseage

#Define Child_id as a factor.
data$child_id <- as.factor(data$newpid)

# Plot outcome of each child as a function of time.
ggplot(data, aes(x = time, y = sqrt_CD4PCT, color = child_id, group = child_id)) +
  geom_line(alpha = 0.8, size = 1) +  # Line plot for individual trends
  labs(
    title = "Outcome of Each Child as a Function of Time",
    x = "Time (Years)",
    y = "Square Root of CD4 Percentage",
    color = "Child ID"
  ) +
  theme_minimal() +
  theme(legend.position = "none") 


#-------------------------------------------------------------------------------
#  ii)
#-------------------------------------------------------------------------------

#Fit the model in R.
model <- lmer(sqrt_CD4PCT ~ time + (1 | child_id), data = data)

# Summarize the model
summary(model)

#The coefficient time = -0.36609 suggests an average annual decline in CD4PCT 
#of -0.36609.


#------------------------------------------------------------------------------
# iii)
#------------------------------------------------------------------------------

#Fit the extended model
extended_model <- lmer(sqrt_CD4PCT ~ time + treatmnt + baseage + (1 | child_id), data = data)
summary(extended_model) 


#-------------------------------------------------------------------------------
# iv) Predict CD4 percentages for each child at a hypothetical next time point
#-------------------------------------------------------------------------------

#n is the hypothetical next time
n <- 1

# Create a new dataset for the hypothetical next time point
# Add one year to the current maximum 'time' value for each child
#We omit the children that only have NA times, as such a value wouldn't make sense to predict.
next_time_data <- data %>%
  group_by(child_id) %>%
  # Filter to remove groups where all 'time' values are NA
  filter(any(!is.na(time))) %>%
  summarize(
    time = max(time, na.rm = TRUE) + n,  # Hypothetical next time point
    treatmnt = first(treatmnt),          # Treatment status remains the same
    baseage = first(baseage)             # Baseline age remains the same
  ) %>%
  ungroup()

# Use the `predict` function with the extended model to generate predictions
next_time_data <- next_time_data %>%
  mutate(
    predicted_sqrt_CD4PCT = predict(extended_model, newdata = ., re.form = ~(1 | child_id))
  )

# Back-transform the square root predictions to the original scale
next_time_data <- next_time_data %>%
  mutate(predicted_CD4PCT = (predicted_sqrt_CD4PCT)^2)

# View the predictions
print(next_time_data)

# Optionally, visualize the predictions
ggplot(next_time_data, aes(x = child_id, y = predicted_CD4PCT, fill = as.factor(treatmnt))) +
  geom_bar(stat = "identity") +
  labs(
    title = "Predicted CD4 Percentages at Hypothetical Next Time Point",
    x = "Child ID",
    y = "Predicted CD4 Percentage",
    fill = "Treatment"  # Legend title
  ) +
  theme_minimal() +
  scale_x_discrete(
    breaks = as.character(seq(5, 250, by = 5))  # Specify every 5th child ID
  ) +
  theme(
    legend.position = "right",  # Show legend for treatment colors
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),  # Smaller text
    axis.title.x = element_text(size = 12),  # Axis title size
    axis.title.y = element_text(size = 12),  # Y-axis title size
    plot.title = element_text(size = 15)    # Title size
  )


#---------------------------------------------------------
# v)
#---------------------------------------------------------


# Create a dataset for six children with meaningful labels
six_children <- data.frame(
  child_id = c("1_T1", "1_T2", "5_T1", "5_T2", "10_T1", "10_T2"),  # Labels with baseage and treatment
  baseage = c(1, 1, 5, 5, 10, 10),   # Base ages
  treatmnt = c(1, 2, 1, 2, 1, 2),    # Treatments
  time = rep(seq(0, 2, by = 0.5), each = 6)  # Time points for predictions (0 to 2 years)
)

# Expand the dataset so each child has predictions at all time points
six_children <- six_children %>%
  group_by(child_id, baseage, treatmnt) %>%
  summarize(time = seq(0, 2, by = 0.5), .groups = "drop") %>%
  ungroup()

# Predict CD4 percentages using the model
six_children <- six_children %>%
  mutate(
    predicted_sqrt_CD4PCT = predict(
      extended_model,
      newdata = six_children,
      re.form = NA,  # Only use fixed effects
      allow.new.levels = TRUE
    ),
    predicted_CD4PCT = (predicted_sqrt_CD4PCT)^2  # Back-transform to the original scale
  )

# Plot the predictions with meaningful labels
ggplot(six_children, aes(x = time, y = predicted_CD4PCT, color = child_id, group = child_id)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Predicted CD4 Percentages for Six Hypothetical Children",
    x = "Time (Years)",
    y = "Predicted CD4 Percentage",
    color = "Baseage_Treatment"  # Legend title reflecting the labels
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 15)
  )



#-------------------------------------------------------
# vi)
#-------------------------------------------------------

# Fit the final model
extended_model2 <- lmer(sqrt_CD4PCT ~ time + treatmnt + baseage + (time | child_id), data = data)

anova(extended_model, extended_model2)  # Likelihood ratio test
anova(model, extended_model2)

# Summarize the model
summary(extended_model2)

#We perform model selection. 
#Therefore we start by fitting a fixed effect model.
# Fit the fixed-effects model
fixed_effect_model <- lm(sqrt_CD4PCT ~ time + treatmnt + baseage, data = data)
anova(fixed_effect_model)
#All fixed effects are significant
ranova(extended_model2)
#Random effect is significant.

#Conclude that we keep the extended_model2.


#Model Diagnostics:

par(mfrow=c(1,2))
plot(extended_model2)


# Check normality of residuals
qqnorm(residuals(extended_model2))
qqline(residuals(extended_model2))


ranef_values <- ranef(extended_model2)
dotplot(extended_model2)


temp<-names(ranef(extended_model2))
temp
qqnorm(unlist(ranef(extended_model2)[[1]]),main=paste(temp[1]),cex.main=1)
lines((-3):3,sd(unlist(ranef(extended_model2)[[1]]))*((-3):3),col="red")

#---------------------------------------------------------------
# vii)
#---------------------------------------------------------------

summary(extended_model2)

VarCorr(extended_model2)




