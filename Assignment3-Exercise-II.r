#####################################
# Exercise II
#####################################

# Load libraries
library(ggplot2)
library(dplyr)
library(lmerTest)
library(lattice)

# Read the data
data <- read_csv("allvar.csv")

#---------------------------------------
# i)
#---------------------------------------

#We transform the variable CD4PCT
data$sqrt_CD4PCT <- sqrt(data$CD4PCT)

#Create a new variable time, defined as difference between time of visit and baseage.
data$time <- data$visage-data$baseage

#Define the children and treatment as factors.
data$child_id <- as.factor(data$newpid)

data$treatmnt <- as.factor(data$treatmnt)

# Graph outcome of each child as a function of time.
# They are color coded based on treatment.
ggplot(data, aes(x = time, y = sqrt_CD4PCT, color = treatmnt, group = child_id)) +
  geom_line(alpha = 0.8, size = 1) +  # Line plot for individual trends
  labs(
    title = "Outcome of Each Child as a Function of Time",
    x = "Time (Years)",
    y = "Square Root of CD4 Percentage",
    color = "Treatment"
  ) +
  theme_minimal() +
  theme(legend.position = "right") 



#------------------------------------------------
# ii)
#------------------------------------------------

#Mathematical expression of model:

#sqrt(CD4PCT)_ij=\alpha + \beta*time_ij + a_i + \epsilon_ij

#Fit the model
model <- lmer(sqrt_CD4PCT ~ time + (1 | child_id), data = data)

summary(model)
#beta = - 0.278

#Time variable: Tells us how much time as past since the initial visit.

#beta*Time: Tells us the change in sqrt_CD4PCT per year since the initial visit.

logLik(model)

#------------------------------------------------
# iii)
#------------------------------------------------

#Mathematical expression when including child-level predictors
#sqrt(CD4PCT)_ij=\alpha + \beta_1*time_ij + \beta_2*treatmnt + \beta_3*baseage 
#               + a_i + \epsilon_ij

extended_model <- lmer(sqrt_CD4PCT ~ time + treatmnt + baseage + (1 | child_id), data = data)
summary(extended_model) 


#------------------------------------------------
# iv)
#------------------------------------------------


# we create a new dataset for the hypothetical next time point
# where we add n years to the current macimum time value for each child.
# We omit the children that only have NA times.
# n is the hypothetical next time we predict
n <- 3

predicted_data <- data %>%
  group_by(child_id) %>%
  
  filter(any(!is.na(time))) %>% #Here we filter the NA 
  summarize(
    time = max(time, na.rm = TRUE) + n,  # Hypothetical next time point
    treatmnt = first(treatmnt),          # Treatment status (same as before)
    baseage = first(baseage)             # Baseline age (same as before)
  ) %>%
  ungroup()

# Generate predictions and add as a new column directly
predicted_data$predicted_sqrt_CD4PCT <- predict(extended_model, newdata = predicted_data, re.form = ~(1 | child_id)
)

# Back-transform the square root predictions to the original scale
predicted_data$predicted_CD4PCT <- (predicted_data$predicted_sqrt_CD4PCT)^2

# View the predictions
print(predicted_data)

# Plot prediction colored coded based on treatment.
ggplot(predicted_data, aes(x = child_id, y = predicted_CD4PCT, fill = as.factor(treatmnt))) +
  geom_bar(stat = "identity") +
  labs(
    title = "Predicted CD4 Percentages at Hypothetical Next Time Point",
    x = "Child ID",
    y = "Predicted CD4 Percentage",
    fill = "Treatment"  # Legend title
  ) +
  theme_minimal() +
  scale_x_discrete(
    breaks = as.character(seq(10, 250, by = 5))  # Specify every 5th child ID
  ) +
  theme(
    legend.position = "right",  # Show legend for treatment colors
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),  # Smaller text
    axis.title.x = element_text(size = 12),  # Axis title size
    axis.title.y = element_text(size = 12),  # Y-axis title size
    plot.title = element_text(size = 15)    # Title size
  )


#--------------------------------------------------------
# v)
#--------------------------------------------------------

# We create a new dataset for two new children that have baseage=4.
# We introduce one child per treatment.
new_children <- data.frame(
  child_id = c("Treatment 1", "Treatment 2"),  # Labels with baseage and treatment
  baseage = c(4, 4),   # Base ages
  treatmnt = c(1, 2),    # Treatments
  time = rep(seq(0, 2, by = 0.25), each = 2)  # Time points for predictions (0 to 2 years)
)


#Make sure variables are defined as factors
new_children$treatmnt <- as.factor(new_children$treatmnt)
new_children$child_id <- as.factor(new_children$child_id)

# Predict square root CD4 percentages
new_children$predicted_sqrt_CD4PCT <- predict(
  extended_model,
  newdata = new_children,
  re.form = NA,  # Only use fixed effects
  allow.new.levels = TRUE
)

# Back-transform the square root predictions to the original scale
new_children$predicted_CD4PCT <- (new_children$predicted_sqrt_CD4PCT)^2

# Plot the predictions with meaningful labels
ggplot(new_children, aes(x = time, y = predicted_CD4PCT, color = treatmnt, group = child_id)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Predicted CD4 Percentages for Six Hypothetical Children",
    x = "Time (Years)",
    y = "Predicted CD4 Percentage",
    color = "Treatment"  # Legend title reflecting the labels
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 15)
  )

#It makes sence that treatment 2 has a higher predicted CD4 Percentage
# because the parameter for treatment is positive and (treatmnt1=0, treatmnt2=1)

#-------------------------------------------------------
# vi)
#-------------------------------------------------------

#Fit a model that also includes varying slopes.

#Mathematical expression is:
#sqrt(CD4PCT)_ij=\alpha + \beta_1*time_ij + \beta_2*treatmnt + \beta_3*baseage 
#               + a_i + b_i + \epsilon_ij


extended_model2 <- lmer(sqrt_CD4PCT ~ time + treatmnt + baseage + (time | child_id), data = data)


#Model selection
anova(extended_model2)
#Treatment effect is insignificant.

#We remove it and fit a model again.
final_model <- lmer(sqrt_CD4PCT ~ time + baseage + (time | child_id), data = data)
anova(final_model)
#Both time and baseage effect is significant.

#Test the random effects.
ranova(final_model)
#Random effects are significant.

#this is thus our final model.


#Model Diagnostics
par(mfrow=c(1,2))
plot(final_model)

# Check normality of residuals
qqnorm(residuals(final_model))
qqline(residuals(final_model))

temp<-names(ranef(final_model))
temp
qqnorm(unlist(ranef(final_model)[[1]]),main=paste(temp[1]),cex.main=1)
lines((-3):3,sd(unlist(ranef(final_model)[[1]]))*((-3):3),col="red")


ranef_values <- ranef(final_model)
dotplot(ranef_values)


#-------------------------------------------
# vii)
#-------------------------------------------

summary(final_model)
