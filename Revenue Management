# Projects-Business-analytics
# Introduction

The main purpose is to illustrate a succinct analysis of how **Ordinary Least Squares (OLS)**, **Ridge Regression (RR)**, **Slab Regression (SR)**, **Stein (St)**, **Diagonal Shrinkage (DSh)**, and **Shrinkage (Sh)** *Multiple Linear Regression (MLR)* estimators could be applied for a real-life dataset

*Data* presented the composite index S&P GSCI (Standard & Poor's Goldman Sachs Commodity Index) from 03/10/2005 to 09/05/2024. S&P GSCI is a widely recognised measure of commodity market performance, tracking a basket of 24 commodities.

\newpage

# Section 1: Data Cleaning

**Data Preparation**

The initial phase of our analysis necessitates a comprehensive preparation of the dataset to ensure uniformity and accessibility of the data. The process begins with loading all necessary R packages that underpin the computational aspects of our analysis.

**Load packages**

```{r load_packages, results='hide', message=FALSE, warning=FALSE}
options(width = 80, results='hide', message=FALSE)
# Clear the environment
rm(list = ls())

# Comprehensive list of all necessary packages
list_pack <- c("MASS", "base", "clusterGeneration", "glmnet", "matrixcalc", "lubridate",
               "riskParityPortfolio", "stats", "lrmest", "caret", "tictoc", "breakfast",
               "quantmod", "PerformanceAnalytics", "dygraphs", "htmlwidgets", "patchwork",
               "webshot", "xts", "readr", "rportfolio", "quarks", "ggplot2", "tidyr",
               "scales", "sandwich", "strucchange", "dplyr", "zoo", "tseries", "reshape2", 
               "mnormt", "expm", "maotai")

# Load all listed packages using lapply, checking for installation
lapply(list_pack, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

# Set options to discourage scientific notation in output
options(scipen = 999)
```

\newpage

**Data Cleaning**

```{r data_cleaning}
options(width = 80)
# Read the CSV file
data <- read_csv("real_data_commodity.csv")
data$Dates <- as.Date(data$Dates, format = "%d/%m/%Y")

# Rename columns to replace spaces with underscores
data <- rename_with(data, ~ gsub(" ", "_", .))

# Check for initial missing values
initial_na <- data %>%
  dplyr::select(SPGSCI_Index, ends_with("Comdty")) %>%
  filter(if_any(everything(), is.na))

if (nrow(initial_na) > 0) {
  print("Initial missing values found:")
  print(initial_na)
} else {
  print("No initial missing values found.")
}

# Identify columns with zero or negative values and filter out
zero_neg_values <- data %>%
  filter(SPGSCI_Index <= 0 | if_any(ends_with("Comdty"), ~ . <= 0))

if (nrow(zero_neg_values) > 0) {
  print("Rows with zero or negative values found:")
  print(zero_neg_values)
}
data_filtered <- data %>%
  filter(SPGSCI_Index > 0) %>%
  filter(if_all(ends_with("Comdty"), ~ . > 0))

# Calculate Log Returns for SP GSCI and each commodity
data_log_returns <- data_filtered %>%
  mutate(Log_Return_SPGSCI = log(SPGSCI_Index / lag(SPGSCI_Index))) %>%
  mutate(across(ends_with("Comdty"), ~ log(. / lag(.)), 
                .names = "Log_Return_{col}")) %>%
  dplyr::select(Dates, Log_Return_SPGSCI, starts_with("Log_Return_"))

# Check for NA rows in the calculated log returns
nan_rows <- data_log_returns %>%
  filter(if_any(starts_with("Log_Return_"), ~ is.na(.)))

if (nrow(nan_rows) > 0) {
  print("Rows with NA values after log return calculation:")
  print(nan_rows)
} else {
  print("No NA values after log return calculation.")
}

data_log_returns <- na.omit(data_log_returns)

# Save the processed data to a new CSV file
write_csv(data_log_returns, "Processed_Data_With_Log_Returns_SPGSCI.csv")
```

\newpage

# Section 2: Change Point Detection for SPGSCI Index

**Change Points Detection**: this section might take around 10 mins to run through the detection

```{r change_points_detection, warning=FALSE}
options(width = 80)
# Convert log returns to an xts object for time series manipulation
returns_spgsci <- xts(x=data_log_returns$Log_Return_SPGSCI, order.by=data_log_returns$Dates)

# Calculate daily portfolio returns without a wealth index
returns_spgsci_xts <- Return.portfolio(R = returns_spgsci, wealth.index = FALSE, rebalance_on = "days")

# Calculate mean returns across days
returns_spgsci_non <- rowMeans(returns_spgsci_xts)

# Calculate wealth index returns for breakpoint analysis
breakpoint <- Return.portfolio(R = returns_spgsci, wealth.index = TRUE)
breakpoint_nonxts <- rowMeans(breakpoint)

# Run breakpoint analysis; 
# Note: this might take around 10 mins to run through the detection
breakpoint_spgsci <- breakpoints(breakpoint ~ c(1:length(breakpoint)), h = 0.1)

# Plot and summarise breakpoints
plot(breakpoint_spgsci)
summary(breakpoint_spgsci)

# Extract breakpoints and convert indices to dates
breakpoint_indices <- breakpoint_spgsci$breakpoints
print(breakpoint_indices)  # Display breakpoint indices

# Convert breakpoint indices to dates
break_dates <- index(returns_spgsci[breakpoint_indices])
break_dates <- as.character(break_dates)
print(break_dates)

# Prepare breakpoints for plotting, including the start and end points of the series
full_breakpoints <- c(1, breakpoint_indices, length(returns_spgsci))
breakpoint_dates <- c(min(data_log_returns$Dates), break_dates, max(data_log_returns$Dates))
```

**Plotting Change Points Detection**

```{r plotting_change_points_detection, warning=FALSE}
options(width = 80)
# Convert data into a data frame
data_frame_1 <- data.frame(Day_Index = 1:length(returns_spgsci_non), 
                           Daily_Return = returns_spgsci_non)

breakpoint_dates[1] <- "2005-10-03"

# Create the plot
plot_spgsci <- ggplot(data_frame_1, aes(x = Day_Index, y = Daily_Return)) +
  geom_bar(stat="identity", fill="darkblue", width = 1.5) +  
  geom_vline(xintercept = breakpoint_indices, color = "orange", 
             linetype = "dashed", linewidth = 0.8) +  
  geom_vline(xintercept = c(1, length(returns_spgsci)), color = "red", 
             linetype = "solid", linewidth = 1) +
  scale_x_continuous(breaks = full_breakpoints, labels = breakpoint_dates) +
  scale_y_continuous(limits = c(-0.1, 0.1)) +  
  theme_light() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                   margin = margin(t = 10), size = 10),  
        axis.title.x = element_text(margin = margin(t = 20, b = 10),size = 12), 
        axis.title.y = element_text(margin = margin(r = 10), size = 12),
        panel.grid.major.y = element_line(color = "gray", linewidth = 0.5, 
                                          linetype = "dotted"),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
        ) +
  labs(title = "SP GSCI Performance (Daily Log Return)",
       x = "Day Index (SP GSCI Year 2005-2024)", 
       y = "Daily Log Return") +
  theme(legend.position = "none")

print(plot_spgsci)

# Save the plot
ggsave("Log_SPGSCI_Performance_Plot.png", plot_spgsci, width = 12, height = 6, dpi = 300)
ggsave("Log_SPGSCI_Performance_Plot.eps", plot_spgsci, width = 12, height = 6, dpi = 300)


# Convert numeric vectors into a data frame for plotting purposes
data_frame_2 <- data.frame(
  Day_Index = 1:length(breakpoint),  # Create an index for each day
  # Numeric conversion may not be necessary if `breakpoint` is already numeric
  Accumulated_Return = as.numeric(breakpoint)  
)

# Initialize ggplot with data 
plot_spgsci_acc <- ggplot(data_frame_2, 
                          aes(x = Day_Index, y = Accumulated_Return)) +
  geom_line(color = "darkblue", size = 0.6) +  
  geom_vline(xintercept = breakpoint_indices, color = "orange", 
             linetype = "dashed", linewidth = 0.8) +  
  geom_vline(xintercept = c(1, length(breakpoint)), color = "red", 
             linetype = "solid", linewidth = 1) +  
  scale_x_continuous(breaks = full_breakpoints, labels = breakpoint_dates) +  
  scale_y_continuous(limits = range(data_frame_2$Accumulated_Return, 
                                    na.rm = TRUE)) +  
  theme_light() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                   margin = margin(t = 10), size = 10),  
        axis.title.x = element_text(margin = margin(t = 20, b = 10), size = 12), 
        axis.title.y = element_text(margin = margin(r = 10), size = 12),  
        panel.grid.major.y = element_line(color = "gray", linewidth = 0.5, 
                                          linetype = "dotted"),  
        panel.grid.major.x = element_blank(),  
        panel.background = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
        ) +
  labs(title = "SP GSCI Accumulated Performance",
       x = "Day Index (SP GSCI Year 2005-2024)", 
       y = "Accumulated Total Return") +
  theme(legend.position = "none")  

print(plot_spgsci_acc)

```

\newpage

# Section 3: Further Data Cleaning

**Further Data Cleaning**

After cleaning the data, we noticed that some commodities have prices that do not change during certain periods. Because these prices stay the same, their variance is very low, and the log return is zero throughout the period. This can cause problems with our regression analysis, such as high multicollinearity, which might lead to incorrect estimates or NA values. To ensure our model is reliable, we will review the data again and remove these commodities with constant prices.

```{r further_data_cleaning}
options(width = 80)
# Initialize a list to store data for each period
periods_data <- list()
# Initialize a vector to keep track of columns that are removed across all periods
all_removed_columns <- c()

# Process each period using a loop. 
# The last period is handled differently to include the final date
for (i in seq_along(breakpoint_dates)[-length(breakpoint_dates)]) {
  if (i < length(breakpoint_dates) - 1) {
    # Filter data for the current period, excluding the endpoint
    period_data <- data_log_returns %>%
      dplyr::filter(Dates >= breakpoint_dates[i] & Dates < breakpoint_dates[i + 1])
  } else {
    # Filter data for the last period, including the endpoint
    period_data <- data_log_returns %>%
      dplyr::filter(Dates >= breakpoint_dates[i] & Dates <= max(data_log_returns$Dates))
  }
  
  # Identify columns with 2 or fewer unique non-NA values, indicating low variance
  zero_var_cols <- sapply(period_data, function(x) length(unique(na.omit(x))) <=2)
  
  # If any columns meet the criteria, remove them from the period data
  if (any(zero_var_cols)) {
    removed_columns <- names(period_data)[zero_var_cols]
    all_removed_columns <- c(all_removed_columns, removed_columns)
    period_data <- period_data[, !zero_var_cols, drop = FALSE]
    
    # Print messages to the console about the removal process
    message(sprintf("Period %d: Removed %d columns with low variance.", 
                    i, sum(zero_var_cols)))
    message(sprintf("Columns removed: %s", paste(removed_columns, collapse = ", ")))
  } else {
    message(sprintf("Period %d: No columns removed.", i))
  }
  
  # Store the cleaned data for the current period
  periods_data[[i]] <- period_data
}

# Identify all unique columns that were removed in any period
columns_to_remove <- unique(all_removed_columns)
print(columns_to_remove)

# Create a new list to store the cleaned data for all periods after removing identified columns
cleaned_periods_data <- list()
for (i in seq_along(breakpoint_dates)[-length(breakpoint_dates)]) {
  if (i < length(breakpoint_dates) - 1) {
    period_data <- data_log_returns %>%
      dplyr::filter(Dates >= breakpoint_dates[i] & Dates < breakpoint_dates[i + 1])
  } else {
    period_data <- data_log_returns %>%
      dplyr::filter(Dates >= breakpoint_dates[i] & Dates <= max(data_log_returns$Dates))
  }
  
  # Remove the identified columns from each period data
  period_data <- period_data %>% dplyr::select(-all_of(columns_to_remove))
  
  # Store the further cleaned data
  cleaned_periods_data[[i]] <- period_data
}

# Access and potentially print out data for each period
# Useful for debugging or specific analyses
data_period_1 <- cleaned_periods_data[[1]]
data_period_2 <- cleaned_periods_data[[2]]
data_period_3 <- cleaned_periods_data[[3]]
data_period_4 <- cleaned_periods_data[[4]]
data_period_5 <- cleaned_periods_data[[5]]
data_period_6 <- cleaned_periods_data[[6]]
data_period_7 <- cleaned_periods_data[[7]]
data_period_8 <- cleaned_periods_data[[8]]

```

\newpage

# Section 4: Period 8 estimates for the Six Estimators

We only focus on **Period 8 estimations** as a proof of concept, which is an interesting example as the composite index S&P GSCI has a volatile behaviour. Recall that Period 8 runs from *06 March 2020* to *30 June 2022.*

**Set X and Y -- we regress Y** (log returns of S&P GSCI) **on X** (log return of 35 commodities) and we aim to identify significant shifts in the commodity markets as represented by the composite index S&P GSCI (Standard & Poor's Goldman Sachs Commodity Index). 

```{r set_X_and_Y}
options(width = 80)
y_8 <- data_period_8[[2]]
X_8 <- as.matrix(data_period_8[-c(1, 2)])
```

**Ordinary Least Squares (OLS)** -- this is the usual estimator

```{r OLS}
options(width = 80)
OLS_model_period_8 <- lm(y_8 ~ X_8)
OLS_beta_coef_8 <- as.vector(OLS_model_period_8$coefficients)
```

**Ridge Regression (RR)** -- this is the **L2** **penalised** version of the usual estimator **based on 10-fold CV** so more computationally expensive than OLS

```{r RR}
options(width = 80)
RR_ost_revised <- function(data_X, data_Y) { 
  xxx.temp <- data_X #input matrix X not \tilde(X)
  yyy.temp <- data_Y 
  glmnet_grid = c(0, 10^seq(-6,2,length = 49))
  # default k=10 fold CV
  glmnet_fit = cv.glmnet(data_X, data_Y, alpha=0, lambda=glmnet_grid) 
  lambda_min_RR_glmnet =  glmnet_fit$lambda.min
  # extract the model coefficients (thetas) corresponding to lambda.min
  theta_RR <- as.vector(coef(glmnet_fit, s="lambda.min"))
  return(theta_RR)
}
RR_model_period_8 <- RR_ost_revised(X_8, y_8)
RR_beta_coef_8 <- as.vector(RR_model_period_8)
```

Input the source code which already contain the functions for the following estimation

```{r source_code}
options(width = 80)
script_path_est_models <- "R_Code_SR_St_DSh_Sh.R"
source(script_path_est_models)
```

**Slab Regression (SR)** -- this is the **Slab** **penalised** version of the usual estimator; **no CV is needed** so computationally equivalent to OLS

```{r SR}
options(width = 80)
SR_model_period_8 <- SR_ost_revised(X_8, y_8)
SR_beta_coef_8 <- as.vector(SR_model_period_8)
```

**Stein (St)** -- this is the **shrinkage** estimator 1; **no CV is needed** so computationally equivalent to OLS

```{r St}
options(width = 80)
St_model_period_8 <- St_ost_revised(X_8, y_8)
St_beta_coef_8 <- as.vector(St_model_period_8)
```

**Diagonal Shrinkage (DSh)** -- this is the **shrinkage** estimator 2; **no CV is needed** so computationally equivalent to OLS

```{r DSh}
options(width = 80)
DSh_model_period_8 <- DSh_ost_revised(X_8, y_8)
DSh_beta_coef_8 <- as.vector(DSh_model_period_8)
```

**Shrinkage (Sh)** -- this is the **shrinkage** estimator 3; **no CV is needed** so computationally equivalent to OLS, though this is not in reality true as its computation requires a non-scalable operation (known as Sylvester equation)

```{r Sh}
options(width = 80)
Sh_model_period_8 <- Sh_ost_revised(X_8, y_8)
Sh_beta_coef_8 <- as.vector(Sh_model_period_8)
```

**Results**

We put all results into a data frame for the ease of comparing these estimates

```{r results}
options(width = 80)
# Create a data frame with each column representing a model's beta coefficients
beta_coefficients_table <- data.frame(
  OLS_beta_coef_8 = OLS_beta_coef_8,
  RR_beta_coef_8 = RR_beta_coef_8,
  SR_beta_coef_8 = SR_beta_coef_8,
  St_beta_coef_8 = St_beta_coef_8,
  DSh_beta_coef_8 = DSh_beta_coef_8,
  Sh_beta_coef_8 = Sh_beta_coef_8
)

# Transpose the data frame to make beta coefficients as rows
beta_coefficients_table_t <- t(beta_coefficients_table)

colnames(beta_coefficients_table_t) <- c("Intercept", colnames(X_8))

# View the updated table
print(beta_coefficients_table_t)
```


```{r }
# Assuming beta_coefficients_table_t is your transposed matrix

# Convert matrix to data frame
beta_df <- as.data.frame(beta_coefficients_table_t)

ols_row_values <- beta_coefficients_table_t["OLS_beta_coef_8", ]

# Step 3: Order the columns of the matrix based on this row's values
sorted_col_indices <- order(ols_row_values, decreasing = TRUE)  # Get indices for sorting high to low

# Step 4: Sort the entire matrix based on these indices
sorted_beta_coefficients_table_t <- beta_coefficients_table_t[, sorted_col_indices]

# Step 5: Check the sorted matrix
colnames_filtered <- colnames(sorted_beta_coefficients_table_t)[colnames(sorted_beta_coefficients_table_t) != "Intercept"]

colnames_transformed <- gsub("Log_Return_(.*)_Comdty", "\\1", colnames_filtered)

# Create a data frame from the column names
colnames_df <- data.frame(Com_name = colnames_transformed)

# Define the filename and path
filename <- "colnames_transformed.csv"

# Write the data frame to a CSV file
write.csv(colnames_df, filename, row.names = FALSE)

```

**Visulisation**

We can observe that the beta of CL1_Comdty, at 0.3888456 in the OLS model, is higher than that of the others. Why is this the case? To explore this, let's visualize the prices of CL1_Comdty and CO1_Comdty along with the S&P GSCI Index to examine their correlations.

```{r visulisation}
visulisation <- data_filtered %>%
  select(Dates, CL1_Comdty, CO1_Comdty) %>%  # Select specific columns
  filter(Dates >= as.Date("2022-06-30") & Dates <= as.Date("2024-05-09"))

start_date <- as.Date("2022-06-30")
end_date <- as.Date("2024-05-09")

ggplot(data = visulisation, aes(x = Dates)) +
  geom_line(aes(y = CL1_Comdty, color = "WTI Crude Oil (CL1_Comdty)")) +  # Plot CL1_Comdty
  geom_line(aes(y = CO1_Comdty, color = "Brent Crude Oil (CO1_Comdty)")) +  # Plot CO1_Comdty
  #geom_line(aes(y = SPGSCI_Index, color = "SPGSCI_Index")) +  # Plot SPGSCI_Index
  labs(title = "Commodity Price Over Time",
       x = "Date",
       y = "Index Value",
       color = "Index Type") +  # Labels for the plot
  theme_minimal() +  # Use a minimal theme
  scale_color_manual(values = c("WTI Crude Oil (CL1_Comdty)" = "blue", 
                                "Brent Crude Oil (CO1_Comdty)" = "green"))+
  geom_vline(aes(xintercept = start_date), color = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = end_date), color = "red", linetype = "dashed")
```

We selected commodities with an OLS beta greater than 5%, including Gasoil Petroleum (QS_Comdty), WTI Crude Oil (CL_Comdty), Wheat (W_Comdty), and Copper (HG_Comdty). We then visualized their log returns compared to the SPGSCI.

```{r visulisation_log_return, fig.width=8, fig.height=12}
visulisation_log_return <- data_period_8 %>%
  select(Dates, Log_Return_SPGSCI, Log_Return_QS1_Comdty, Log_Return_CL1_Comdty,
         Log_Return_HG1_Comdty, Log_Return_W_1_Comdty)

visulisation_log_return_long <- visulisation_log_return %>%
  pivot_longer(
    cols = -Dates,  # Exclude the Dates column from the reshaping
    names_to = "Commodity",  # New column for commodity names
    values_to = "Log_Return"  # New column for log return values
  )


# Reorder the levels of the Commodity factor, if necessary, to control the order of the panels
visulisation_log_return_long$Commodity <- factor(visulisation_log_return_long$Commodity, levels = c("Log_Return_SPGSCI",
                                                          "Log_Return_QS1_Comdty",
                                                          "Log_Return_CL1_Comdty",
                                                          "Log_Return_W_1_Comdty",
                                                          "Log_Return_HG1_Comdty" 
                                                          ),
                                                 labels = c("SP GSCI Over Time",
                                                            "Gasoil Petroleum Over Time (QS1 Comdty)",
                                                            "WTI Crude Oil Over Time (CL1 Comdty)",
                                                            "Wheat Over Time (W1 Comdty)",
                                                            "Copper Over Time (HG1 Comdty)"
                                                            ))


custom_colors <- c("SP GSCI Over Time" = "red", "Gasoil Petroleum Over Time (QS1 Comdty)" = "blue",
                   "WTI Crude Oil Over Time (CL1 Comdty)" = "green", "Copper Over Time (HG1 Comdty)" = "orange", "Wheat Over Time (W1 Comdty)" = "brown")

# Generate the plot with customized facet labels
log_return_plot <- ggplot(visulisation_log_return_long, aes(x = Dates, y = Log_Return, color = Commodity)) +
  geom_line() +
  scale_color_manual(values = custom_colors) +
  facet_wrap(~Commodity, scales = "free", ncol = 1) +
  theme_minimal() +
  labs(title = "Log Returns of Commodities Over Time",
       x = "Date",
       y = "Log Return") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=6),
        strip.background = element_blank(),  # Removes background from facet labels
        strip.text.x = element_text(face = "bold"))  # Bolds the facet labels

# Print the plot
print(log_return_plot)


```
