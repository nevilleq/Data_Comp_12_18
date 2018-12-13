Forward, Backward Subset Selection - MSE
================
Quinton Neville
12/8/2018

Load, clean, manipulate, and tidy the data
==========================================

``` r
# Import data
cancer_raw = readr::read_csv("./Data/Cancer_Registry.csv") %>% 
  janitor::clean_names() 

#dim(cancer_raw)
#head(cancer_raw)

# Check NA values for each column
#n_NA = sapply(cancer_raw[1:34], function(x) sum(length(which(is.na(x)))))
#n_NA

# Check the percentage of NA values for each column
#percentage_NA = sapply(cancer_raw[1:34], function(x) sum(length(which(is.na(x)))) / nrow(cancer_raw))
#percentage_NA %>% data.frame()

#Pulling quartiles for study_per_cap categorical manipulation
study.quart <- with(cancer_raw, study_per_cap[study_per_cap > 0]) %>%
  quantile(., probs = c(0.25, 0.5, 0.75))

#To add a variable for region, state has too many factors for CV and data is not large enough for 51 fct var to be useful
#From data(state) built in R, going to left_join with cancer.df
#District of Columbia is missing, so I'm just going to impute it with the same factor as Maryland, 2
data(state)
state.df <- cbind(state.name, state.region) %>% as.tibble() %>% rename(state = state.name, region = state.region)

#Variable Manipulation
cancer.df <- cancer_raw %>%    #Remove Rows with > 20% missing 
  dplyr::select(-pct_some_col18_24) %>%  #Remove for too many missing
  mutate(
    pct_non_white = pct_black + pct_asian + pct_other_race, #Creating white, non-white percentages variables
    state = str_split_fixed(geography, ", ", 2)[ ,2] %>% as.factor(), #pulling state variable and casting as factor, possible region?
    binned_inc_lb = str_split_fixed(binned_inc, ", ", 2)[ ,1] %>% parse_number(), #pulling numeric lower bound
    binned_inc_ub = str_split_fixed(binned_inc, ", ", 2)[ ,2] %>% parse_number(), #pulling numeric upper bound
    binned_inc_point = (binned_inc_lb + binned_inc_ub)/2, #computing point estimate from ub,lb (intervals symmetric)
    study_quantile = ifelse(study_per_cap == 0, "None", 
                           ifelse(study_per_cap > 0 & study_per_cap <= study.quart[1], "Low", 
                                  ifelse(study_per_cap > study.quart[1] & study_per_cap <= study.quart[2], "Moderate", 
                                         ifelse(study_per_cap > study.quart[2] & study_per_cap <= study.quart[3], "High", 
                                                "Very High")))),
    study_quantile = as.factor(study_quantile) %>% fct_relevel(., "None", "Low", "Moderate", "High", "Very High"),
    avg_deaths_yr_pop = avg_deaths_per_year/pop_est2015,  #incorporate two vars into one (multicollinearity)
    avg_ann_count_pop = avg_ann_count/pop_est2015 #incorporate two vars into one (multicollinearity)
  ) %>%
  left_join(., state.df) %>%
  mutate(region = ifelse(is.na(region), 2, region),
         region = as.factor(region)) %>%
  dplyr::select(-c(binned_inc, geography, state, study_per_cap))

######TO avoid singular matrixes for our model fits, we need remove highly correlated and/or direct linear combos
###### Also lets make the response index 1, makes everything easier

cancer.df <- cancer.df %>%
  dplyr::select(-c(avg_deaths_per_year, avg_ann_count, pct_black, pct_asian, pct_other_race)) %>%
  dplyr::select(target_death_rate, study_quantile, region, everything())


#mostly need to decide between white/non-white or all race % AND binned_inc_lb, ub or median income
```

Imputing Values with less than 20% missing (two variables)
==========================================================

-   pct\_employed16\_over ~ 4%
-   pct\_private\_coverage\_alone ~ 20%

``` r
#Impute those missing less than 20%
#1. pct_employed16_over
#2. pct_private_coverage_alone

#Set up appropriate test and train for pct_employed16_over (removing other missing % variable and response (target death))
train.df <- cancer.df %>% dplyr::select(-c(pct_private_coverage_alone, target_death_rate)) %>% filter(!is.na(pct_employed16_over))
test.df <- cancer.df %>% dplyr::select(-c(pct_private_coverage_alone, target_death_rate)) %>% filter(is.na(pct_employed16_over))

#Function for imputation (after correct test, train set up), charstring must literally be the character of impute variable i.e. "var1"
impute.lasso <- function(train.df, test.df, charstring){

  if ((charstring %in% names(train.df))) {
    
#pull variable index
index <- which(names(train.df) == charstring)
  
#Set up Matrices
#Create Design Matrix Train
X <- train.df[ , -index] %>%
  names() %>% 
  paste(., collapse = "+") %>%    
  paste("~ ", .) %>%
  formula() %>%
  model.matrix(.,train.df)
  
#Create Design Matrix Test
X1 <- test.df[, -index] %>%
  names() %>% 
  paste(., collapse = "+") %>%    
  paste("~ ", .) %>%
  formula() %>%
  model.matrix(., test.df)

#Remove Intercept  
X <- X[,-1]
X1 <- X1[,-1]

#Create Response vector (as matrix)
Y <- train.df[, index] %>% as.matrix()

#Optimize lambda
lambda.grid <- 10^seq(-3,1,length = 100)

#CV n = 10
cv.lasso <- cv.glmnet(X, Y, alpha = 1, intercept = TRUE, lambda = lambda.grid, family = "gaussian")

#Grab optimal lambda
opt.lambda.lasso <- cv.lasso$lambda.min

#Run model
unemploy.lasso <- glmnet(X, Y, alpha = 1, intercept = TRUE, lambda = opt.lambda.lasso, family = "gaussian")

#Return predictions
predict(unemploy.lasso, newx = X1)
  }else{
    stop("Error: Incorrect variable name")
  }
}

#Impute employed16_over_preds (first since it has less missing data ~4%)
employed16_over_preds <- impute.lasso(train.df = train.df, test.df, "pct_employed16_over")

#Set up appropriate test and train
train.df <- cancer.df %>% dplyr::select(-c(pct_employed16_over, target_death_rate)) %>% filter(!is.na(pct_private_coverage_alone))
test.df <- cancer.df %>% dplyr::select(-c(pct_employed16_over, target_death_rate)) %>% filter(is.na(pct_private_coverage_alone))

#Impute pct_private_coverage_alone (second since it has more missing data ~20%)
pct_private_coverage_alone_preds <- impute.lasso(train.df = train.df, test.df, "pct_private_coverage_alone")

#Replace Imputed values
cancer.df <- cancer.df %>%
  mutate(imp_pct_employed16_over = ifelse(is.na(pct_employed16_over),
                                          employed16_over_preds, pct_employed16_over),
         imp_pct_private_coverage_alone = ifelse(is.na(pct_private_coverage_alone),
                                          pct_private_coverage_alone_preds, pct_private_coverage_alone)
        )

#Looks good, so we will replace imputed variables in our final data set
cancer.df <- cancer.df %>%
  dplyr::select(-c(pct_employed16_over, pct_private_coverage_alone))
```

See generally what the optimal number of preds in an lm() ought to be (CP, adj*R*<sup>2</sup>, BIC (similar to AIC)).

``` r
#Summary of models for each size (one model per size)
#Set max number of predictors
nvmax <- ncol(cancer.df) - 1
reg.subsets <- cancer.df %>% 
  dplyr::select(-binned_inc_point) %>% #Remove linear dependency (average of lower, upper bound)
  regsubsets(target_death_rate ~ ., data = ., really.big = FALSE, nvmax = nvmax)
rs <- summary(reg.subsets)

# Plots of Cp and Adj-R2 as functions of parameters
r.df <- tibble(
  preds = 1:nvmax,
  cp = rs$cp,
  adjr2 = rs$adjr2,
  bic = rs$bic,
  step = 1:nvmax
)

#Grab max r^2 pred number
max.r2 <- with(r.df, which.max(adjr2))

#Grab min bic pred number
min.bic <- with(r.df, which.min(bic))

cp.plot <- r.df %>% ggplot(aes(x = preds, y = cp, size = cp)) +
  geom_point(colour = "purple") +
  geom_line(alpha = 0.5, colour = "purple", size = 1.25) +
 # geom_point(aes(x = preds, step),color = "black",size = 0.5) +
  geom_line(aes(x = preds, step), size = 1, linetype = 2, color = "red") +
  labs(
    x = "Number of Preds",
    y = "CP Criterion",
    title = "Optimal Number of Preds, CP"
  ) + 
  ylim(c(0, 300)) + 
  theme(legend.position = "none") 

adjr2.plot <- r.df %>% ggplot(aes(x = preds, y = adjr2, size = 1 - adjr2)) +
  geom_point(colour = "purple") +
  geom_line(alpha = 0.5, colour = "purple", size = 1.25) +
  geom_vline(xintercept = max.r2, size = 1, linetype = 2, color = "red") +
  labs(
    x = "Number of Preds",
    y = "Adjusted R^2",
    title = "Optimal Number of Preds, Adj.R^2"
  ) + theme(legend.position = "none")

bic.plot <- r.df %>% ggplot(aes(x = preds, y = bic, size = bic)) +
  geom_point(colour = "purple") +
  geom_line(alpha = 0.5, colour = "purple", size = 1.25) +
  geom_vline(xintercept = min.bic, size = 1, linetype = 2, color = "red") +
  labs(
    x = "Number of Preds",
    y = "Bayesian Information Criterion",
    title = "Optimal Number of Preds, BIC"
  ) + theme(legend.position = "none")

(cp.plot + adjr2.plot) / bic.plot
```

    ## Warning: Removed 8 rows containing missing values (geom_point).

    ## Warning: Removed 8 rows containing missing values (geom_path).

<img src="Q_subset_selection_alg_files/figure-markdown_github/unnamed-chunk-5-1.png" width="90%" style="display: block; margin: auto;" />

Based on the plots above, *C**P* criterion in the upper left with *p* ≤ *C**P* constraint in red, that somewhere around a 25-27 predictor model ought to be optimal. With respect to adjusted *R*<sup>2</sup>, it appears as though we reach a converging maximum starting around 20 and negligible increase after, where additional predictors have diminishing marginal return. Lastly, BIC is more conservative (penalizing more strongly for more predictors) and seems to suggest between a 15-20 predictor model (closer to 20). This may inform our subset selection criterion.

Feature Selection Lasso
=======================

Inputes are data (cancer.df), response (default target\_death\_rate), lambda (penalty, default = 1, higher penalty removes more coefficients), print(logical whether or not to print the results)

``` r
lasso.feature.removal <- function(data = cancer.df, response = "target_death_rate", lambda = 1, print = TRUE) {
  
  #Snag the index for matrices
  index <- which(names(data) == response)
  
  #Response 
  Y <- as.matrix(data[, index])
  
  #Design matrix
  X <- data[ ,-index] %>%
    names() %>% 
    paste(., collapse = "+") %>%    
    paste("~ ", .) %>%
    formula() %>%
    model.matrix(.,data)
  X <- X[,-1]  
  #Fit Model
  mod.lasso <- glmnet(X,Y,alpha = 1,intercept = T,lambda = lambda, family = "gaussian")
  #Print if true
  if (isTRUE(print)) {
  coef.lasso   <- coef(mod.lasso)[,1]
  remove.index <- which(coef.lasso == 0)
  keep.index   <- which(coef.lasso != 0)
                        
    print(sprintf("Lambda = %f : Lasso method removed %i variables:", lambda, length(remove.index)))
    print(paste(names(remove.index)))
    print(sprintf("Lambda = %f : Lasso method selected %i variables:",lambda, length(keep.index)))
    print(paste(names(keep.index)))
    cat("\n")
  }
 # return(names(which(coef.lasso == 0)))
}

map(.x = c(0.01, .1, 0.5, 1, 2), ~lasso.feature.removal(cancer.df, lambda = .x))
```

    ## [1] "Lambda = 0.010000 : Lasso method removed 0 variables:"
    ## character(0)
    ## [1] "Lambda = 0.010000 : Lasso method selected 38 variables:"
    ##  [1] "(Intercept)"                    "study_quantileLow"             
    ##  [3] "study_quantileModerate"         "study_quantileHigh"            
    ##  [5] "study_quantileVery High"        "region2"                       
    ##  [7] "region3"                        "region4"                       
    ##  [9] "incidence_rate"                 "med_income"                    
    ## [11] "pop_est2015"                    "poverty_percent"               
    ## [13] "median_age"                     "median_age_male"               
    ## [15] "median_age_female"              "avg_household_size"            
    ## [17] "percent_married"                "pct_no_hs18_24"                
    ## [19] "pct_hs18_24"                    "pct_bach_deg18_24"             
    ## [21] "pct_hs25_over"                  "pct_bach_deg25_over"           
    ## [23] "pct_unemployed16_over"          "pct_private_coverage"          
    ## [25] "pct_emp_priv_coverage"          "pct_public_coverage"           
    ## [27] "pct_public_coverage_alone"      "pct_white"                     
    ## [29] "pct_married_households"         "birth_rate"                    
    ## [31] "pct_non_white"                  "binned_inc_lb"                 
    ## [33] "binned_inc_ub"                  "binned_inc_point"              
    ## [35] "avg_deaths_yr_pop"              "avg_ann_count_pop"             
    ## [37] "imp_pct_employed16_over"        "imp_pct_private_coverage_alone"
    ## 
    ## [1] "Lambda = 0.100000 : Lasso method removed 6 variables:"
    ## [1] "region3"           "poverty_percent"   "median_age"       
    ## [4] "pct_bach_deg18_24" "binned_inc_lb"     "binned_inc_point" 
    ## [1] "Lambda = 0.100000 : Lasso method selected 32 variables:"
    ##  [1] "(Intercept)"                    "study_quantileLow"             
    ##  [3] "study_quantileModerate"         "study_quantileHigh"            
    ##  [5] "study_quantileVery High"        "region2"                       
    ##  [7] "region4"                        "incidence_rate"                
    ##  [9] "med_income"                     "pop_est2015"                   
    ## [11] "median_age_male"                "median_age_female"             
    ## [13] "avg_household_size"             "percent_married"               
    ## [15] "pct_no_hs18_24"                 "pct_hs18_24"                   
    ## [17] "pct_hs25_over"                  "pct_bach_deg25_over"           
    ## [19] "pct_unemployed16_over"          "pct_private_coverage"          
    ## [21] "pct_emp_priv_coverage"          "pct_public_coverage"           
    ## [23] "pct_public_coverage_alone"      "pct_white"                     
    ## [25] "pct_married_households"         "birth_rate"                    
    ## [27] "pct_non_white"                  "binned_inc_ub"                 
    ## [29] "avg_deaths_yr_pop"              "avg_ann_count_pop"             
    ## [31] "imp_pct_employed16_over"        "imp_pct_private_coverage_alone"
    ## 
    ## [1] "Lambda = 0.500000 : Lasso method removed 18 variables:"
    ##  [1] "study_quantileLow"              "study_quantileModerate"        
    ##  [3] "study_quantileHigh"             "study_quantileVery High"       
    ##  [5] "region3"                        "region4"                       
    ##  [7] "med_income"                     "pop_est2015"                   
    ##  [9] "poverty_percent"                "median_age"                    
    ## [11] "pct_bach_deg18_24"              "pct_public_coverage"           
    ## [13] "pct_married_households"         "pct_non_white"                 
    ## [15] "binned_inc_lb"                  "binned_inc_point"              
    ## [17] "imp_pct_employed16_over"        "imp_pct_private_coverage_alone"
    ## [1] "Lambda = 0.500000 : Lasso method selected 20 variables:"
    ##  [1] "(Intercept)"               "region2"                  
    ##  [3] "incidence_rate"            "median_age_male"          
    ##  [5] "median_age_female"         "avg_household_size"       
    ##  [7] "percent_married"           "pct_no_hs18_24"           
    ##  [9] "pct_hs18_24"               "pct_hs25_over"            
    ## [11] "pct_bach_deg25_over"       "pct_unemployed16_over"    
    ## [13] "pct_private_coverage"      "pct_emp_priv_coverage"    
    ## [15] "pct_public_coverage_alone" "pct_white"                
    ## [17] "birth_rate"                "binned_inc_ub"            
    ## [19] "avg_deaths_yr_pop"         "avg_ann_count_pop"        
    ## 
    ## [1] "Lambda = 1.000000 : Lasso method removed 22 variables:"
    ##  [1] "study_quantileLow"              "study_quantileModerate"        
    ##  [3] "study_quantileHigh"             "study_quantileVery High"       
    ##  [5] "region3"                        "region4"                       
    ##  [7] "med_income"                     "pop_est2015"                   
    ##  [9] "poverty_percent"                "median_age"                    
    ## [11] "avg_household_size"             "pct_no_hs18_24"                
    ## [13] "pct_bach_deg18_24"              "pct_private_coverage"          
    ## [15] "pct_public_coverage"            "pct_white"                     
    ## [17] "pct_married_households"         "pct_non_white"                 
    ## [19] "binned_inc_lb"                  "binned_inc_point"              
    ## [21] "imp_pct_employed16_over"        "imp_pct_private_coverage_alone"
    ## [1] "Lambda = 1.000000 : Lasso method selected 16 variables:"
    ##  [1] "(Intercept)"               "region2"                  
    ##  [3] "incidence_rate"            "median_age_male"          
    ##  [5] "median_age_female"         "percent_married"          
    ##  [7] "pct_hs18_24"               "pct_hs25_over"            
    ##  [9] "pct_bach_deg25_over"       "pct_unemployed16_over"    
    ## [11] "pct_emp_priv_coverage"     "pct_public_coverage_alone"
    ## [13] "birth_rate"                "binned_inc_ub"            
    ## [15] "avg_deaths_yr_pop"         "avg_ann_count_pop"        
    ## 
    ## [1] "Lambda = 2.000000 : Lasso method removed 26 variables:"
    ##  [1] "study_quantileLow"              "study_quantileModerate"        
    ##  [3] "study_quantileHigh"             "study_quantileVery High"       
    ##  [5] "region3"                        "region4"                       
    ##  [7] "med_income"                     "pop_est2015"                   
    ##  [9] "poverty_percent"                "median_age"                    
    ## [11] "avg_household_size"             "pct_no_hs18_24"                
    ## [13] "pct_bach_deg18_24"              "pct_private_coverage"          
    ## [15] "pct_emp_priv_coverage"          "pct_public_coverage"           
    ## [17] "pct_white"                      "pct_married_households"        
    ## [19] "birth_rate"                     "pct_non_white"                 
    ## [21] "binned_inc_lb"                  "binned_inc_ub"                 
    ## [23] "binned_inc_point"               "avg_ann_count_pop"             
    ## [25] "imp_pct_employed16_over"        "imp_pct_private_coverage_alone"
    ## [1] "Lambda = 2.000000 : Lasso method selected 12 variables:"
    ##  [1] "(Intercept)"               "region2"                  
    ##  [3] "incidence_rate"            "median_age_male"          
    ##  [5] "median_age_female"         "percent_married"          
    ##  [7] "pct_hs18_24"               "pct_hs25_over"            
    ##  [9] "pct_bach_deg25_over"       "pct_unemployed16_over"    
    ## [11] "pct_public_coverage_alone" "avg_deaths_yr_pop"

    ## [[1]]
    ## NULL
    ## 
    ## [[2]]
    ## NULL
    ## 
    ## [[3]]
    ## NULL
    ## 
    ## [[4]]
    ## NULL
    ## 
    ## [[5]]
    ## NULL

K-fold cross validation funciton for subset selecion (MSE, AIC, BIC criterion)
==============================================================================

``` r
sampleSize <- nrow(cancer.df)

mseCV <- function(data.df, kfolds = 10){
  folds <- sample(1:kfolds, sampleSize, rep = T)
  mse <- rep(0, kfolds)
  for (k in 1:kfolds) {
    train.df <- data.df[folds != k,]
    test.df <- data.df[folds == k,]
    
    lm.mod <- lm(target_death_rate ~ ., data = train.df)
    preds <- predict(lm.mod, newdata = test.df)
    mse[k] <- with(test.df,mean((target_death_rate - preds)^2))
  }
  mean(mse)
}

aicCV <- function(data.df, kfolds = 10){
  folds <- sample(1:kfolds, sampleSize, rep = T)
  aic <- rep(0, kfolds)
  for (k in 1:kfolds) {
    train.df <- data.df[folds != k,]
    test.df <- data.df[folds == k,]
    
    lm.mod <- lm(target_death_rate ~ ., data = train.df)
    aic[k] <- AIC(lm.mod)
  }
  mean(aic)
}

bicCV <- function(data.df, kfolds = 10){
  folds <- sample(1:kfolds, sampleSize, rep = T)
  bic <- rep(0, kfolds)
  for (k in 1:kfolds) {
    train.df <- data.df[folds != k,]
    test.df <- data.df[folds == k,]
    
    lm.mod <- lm(target_death_rate ~ ., data = train.df)
    bic[k] <- BIC(lm.mod)
  }
  mean(bic)
}
```

1. Models
=========

Scale continuous predictors

``` r
cancer.df <- bind_cols(
    cancer.df %>%
    dplyr::select(c(target_death_rate:region)),
    cancer.df %>%
    dplyr::select(-c(target_death_rate:region)) %>%
    scale() %>% as.tibble()
)
```

Forward Subset Selection with k-Fold CV and MSE criterion
---------------------------------------------------------

Function to run the subset algorithm, it will output the variable selection process so you visualize how it works. Comment out the set seed to get true variability in the subsets, only leave it in for reproducibility.

-   The function `f.subset.mse(sub.cancer.df, maxPreds = 10, nfolds = 5)` takes in a data frame, max number of predictors you want the algorithm to run through (can't be larger than the number of predictors in the data frame), and the number of folds for the cross validation.

##### Note mustreorder data frame so response is the first column of the data frame (see above code chunk)

##### Also have to take out state variable, too many factor levels and CV process can fail (re add later to see if it's useful)

-   The way the algorithm works is

1.  Starts with an empty (null) set of current predictors in the data frame input
2.  Starts with a full set of available predictors
3.  Loops through all available preds and adds one at a time, splits the data into kfolds test/train (CV), fits lm() with current preds, stores MSE
4.  Selects the predictor whose CV MSE was smallest
5.  Adds 'best' predictor to current predictors, removes it from available predictors
6.  Stores that predictor set as 'best' and records the indices and CV MSE in a matrix
7.  Repeats 3.- 5. until you have met the maximum number of predictors you specified
8.  Returns a list with all 'best' predictor sets at each iteration, matrix of results, and prespecified max preds
9.  Last few lines of code output the optimal (minimum MSE) predictor set

*notes - CV is set to 5 for speed, changing to 10 will take ~1.5 times as long to run* *notes - the print() lines within the function are so you can visualize how the algorithm is working, they can be commented out to reduce clutter when you want to run the selection function iteratively (avoids reams of output)*

``` r
#Forward Subset#
#set.seed(4) #Comment this out if you want to really get a different subset each time

#Function
#inputs:
#full data frame (whatever that ends up being)
#maximum number you want to look at subsets up to (no larger than p)
#Number of folds for the cross validation step, 10 is a little slow but better for final paper
#criterion to minimize -  either "mse", or "aic". Must put those in just as they appear with "

f.subset <- function(sub.cancer.df, maxPreds = 10, nfolds = 5, criterion = "mse") {

#number of possible predictors
## Allthe predictors (their indices).
allPreds <- 1:ncol(sub.cancer.df)
allPreds <- allPreds[-1] #substract the response

#current set of preds (starts empty), and available
currPreds <- c()
availPreds <- setdiff(allPreds,currPreds)
#Record the min errors
minError <- c()
#Initalize pred list and mse result matrix (just for storage)
pred.list <- list()
result.mat <- matrix(nrow = ncol(cancer.df) - 1, ncol = 2)
#Initalize iteration
i <- 1
#Forward selection loop
while (length(currPreds) < maxPreds) {
  ##add predictor which decreases MSE (as determined by CV or
  ##Bootstrapping)
  ## The MSE/AIC computed as we add each of the available predictors
  allError <- c()
  for (id in availPreds) {
    data.df <- sub.cancer.df[,c(currPreds,id,1)] #Only the preds we want in the df
    if (criterion == "mse") {
      error <- mseCV(data.df, nfolds)  #cross validate mse for the model with the added pred
    } else if (criterion == "aic") {
      error <- aicCV(data.df, nfolds) #same but with AIC
    } else if (criterion == "bic") {
      error <- bicCV(data.df, nfolds) #same bic
    } else {
      print("Invalid criterion")
      stop("Wrong criterion input")
    }
    allError <- c(allError,error)  #Collect all the errors for every variable added
  }
  ##Find the min
  id <- which.min(allError)
  ##get the best predictor and lowest MSE/AIC
  bestPred <- availPreds[id]
  bestError <- min(allError)
  ##Add these into the collection
  currPreds <- c(currPreds,bestPred)
  minError <- c(minError,bestError)
  ##Take out the pred you just added from consideration
  availPreds <- setdiff(allPreds,currPreds)
  ## Print stuff out for debugging and attention-grabbing
  print(sprintf("Iteration: %i Predictor Added: %s %s Value: %s",i, names(sub.cancer.df[,bestPred]), criterion, bestError))
  print(currPreds)
  result.mat[i,] <- c(bestPred,bestError)         #You can also comment out all print() later, it's just to watch the algorithm work
  pred.list[[i]] <- currPreds
  i <- i + 1
    }
  return(list(pred.list = pred.list, result.mat = result.mat, maxPreds = maxPreds)) #returns list of preds, results, and max num preds
}

#Run Subset, call function, output is a list with predictor sets, reslut matrix and maxPreds
f.mse.list <- f.subset(cancer.df, maxPreds = ncol(cancer.df) - 1, nfolds = 5, criterion = "mse")
```

    ## [1] "Iteration: 1 Predictor Added: avg_deaths_yr_pop mse Value: 527.109955924533"
    ## [1] 30
    ## [1] "Iteration: 2 Predictor Added: median_age_female mse Value: 284.351419300116"
    ## [1] 30 10
    ## [1] "Iteration: 3 Predictor Added: region mse Value: 245.589092436763"
    ## [1] 30 10  3
    ## [1] "Iteration: 4 Predictor Added: incidence_rate mse Value: 217.473871372474"
    ## [1] 30 10  3  4
    ## [1] "Iteration: 5 Predictor Added: pct_unemployed16_over mse Value: 201.582565778764"
    ## [1] 30 10  3  4 18
    ## [1] "Iteration: 6 Predictor Added: pct_hs18_24 mse Value: 190.359016634293"
    ## [1] 30 10  3  4 18 14
    ## [1] "Iteration: 7 Predictor Added: pct_emp_priv_coverage mse Value: 182.467718240032"
    ## [1] 30 10  3  4 18 14 20
    ## [1] "Iteration: 8 Predictor Added: pct_private_coverage mse Value: 174.221251373531"
    ## [1] 30 10  3  4 18 14 20 19
    ## [1] "Iteration: 9 Predictor Added: binned_inc_ub mse Value: 170.846732455674"
    ## [1] 30 10  3  4 18 14 20 19 28
    ## [1] "Iteration: 10 Predictor Added: median_age_male mse Value: 168.561746710058"
    ##  [1] 30 10  3  4 18 14 20 19 28  9
    ## [1] "Iteration: 11 Predictor Added: pct_hs25_over mse Value: 165.987124167932"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16
    ## [1] "Iteration: 12 Predictor Added: avg_ann_count_pop mse Value: 165.678929278739"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31
    ## [1] "Iteration: 13 Predictor Added: pct_public_coverage mse Value: 163.929611809892"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21
    ## [1] "Iteration: 14 Predictor Added: pct_public_coverage_alone mse Value: 162.269843629525"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22
    ## [1] "Iteration: 15 Predictor Added: imp_pct_employed16_over mse Value: 159.617704106819"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32
    ## [1] "Iteration: 16 Predictor Added: median_age mse Value: 158.958783023713"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8
    ## [1] "Iteration: 17 Predictor Added: pct_no_hs18_24 mse Value: 159.176447894592"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13
    ## [1] "Iteration: 18 Predictor Added: percent_married mse Value: 158.901090588193"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12
    ## [1] "Iteration: 19 Predictor Added: med_income mse Value: 159.140420184157"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5
    ## [1] "Iteration: 20 Predictor Added: binned_inc_lb mse Value: 158.383537481541"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27
    ## [1] "Iteration: 21 Predictor Added: imp_pct_private_coverage_alone mse Value: 158.300808763086"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33
    ## [1] "Iteration: 22 Predictor Added: birth_rate mse Value: 158.64876490798"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25
    ## [1] "Iteration: 23 Predictor Added: pct_bach_deg25_over mse Value: 158.427590177624"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [1] "Iteration: 24 Predictor Added: pct_married_households mse Value: 157.103533223164"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24
    ## [1] "Iteration: 25 Predictor Added: pct_bach_deg18_24 mse Value: 158.478131276252"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24 15
    ## [1] "Iteration: 26 Predictor Added: pop_est2015 mse Value: 158.264773859822"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24 15  6
    ## [1] "Iteration: 27 Predictor Added: study_quantile mse Value: 158.847796156561"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24 15  6  2
    ## [1] "Iteration: 28 Predictor Added: avg_household_size mse Value: 158.399711676517"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24 15  6  2 11
    ## [1] "Iteration: 29 Predictor Added: pct_non_white mse Value: 158.42849203918"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24 15  6  2 11 26
    ## [1] "Iteration: 30 Predictor Added: binned_inc_point mse Value: 158.76913039041"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24 15  6  2 11 26 29
    ## [1] "Iteration: 31 Predictor Added: pct_white mse Value: 157.663458935464"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24 15  6  2 11 26 29 23
    ## [1] "Iteration: 32 Predictor Added: poverty_percent mse Value: 158.728603186935"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  8 13 12  5 27 33 25 17
    ## [24] 24 15  6  2 11 26 29 23  7

``` r
f.aic.list <- f.subset(cancer.df, maxPreds = ncol(cancer.df) - 1, nfolds = 5, criterion = "aic")
```

    ## [1] "Iteration: 1 Predictor Added: avg_deaths_yr_pop aic Value: 22198.3407562201"
    ## [1] 30
    ## [1] "Iteration: 2 Predictor Added: median_age_female aic Value: 20684.9420954273"
    ## [1] 30 10
    ## [1] "Iteration: 3 Predictor Added: region aic Value: 20327.9257312174"
    ## [1] 30 10  3
    ## [1] "Iteration: 4 Predictor Added: incidence_rate aic Value: 20038.4854015885"
    ## [1] 30 10  3  4
    ## [1] "Iteration: 5 Predictor Added: pct_unemployed16_over aic Value: 19846.488243267"
    ## [1] 30 10  3  4 18
    ## [1] "Iteration: 6 Predictor Added: pct_hs18_24 aic Value: 19698.3346730745"
    ## [1] 30 10  3  4 18 14
    ## [1] "Iteration: 7 Predictor Added: pct_emp_priv_coverage aic Value: 19593.3096724274"
    ## [1] 30 10  3  4 18 14 20
    ## [1] "Iteration: 8 Predictor Added: pct_private_coverage aic Value: 19495.2064133987"
    ## [1] 30 10  3  4 18 14 20 19
    ## [1] "Iteration: 9 Predictor Added: binned_inc_ub aic Value: 19445.3468115051"
    ## [1] 30 10  3  4 18 14 20 19 28
    ## [1] "Iteration: 10 Predictor Added: median_age_male aic Value: 19408.153934971"
    ##  [1] 30 10  3  4 18 14 20 19 28  9
    ## [1] "Iteration: 11 Predictor Added: pct_hs25_over aic Value: 19378.579829986"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16
    ## [1] "Iteration: 12 Predictor Added: avg_ann_count_pop aic Value: 19352.6123558114"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31
    ## [1] "Iteration: 13 Predictor Added: pct_public_coverage aic Value: 19341.5338442161"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21
    ## [1] "Iteration: 14 Predictor Added: pct_public_coverage_alone aic Value: 19295.8578747603"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22
    ## [1] "Iteration: 15 Predictor Added: imp_pct_employed16_over aic Value: 19267.5143001291"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32
    ## [1] "Iteration: 16 Predictor Added: birth_rate aic Value: 19260.0722651453"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25
    ## [1] "Iteration: 17 Predictor Added: med_income aic Value: 19254.8189449111"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5
    ## [1] "Iteration: 18 Predictor Added: pct_bach_deg25_over aic Value: 19250.6131196436"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17
    ## [1] "Iteration: 19 Predictor Added: percent_married aic Value: 19245.1258802374"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12
    ## [1] "Iteration: 20 Predictor Added: imp_pct_private_coverage_alone aic Value: 19238.1718249234"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33
    ## [1] "Iteration: 21 Predictor Added: pop_est2015 aic Value: 19232.8664696153"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6
    ## [1] "Iteration: 22 Predictor Added: pct_non_white aic Value: 19233.189293977"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26
    ## [1] "Iteration: 23 Predictor Added: pct_white aic Value: 19215.7884199275"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [1] "Iteration: 24 Predictor Added: pct_no_hs18_24 aic Value: 19215.1012460789"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13
    ## [1] "Iteration: 25 Predictor Added: pct_bach_deg18_24 aic Value: 19214.5479816211"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13 15
    ## [1] "Iteration: 26 Predictor Added: pct_married_households aic Value: 19215.6288030833"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13 15 24
    ## [1] "Iteration: 27 Predictor Added: median_age aic Value: 19213.2528664267"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13 15 24  8
    ## [1] "Iteration: 28 Predictor Added: poverty_percent aic Value: 19216.9827864039"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13 15 24  8  7
    ## [1] "Iteration: 29 Predictor Added: binned_inc_lb aic Value: 19216.4912556002"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13 15 24  8  7 27
    ## [1] "Iteration: 30 Predictor Added: binned_inc_point aic Value: 19219.403509208"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13 15 24  8  7 27 29
    ## [1] "Iteration: 31 Predictor Added: avg_household_size aic Value: 19222.2515913101"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13 15 24  8  7 27 29 11
    ## [1] "Iteration: 32 Predictor Added: study_quantile aic Value: 19225.3088314616"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32 25  5 17 12 33  6 26 23
    ## [24] 13 15 24  8  7 27 29 11  2

``` r
f.bic.list <- f.subset(cancer.df, maxPreds = ncol(cancer.df) - 1, nfolds = 5, criterion = "bic")
```

    ## [1] "Iteration: 1 Predictor Added: avg_deaths_yr_pop bic Value: 22215.0033763571"
    ## [1] 30
    ## [1] "Iteration: 2 Predictor Added: median_age_female bic Value: 20707.8170403437"
    ## [1] 30 10
    ## [1] "Iteration: 3 Predictor Added: region bic Value: 20369.8067049957"
    ## [1] 30 10  3
    ## [1] "Iteration: 4 Predictor Added: incidence_rate bic Value: 20084.4463875334"
    ## [1] 30 10  3  4
    ## [1] "Iteration: 5 Predictor Added: pct_unemployed16_over bic Value: 19899.7166732842"
    ## [1] 30 10  3  4 18
    ## [1] "Iteration: 6 Predictor Added: pct_hs18_24 bic Value: 19759.0367010144"
    ## [1] 30 10  3  4 18 14
    ## [1] "Iteration: 7 Predictor Added: pct_emp_priv_coverage bic Value: 19656.9506398072"
    ## [1] 30 10  3  4 18 14 20
    ## [1] "Iteration: 8 Predictor Added: pct_private_coverage bic Value: 19562.6386938255"
    ## [1] 30 10  3  4 18 14 20 19
    ## [1] "Iteration: 9 Predictor Added: binned_inc_ub bic Value: 19518.4890122503"
    ## [1] 30 10  3  4 18 14 20 19 28
    ## [1] "Iteration: 10 Predictor Added: median_age_male bic Value: 19489.2515604482"
    ##  [1] 30 10  3  4 18 14 20 19 28  9
    ## [1] "Iteration: 11 Predictor Added: pct_hs25_over bic Value: 19465.1823400069"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16
    ## [1] "Iteration: 12 Predictor Added: avg_ann_count_pop bic Value: 19445.6206705907"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31
    ## [1] "Iteration: 13 Predictor Added: pct_public_coverage bic Value: 19440.2066653677"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21
    ## [1] "Iteration: 14 Predictor Added: pct_public_coverage_alone bic Value: 19396.9151522338"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22
    ## [1] "Iteration: 15 Predictor Added: imp_pct_employed16_over bic Value: 19378.7308346852"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32
    ## [1] "Iteration: 16 Predictor Added: med_income bic Value: 19374.6865453667"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5
    ## [1] "Iteration: 17 Predictor Added: birth_rate bic Value: 19374.7604777838"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25
    ## [1] "Iteration: 18 Predictor Added: percent_married bic Value: 19375.9883591105"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12
    ## [1] "Iteration: 19 Predictor Added: pct_bach_deg25_over bic Value: 19376.856207043"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17
    ## [1] "Iteration: 20 Predictor Added: pct_non_white bic Value: 19375.9386884897"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26
    ## [1] "Iteration: 21 Predictor Added: pct_white bic Value: 19364.3815218732"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23
    ## [1] "Iteration: 22 Predictor Added: imp_pct_private_coverage_alone bic Value: 19368.0424138424"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33
    ## [1] "Iteration: 23 Predictor Added: pop_est2015 bic Value: 19371.8294688709"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [1] "Iteration: 24 Predictor Added: pct_no_hs18_24 bic Value: 19375.9081460389"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13
    ## [1] "Iteration: 25 Predictor Added: median_age bic Value: 19381.7226294224"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13  8
    ## [1] "Iteration: 26 Predictor Added: pct_married_households bic Value: 19390.4305658794"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13  8 24
    ## [1] "Iteration: 27 Predictor Added: poverty_percent bic Value: 19394.5660859162"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13  8 24  7
    ## [1] "Iteration: 28 Predictor Added: avg_household_size bic Value: 19402.3168714622"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13  8 24  7 11
    ## [1] "Iteration: 29 Predictor Added: pct_bach_deg18_24 bic Value: 19406.0481583823"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13  8 24  7 11 15
    ## [1] "Iteration: 30 Predictor Added: binned_inc_lb bic Value: 19415.4249173136"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13  8 24  7 11 15 27
    ## [1] "Iteration: 31 Predictor Added: binned_inc_point bic Value: 19414.9057621687"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13  8 24  7 11 15 27 29
    ## [1] "Iteration: 32 Predictor Added: study_quantile bic Value: 19442.8979433062"
    ##  [1] 30 10  3  4 18 14 20 19 28  9 16 31 21 22 32  5 25 12 17 26 23 33  6
    ## [24] 13  8 24  7 11 15 27 29  2

``` r
#Show the 'best' final selection with minimal MSE (takes in a list object from the previous function, with criterion (MSE/AIC))
present.fs.result <- function(f.result.list, criterion) {
lm.fs.preds <- with(f.result.list, pred.list[[which.min(result.mat[,2])]]) #Pick out indices from best model
fs.mse <- with(f.result.list, result.mat[which.min(result.mat[,2]), 2])
print(sprintf("The best predictor set of %s predictors, out of a max of %s, (%s = %s)",
              length(lm.fs.preds), f.result.list$maxPreds, criterion, round(fs.mse, 3)))
print(names(cancer.df[,c(lm.fs.preds)]))
}

#Show the final results and subset selections, with CV criterion stats
present.fs.result(f.mse.list, "MSE")
```

    ## [1] "The best predictor set of 24 predictors, out of a max of 32, (MSE = 157.104)"
    ##  [1] "avg_deaths_yr_pop"              "median_age_female"             
    ##  [3] "region"                         "incidence_rate"                
    ##  [5] "pct_unemployed16_over"          "pct_hs18_24"                   
    ##  [7] "pct_emp_priv_coverage"          "pct_private_coverage"          
    ##  [9] "binned_inc_ub"                  "median_age_male"               
    ## [11] "pct_hs25_over"                  "avg_ann_count_pop"             
    ## [13] "pct_public_coverage"            "pct_public_coverage_alone"     
    ## [15] "imp_pct_employed16_over"        "median_age"                    
    ## [17] "pct_no_hs18_24"                 "percent_married"               
    ## [19] "med_income"                     "binned_inc_lb"                 
    ## [21] "imp_pct_private_coverage_alone" "birth_rate"                    
    ## [23] "pct_bach_deg25_over"            "pct_married_households"

``` r
present.fs.result(f.aic.list, "AIC")
```

    ## [1] "The best predictor set of 27 predictors, out of a max of 32, (AIC = 19213.253)"
    ##  [1] "avg_deaths_yr_pop"              "median_age_female"             
    ##  [3] "region"                         "incidence_rate"                
    ##  [5] "pct_unemployed16_over"          "pct_hs18_24"                   
    ##  [7] "pct_emp_priv_coverage"          "pct_private_coverage"          
    ##  [9] "binned_inc_ub"                  "median_age_male"               
    ## [11] "pct_hs25_over"                  "avg_ann_count_pop"             
    ## [13] "pct_public_coverage"            "pct_public_coverage_alone"     
    ## [15] "imp_pct_employed16_over"        "birth_rate"                    
    ## [17] "med_income"                     "pct_bach_deg25_over"           
    ## [19] "percent_married"                "imp_pct_private_coverage_alone"
    ## [21] "pop_est2015"                    "pct_non_white"                 
    ## [23] "pct_white"                      "pct_no_hs18_24"                
    ## [25] "pct_bach_deg18_24"              "pct_married_households"        
    ## [27] "median_age"

``` r
present.fs.result(f.bic.list, "BIC")
```

    ## [1] "The best predictor set of 21 predictors, out of a max of 32, (BIC = 19364.382)"
    ##  [1] "avg_deaths_yr_pop"         "median_age_female"        
    ##  [3] "region"                    "incidence_rate"           
    ##  [5] "pct_unemployed16_over"     "pct_hs18_24"              
    ##  [7] "pct_emp_priv_coverage"     "pct_private_coverage"     
    ##  [9] "binned_inc_ub"             "median_age_male"          
    ## [11] "pct_hs25_over"             "avg_ann_count_pop"        
    ## [13] "pct_public_coverage"       "pct_public_coverage_alone"
    ## [15] "imp_pct_employed16_over"   "med_income"               
    ## [17] "birth_rate"                "percent_married"          
    ## [19] "pct_bach_deg25_over"       "pct_non_white"            
    ## [21] "pct_white"

``` r
#Pull the indices for the selected models, for later comparison
fs.mse.preds <- with(f.mse.list, pred.list[[which.min(result.mat[,2])]])
fs.aic.preds <- with(f.aic.list, pred.list[[which.min(result.mat[,2])]])
fs.bic.preds <- with(f.bic.list,  pred.list[[which.min(result.mat[,2])]])
```

If you want to repeat the process to get a feel for what you think might actually be the best of the 'best' subset selections (they very slightly from iteration to iteration by cross validation). Pick the number of iterations `len` and let the chunk run, you can see the process work (reccomend commenting out the print() calls in the functions above to reduce all the extra output), then at the end it will print all the subsets selected and you can sift through which variables you want for a final subset model. See which ones get selected most often basically.

``` r
#Repeat if you want to
#Number of repetitions
len <- 10
f.result.list <- list()

#Repeat Subset Selection algorithm
for (i in 1:len) {
  f.result.list[[i]] <- f.subset(sub.cancer.df, maxPreds = 10, nfolds = 5, criterion = "mse")
}

#View the results
for (i in 1:len) {
present.fs.result(f.result.list[[i]], "MSE")
}

#See which vars seem to be selected most often, make a subjective call what subset you like best.

#fs.lm <- lm(target_death_rate ~ ., data = cancer.df[,c(lm.fs.preds,1)])
```

Repeat the process for Backwards Subset
=======================================

NOTE this will take much longer to run
--------------------------------------

##### Again must reorder data frame so response is the first column of the data frame (see above code chunk)

##### Again state variable removed, too many factor levels and CV process can fail (re add later to see if it's useful)

-   The way the algorithm works is

1.  Starts with an exhaustive set of current predictors in the data frame input (all preds in the data frame)
2.  Starts with a empty set of predictors removed
3.  Loops through all current preds and removes one at a time, splits the data into kfolds test/train (CV), fits lm() with current preds, stores MSE
4.  Selects the set whose CV MSE was smallest
5.  The optimal set had one of the predictors removed, selects that predictor as 'worst; removes it from current predictors
6.  Stores that predictor set as worst and records the indices and CV MSE in a matrix for the model without that pred
7.  Repeats 3.- 5. until you have met the minimum number of predictors you specified
8.  Returns a list with all 'best' predictor sets at each iteration, matrix of results, and prespecified max preds
9.  Last few lines of code output the optimal (minimum MSE) predictor set

*notes - CV is set to 5 for speed, changing to 10 will take ~1.5 times as long to run* *notes - the print() lines within the function are so you can visualize how the algorithm is working, they can be commented out to reduce clutter when you want to run the selection function iteratively (avoids reams of output)*

``` r
#Backward Selection
#set.seed(4444)  #Set seed for reproducibility
###############################################

#Num all possible preds
maxPreds <- ncol(cancer.df - 1)

#Function
#Inputs:
#full data set (whatever that ends up being, ideally non-singular)
#minimum number of preds you want the algorithm to reduce subsets to
#Number of folds for CV step
#Criterion to minimize: "mse" or "aic"

b.subset <- function(sub.cancer.df, minPreds = 1, nfolds = 5, criterion = "mse"){

## Allthe predictors (their indices).
allPreds <- 1:ncol(sub.cancer.df)
allPreds <- allPreds[-1] #take out response

#current set of preds (starts full), and available
currPreds <- allPreds
availPreds <- allPreds
#Record the min errors
minError <- c()
#Set up storage objects
pred.list <- list()
result.mat <- matrix(nrow = ncol(sub.cancer.df) - 1, ncol = 2)
#initialize iteration
i <- 1
#Forward selection loop
while (length(currPreds) >= minPreds) {
  ##remove predictor which has CV lm with lowest MSE/AIC (as determined by CV)
  ## The MSE/AIC computed as we remove of the available predictors
  allError <- c()
  for (id in availPreds) {
    data.df <- sub.cancer.df[,c(currPreds[-id],1)]  #Take data frame without the predictor were removing
    if (criterion == "mse") {
      error <- mseCV(data.df, nfolds)  #calculate CV mse
    } else if (criterion == "aic") {
      error <- aicCV(data.df, nfolds)  #calculate CV aic
    } else if (criterion == "bic") {
      error <- bicCV(data.df, nfolds)  # CV bic
    } else {
      print("Invalid criterion")
      stop("Wrong criterion input")
    }
    allError <- c(allError,error)  #Collect all our criterions (error is misleading, can be mse or aic)
  }
  ##Find the min
  id <- which.min(allError)
  ##get the worst predictor and MSE
  worstPred <- availPreds[id]
  bestError <- min(allError)
  ##Remove these from the collection
  currPreds <- currPreds[-id]
  minError <- c(minError, bestError)
  availPreds <- currPreds
  ## Print stuff out for debugging and attention-grabbing
  print(sprintf("Predictor Removed: %s  %s Value: %s",names(sub.cancer.df[,worstPred]), criterion, bestError))
  print(currPreds)
  result.mat[i,] <- c(worstPred,bestError)   #All print() can be commented out at any time, just to visualize alg. process
  pred.list[[i]] <- currPreds
  
  i <- i + 1
  }
  list(pred.list = pred.list, result.mat = result.mat)
}

#Call functions, output is a list of preds and a result matrix of criterion and predictor removed
b.mse.list <- b.subset(cancer.df, minPreds = 10, nfolds = 5, criterion = "mse")
```

    ## [1] "Predictor Removed: poverty_percent  mse Value: 157.446998590786"
    ##  [1]  2  3  4  5  6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [24] 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: pct_public_coverage_alone  mse Value: 157.048286112719"
    ##  [1]  2  3  4  5  6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26
    ## [24] 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: med_income  mse Value: 161.107162277694"
    ##  [1]  2  3  4  6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26 27
    ## [24] 28 29 30 31 32 33
    ## [1] "Predictor Removed: binned_inc_ub  mse Value: 161.221538609294"
    ##  [1]  2  3  4  6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26 27
    ## [24] 29 30 31 32 33
    ## [1] "Predictor Removed: pct_no_hs18_24  mse Value: 161.384376758388"
    ##  [1]  2  3  4  6  8  9 10 11 12 14 15 16 17 18 19 20 21 23 24 25 26 27 29
    ## [24] 30 31 32 33
    ## [1] "Predictor Removed: avg_deaths_yr_pop  mse Value: 161.330171553721"
    ##  [1]  2  3  4  6  8  9 10 11 12 14 15 16 17 18 19 20 21 23 24 25 26 27 29
    ## [24] 31 32 33
    ## [1] "Predictor Removed: pct_married_households  mse Value: 357.573620331749"
    ##  [1]  2  3  4  6  8  9 10 11 12 14 15 16 17 18 19 20 21 23 25 26 27 29 31
    ## [24] 32 33
    ## [1] "Predictor Removed: median_age  mse Value: 363.994807855301"
    ##  [1]  2  3  4  6  9 10 11 12 14 15 16 17 18 19 20 21 23 25 26 27 29 31 32
    ## [24] 33
    ## [1] "Predictor Removed: imp_pct_private_coverage_alone  mse Value: 365.770067904924"
    ##  [1]  2  3  4  6  9 10 11 12 14 15 16 17 18 19 20 21 23 25 26 27 29 31 32
    ## [1] "Predictor Removed: avg_ann_count_pop  mse Value: 366.029778073094"
    ##  [1]  2  3  4  6  9 10 11 12 14 15 16 17 18 19 20 21 23 25 26 27 29 32
    ## [1] "Predictor Removed: birth_rate  mse Value: 365.721982066546"
    ##  [1]  2  3  4  6  9 10 11 12 14 15 16 17 18 19 20 21 23 26 27 29 32
    ## [1] "Predictor Removed: binned_inc_lb  mse Value: 365.138857067602"
    ##  [1]  2  3  4  6  9 10 11 12 14 15 16 17 18 19 20 21 23 26 29 32
    ## [1] "Predictor Removed: median_age_female  mse Value: 370.789712408652"
    ##  [1]  2  3  4  6  9 11 12 14 15 16 17 18 19 20 21 23 26 29 32
    ## [1] "Predictor Removed: pct_public_coverage  mse Value: 369.841008008041"
    ##  [1]  2  3  4  6  9 11 12 14 15 16 17 18 19 20 23 26 29 32
    ## [1] "Predictor Removed: median_age_male  mse Value: 371.322455351234"
    ##  [1]  2  3  4  6 11 12 14 15 16 17 18 19 20 23 26 29 32
    ## [1] "Predictor Removed: pct_hs18_24  mse Value: 371.866926485572"
    ##  [1]  2  3  4  6 11 12 15 16 17 18 19 20 23 26 29 32
    ## [1] "Predictor Removed: pct_non_white  mse Value: 373.484313865993"
    ##  [1]  2  3  4  6 11 12 15 16 17 18 19 20 23 29 32
    ## [1] "Predictor Removed: pct_hs25_over  mse Value: 377.795544754047"
    ##  [1]  2  3  4  6 11 12 15 17 18 19 20 23 29 32
    ## [1] "Predictor Removed: binned_inc_point  mse Value: 382.794831329532"
    ##  [1]  2  3  4  6 11 12 15 17 18 19 20 23 32
    ## [1] "Predictor Removed: pct_bach_deg18_24  mse Value: 384.071333294169"
    ##  [1]  2  3  4  6 11 12 17 18 19 20 23 32
    ## [1] "Predictor Removed: avg_household_size  mse Value: 383.489086975768"
    ##  [1]  2  3  4  6 12 17 18 19 20 23 32
    ## [1] "Predictor Removed: percent_married  mse Value: 383.950472836303"
    ##  [1]  2  3  4  6 17 18 19 20 23 32
    ## [1] "Predictor Removed: pct_emp_priv_coverage  mse Value: 384.417796504025"
    ## [1]  2  3  4  6 17 18 19 23 32

``` r
b.aic.list <- b.subset(cancer.df, minPreds = 10, nfolds = 5, criterion = "aic")
```

    ## [1] "Predictor Removed: median_age_female  aic Value: 19219.210782943"
    ##  [1]  2  3  4  5  6  7  8  9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [24] 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: pct_no_hs18_24  aic Value: 19399.9825323619"
    ##  [1]  2  3  4  5  6  7  8  9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26
    ## [24] 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: avg_ann_count_pop  aic Value: 19400.3807393835"
    ##  [1]  2  3  4  5  6  7  8  9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26
    ## [24] 27 28 29 30 32 33
    ## [1] "Predictor Removed: poverty_percent  aic Value: 19427.6913347221"
    ##  [1]  2  3  4  5  6  8  9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    ## [24] 28 29 30 32 33
    ## [1] "Predictor Removed: med_income  aic Value: 19428.5699521429"
    ##  [1]  2  3  4  6  8  9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
    ## [24] 29 30 32 33
    ## [1] "Predictor Removed: pct_white  aic Value: 19430.4610370696"
    ##  [1]  2  3  4  6  8  9 11 12 14 15 16 17 18 19 20 21 22 24 25 26 27 28 29
    ## [24] 30 32 33
    ## [1] "Predictor Removed: imp_pct_private_coverage_alone  aic Value: 19455.4925463237"
    ##  [1]  2  3  4  6  8  9 11 12 14 15 16 17 18 19 20 21 22 24 25 26 27 28 29
    ## [24] 30 32
    ## [1] "Predictor Removed: pct_public_coverage  aic Value: 19457.9295607132"
    ##  [1]  2  3  4  6  8  9 11 12 14 15 16 17 18 19 20 22 24 25 26 27 28 29 30
    ## [24] 32
    ## [1] "Predictor Removed: binned_inc_point  aic Value: 19604.5192307353"
    ##  [1]  2  3  4  6  8  9 11 12 14 15 16 17 18 19 20 22 24 25 26 27 28 30 32
    ## [1] "Predictor Removed: pct_hs25_over  aic Value: 19603.8351407264"
    ##  [1]  2  3  4  6  8  9 11 12 14 15 17 18 19 20 22 24 25 26 27 28 30 32
    ## [1] "Predictor Removed: pct_non_white  aic Value: 19610.8597610751"
    ##  [1]  2  3  4  6  8  9 11 12 14 15 17 18 19 20 22 24 25 27 28 30 32
    ## [1] "Predictor Removed: binned_inc_lb  aic Value: 19623.6897673599"
    ##  [1]  2  3  4  6  8  9 11 12 14 15 17 18 19 20 22 24 25 28 30 32
    ## [1] "Predictor Removed: pct_bach_deg18_24  aic Value: 19623.5664220443"
    ##  [1]  2  3  4  6  8  9 11 12 14 17 18 19 20 22 24 25 28 30 32
    ## [1] "Predictor Removed: pct_hs18_24  aic Value: 19625.7184577284"
    ##  [1]  2  3  4  6  8  9 11 12 17 18 19 20 22 24 25 28 30 32
    ## [1] "Predictor Removed: avg_deaths_yr_pop  aic Value: 19698.3333552655"
    ##  [1]  2  3  4  6  8  9 11 12 17 18 19 20 22 24 25 28 32
    ## [1] "Predictor Removed: pct_married_households  aic Value: 21365.2982159228"
    ##  [1]  2  3  4  6  8  9 11 12 17 18 19 20 22 25 28 32
    ## [1] "Predictor Removed: pop_est2015  aic Value: 21411.8217913938"
    ##  [1]  2  3  4  8  9 11 12 17 18 19 20 22 25 28 32
    ## [1] "Predictor Removed: pct_private_coverage  aic Value: 21414.3442922993"
    ##  [1]  2  3  4  8  9 11 12 17 18 20 22 25 28 32
    ## [1] "Predictor Removed: birth_rate  aic Value: 21420.4457726375"
    ##  [1]  2  3  4  8  9 11 12 17 18 20 22 28 32
    ## [1] "Predictor Removed: incidence_rate  aic Value: 21426.9477662911"
    ##  [1]  2  3  8  9 11 12 17 18 20 22 28 32
    ## [1] "Predictor Removed: pct_public_coverage_alone  aic Value: 21995.6394867955"
    ##  [1]  2  3  8  9 11 12 17 18 20 28 32
    ## [1] "Predictor Removed: binned_inc_ub  aic Value: 22055.0373646202"
    ##  [1]  2  3  8  9 11 12 17 18 20 32
    ## [1] "Predictor Removed: region  aic Value: 22053.8138549103"
    ## [1]  2  8  9 11 12 17 18 20 32

``` r
b.bic.list <- b.subset(cancer.df, minPreds = 10, nfolds = 5, criterion = "bic")
```

    ## [1] "Predictor Removed: poverty_percent  bic Value: 19437.2525146574"
    ##  [1]  2  3  4  5  6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [24] 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: median_age_male  bic Value: 19427.4708869411"
    ##  [1]  2  3  4  5  6  8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
    ## [24] 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: percent_married  bic Value: 19437.6561134773"
    ##  [1]  2  3  4  5  6  8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    ## [24] 28 29 30 31 32 33
    ## [1] "Predictor Removed: median_age  bic Value: 19444.3388439238"
    ##  [1]  2  3  4  5  6 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
    ## [24] 29 30 31 32 33
    ## [1] "Predictor Removed: median_age_female  bic Value: 19436.0668526705"
    ##  [1]  2  3  4  5  6 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
    ## [24] 30 31 32 33
    ## [1] "Predictor Removed: incidence_rate  bic Value: 19982.730627855"
    ##  [1]  2  3  5  6 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
    ## [24] 31 32 33
    ## [1] "Predictor Removed: region  bic Value: 20294.5703800557"
    ##  [1]  2  5  6 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
    ## [24] 32 33
    ## [1] "Predictor Removed: study_quantile  bic Value: 20383.8391388234"
    ##  [1]  5  6 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
    ## [24] 33
    ## [1] "Predictor Removed: pct_married_households  bic Value: 20366.5796232144"
    ##  [1]  5  6 11 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: pop_est2015  bic Value: 20371.4050608852"
    ##  [1]  5 11 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: pct_hs18_24  bic Value: 20372.072572546"
    ##  [1]  5 11 13 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: med_income  bic Value: 20375.0892748918"
    ##  [1] 11 13 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: pct_emp_priv_coverage  bic Value: 20374.3621747442"
    ##  [1] 11 13 15 16 17 18 19 21 22 23 25 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: pct_private_coverage  bic Value: 20394.1077629478"
    ##  [1] 11 13 15 16 17 18 21 22 23 25 26 27 28 29 30 31 32 33
    ## [1] "Predictor Removed: avg_deaths_yr_pop  bic Value: 20413.7719808343"
    ##  [1] 11 13 15 16 17 18 21 22 23 25 26 27 28 29 31 32 33
    ## [1] "Predictor Removed: pct_bach_deg18_24  bic Value: 22172.9779233212"
    ##  [1] 11 13 16 17 18 21 22 23 25 26 27 28 29 31 32 33
    ## [1] "Predictor Removed: imp_pct_private_coverage_alone  bic Value: 22172.6740151859"
    ##  [1] 11 13 16 17 18 21 22 23 25 26 27 28 29 31 32
    ## [1] "Predictor Removed: pct_public_coverage  bic Value: 22177.2968238244"
    ##  [1] 11 13 16 17 18 22 23 25 26 27 28 29 31 32
    ## [1] "Predictor Removed: pct_no_hs18_24  bic Value: 22181.8578937501"
    ##  [1] 11 16 17 18 22 23 25 26 27 28 29 31 32
    ## [1] "Predictor Removed: pct_unemployed16_over  bic Value: 22229.1282159039"
    ##  [1] 11 16 17 22 23 25 26 27 28 29 31 32
    ## [1] "Predictor Removed: avg_household_size  bic Value: 22260.4927199041"
    ##  [1] 16 17 22 23 25 26 27 28 29 31 32
    ## [1] "Predictor Removed: binned_inc_ub  bic Value: 22263.044875725"
    ##  [1] 16 17 22 23 25 26 27 29 31 32
    ## [1] "Predictor Removed: imp_pct_employed16_over  bic Value: 22262.0817884106"
    ## [1] 16 17 22 23 25 26 27 29 31

``` r
#Visualize function takes in a list object from the last function
present.bs.result <- function(result.list, criterion) {
lm.bs.preds <- with(result.list, pred.list[[which.min(result.mat[,2])]])
lm.bs.mse <- with(result.list, result.mat[which.min(result.mat[,2]), 2])
print(sprintf("The best predictor set of %s predictors, out of a max of %s, (%s = %s)",
              length(lm.bs.preds), maxPreds, criterion, round(lm.bs.mse, 3)))
names(cancer.df[,c(lm.bs.preds)])
}

#Present the final subsets which were selected
present.bs.result(b.mse.list, "MSE")
```

    ## [1] "The best predictor set of 30 predictors, out of a max of 33, (MSE = 157.048)"

    ##  [1] "study_quantile"                 "region"                        
    ##  [3] "incidence_rate"                 "med_income"                    
    ##  [5] "pop_est2015"                    "median_age"                    
    ##  [7] "median_age_male"                "median_age_female"             
    ##  [9] "avg_household_size"             "percent_married"               
    ## [11] "pct_no_hs18_24"                 "pct_hs18_24"                   
    ## [13] "pct_bach_deg18_24"              "pct_hs25_over"                 
    ## [15] "pct_bach_deg25_over"            "pct_unemployed16_over"         
    ## [17] "pct_private_coverage"           "pct_emp_priv_coverage"         
    ## [19] "pct_public_coverage"            "pct_white"                     
    ## [21] "pct_married_households"         "birth_rate"                    
    ## [23] "pct_non_white"                  "binned_inc_lb"                 
    ## [25] "binned_inc_ub"                  "binned_inc_point"              
    ## [27] "avg_deaths_yr_pop"              "avg_ann_count_pop"             
    ## [29] "imp_pct_employed16_over"        "imp_pct_private_coverage_alone"

``` r
present.bs.result(b.aic.list, "AIC")
```

    ## [1] "The best predictor set of 31 predictors, out of a max of 33, (AIC = 19219.211)"

    ##  [1] "study_quantile"                 "region"                        
    ##  [3] "incidence_rate"                 "med_income"                    
    ##  [5] "pop_est2015"                    "poverty_percent"               
    ##  [7] "median_age"                     "median_age_male"               
    ##  [9] "avg_household_size"             "percent_married"               
    ## [11] "pct_no_hs18_24"                 "pct_hs18_24"                   
    ## [13] "pct_bach_deg18_24"              "pct_hs25_over"                 
    ## [15] "pct_bach_deg25_over"            "pct_unemployed16_over"         
    ## [17] "pct_private_coverage"           "pct_emp_priv_coverage"         
    ## [19] "pct_public_coverage"            "pct_public_coverage_alone"     
    ## [21] "pct_white"                      "pct_married_households"        
    ## [23] "birth_rate"                     "pct_non_white"                 
    ## [25] "binned_inc_lb"                  "binned_inc_ub"                 
    ## [27] "binned_inc_point"               "avg_deaths_yr_pop"             
    ## [29] "avg_ann_count_pop"              "imp_pct_employed16_over"       
    ## [31] "imp_pct_private_coverage_alone"

``` r
present.bs.result(b.bic.list, "BIC")
```

    ## [1] "The best predictor set of 30 predictors, out of a max of 33, (BIC = 19427.471)"

    ##  [1] "study_quantile"                 "region"                        
    ##  [3] "incidence_rate"                 "med_income"                    
    ##  [5] "pop_est2015"                    "median_age"                    
    ##  [7] "median_age_female"              "avg_household_size"            
    ##  [9] "percent_married"                "pct_no_hs18_24"                
    ## [11] "pct_hs18_24"                    "pct_bach_deg18_24"             
    ## [13] "pct_hs25_over"                  "pct_bach_deg25_over"           
    ## [15] "pct_unemployed16_over"          "pct_private_coverage"          
    ## [17] "pct_emp_priv_coverage"          "pct_public_coverage"           
    ## [19] "pct_public_coverage_alone"      "pct_white"                     
    ## [21] "pct_married_households"         "birth_rate"                    
    ## [23] "pct_non_white"                  "binned_inc_lb"                 
    ## [25] "binned_inc_ub"                  "binned_inc_point"              
    ## [27] "avg_deaths_yr_pop"              "avg_ann_count_pop"             
    ## [29] "imp_pct_employed16_over"        "imp_pct_private_coverage_alone"

``` r
#Pull the indices for the selected models, for later comparison
bs.mse.preds <- with(b.mse.list, pred.list[[which.min(result.mat[,2])]])
bs.aic.preds <- with(b.aic.list, pred.list[[which.min(result.mat[,2])]])
bs.bic.preds <- with(b.bic.list,  pred.list[[which.min(result.mat[,2])]])
```

In general, backwards selection selects larger models than forawrds. I would suggest only using forward selection in this case to choose a smaller model. Or stepwise based on some other criterion.

Repeat to get a good feel for which preds get selected most often (note this will take a long time if the original p-set you are using is large, wouldnt use the full set for these iterations). This can be useful to get a feel for what you think might actually be the best of the 'best' subset selections (they very slightly from iteration to iteration by cross validation). Pick the number of iterations `len` and let the chunk run, you can see the process work (reccomend commenting out the print() calls in the functions above to reduce all the extra output), then at the end it will print all the subsets selected and you can sift through which variables you want for a final subset model. See which ones get selected most often basically.

Ridge and Lasso Functions for Comparison
========================================

``` r
#Initialize Sample Size
sampleSize <- nrow(cancer.df)

#Initialize lambda grid
lambda.grid <- 10^seq(-3,2,length = 200)

#Lambda optimization Ridge
ridge.opt <- function(data, grid = lambda.grid, response = "target_death_rate"){
    #Snag the index for the matrices  
  index <- which(names(data) == response)
  
  #Response 
  Y <- as.matrix(data[ ,index])
  
  #Design matrix
  X <- data[ ,-index] %>%
    names() %>% 
    paste(., collapse = "+") %>%    
    paste("~ ", .) %>%
    formula() %>%
    model.matrix(., data)
  X <- X[,-1]
  
  #Lambda optimize
  cv.ridge <- cv.glmnet(X,Y,alpha = 0,intercept = T, lambda = grid, family = "gaussian")
  cv.ridge$lambda.min
}

lambda.opt.ridge <- ridge.opt(cancer.df, grid = lambda.grid)

#Lambda optimization Lasso
lasso.opt <- function(data, grid = lambda.grid, response = "target_death_rate"){  
  #Snag the index for matrices
  index <- which(names(data) == response)
  
  #Response 
  Y <- as.matrix(data[, index])
  
  #Design matrix
  X <- data[ ,-index] %>%
    names() %>% 
    paste(., collapse = "+") %>%    
    paste("~ ", .) %>%
    formula() %>%
    model.matrix(.,data)
  X <- X[,-1]  
  
  #Optimize Lambda
  cv.lasso <- cv.glmnet(X,Y,alpha = 1,intercept = T,lambda = grid, family = "gaussian")
  cv.lasso$lambda.min
}

lambda.opt.lasso <- lasso.opt(cancer.df, grid = lambda.grid)

ridge <- function(data, response = "target_death_rate") {

  #Snag the index for the matrices  
  index <- which(names(data) == response)
  
  #Response 
  Y <- as.matrix(data[ ,index])
  
  #Design matrix
  X <- data[ ,-index] %>%
    names() %>% 
    paste(., collapse = "+") %>%    
    paste("~ ", .) %>%
    formula() %>%
    model.matrix(., data)
  X <- X[,-1]
  
  #Lambda optimize
  cv.ridge <- cv.glmnet(X,Y,alpha = 0,intercept = T, lambda = lambda.grid, family = "gaussian")
  lambda.opt.ridge <- cv.ridge$lambda.min
  
  #Return model
  return(glmnet(X,Y,alpha = 0,intercept = T,lambda = lambda.opt.ridge, family = "gaussian"))
}

lasso <- function(data, response = "target_death_rate", print = FALSE) {
  
  #Snag the index for matrices
  index <- which(names(data) == response)
  
  #Response 
  Y <- as.matrix(data[, index])
  
  #Design matrix
  X <- data[ ,-index] %>%
    names() %>% 
    paste(., collapse = "+") %>%    
    paste("~ ", .) %>%
    formula() %>%
    model.matrix(.,data)
  X <- X[,-1]  
  
  #Optimize Lambda
  cv.lasso <- cv.glmnet(X,Y,alpha = 1,intercept = T,lambda = lambda.grid, family = "gaussian")
  lambda.opt.lasso <- cv.lasso$lambda.min
  mod.lasso <- glmnet(X,Y,alpha = 1,intercept = T,lambda = lambda.opt.lasso,family = "gaussian")
  #Print if true
  if (isTRUE(print)) {
  coef.lasso <- coef(mod.lasso)[,1]
    print("The Lasso method selected the variables:")
    print(paste(names(which(coef.lasso != 0))))
  }
  return(mod.lasso)
}

glmnet.mse <- function(model, data, response = "target_death_rate") {
  #Snag the index for matrices
  index <- which(names(data) == response)
  #Design matrix
  X1 <- data[ ,-index] %>%
    names() %>% 
    paste(., collapse = "+") %>%    
    paste("~ ", .) %>%
    formula() %>%
    model.matrix(., data)
  X1 <- X1[,-1] 
  preds <- predict(model, newx = X1, type = "response")
  return(with(data, mean((target_death_rate - preds)^2)))
}

lassoCV <- function(data.df = cancer.df, kfolds = 10, response = "target_death_rate"){
  index <- which(names(data.df) == response) #Grab response index for matrices
  folds <- sample(1:kfolds, sampleSize, rep = T) #Randomly sample indeices for kfolds
  error <- rep(0, kfolds) #initialize error vector
  for (k in 1:kfolds) {
    train.df <- data.df[folds != k,] #Split into test train
    test.df <- data.df[folds == k,]
    X1 <- test.df[ , -index] %>%
          names() %>% 
          paste(., collapse = "+") %>%    
          paste("~ ", .) %>%
          formula() %>%
          model.matrix(.,test.df)
    X1 <- X1[,-1]
     
    lasso.mod <- lasso(train.df, response)
       preds <- predict(lasso.mod, newx = X1, type = "response")
    error[k] <- with(test.df, mean((target_death_rate - preds)^2))
  }
  mean(error) #return kfold mean mse
}

ridgeCV <- function(data.df, kfolds = 10, response){
  index <- which(names(data.df) == response) #Grab response index for matrices
  folds <- sample(1:kfolds, sampleSize, rep = T) #Randomly sample indeices for kfolds
  error <- rep(0, kfolds) #initialize error vector
  for (k in 1:kfolds) {
    train.df <- data.df[folds != k,] #Split into test train
     test.df <- data.df[folds == k,]
    #Design test matrix
  X1 <- test.df[ , -index] %>%
  names() %>% 
  paste(., collapse = "+") %>%    
  paste("~ ", .) %>%
  formula() %>%
  model.matrix(.,test.df)
  #Remove Intercept
   X1 <- X1[ ,-1] 
  #Fit model, predict on kth fold, record mse
   ridge.mod <- ridge(train.df, response)
       preds <- predict(ridge.mod, newx = X1, type = "response")
    error[k] <- with(test.df, mean((target_death_rate - preds)^2))
  }
  mean(error) #return kfold mean mse
}

lasso1 <- cancer.df %>% lasso(., print = TRUE)
```

    ## [1] "The Lasso method selected the variables:"
    ##  [1] "(Intercept)"                    "study_quantileLow"             
    ##  [3] "study_quantileModerate"         "study_quantileHigh"            
    ##  [5] "study_quantileVery High"        "region2"                       
    ##  [7] "region3"                        "region4"                       
    ##  [9] "incidence_rate"                 "med_income"                    
    ## [11] "pop_est2015"                    "poverty_percent"               
    ## [13] "median_age"                     "median_age_male"               
    ## [15] "median_age_female"              "avg_household_size"            
    ## [17] "percent_married"                "pct_no_hs18_24"                
    ## [19] "pct_hs18_24"                    "pct_bach_deg18_24"             
    ## [21] "pct_hs25_over"                  "pct_bach_deg25_over"           
    ## [23] "pct_unemployed16_over"          "pct_private_coverage"          
    ## [25] "pct_emp_priv_coverage"          "pct_public_coverage"           
    ## [27] "pct_public_coverage_alone"      "pct_white"                     
    ## [29] "pct_married_households"         "birth_rate"                    
    ## [31] "pct_non_white"                  "binned_inc_ub"                 
    ## [33] "binned_inc_point"               "avg_deaths_yr_pop"             
    ## [35] "avg_ann_count_pop"              "imp_pct_employed16_over"       
    ## [37] "imp_pct_private_coverage_alone"

``` r
coef(lasso1)
```

    ## 38 x 1 sparse Matrix of class "dgCMatrix"
    ##                                          s0
    ## (Intercept)                    173.42792657
    ## study_quantileLow                1.24562764
    ## study_quantileModerate          -0.22753729
    ## study_quantileHigh              -0.80262350
    ## study_quantileVery High          0.89436343
    ## region2                          9.49796761
    ## region3                          1.57935943
    ## region4                          2.11486849
    ## incidence_rate                   4.68388960
    ## med_income                       1.89548896
    ## pop_est2015                     -0.33713730
    ## poverty_percent                  0.04389737
    ## median_age                       0.01337026
    ## median_age_male                 -3.19517566
    ## median_age_female              -11.74035772
    ## avg_household_size               0.07672839
    ## percent_married                 -2.14210316
    ## pct_no_hs18_24                   0.57694662
    ## pct_hs18_24                      2.50640143
    ## pct_bach_deg18_24                0.09313643
    ## pct_hs25_over                    1.02295940
    ## pct_bach_deg25_over             -1.82411775
    ## pct_unemployed16_over            3.67887331
    ## pct_private_coverage            -1.12674886
    ## pct_emp_priv_coverage            3.37072041
    ## pct_public_coverage            -10.76114146
    ## pct_public_coverage_alone        9.41557963
    ## pct_white                       -2.74324632
    ## pct_married_households           0.92404581
    ## birth_rate                      -0.78617512
    ## pct_non_white                   -2.91479668
    ## binned_inc_lb                    .         
    ## binned_inc_ub                    1.87719723
    ## binned_inc_point                 0.03773051
    ## avg_deaths_yr_pop               29.51533089
    ## avg_ann_count_pop               -1.38388588
    ## imp_pct_employed16_over         -1.74148965
    ## imp_pct_private_coverage_alone   0.79775501

``` r
#test - lasso, ridge - all good
#ridge(data = cancer.df)
#lasso(data = cancer.df, print = TRUE)

#test - CV - all good
lassoCV(cancer.df, 10, "target_death_rate")
```

    ## [1] 158.4419

``` r
ridgeCV(cancer.df, 10, "target_death_rate")
```

    ## [1] 157.8064

Model Comparison with Iterative 10-Fold CV
==========================================

Takes a bit to run.

``` r
#Number of CV iterations
nrep <- 100

#make sure target_death_rate is the first index in the data set (just way easier that way) District of columbia is a problem
#State is a problem, need to make a region variable
library(modelr)
cv.df <- cancer.df %>%
  mutate(median_age_all_gender = (median_age_male + median_age_female)/2) %>% 
  mutate(south = ifelse(region == "South", 1, 0)) %>%
         crossv_mc(., nrep) %>%
  mutate(train = map(train, as.tibble),
         test = map(test, as.tibble)) %>%
  mutate(fs.mse.lm = map(train, ~lm(target_death_rate ~ ., data = .x[, c(fs.mse.preds, 1)])),
         fs.aic.lm = map(train, ~lm(target_death_rate ~ ., data = .x[, c(fs.aic.preds, 1)])),
         fs.bic.lm = map(train, ~lm(target_death_rate ~ ., data = .x[, c(fs.bic.preds, 1)])),
         bs.mse.lm = map(train, ~lm(target_death_rate ~ ., data = .x[, c(bs.mse.preds, 1)])),
         bs.aic.lm = map(train, ~lm(target_death_rate ~ ., data = .x[, c(bs.aic.preds, 1)])),
         bs.bic.lm = map(train, ~lm(target_death_rate ~ ., data = .x[, c(bs.bic.preds, 1)])),
         lasso.lm  = map(train, ~lasso(data = .x)),
         ridge.lm  = map(train, ~ridge(data = .x)),
         auto.back  = map(train, ~lm(target_death_rate ~ pct_white + pct_hs25_over + 
          imp_pct_employed16_over + pct_private_coverage + pct_public_coverage + 
          avg_ann_count_pop + avg_deaths_yr_pop + med_income + study_quantile + 
          median_age_all_gender + south, data = .x)),
         auto.front  = map(train, ~lm(target_death_rate ~ pct_white + pct_hs25_over +
          imp_pct_employed16_over + pct_private_coverage + pct_public_coverage + 
          avg_ann_count_pop + avg_deaths_yr_pop + med_income + study_quantile +
          median_age_all_gender + south + avg_household_size, data = .x))) %>%
  mutate(fs.mse.mse = map2_dbl(fs.mse.lm, test, ~mse(model = .x, data = .y)),
         fs.aic.mse = map2_dbl(fs.aic.lm, test, ~mse(model = .x, data = .y)),
         fs.bic.mse = map2_dbl(fs.bic.lm, test, ~mse(model = .x, data = .y)),
         bs.mse.mse = map2_dbl(bs.mse.lm, test, ~mse(model = .x, data = .y)),
         bs.aic.mse = map2_dbl(bs.aic.lm, test, ~mse(model = .x, data = .y)),
         bs.bic.mse = map2_dbl(bs.bic.lm, test, ~mse(model = .x, data = .y)),
         lasso.mse  = map2_dbl( lasso.lm, test, ~glmnet.mse(model = .x, data = .y)),
         ridge.mse  = map2_dbl( lasso.lm, test, ~glmnet.mse(model = .x, data = .y)),
         auto.front.mse = map2_dbl(auto.back, test, ~mse(model = .x, data = .y)),
         auto.back.mse = map2_dbl(auto.front, test, ~mse(model = .x, data = .y)))
```

Visualize Model Comparison
==========================

``` r
#Violin box plots
violin.mse <- cv.df %>% 
  dplyr::select(ends_with(".mse")) %>% 
  gather(key = model, value = mse) %>% 
  mutate(model = str_replace(model, "mse", ""),
         model = fct_reorder(model, mse, .desc = FALSE, .fun = mean)) %>% 
  ggplot(aes(x = model, y = mse)) + 
  geom_violin(aes(fill = model), trim = FALSE, alpha = 0.3) + 
  geom_boxplot(width = 0.25) +
  labs(
    y = "CV Mean Sqaure Error",
    x = "Model",
    title = sprintf("Model MSE Comparison by 10-Fold CV: %i Iterations", nrep)
  ) +
  viridis::scale_fill_viridis(
    option = "magma",
    name = "MSE",
    begin = 1,
    end = 0,
    discrete = TRUE) 

density.mse <- cv.df %>%
   dplyr::select(ends_with(".mse")) %>% 
  gather(key = model, value = mse) %>% 
  mutate(model = str_replace(model, ".mse", ""),
         model = fct_reorder(model, mse, .desc = FALSE, .fun = mean)) %>% 
  ggplot(aes(x = mse, colour = model, fill = model)) + 
  geom_density(position = "stack", alpha = 0.3) + 
  #geom_boxplot(width = 0.25) +
  labs(
    y = "CV Mean Sqaure Error",
    x = "Model",
    title = sprintf("Model MSE Comparison by 10-Fold CV: %i Iterations", nrep)
  ) +
  viridis::scale_fill_viridis(
    option = "magma",
    name = "MSE",
    begin = 1,
    end = 0,
    discrete = TRUE) +
    viridis::scale_colour_viridis(
    option = "magma",
    name = "MSE",
    begin = 1,
    end = 0,
    discrete = TRUE)
```

<img src="Q_subset_selection_alg_files/figure-markdown_github/unnamed-chunk-15-1.png" width="90%" style="display: block; margin: auto;" /><img src="Q_subset_selection_alg_files/figure-markdown_github/unnamed-chunk-15-2.png" width="90%" style="display: block; margin: auto;" />

*This will take a long time to run*

``` r
#Repeat if you want to
#Number of repetitions
len <- 10
b.result.list <- list()

#Repeat Subset Selection algorithm
for (i in 1:len) {
  b.result.list[[i]] <- b.subset(sub.cancer.df, minPreds = 30, nfolds = 5, "mse")
}

#View the results
for (i in 1:len) {
present.bs.result(b.result.list[[i]], "MSE")
}

#See which vars seem to be selected most often, make a subjective call what subset you like best.

#fs.lm <- lm(target_death_rate ~ ., data = cancer.df[,c(lm.fs.preds,1)])
```
