# data wrangling to connect environmental and MRI data
# 
# Max Korbmacher, July 2025
#
# ------------------------------------------------- #
# -----------------Contents------------------------ #
# ------------------------------------------------- #
# 1. Prep------------------------------------------ #
# 2. Global Correlates----------------------------- #
# 2.1 Cross-sectional------------------------------ #
# 2.2 Longitudinal--------------------------------- #
# 3. Regional Correlates--------------------------- #
# ------------------------------------------------- #
# ------------------------------------------------- #
#
#
#
# 1.Prep-------------------------------------------

# 1.1 Packages-------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, dplyr, ggplot2, reshape2,lme4,ggpubr,
               psych, lmerTest, performance, broom, tidyr,
               gtsummary,ggrepel
               )

# 1.2 Load data------------------------------------

# Specify the data path
savepath = "/cluster/projects/p33/users/maxk/UKB/environment/results/"
datapath = "/cluster/projects/p33/users/maxk/UKB/environment/data/"

cross_env = read.csv(paste(datapath,"final_cross.csv",sep=""))
long_env = read.csv(paste(datapath,"final_long.csv",sep=""))

# 1.3 Demographics ------
# remove 47 Bristol subs
cross_env = cross_env %>% filter(site != "Bristol")

# define outcomes and predictors (necessary for analyses)
outcomes = cross_env %>% select(FA_Mean, MD_Mean, RD_Mean, AD_Mean,
                                #                     "v_intra_Mean" ,          "v_extra_Mean"        ,   "v_csf_Mean"     ,        "micro_Rd_Mean"   ,       "micro_FA_Mean",         
                                #                     "micro_Ax_Mean"      ,    "micro_ADC_Mean"      ,   "Dax_intra_Mean"     ,    "Dax_extra_Mean",
                                CortexVol,CorticalWhiteMatterVol) %>% names


predictors=cross_env %>% dplyr::select(ends_with("at_age"),"RDS1", "N1", 
                                       "ANX", "ADHD", "ASD", "BIP", "MDD", "OCD", "SCZ", "AD") %>% names
# make summary table
pgs_cols = c("ANX", "ADHD", "ASD", "BIP", "MDD", "OCD", "SCZ", "AD")
for (col in pgs_cols) {
  new_col <- paste0(col, "_Z")
  cross_env[[new_col]] <- as.vector(scale(cross_env[[col]]))
}
Zs =  c("ANX_Z", "ADHD_Z", "ASD_Z", "BIP_Z", "MDD_Z", "OCD_Z", "SCZ_Z", "AD_Z")
tbl = cross_env %>% select(outcomes, Zs, RDS1, N1, age, sex,SurfaceHoles, EstimatedTotalIntraCranialVol, income, site) %>% tbl_summary(by=site)
tbl |> as_gt() |> gt::gtsave(filename = paste(savepath,"Demographics.html",sep="")) # use extensions .png, .html, .docx, .rtf, .tex, .ltx
tbl |> as_gt() |> gt::gtsave(filename = paste(savepath,"Demographics.tex",sep="")) # use extensions .png, .html, .docx, .rtf, .tex, .ltx
#
#
# outcomes = long_env %>% select(ends_with("_Mean"),CortexVol,CorticalWhiteMatterVol,BrainAge) %>% names
demo = read.csv(paste(datapath,"demo/environment.csv",sep=""))
demo = rbind(demo %>% select(eid,X738.2.0) %>% rename(income = X738.2.0),demo %>% select(eid, X738.3.0)%>%rename(income=X738.3.0))
long_env = merge(long_env,demo,by="eid")

# The same for long data
long_env=long_env %>% filter(site != "Bristol")
time2 = long_env %>% filter(session>1)
# define outcomes
outcomes = c(outcomes,"BrainAge")
predictors=c("si10","d2m","t2m","u10","v10","sp","avg_slhtf","avg_snswrf","avg_sdlwrf","avg_lsprate","uvb","avg_sdirswrf","avg_sduvrf","avg_snlwrf","tp","RDS", "N", 
             "ANX", "ADHD", "ASD", "BIP", "MDD", "OCD", "SCZ", "AD")
# 
# make summary table
for (col in pgs_cols) {
  new_col <- paste0(col, "_Z")
  time2[[new_col]] <- as.vector(scale(time2[[col]]))
}

tbl = time2 %>% select(outcomes, Zs, RDS, N, age, sex,SurfaceHoles, EstimatedTotalIntraCranialVol, income, site) %>% tbl_summary(by=site)
tbl |> as_gt() |> gt::gtsave(filename = paste(savepath,"Demographics_long.html",sep="")) # use extensions .png, .html, .docx, .rtf, .tex, .ltx
tbl |> as_gt() |> gt::gtsave(filename = paste(savepath,"Demographics-long.tex",sep="")) # use extensions .png, .html, .docx, .rtf, .tex, .ltx

# 1.4 Reduce data: PCA FOR LONGITUDINAL ANALYSES!----------------------------------
print("#")
print("#")
print("#")
print("#")
print("#")
print("#")
print("#")
print("#")
print("#")
print("#")
print("#")
print("THIS SECTION IS BEING SKIPPED AT THE MOMENT!!!!")
print("#")
print("#")
print("#")
print("#")
print("#")
print("#")
print("#")

if (FALSE){
  print("#####################")
  print("#####################")
  print("First: PGRS.")
  print("PGRS are not strongly related.")
  cor(pgrs[!pgrs$FID %in% long_env$eid,]%>%select(-X,-FID))
  print("PCA not necessary.")
  print("#####################")
  print("Second: Mental Health.")
  print("Only two scores. Move on,")
  print("#####################")
  print("Third: brain stats.")
  print("Quite strongly correlated. At least for multiple metrics")
  cor(dMRI_cross[!dMRI_cross$eid %in% long_env$eid,]%>%select(ends_with("Mean")))
  
  dMRI_pca_model_train <- prcomp(dMRI_cross[!dMRI_cross$eid %in% long_env$eid,]%>%select(-eid,-sex,-site,-age,-X), scale = TRUE, center = TRUE)
  
  T1w_pca_model_train <- prcomp(T1w[!T1w$eid %in% long_env$eid,]%>%
                                  select(-eid, 
                                         -Left.WM.hypointensities, -Right.WM.hypointensities,
                                         -Left.non.WM.hypointensities,-Right.non.WM.hypointensities), 
                                scale = TRUE, center = TRUE)
  print("check results")
  summary(T1w_pca_model_train)
  plot(T1w_pca_model_train, type = "l", main = "Scree Plot of Principal Components")
  print("Two components for T1w data!")
  
  summary(dMRI_pca_model_train)
  plot(dMRI_pca_model_train, type = "l", main = "Scree Plot of Principal Components")
  print("We should keep 3 components for dMRI data.")
  
  print("##########")
  print("Finally, let's try and combine dMRI and T1w data.")
  MRI = merge(dMRI_cross,T1w, by="eid")
  MRI_pca_model_train <- prcomp(MRI[!MRI$eid %in% long_env$eid,]%>%select(-eid,-sex,
                                                                          -Left.WM.hypointensities, -Right.WM.hypointensities,
                                                                          -Left.non.WM.hypointensities,-Right.non.WM.hypointensities,
                                                                          -site,-age,-X), scale = TRUE, center = TRUE)
  plot(MRI_pca_model_train, type = "l", main = "Scree Plot of Principal Components")
  print("3 components required.")
  chosen_num_components = 3
  print("##########")
  print("##########")
  print("##########")
  
  # Apply function to see top loading vars
  # This assumes hat pca_model_train is already fitted and chosen_num_components is set
  # from previous steps.
  
  # --- 1. Prepare Loadings for Analysis ---
  loadings_df_dMRI <- as.data.frame(dMRI_pca_model_train$rotation) %>%
    tibble::rownames_to_column(var = "Brain_Variable") %>%
    select(Brain_Variable, paste0("PC", 1:chosen_num_components))
  loadings_df_T1w <- as.data.frame(T1w_pca_model_train$rotation) %>%
    tibble::rownames_to_column(var = "Brain_Variable") %>%
    select(Brain_Variable, paste0("PC", 1:chosen_num_components))
  loadings_df_MRI <- as.data.frame(MRI_pca_model_train$rotation) %>%
    tibble::rownames_to_column(var = "Brain_Variable") %>%
    select(Brain_Variable, paste0("PC", 1:chosen_num_components))
  
  # --- 2. Examine Top Loadings for Each Component Systematically ---
  top_loadings = function(loadings_df){
    # We'll set a threshold for "significant" loadings
    loading_threshold <- 0.1 # Adjust this based on your data and desired stringency
    # Common values are 0.1, 0.2, or even 0.3 for very large datasets
    
    cat("\n--- Detailed Loadings for Chosen Components (Absolute Value >", loading_threshold, ") ---\n")
    for (i in 1:chosen_num_components) {
      pc_name <- paste0("PC", i)
      cat(paste0("\nComponent: ", pc_name, "\n"))
      
      # Filter for variables with absolute loading above the threshold for the current PC
      component_loadings <- loadings_df %>%
        select(Brain_Variable, !!sym(pc_name)) %>% # Select current PC's loadings
        filter(abs(.[[pc_name]]) >= loading_threshold) %>% # Filter by absolute threshold
        arrange(desc(abs(.[[pc_name]]))) # Sort by absolute loading magnitude
      
      if (nrow(component_loadings) > 0) {
        # Print a limited number of top variables, or all if less than 20
        print_limit <- min(20, nrow(component_loadings)) # Print top 20 or fewer
        
        cat(paste0("  Variables with |loading| >= ", loading_threshold, " (Top ", print_limit, " shown):\n"))
        for (j in 1:print_limit) {
          cat(sprintf("    - %-25s: %8.3f\n",
                      component_loadings$Brain_Variable[j],
                      component_loadings[[pc_name]][j]))
        }
      } else {
        cat("  No variables met the loading threshold for this component.\n")
      }
      #cat(paste0("  -> Tentative Name for ", pc_name, ": [Your Interpretation Here]\n"))
    }
  }
  top_loadings(loadings_df_dMRI)# >> First component = WM structural integrity, Second = Microstructural integrity, Third = extracellular diffusivity
  top_loadings(loadings_df_T1w) # >> First component = brain tissue, Second component = ventricles and CSF
  top_loadings(loadings_df_MRI) # These components do not seem to be well explained
  
  
  
  
  
  
  
  
  
  
  
  
  # next step: check whether there are some anatomical landmarks.
  
}
























# --- 3. (Conceptual) Grouping by Anatomical Regions ---
# This part is highly conceptual as it depends on how your 1865 variables are named
# or if you have a separate mapping file.

# Assuming you have a mapping where each Brain_Mxxx corresponds to a Region
# Example mapping (you'd need to create this for all 1865 vars)
# brain_region_map <- data.frame(
#   Brain_Variable = paste0("Brain_M", 1:1865),
#   Region = sample(c("Frontal_Lobe", "Temporal_Lobe", "Parietal_Lobe",
#                     "Occipital_Lobe", "Cerebellum", "Subcortical"),
#                   1865, replace = TRUE)
# )
#
# # Join loadings with region map
# loadings_with_regions <- loadings_df %>%
#   left_join(brain_region_map, by = "Brain_Variable")
#
# # Calculate average/sum of loadings by region for each component
# cat("\n--- Average Absolute Loadings by Region ---\n")
# for (i in 1:chosen_num_components) {
#   pc_name <- paste0("PC", i)
#   avg_abs_loadings_by_region <- loadings_with_regions %>%
#     group_by(Region) %>%
#     summarise(Avg_Abs_Loading = mean(abs(!!sym(pc_name)))) %>%
#     arrange(desc(Avg_Abs_Loading))
#
#   cat(paste0("\nComponent: ", pc_name, "\n"))
#   print(avg_abs_loadings_by_region)
# }




# 2. Correlates------------------------------------
# 2.1 Cross-sectional------------------------------
# 2.1.1 Simple correlation matrix
# cormat = cor(Filter(is.numeric, long_env %>% filter(session==1) %>%
#                       select(-eid,-age,-birthyear,-birthmonth, -session, -sex)),use = "pairwise.complete.obs")
# melted_cormat = melt(cormat)
# plot1 = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile()+  scale_fill_gradient2(low = "blue",mid = "white", high = "red") +
#   theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# cormat2 = cor(Filter(is.numeric, long_env %>% filter(session==2) %>%
#                       select(-eid,-age,-birthyear,-birthmonth, -session, -sex)),use = "pairwise.complete.obs")
# melted_cormat = melt(cormat)
# plot2 = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile()+  scale_fill_gradient2(low = "blue",mid = "white", high = "red") +
#   theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# ggarrange(plot1, plot2)
# 
# 2.1.1 Linear models accounting for covariates and testing for income as a strong predictor ----
results <- list()

#outcomes = cross_env %>% select(ends_with("_Mean"),CortexVol,CorticalWhiteMatterVol) %>% names

# Standardize ONLY numeric predictors/outcomes (safe for factor predictors)
cross_scaled <- cross_env %>% mutate(across(all_of(c(predictors, outcomes)), ~ if (is.numeric(.x)) as.numeric(scale(.x)) else .x))
cross_scaled$site = cross_env$site
results <- list()
lrt_results <- list()

for (outcome in outcomes) {
  for (predictor in predictors) {
    if (predictor == outcome) next
    
    vars_needed <- c(outcome, predictor,
                     "sex", "age", "SurfaceHoles", "EstimatedTotalIntraCranialVol", "income", "site")
    
    dat_sub <- cross_env %>%
      select(all_of(vars_needed)) %>%
      filter(complete.cases(.)) %>%
      mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))
    
    if (nrow(dat_sub) < 10) next
    if (length(unique(dat_sub$income)) < 2) next
    
    # ------------------------
    # Models
    # ------------------------
    full_f <- paste0(outcome, " ~ ", predictor,
                     " + sex + age + SurfaceHoles + EstimatedTotalIntraCranialVol + income + site")
    model_full <- lm(as.formula(full_f), data = dat_sub)
    
    reduced_income_f <- paste0(outcome, " ~ ", predictor,
                               " + sex + age + SurfaceHoles + EstimatedTotalIntraCranialVol + site")
    model_reduced_income <- lm(as.formula(reduced_income_f), data = dat_sub)
    
    td <- tidy(model_full, conf.int = TRUE)
    
    # ------------------------
    # Predictor
    # ------------------------
    pred_exact <- td %>% filter(term == predictor)
    if (nrow(pred_exact) == 1) {
      std_beta_predictor <- pred_exact$estimate
      conf_low_predictor <- pred_exact$conf.low
      conf_high_predictor <- pred_exact$conf.high
      p_value_predictor <- pred_exact$p.value
      predictor_omnibus_p <- NA_real_
    } else {
      reduced_no_pred_f <- paste0(outcome, " ~ sex + age + SurfaceHoles + EstimatedTotalIntraCranialVol + income + site")
      model_reduced_no_pred <- lm(as.formula(reduced_no_pred_f), data = dat_sub)
      an_pred <- anova(model_reduced_no_pred, model_full)
      predictor_omnibus_p <- an_pred$`Pr(>F)`[2]
      
      std_beta_predictor <- NA_real_
      conf_low_predictor <- NA_real_
      conf_high_predictor <- NA_real_
      p_value_predictor <- NA_real_
    }
    
    # ------------------------
    # Likelihood ratio test (income)
    # ------------------------
    an_income <- anova(model_reduced_income, model_full)
    
    r2_full <- summary(model_full)$r.squared
    r2_reduced <- summary(model_reduced_income)$r.squared
    delta_r2 <- r2_full - r2_reduced
    
    lrt_results[[paste0(outcome, "_by_", predictor)]] <- data.frame(
      outcome   = outcome,
      predictor = predictor,
      n_obs     = nrow(dat_sub),
      df        = an_income$Df[2],
      f_value   = an_income$F[2],
      p_value   = an_income$`Pr(>F)`[2],
      r2_full   = r2_full,
      r2_reduced = r2_reduced,
      delta_r2   = delta_r2,
      stringsAsFactors = FALSE
    )
    
    # ------------------------
    # Store predictor + income together
    # ------------------------
    # results[[paste0(outcome, "_by_", predictor)]] <- income_results %>%
    #   mutate(
    #     std_beta_predictor = std_beta_predictor,
    #     conf_low_predictor = conf_low_predictor,
    #     conf_high_predictor = conf_high_predictor,
    #     p_value_predictor = p_value_predictor,
    #     predictor_omnibus_p = predictor_omnibus_p
    #   )
    # Only store in results if predictor is NOT income
    if (!grepl("income", predictor, ignore.case = TRUE)) {
      results[[paste0(outcome, "_", predictor)]] <- data.frame(
        outcome = outcome,
        predictor = predictor,
        beta = std_beta_predictor,
        conf_low <- conf_low_predictor,
        conf_high <- conf_high_predictor,
        p_value <- p_value_predictor
      )
    }
  }
}
cross_res <- bind_rows(results)
#
#
#
#
#
########### OLD
# 
# # Loop over outcome and predictor pairs
# # Create a standardized version of the data
# cross_scaled <- cross_env %>%
#   mutate(across(all_of(c(predictors, outcomes)), scale))
# # Loop through all (outcome, predictor) pairs
# for (outcome in outcomes) {
#   for (predictor in predictors) {
#     if (predictor != outcome) {
#       
#       # Fit linear model
#       formula_str <- paste0(outcome, " ~ ", predictor," + sex + age + SurfaceHoles + EstimatedTotalIntraCranialVol + income + site")
#       model <- lm(as.formula(formula_str), data = cross_scaled)
#       
#       # Tidy output and extract relevant row
#       coef_df <- tidy(model, conf.int = TRUE) %>%
#         filter(term == predictor)
#       
#       # Store result
#       results[[paste0(outcome, "_by_", predictor)]] <- data.frame(
#         outcome   = outcome,
#         predictor = predictor,
#         std_beta  = coef_df$estimate,
#         conf_low  = coef_df$conf.low,
#         conf_high = coef_df$conf.high,
#         p_value   = coef_df$p.value
#       )
#     }
#   }
# }

# corrected p-vals

cross_res$p.correct = p.adjust(cross_res$p_value....p_value_predictor,method="BH")
#cross_res = cross_res %>% filter(p.correct<.05)

# summarise effects by predictor
cross_res %>% group_by(predictor) %>% summarize(M = mean(abs(beta)), SD = sd(abs(beta)), Md = median(abs(beta)), MAD = mad(abs(beta)))

cross_res <- cross_res %>%
  mutate(sig_label = ifelse(p.correct < 0.05, "*", ""))

# label data 
cross_res$outcome = factor(cross_res$outcome)
cross_res$predictor = factor(cross_res$predictor)
# levels(cross_res$outcome) = c("DTI-AD","DKI-AK","WMTI-AWF","CGMV","WMV",
#                              "BRIA-DAXextra","BRIA-DAXintra","DTI-FA","DTI-MD", 
#                              "BRIA-microADC","BRIA-microAX","BRIA-microFA","BRIA-microRD", 
#                              "DKI-MK", "WMTI-radEAD", "DTI-RD", "DKI-RK", 
#                              "SMT-FA", "SMT-long", "SMTmc-Diff", "SMTmc-extraMD", 
#                              "SMTmc-extratrans", "SMTmc-intra", "SMT-MD", "SMT-trans", 
#                              "BRIA-vCSF", "BRIA-vExtra", "BRIA-vIntra") 
levels(cross_res$outcome) = c("DTI-AD","CGMV","WMV",
                              #"BRIA-DAXextra","BRIA-DAXintra",
                              "DTI-FA","DTI-MD", 
                              #"BRIA-microADC","BRIA-microAX","BRIA-microFA","BRIA-microRD", 
                              "DTI-RD") 
                              #"BRIA-vCSF", "BRIA-vExtra", "BRIA-vIntra") 
# CGMV = Cortical grey matter volume, WMV = Cortical white matter volume
levels(cross_res$predictor) = c("PGRS: Alzheimer's Disease",
                                "PGRS: Attention-Deficit/Hyperactivity Disorder",
                                "PGRS: Anxiety Disorder",
                                "PGRS: Autism Spectrum Disorder",
                                "Weather: Time-mean large-scale precipitation rate", 
                               "Weather: Time-mean surface direct short-wave radiation flux",
                               "Weather: Time-mean surface downward long-wave radiation flux",
                               "Weather: Time-mean surface downward UV radiation flux",
                               "Weather: Time-mean surface latent heat flux",
                               "Weather: Time-mean surface net long-wave radiation flux",
                               "Weather: Time-mean surface net short-wave radiation flux",
                               
                               "PGRS: Bipolar Disorder",
                               
                               "Weather: Surface net short-wave radiation flux",
                               
                               "PGRS: Major Depressive Disorder",
                               
                               "Self-Reports: Neuroticism",
                               
                               "PGRS: Obsessive-Compulsive Disorder",
                               
                               "Self-Reports: Recent Depressive Symptoms",
                               
                               "PGRS: Schizophrenia",
                               
                               "Weather: 10 metre wind speed",
                               "Weather: Surface pressure",
                               "Weather: 2 metre temperature",
                               "Weather: Total precipitation",
                               "Weather: 10 metre U wind component",
                               "Weather: Surface downward UV radiation",
                               "Weather: 10 metre V wind component"
)


# plot
p_cross = ggplot(cross_res, aes(x = outcome, y = predictor, fill = beta)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = sig_label), size = 6, color = "black") +
  scale_fill_gradient2(
    low = "#a6bddb", mid = "white", high = "#fbb4b9",  # pastel blue to white to pastel red
    midpoint = 0,
    name = "Std. Beta",
    limits = c(min(cross_res$beta), max(cross_res$beta))  # consistent scale
  ) +
  #scale_x_discrete(labels = outcome_labels) +
  #scale_y_discrete(labels = predictor_labels) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Outcome",
    y = "Predictor"
    #,title = "Standardized Associations between Climate and Brain Variables"
  ) +
  xlim(sort(levels(cross_res$outcome))) + ylim(sort(levels(cross_res$predictor)))

# save figure and table
ggsave(filename = "/tsd/p33/home/p33-maxk/p1.pdf",plot = p_cross, width = 10, height = 8)
ggsave(filename = paste(savepath,"Associations.pdf",sep=""),plot = p_cross, width = 10, height = 8)

#

write.csv(file = paste(savepath,"correlates.csv",sep=""),cross_res)
#
#
# Show some descriptives of the correlational analysis
filter(cross_res, grepl("Weather",predictor)) %>% summarise(Md = median(abs(beta)), MAD = mad(abs(beta)))
filter(cross_res, grepl("PGRS",predictor)) %>% summarise(Md = median(abs(beta)), MAD = mad(abs(beta)))
filter(cross_res, grepl("Self-Reports",predictor)) %>% summarise(Md = median(abs(beta)), MAD = mad(abs(beta)))
#
#
# 2.1.2 LRTs---------------------------------
# Goal here is to compare a comparison model containing each predictor of interest with the respective
## re-define predictors
predictors=cross_env %>% dplyr::select(ends_with("at_age"),"RDS1", "N1", 
                                       "ANX", "ADHD", "ASD", "BIP", "MDD", "OCD", "SCZ", "AD") %>% names

## Base names of weather predictors
weather_base <- c("si10","d2m","t2m","u10","v10","sp","avg_slhtf","avg_snswrf",
                  "avg_sdlwrf","avg_lsprate","uvb","avg_sdirswrf","avg_sduvrf",
                  "avg_snlwrf","tp")

## Add the "_at_age" suffix
weather_vars <- paste0(weather_base, "_at_age")

## Everything else from predictors is "other"
other_vars <- setdiff(predictors, weather_vars)

lrt_weather <- list()

## loop trough the list comparing models with or without added PGRS or self-reported scores.
for (outcome in outcomes) {
  for (weather in weather_vars) {
    for (other in other_vars) {
      
      vars_needed <- c(outcome, weather, other,
                       "sex", "age", "SurfaceHoles", "EstimatedTotalIntraCranialVol", "income")
      
      dat_sub <- cross_env %>%
        select(all_of(vars_needed)) %>%
        filter(complete.cases(.)) %>%
        mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))
      
      if (nrow(dat_sub) < 10) next
      if (length(unique(dat_sub$income)) < 2) next
      
      # ------------------------
      # Models
      # ------------------------
      base_f <- paste0(outcome, " ~ ", weather,
                       " + sex + age + SurfaceHoles + EstimatedTotalIntraCranialVol + income")
      extended_f <- paste0(base_f, " + ", other)
      
      model_base <- lm(as.formula(base_f), data = dat_sub)
      model_ext  <- lm(as.formula(extended_f), data = dat_sub)
      
      # ------------------------
      # Likelihood ratio test
      # ------------------------
      an_lrt <- anova(model_base, model_ext)
      
      r2_base <- summary(model_base)$r.squared
      r2_ext  <- summary(model_ext)$r.squared
      delta_r2 <- r2_ext - r2_base
      
      lrt_weather[[paste(outcome, weather, other, sep = "_")]] <- data.frame(
        outcome   = outcome,
        weather   = weather,
        added_var = other,
        n_obs     = nrow(dat_sub),
        df        = an_lrt$Df[2],
        f_value   = an_lrt$F[2],
        p_value   = an_lrt$`Pr(>F)`[2],
        r2_base   = r2_base,
        r2_ext    = r2_ext,
        delta_r2  = delta_r2,
        stringsAsFactors = FALSE
      )
    }
  }
}

final_weather_lrt <- bind_rows(lrt_weather)
final_weather_lrt$p_corrected = p.adjust(final_weather_lrt$p_value, method = "BH")

#
# write LRT table
write.csv(file = "/cluster/projects/p33/users/maxk/UKB/environment/results/LRT.csv",final_weather_lrt)
#
#
# sum up the added variance explained of variables of interest
final_weather_lrt %>% group_by(added_var) %>% summarize(M = mean(abs(delta_r2)), SD = sd(abs(delta_r2)), Md = median(abs(delta_r2)), MAD = mad(abs(delta_r2)))
final_weather_lrt %>% filter(p_corrected < 0.05)

# # 2.2 Longitudinal---------------------------------
# Split into two dataframes and correlate
# Split by session
df_t1 <- long_env %>% filter(session == 1)  # or "T1", "baseline", etc.
df_t2 <- long_env %>% filter(session == 2)  # or "T2", "followup", etc.

# Select only numeric columns
numeric_vars <- long_env %>% 
  select(where(is.numeric)) %>% 
  select(-session, -eid) %>%  # Remove ID and session
  names()

# Correlate numeric variables between sessions
correlations <- sapply(numeric_vars, function(var) {
  # Match by ID
  merged <- inner_join(
    df_t1 %>% select(eid, !!sym(var)) %>% rename(var_t1 = !!sym(var)),
    df_t2 %>% select(eid, !!sym(var)) %>% rename(var_t2 = !!sym(var)),
    by = "eid"
  )
  
  cor(merged$var_t1, merged$var_t2, use = "complete.obs")
})

# View results
print("Correlations between T1 and T2:")
print(round(correlations, 3))




# Loop over outcome and predictor combinations
associations =
  function(data_frame){
    results <- list()
    for (outcome in outcomes) {
      for (predictor in predictors) {
        if (predictor != outcome) {
          # Standardize predictor and outcome within each time point
          df0 = data_frame %>% dplyr::select(sex,age)
          df_scaled = data_frame %>% dplyr::select(-sex,-age)
          df_scaled <- df_scaled %>%
            group_by(session) %>%
            mutate(
              !!outcome := scale(.data[[outcome]])[,1],
              !!predictor := scale(.data[[predictor]])[,1]
            ) %>%
            ungroup()
          df_scaled = cbind(df_scaled,df0)
          # Fit the mixed model
          formula_str <- paste0(outcome, " ~ ", predictor, " + sex + age +  SurfaceHoles + EstimatedTotalIntraCranialVol + income + (1|eid)")
          model <- lmer(as.formula(formula_str), data = df_scaled)

          # Extract summary and CI
          coef_summary <- summary(model)$coefficients
          ci <- suppressMessages(confint(model, parm = predictor, method = "Wald"))

          # Extract results
          std_beta <- coef_summary[predictor, "Estimate"]
          p_val    <- coef_summary[predictor, "Pr(>|t|)"]

          # Store results
          results[[paste0(outcome, "_by_", predictor)]] <- data.frame(
            outcome   = outcome,
            predictor = predictor,
            std_beta  = std_beta,
            conf_low  = ci[1],
            conf_high = ci[2],
            p_value   = p_val
          )
        }
      }
    }
    return(results)
  }
res = associations(long_env)
# Combine all results into one dataframe
final_df <- bind_rows(res)

# View result
print(final_df)
# corrected p-vals
final_df$p.correct = p.adjust(final_df$p_value,method="BH")
#final_df = final_df %>% filter(p.correct<.05)

# check which predictors
levels(factor(final_df$predictor))

# summarise effects by predictor
final_df %>% group_by(predictor) %>% summarize(M = mean(abs(std_beta)), SD = sd(abs(std_beta)), Md = median(abs(std_beta)), MAD = mad(abs(std_beta)))


# largest effects detected for avg_lsprate, d2m, sp, t2m, tp
#
# avg_lsprate = Time-mean large-scale precipitation rate
# d2m = 2 metre dewpoint temperature
# sp = surface pressure
# t2m = 2 metre temperature
# tp = total precipation

# show only largest effects
final_df %>% filter(abs(std_beta) > 0.05) %>% arrange(predictor)

final_df <- final_df %>%
  mutate(sig_label = ifelse(p.correct < 0.05, "*", ""))

final_df$outcome = factor(final_df$outcome)
final_df$predictor = factor(final_df$predictor)
levels(final_df$outcome) = c("DTI-AD","Brain Age","CGMV","WMV",
                              #"BRIA-DAXextra","BRIA-DAXintra",
                              "DTI-FA","DTI-MD", 
                              #"BRIA-microADC","BRIA-microAX","BRIA-microFA","BRIA-microRD", 
                              "DTI-RD") 
#"BRIA-vCSF", "BRIA-vExtra", "BRIA-vIntra") 
# CGMV = Cortical grey matter volume, WMV = Cortical white matter volume
levels(final_df$predictor) = c("PGRS: Alzheimer's Disease",
                                "PGRS: Attention-Deficit/Hyperactivity Disorder",
                                "PGRS: Anxiety Disorder",
                                "PGRS: Autism Spectrum Disorder",
                                "Weather: Time-mean large-scale precipitation rate", 
                                "Weather: Time-mean surface direct short-wave radiation flux",
                                "Weather: Time-mean surface downward long-wave radiation flux",
                                "Weather: Time-mean surface downward UV radiation flux",
                                "Weather: Time-mean surface latent heat flux",
                                "Weather: Time-mean surface net long-wave radiation flux",
                                "Weather: Time-mean surface net short-wave radiation flux",
                                
                                "PGRS: Bipolar Disorder",
                                
                                "Weather: Surface net short-wave radiation flux",
                                
                                "PGRS: Major Depressive Disorder",
                                
                                "Self-Reports: Neuroticism",
                                
                                "PGRS: Obsessive-Compulsive Disorder",
                                
                                "Self-Reports: Recent Depressive Symptoms",
                                
                                "PGRS: Schizophrenia",
                                
                                "Weather: 10 metre wind speed",
                                "Weather: Surface pressure",
                                "Weather: 2 metre temperature",
                                "Weather: Total precipitation",
                                "Weather: 10 metre U wind component",
                                "Weather: Surface downward UV radiation",
                                "Weather: 10 metre V wind component"
)

# Tile plot
p2 = ggplot(final_df, aes(x = outcome, y = predictor, fill = std_beta)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = sig_label), size = 6, color = "black") +
  scale_fill_gradient2(
    low = "#a6bddb", mid = "white", high = "#fbb4b9",  # pastel blue to white to pastel red
    midpoint = 0,
    name = "Std. Beta",
    limits = c(min(final_df$std_beta), max(final_df$std_beta))  # consistent scale
  ) +
  #scale_x_discrete(labels = outcome_labels) +
  #scale_y_discrete(labels = predictor_labels) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Outcome",
    y = "Predictor"
    #title = "Standardized Associations between Climate and Brain Variables"
  ) + xlim(sort(levels(final_df$outcome))) + ylim(sort(levels(final_df$predictor)))
ggsave(filename = "/tsd/p33/home/p33-maxk/p1.pdf",plot = p1, width = 10, height = 8)
ggsave(filename = paste(savepath,"Longitudinal.pdf",sep=""),plot = p1, width = 10, height = 8)

write.csv(file = paste(savepath,"correlates_long.csv",sep=""),final_df)


# 3. Regional Correlates ------------
# We use tracts! And we use GM regions
outcomes = cross_env %>% select(FA_Mean, MD_Mean, RD_Mean, AD_Mean,
                                CortexVol,CorticalWhiteMatterVol) %>% names
T1w = read.csv("/cluster/projects/p33/users/maxk/UKB/data/T1w_50k/merged/T1_data.csv")
T1w = T1w %>% filter(euler_mean < mean(euler_mean)+3*sd(euler_mean)) %>% select(-Sex,-Age,-Scanner)
cross_env  = cross_env %>% select(!starts_with("Left")) %>% select(!starts_with("Right"))
regional_df = merge(cross_env,T1w, by="eid")
tracts = regional_df %>% select(ends_with("ATRL"),ends_with("ATRR"),
                                ends_with("CSTL"),ends_with("CSTR"),
                                ends_with("CINGL"),ends_with("CINGR"),
                                ends_with("CGL"),ends_with("CGR"),
                                ends_with("IFOFL"),ends_with("IFOFR"),
                                ends_with("ILFL"),ends_with("ILFR"),
                                ends_with("SLFL"),ends_with("SLFR"),
                                ends_with("UFL"),ends_with("UFR"),
                                ends_with("SLFTL"),ends_with("SLTFR"),
                                ends_with("FMIN"),ends_with("FMAJ")) %>%
  select(starts_with("AD_"), starts_with("MD_"), starts_with("RD_"), starts_with("FA_")) %>%
  names()
volumes = regional_df %>% select(starts_with("lh_"), starts_with("rh_"), 
                                 starts_with("Right"), starts_with("Left")) %>%
  select(!ends_with("area")) %>% select(!ends_with("thickness")) %>%
  names
# that gives us DTI scalars for tracts and volumes for DK atlas
outcomes = c(volumes,tracts)
# We define again predictors
predictors=cross_env %>% dplyr::select(ends_with("at_age"),"RDS1", "N1", 
                                       "ANX", "ADHD", "ASD", "BIP", "MDD", "OCD", "SCZ", "AD") %>% names
#
#
#
#
# Run models
cross_scaled <- regional_df %>% mutate(across(all_of(c(predictors, outcomes)), ~ if (is.numeric(.x)) as.numeric(scale(.x)) else .x))
cross_scaled$site = regional_df$site
results <- list()
lrt_results <- list()

for (outcome in outcomes) {
  for (predictor in predictors) {
    if (predictor == outcome) next
    
    vars_needed <- c(outcome, predictor,
                     "sex", "age", "SurfaceHoles", "EstimatedTotalIntraCranialVol.x", "income", "site")
    
    dat_sub <- regional_df %>%
      select(all_of(vars_needed)) %>%
      filter(complete.cases(.)) %>%
      mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))
    
    if (nrow(dat_sub) < 10) next
    if (length(unique(dat_sub$income)) < 2) next
    
    # ------------------------
    # Models
    # ------------------------
    full_f <- paste0(outcome, " ~ ", predictor,
                     " + sex + age + SurfaceHoles + EstimatedTotalIntraCranialVol.x + income + site")
    model_full <- lm(as.formula(full_f), data = dat_sub)
    
    reduced_income_f <- paste0(outcome, " ~ ", predictor,
                               " + sex + age + SurfaceHoles + EstimatedTotalIntraCranialVol.x + site")
    model_reduced_income <- lm(as.formula(reduced_income_f), data = dat_sub)
    
    td <- tidy(model_full, conf.int = TRUE)
    
    # ------------------------
    # Predictor
    # ------------------------
    pred_exact <- td %>% filter(term == predictor)
    if (nrow(pred_exact) == 1) {
      std_beta_predictor <- pred_exact$estimate
      conf_low_predictor <- pred_exact$conf.low
      conf_high_predictor <- pred_exact$conf.high
      p_value_predictor <- pred_exact$p.value
      predictor_omnibus_p <- NA_real_
    } else {
      reduced_no_pred_f <- paste0(outcome, " ~ sex + age + SurfaceHoles + EstimatedTotalIntraCranialVol.x + income + site")
      model_reduced_no_pred <- lm(as.formula(reduced_no_pred_f), data = dat_sub)
      an_pred <- anova(model_reduced_no_pred, model_full)
      predictor_omnibus_p <- an_pred$`Pr(>F)`[2]
      
      std_beta_predictor <- NA_real_
      conf_low_predictor <- NA_real_
      conf_high_predictor <- NA_real_
      p_value_predictor <- NA_real_
    }
    
    # ------------------------
    # Likelihood ratio test (income)
    # ------------------------
    an_income <- anova(model_reduced_income, model_full)
    
    r2_full <- summary(model_full)$r.squared
    r2_reduced <- summary(model_reduced_income)$r.squared
    delta_r2 <- r2_full - r2_reduced
    
    lrt_results[[paste0(outcome, "_by_", predictor)]] <- data.frame(
      outcome   = outcome,
      predictor = predictor,
      n_obs     = nrow(dat_sub),
      df        = an_income$Df[2],
      f_value   = an_income$F[2],
      p_value   = an_income$`Pr(>F)`[2],
      r2_full   = r2_full,
      r2_reduced = r2_reduced,
      delta_r2   = delta_r2,
      stringsAsFactors = FALSE
    )
    
    # ------------------------
    # Store predictor + income together
    # ------------------------
    # results[[paste0(outcome, "_by_", predictor)]] <- income_results %>%
    #   mutate(
    #     std_beta_predictor = std_beta_predictor,
    #     conf_low_predictor = conf_low_predictor,
    #     conf_high_predictor = conf_high_predictor,
    #     p_value_predictor = p_value_predictor,
    #     predictor_omnibus_p = predictor_omnibus_p
    #   )
    # Only store in results if predictor is NOT income
    if (!grepl("income", predictor, ignore.case = TRUE)) {
      results[[paste0(outcome, "_", predictor)]] <- data.frame(
        outcome = outcome,
        predictor = predictor,
        beta = std_beta_predictor,
        conf_low <- conf_low_predictor,
        conf_high <- conf_high_predictor,
        p_value <- p_value_predictor
      )
    }
  }
}
cross_res <- bind_rows(results)
#
#
#
#
# plotting many to many relationships is difficult.
# it is possible to select a single predictors such as the precipitation, which shows some of the strongest associations in global / whole brain associations
#
#
# first we save all results
write.csv(file = paste(savepath,"correlates_regional.csv",sep=""),cross_res)
#
#
# now, we filter the data for total precipitation within the last month
cross_res = cross_res %>% filter(predictor=="tp_at_age")

# correct the p-val
cross_res$p.correct = p.adjust(cross_res$p_value....p_value_predictor,method="BH")
#cross_res = cross_res %>% filter(p.correct<.05)
#
# give a star label
cross_res <- cross_res %>%
  mutate(sig_label = ifelse(p.correct < 0.05, "*", ""))

# log p vals
cross_res$p.correct = -log10(cross_res$p.correct)
# add a column of directionality for beta values / slopes
cross_res$Slope <- "No relation (p > 0.05)"
cross_res$Slope[cross_res$beta > 0 & cross_res$p.correct > -log10(0.05)] <- "Positive relation"
cross_res$Slope[cross_res$beta < 0 & cross_res$p.correct > -log10(0.05)] <- "Negativel relation"
#
# edit names as well
cross_res$defeature = ifelse(cross_res$p.correct>10,cross_res$outcome,NA)
cross_res$defeature = gsub("_"," ", cross_res$defeature)
cross_res$defeature = gsub("lh","left", cross_res$defeature)
cross_res$defeature = gsub("rh","right", cross_res$defeature)
# plot
volcano1 = ggplot(data=cross_res, aes(x=beta, y=p.correct,col = Slope, label=defeature)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(data = subset(cross_res, beta < -0.01),colour='black', nudge_x = -0.05, direction = "y", segment.size = 0.1, xlim = c(-0.15,-0.35))+ #c(-0.3,-0.6))+
  geom_text_repel(data = subset(cross_res, beta > 0.01),colour='black', nudge_x = 0.05, direction = "y",segment.size = 0.1, xlim = c(0.15,0.35))+
  scale_color_manual(values=c("#0072B2", "#999999","#D55E00", "#56B4E9", "#E69F00"))+
  #scale_color_manual(values=c("#0072B2","#D55E00", "#999999", "#56B4E9", "#E69F00"))+
  xlab("Corrected standardized effect of total precipitation")+ylab("-log10(FDR-corrected p)")+
  xlim(-0.25,0.25)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
#ggsave(filename = "/tsd/p33/home/p33-maxk/p2.pdf",plot = volcano1, width = 9, height = 7)
ggsave(filename = paste(savepath,"Volcano_Total_Precipitation.pdf",sep=""),plot = volcano1, width = 9, height = 7)
#
#
#
#
#
#
#
#
#
#
#
# Do the same for 10-meter wind speed = si10_at_age
cross_res <- bind_rows(results)
#
# now, we filter the data for total precipitation within the last month
cross_res = cross_res %>% filter(predictor=="si10_at_age")

# correct the p-val
cross_res$p.correct = p.adjust(cross_res$p_value....p_value_predictor,method="BH")
#cross_res = cross_res %>% filter(p.correct<.05)
#
# give a star label
cross_res <- cross_res %>%
  mutate(sig_label = ifelse(p.correct < 0.05, "*", ""))

# log p vals
cross_res$p.correct = -log10(cross_res$p.correct)
# add a column of directionality for beta values / slopes
cross_res$Slope <- "No relation (p > 0.05)"
cross_res$Slope[cross_res$beta > 0 & cross_res$p.correct > -log10(0.05)] <- "Positive relation"
cross_res$Slope[cross_res$beta < 0 & cross_res$p.correct > -log10(0.05)] <- "Negativel relation"
#
# edit names as well
cross_res$defeature = ifelse(cross_res$p.correct>20,cross_res$outcome,NA)
cross_res$defeature = gsub("_"," ", cross_res$defeature)
cross_res$defeature = gsub("lh","left", cross_res$defeature)
cross_res$defeature = gsub("rh","right", cross_res$defeature)
# plot
volcano2 = ggplot(data=cross_res, aes(x=beta, y=p.correct,col = Slope, label=defeature)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(data = subset(cross_res, beta < -0.01),colour='black', nudge_x = -0.05, direction = "y", segment.size = 0.1, xlim = c(-0.15,-0.35))+ #c(-0.3,-0.6))+
  geom_text_repel(data = subset(cross_res, beta > 0.01),colour='black', nudge_x = 0.05, direction = "y",segment.size = 0.1, xlim = c(0.15,0.35))+
  scale_color_manual(values=c("#0072B2", "#999999","#D55E00", "#56B4E9", "#E69F00"))+
  #scale_color_manual(values=c("#0072B2","#D55E00", "#999999", "#56B4E9", "#E69F00"))+
  xlab("Corrected standardized effect of 10m Wind Speed")+ylab("-log10(FDR-corrected p)")+
  xlim(-0.25,0.25)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
#ggsave(filename = "/tsd/p33/home/p33-maxk/p3.pdf",plot = volcano2, width = 10, height = 7)
ggsave(filename = paste(savepath,"Volcano_Wind_Speed.pdf",sep=""),plot = volcano2, width = 10, height = 7)


