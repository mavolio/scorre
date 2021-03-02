################################################################################
##  modified from: fit_cumulative_gams.R; Author: Andrew Tredennick (atredenn@gmail.com); Date created: March 21, 2018
################################################################################

# NOTES:
#  (1) for deviance and AIC, negative deltas indicate better models
#  (2) the p-value then says whether the deviance difference is significant
#  (3) some p-values will be NA -- this is OK and indicates the full the model
#      is DEFINITELY NOT BETTER than the null model. So, think of NA as p>0.05.

####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
#library(ggthemes)
library(mgcv)


#laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data')


#functions
###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

#theme set
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###read in data
trt <- read.csv('CoRRE_raw_abundance_Nov2019.csv')%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id)%>%
  unique()%>%
  left_join(read.csv('CoRRE_project_summary.csv'))%>%
  left_join(read.csv('CoRRE_treatment_summary.csv'))

pDiv <- read.csv('CoRRE_pd_metrics.csv')%>%
  separate(plot_id2, c('site_code', 'project_name', 'community_type', 'treatment_year', 'plot_id'), sep='::')%>%
  mutate(treatment_year=as.integer(treatment_year))%>%
  left_join(trt)%>%
  filter(treatment_year>0)%>%
  filter(!is.na(trt_type)) #NGBER has a trt type of unknown, fix for final analysis

#subset down to just the treatment types that have been replicated across many experiments
dat <- pDiv%>%
  filter(trt_type %in% c('CO2', 'drought', 'irr', 'mult_nut', 'N', 'P', 'N*P', 'temp', 'precip_vari'))%>%
  select(site_code, project_name, community_type)%>%
  unique()%>%
  mutate(keep=1)%>%
  left_join(pDiv)%>%
  filter(trt_type %in% c('control', 'CO2', 'drought', 'irr', 'mult_nutrient', 'N', 'P', 'N*P', 'temp', 'precip_vari'))



####
####  SET WORKING DIRECTORIES AND FILENAMES ------------------------------------
####
work_dir  <- "C:\\Users\\lapie\\Desktop\\R files laptop\\scorre\\Phylogenetic Analyses\\" # change as needed
data_dir  <- "C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\Div\\sDiv_sCoRRE_shared\\CoRRE data\\"
results_dir <- "C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\phylogeny and continuous traits group\\2021results\\"
data_file <- "CoRRE_pd_metrics.csv"
# setwd(work_dir)

####
####  DEFINE MODEL FITTING FUNCTION --------------------------------------------
####

fit_compare_gams <- function(df, response, diff_type = "last_year"){
  # Fits two GAMs and compares them with AIC and LLR
  #
  # Args:
  #  data: a dataframe with necessary columns for fitting the GAMMs
  #  response: name of the response variable, must be a column in the dataframe
  #  diff_type: type of treatment-control difference to return. options are one
  #   of "all_years" (mean of all diffs), "year_five" (diff at year 5), and
  #   "last_year" (diff at final year of data). Note that the option "year_five"
  #   requires the dataset to have at least five years of data. "last_year" is
  #   default.
  #
  # Returns:
  #  A tibble with LLR delta deviance, LLR p-value, and delta AIC
  
  # Check that there aren't too many NAs and skip modeling if fraction > 0.5
  y <- df[ , response]
  num_nas <- length(which(is.na(y)))
  fraction_nas <- num_nas/nrow(y)
  
  if(fraction_nas >= 0.5){
    return(
      tibble(
        response_var = response,
        p_value = -9999,
        delta_deviance = NA,
        delta_aic = NA,
        diff = NA,
        diff_se = NA,
        diff_lower = NA,
        diff_upper = NA,
        diff_treatment_year = NA
      )
    )
  }
  
  if(fraction_nas <= 0.5){
    test_formula <- as.formula(
      paste(response, 
            "~ treatment + s(treatment_year, by = treatment, k = (num_years-1)) + 
            s(plot_id, bs='re')"
      )
      )
    
    null_formula <- as.formula(
      paste(response, 
            "~ s(treatment_year, k = (num_years-1)) + 
            s(plot_id, bs='re')"
      )
      )
    
    gam_test <- gam(
      test_formula, 
      data = df,
      method = "REML"
    )
    
    gam_null <- gam(
      null_formula,
      data = df, 
      method = "REML"
    )
    
    # LLR tests
    pvalue <- anova(gam_null, gam_test, test="Chisq")$`Pr(>Chi)`[2]
    dev <- anova(gam_null, gam_test, test="Chisq")$`Resid. Dev`
    delta_div <- diff(dev)  # full - null
    
    # AIC tests
    aics <- AIC(gam_null, gam_test)$AIC
    delta_aic <- diff(aics)  # full - null
    
    # Simulate predictions from the interaction model
    min_year <- min(df$treatment_year)
    max_year <- max(df$treatment_year)
    pdat <- expand.grid(
      treatment_year = seq(min_year, max_year, by = 1),
      treatment = unique(df$treatment)
    )
    pdat$plot_id <- unique(df$plot_id)[1]
    
    control_name <- filter(df, plot_mani == 0) %>% pull(treatment) %>% unique()
    treat_name <- filter(df, plot_mani != 0) %>% pull(treatment) %>% unique()
    tmp_diffs <- smooth_diff(model = gam_test, newdata = pdat, 
                             f1 = treat_name, 
                             f2 = control_name, 
                             var = "treatment", alpha = 0.05, 
                             unconditional = FALSE)
    
    if(diff_type == "all_years"){
      outdiff <- as.data.frame(t(colMeans(tmp_diffs[c("diff","se","upper","lower")])))
      outdiff$treatment_year <- NA
    }
    
    if(diff_type == "mid_year"){
      # find median index, rounds up
      # mid_year <- floor(0.5 + median(1:nrow(tmp_diffs)))
      mid_year <- 5
      outdiff <- tmp_diffs[mid_year,] 
    }
    
    if(diff_type == "last_year"){
      outdiff <- tail(tmp_diffs, 1)
    }
    
    return(
      tibble(
        response_var = response,
        p_value = pvalue,
        delta_deviance = delta_div,
        delta_aic = delta_aic,
        diff = outdiff$diff,
        diff_se = outdiff$se,
        diff_lower = outdiff$lower,
        diff_upper = outdiff$upper,
        diff_treatment_year = outdiff$treatment_year
      )
    )
  }
  
}  # end of model fit and comparison function


####
####  DEFINE FUNCTION FOR CALCULATING SMOOTHER DIFFS ---------------------------
####

smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  # Calculates the difference between two fitted smooths from
  # a GAM model with a 'by' interaction term. This function written
  # by Gavin Simpson:
  # (https://www.fromthebottomoftheheap.net/2017/10/10/difference-splines-i/)
  # Modified slightly by Andrew Tredennick.
  #
  # Args:
  #   model: A fitted GAM model object. Must include an interaction smooth
  #    term via 's(by = xterm)'
  #  newdata: A dataframe with specific predictors for the model.
  #  f1: First factor variable for comparison (control id)
  #  f2: Second factor variable for comparison (treatment id)
  #  var: Name of the 'by' variable (e.g., column name that contains f1 and f2)
  #  alpha: Significance level. Default is 0.05.
  #  unconditional: Correct for uncertainty in paramter covariance (TRUE) or 
  #   or not (FALSE). Default is FALSE.
  #
  # Returns:
  #  A data frame with the difference, se, and 95% CIs for each of the time
  #  points in newdata.
  
  xp <- predict(model, newdata = newdata, type = 'lpmatrix', 
                exclude = "s(plot_id)")  # make prediction for avareage plot
  
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  
  # difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  
  # zero out cols of X related to splines for other factors (redundant here)
  X[, ! (c1 | c2)] <- 0
  
  # zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  
  data.frame(treatment_year = unique(newdata$treatment_year),
             pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}  # end of smooth_diff function



####
####  DEFINE FUNCTION TO FILL TIBBLE WHEN YEARS < 4 ----------------------------
####

fill_empties <- function(...){
  return(
    tibble(
      response_var = c(
        "pd_abs", 
        "mpd_abs", 
        "mntd_abs"
      ),
      p_value = NA,
      delta_deviance = NA,
      delta_aic = NA,
      diff = NA,
      diff_se = NA,
      diff_lower = NA,
      diff_upper = NA,
      diff_treatment_year = NA
    )
  )
}



####
####  READ IN DATA AND CALCULATE CUMULATIVE CHANGE -----------------------------
####
change_metrics <- as_tibble(dat) %>%
  mutate(treatment=as.factor(treatment)) %>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))

##  Calculate cumulative sums of each metric (from Kevin)
change_cumsum <- change_metrics %>%
  filter(treatment_year!=0)%>%
  filter(experiment_length>3)%>%
  mutate(treatment_year=as.numeric(treatment_year)) %>%
  arrange(site_project_comm, plot_id, treatment_year) %>%
  group_by(site_project_comm, treatment, plot_id) %>%
  mutate(pd_abs = abs(pd.ses)) %>%
  mutate(mpd_abs = abs(mpd.ses)) %>%
  mutate(mntd_abs = abs(mntd.ses)) %>%
  mutate_at(vars(pd_abs, 
                 mpd_abs, 
                 mntd_abs), 
            list(cumsum)) %>%
  ungroup()%>%
  mutate(control = ifelse(plot_mani==0,"control","treatment")) %>%
  select(site_project_comm, plot_id, treatment_year, plot_mani, control, treatment, pd_abs, mpd_abs, mntd_abs)



####
####  LOOP OVER SITE_PROJECT_COMMS AND COMPARE CONTROLS VS. TREATMENTS ---------
####
all_sites <- unique(change_cumsum$site_project_comm)
all_comparisons <- {} # empty object for storage
diff_type <- "last_year"

for(do_site in all_sites){
  site_data <- filter(change_cumsum, site_project_comm == do_site)
  site_controls <- filter(site_data, plot_mani == 0)
  site_treatments <- filter(site_data, plot_mani != 0)
  all_treatments <- unique(site_treatments$treatment)
  
  for(do_treatment in all_treatments){
    treatment_data <- filter(site_treatments, treatment == do_treatment)
    model_data <- rbind(site_controls, treatment_data)
    num_years <- length(unique(model_data$treatment_year))
    
    # Skip data with less than three pairs of years
    if(num_years < 4){
      tmp_out <- fill_empties() %>%
        mutate(
          site_proj_comm = do_site,
          treatment = do_treatment
        ) %>%
        dplyr::select(
          site_proj_comm,
          treatment,
          response_var,
          p_value,
          delta_deviance,
          delta_aic
        )
    }
    
    # Compare models for data with more than four pairs of years
    if(num_years > 3){
      
      # PD
      pd_test <- fit_compare_gams(
        df = model_data,
        response = "pd_abs",
        diff_type = "last_year"
      )
      
      # MPD
      mpd_test <- fit_compare_gams(
        df = model_data,
        response = "mpd_abs",
        diff_type
      )
      
      # MNTD
      mntd_test <- fit_compare_gams(
        df = model_data,
        response = "mntd_abs",
        diff_type
      )
      
      tmp_out <- bind_rows(
        pd_test,
        mpd_test,
        mntd_test
      ) %>%
        mutate(
          site_proj_comm = do_site,
          treatment = do_treatment
        ) %>%
        dplyr::select(
          site_proj_comm,
          treatment,
          response_var,
          p_value,
          delta_deviance,
          delta_aic,
          diff,
          diff_se,
          diff_lower,
          diff_upper,
          diff_treatment_year
        )
      
    } # end if for num_years
    
    all_comparisons <- all_comparisons %>%
      bind_rows(tmp_out)
    
  } # end treatment loop
  
  print(paste("Done with site:", do_site))
  
} # end site loop



####
####  SAVE DELTA_AIC TABLE -----------------------------------------------------
####
save_comparisons <- all_comparisons %>%
  filter(is.na(delta_deviance) == FALSE) %>%  # remove sites we don't model
  mutate(
    sig_diff_cntrl_trt = ifelse(
      p_value <= 0.05 & sign(delta_deviance) == -1,
      "yes",
      "no"
    )
  ) %>%
  mutate(
    sig_diff_cntrl_trt = ifelse(is.na(sig_diff_cntrl_trt) == TRUE, "no", sig_diff_cntrl_trt)
  )

outfile <- paste0("gam_comparison_table_", diff_type, ".csv")
write_csv(
  x = save_comparisons,
  path = paste0(results_dir, outfile)
)


# ggplot(save_comparisons, aes(x = final_diff, fill = sig_diff_cntrl_trt)) +
#   geom_histogram() +
#   geom_vline(aes(xintercept = 0)) +
#   facet_wrap(~response_var, scales = "free")


####
####  TALLY THE RESULTS --------------------------------------------------------
####
# gam_results <- read_csv(paste0(results_dir, "gam_comparison_table.csv"))
# 
# sig_tally <- gam_results %>%
#   group_by(response_var) %>%
#   summarise(
#     num_sig = length(which(sig_diff_cntrl_trt == "yes")),
#     num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
#   ) %>%
#   gather(key = sig, value = value, -response_var) %>%
#   mutate(
#     response_var = ifelse(response_var == "evenness_change_abs", "Evenness", response_var),
#     response_var = ifelse(response_var == "gains", "Species gains", response_var),
#     response_var = ifelse(response_var == "losses", "Species losses", response_var),
#     response_var = ifelse(response_var == "rank_change", "Rank change", response_var),
#     response_var = ifelse(response_var == "richness_change_abs", "Richness", response_var),
#     response_var = ifelse(response_var == "composition_change", "Composition", response_var)
#   )
# 
# # Of them all, how many had treatment-time interactions for composition change?
# percent_with_compositional_change <- sig_tally %>%
#   filter(response_var == "Composition") %>%
#   spread(sig, value) %>%
#   mutate(
#     percent_sig = num_sig/(num_nonsig + num_sig)
#   )
# 
# # Of the significant treatment-time interactions, tally other aspects
# sites_with_trt_time_itxn <- gam_results %>%
#   filter(response_var == "composition_change") %>%
#   filter(sig_diff_cntrl_trt == "yes") %>%
#   dplyr::select(site_proj_comm, treatment)
# 
# aspects_of_sigones <- gam_results %>%
#   right_join(sites_with_trt_time_itxn, by = c("site_proj_comm", "treatment")) %>%
#   filter(response_var != "composition_change")
# 
# sig_tally_changers <- aspects_of_sigones %>%
#   group_by(response_var) %>%
#   summarise(
#     num_sig = length(which(sig_diff_cntrl_trt == "yes")),
#     num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
#   ) %>%
#   gather(key = sig, value = value, -response_var) %>%
#   mutate(
#     response_var = ifelse(response_var == "evenness_change_abs", "Evenness", response_var),
#     response_var = ifelse(response_var == "gains", "Species gains", response_var),
#     response_var = ifelse(response_var == "losses", "Species losses", response_var),
#     response_var = ifelse(response_var == "rank_change", "Rank change", response_var),
#     response_var = ifelse(response_var == "richness_change_abs", "Richness", response_var)
#   ) %>%
#   spread(sig, value) %>%
#   mutate(
#     proportion_nonsig = num_nonsig / (num_nonsig + num_sig),
#     proportion_sig = num_sig / (num_nonsig + num_sig)
#   ) %>%
#   dplyr::select(-num_nonsig, -num_sig) %>%
#   gather(key = sig, value = value, -response_var)
# 
# ggplot(sig_tally_changers, aes(x = response_var, y = value, fill = sig)) +
#   geom_col(width = 0.7) +
#   geom_hline(aes(yintercept = 0.5)) +
#   coord_flip() +
#   theme_minimal() +
#   scale_fill_manual(values = c("grey25", "grey75"), name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
#   # scale_fill_brewer(type = "bw", name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
#   labs(x = "Aspect of community change", y = "Proportion of communities") +
#   theme(legend.position = "top")
# 
# ggsave(
#   filename = paste0(results_dir, "cumulative_changers_summary.png"),
#   width = 6,
#   height = 4,
#   units = "in"
# )



####
####  TEST FIT USING JUST PPLOTS -----------------------------------------------
####
# #  Subset controls and one treatment for pplots
# do_site <- "CDR_e001_D"
# do_treatment <- 6
# pplots_cumsum  <- filter(change_cumsum, site_project_comm == do_site)
# pplot_controls <- filter(pplots_cumsum, treatment == 9 | treatment == do_treatment)
# 
# ##  Fit nested GAMs
# ### TODO: CHECK THESE MODEL STRUCTURES!!! Email Gavin?
# gam_test <- gam(evenness_change_abs ~ s(treatment_year, treatment, bs = "fs") +
#                   s(plot_id, bs="re"), data = pplot_controls)
# gam_null <- gam(evenness_change_abs ~ s(treatment_year) + s(plot_id, bs="re"),
#                 data = pplot_controls)
# 
# ##  Compare nested models
# AIC(gam_test, gam_null) # less conservative
# BIC(gam_test, gam_null) # more conservative
# deltaAIC <- diff(AIC(gam_test, gam_null)$AIC)
# 
# ##  Make example plot of predictions from the complex model
# newdata <- data.frame(treatment_year = rep(0:11,2),
#                       treatment = rep(c("N1P0","N2P0"), each = 12),
#                       plot_id = rep(factor(27),times = 12*2))
# gam_preds <- predict(gam_test,
#                      newdata = newdata,
#                      type = "response",
#                      exclude = "s(plot_id)")
# plot_preds <- data.frame(predictions = gam_preds,
#                          treatment_year = rep(0:11,2),
#                          treatment = rep(c("N1P0","N2P0"), each = 12))
# plot_data <- pplot_controls %>%
#   ungroup() %>%
#   select(treatment_year, rank_change, treatment)
# 
# ggplot()+
#   geom_point(data = plot_data,
#              aes(x = treatment_year, y = rank_change, color = treatment))+
#   geom_line(data = plot_preds,
#             aes(x = treatment_year, y = predictions, color = treatment),
#             size = 1.2)+
#   scale_color_brewer(type = "qual", name = "Treatment")+
#   xlab("Treatment year")+
#   ylab("Rank change")+
#   ggtitle("KNZ PPlots Example", subtitle = "Control (N1P0) vs. N2P0")
# ggsave(filename = paste0(data_dir,"figures/pplot_gam_example.pdf"),
#        height = 4,
#        width = 5,
#        units = "in")


##making pplots examples for the talk
# theme_set(theme_bw(12))
# pplots<-change_cumsum%>%
#   filter(site_project_comm=="KNZ_pplots_0")%>%
#   filter(treatment == "N2P0"|treatment=="N2P3"|treatment=="N1P0")
# 
# ggplot(data = subset(pplots, treatment!="N2P0"),
#        aes(x = treatment_year, y = richness_change_abs, color = treatment))+
#   geom_point()+
#   geom_smooth(method = "lm", aes(color = treatment), se = F)+
#   scale_color_brewer(type = "qual", name = "Treatment", labels= c("Control", "N+P"))+
#   xlab("Treatment year")+
#   ylab("Richness change")
# 
# ggplot(data = subset(pplots, treatment!="N2P3"),
#        aes(x = treatment_year, y = richness_change_abs, color = treatment))+
#   geom_point()+
#   geom_smooth(method = "lm", aes(color = treatment), se = F)+
#   scale_color_brewer(type = "qual", name = "Treatment", labels= c("Control", "N"))+
#   xlab("Treatment year")+
#   ylab("Richness change")


