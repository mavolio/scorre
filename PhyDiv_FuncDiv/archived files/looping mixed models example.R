### Looping mixed models


### TB 2019 by time period
time_period_vec_tb19 <- levels(factor(subset(smoist_data_all, Site=="TB" & Year==2019)$Time_Period))
smoist_model_tb_2019_byperiod_master <- {}

for(TIME in 1:length(time_period_vec_tb19)){
  data_temp <- smoist_data_all %>%
    filter(Site=="TB" & Year==2019 & Time_Period==time_period_vec_tb19[TIME])
  model_temp <- lme(Soil_Moisture ~ Drought
                    , data=data_temp
                    , random = ~1 |Block/Paddock/Plot
                    , na.action = na.omit)
  model_out_temp <- data.frame(Site="TB", 
                               Year="2019", 
                               Time_Period=time_period_vec_tb19[TIME],
                               Model_term=row.names(anova.lme(model_temp, type="marginal")[2,]),
                               anova.lme(model_temp, type="marginal")[2,])
  smoist_model_tb_2019_byperiod_master <- rbind(smoist_model_tb_2019_byperiod_master, model_out_temp)
  rm(data_temp, model_temp, model_out_temp)
}
