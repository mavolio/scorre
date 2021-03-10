library(lme4)


model_invasion<- glmer(
  presence ~ 
    NPP+ ddeg + seasonality + footprint +  mdwnc +
    trait_SLAs:NPP + trait_SLAs:ddeg + trait_SLAs:seasonality + trait_SLAs:footprint +trait_SLAs:mdwnc +
    trait_PHs:NPP + trait_PHs:ddeg +  trait_PHs:seasonality + trait_PHs:footprint + trait_PHs:mdwnc + 
    trait_SMs:NPP + trait_SMs:ddeg +trait_SMs:seasonality +  trait_SMs:footprint +trait_SMs:mdwnc +
    (1 + NPP +  seasonality + footprint + ddeg + mdwnc | invasive.species),
  family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10241)), data=invasion
)


library(lme4)
a1 <- lmer(ave_diff ~ trt_type2 + (1|species_matched), data = species.data)
summary(a1)
