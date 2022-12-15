######
#####
##### code to study how species traits affect response to GCDs
##### code by M.Avolio with help from Adam Clark and Tamara Munkemuller
#### created summer 2021, updated with new trait dataset June 20, 2022
#### updated 14.12.2022 by A. Clark
######
rm(list=ls())
setwd("~/Dropbox/tmp/desktop_dropbox/tmp/sCORRE_analysis/")

library(tidyverse)
library(lme4)
library(emmeans)
require(RColorBrewer)
#library(relaimpo)

theme_set(theme_bw(12))

scale_fun = function(x) {
  x.new = log(x)
  x.new = x-mean(x,na.rm=T)
  x.new = x/sd(x, na.rm=T)
  x.new
}

### Trait data

# #barGraphStats(data=, variable="", byFactorNames=c(""))
# barGraphStats <- function(data, variable, byFactorNames) {
#   count <- length(byFactorNames)
#   N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
#   names(N)[1:count] <- byFactorNames
#   names(N) <- sub("^x$", "N", names(N))
#   mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
#   names(mean)[1:count] <- byFactorNames
#   names(mean) <- sub("^x$", "mean", names(mean))
#   sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
#   names(sd)[1:count] <- byFactorNames
#   names(sd) <- sub("^x$", "sd", names(sd))
#   preSummaryStats <- merge(N, mean, by=byFactorNames)
#   finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
#   finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
#   return(finalSummaryStats)
# }  

#read in data
#contTraits1 <- read.csv('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\Final TRY Traits\\Imputed Continuous_Traits\\data to play with\\imputed_continuous_20220620.csv')
contTraits1 <- read.csv('~/Dropbox/SharedFolders/sDiv_sCoRRE_shared/CoRRE data/trait data/Final TRY Traits/Imputed Continuous_Traits/data to play with/imputed_continuous_20220620.csv')
contTraits<-contTraits1%>%
  select(-X.1, -X, -family, -genus, -observation)%>%
  select(species_matched, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass, leaf_C.N) %>% 
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()  %>% 
  filter(seed_dry_mass<30, plant_height_vegetative<10, rooting_depth<3, SLA<75)

# # ##code to drop outliers
# contLong<-contTraits %>%
#   pivot_longer(LDMC:seed_dry_mass, names_to = "trait", values_to = "value")
# ggplot(data=contLong, aes(x=value))+
#   geom_histogram()+
#   facet_wrap(~trait, scales = "free")

#why do I need to scale? Kim, Kevin and I discussed this on 8/10/22 and decided I did not have to do this.
# traitsScaled <- contTraits%>%
#   mutate_at(vars(LDMC:seed_dry_mass), scale)



#catTraits <- read.csv('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\Final TRY Traits\\sCoRRE categorical trait data_final_20211209.csv') %>% 
catTraits <- read.csv('~/Dropbox/SharedFolders/sDiv_sCoRRE_shared/CoRRE data/trait data/Final TRY Traits/sCoRRE categorical trait data_final_20211209.csv')
catTraits <- catTraits %>% 
  select(species_matched, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation) %>% 
  filter(growth_form!="moss",species_matched!="", growth_form!="lycophyte") %>% 
  mutate(mycorrhizal=ifelse(mycorrhizal_type=="none", 'no', ifelse(mycorrhizal_type=="uncertain", "unk", "yes"))) %>% 
  select(-mycorrhizal_type) %>% 
  mutate(photo_path=ifelse(photosynthetic_pathway=="possible C4"|photosynthetic_pathway=="possible C4/CAM", "C4", ifelse(photosynthetic_pathway=="possible CAM", "CAM",photosynthetic_pathway))) %>% 
  select(-photosynthetic_pathway)

## simplify groupings
catTraits$photo_path_simple = catTraits$photo_path
#catTraits$photo_path_simple[catTraits$photo_path_simple%in%c("C4", "CAM")] = "C4.CAM"
catTraits$photo_path_simple[catTraits$photo_path_simple%in%c("parasitic", "possible hybrid")] = NA
catTraits$photo_path = catTraits$photo_path_simple

catTraits$lifespan_simple = catTraits$lifespan
catTraits$lifespan_simple[catTraits$lifespan_simple%in%c("annual", "biennial")] = "ann.bien"
catTraits$lifespan_simple[catTraits$lifespan_simple%in%c("uncertain")] = NA
catTraits$lifespan = catTraits$lifespan_simple

pairs(contTraits[,2:6])


# Read in dci diff
# dcidiff<-read.csv("C:/Users/megha/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_newtrts.csv")
#dcidiff_models<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodelsNov22.csv")
dcidiff_models<-read.csv("~/Dropbox/SharedFolders/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodelsNov22.csv")

test<-dcidiff_models %>% 
  filter(trt_type2=="all mult") %>% 
  select(species_matched) %>% 
  unique()

alldat_cont<-dcidiff_models%>%
  right_join(contTraits)%>%
  gather(LDMC:leaf_C.N, key="trait", value="value")%>%
  na.omit()

#Making a graph of all the models I am goign to run now, just to see if there are patterns.
#ggplot(data=alldat_cont, aes(x=value, y=diff))+
#  geom_point()+
#  geom_smooth(method="lm")+
#  facet_grid(trt_type2~trait, scales="free")
###
#making mixed models - run through all traits.

##-1 says don't give me an overall intercept. just changing the output. I want an intercept per level of each trt_type2
##fit a seperate slope for the value term for each of trt_type levels
## by using trt|species - come up with a value for each species, but do that for each trt type category
#(trt_type|species matches) = get a separate estimate for every species in each trt type
#(1|species matches) = each species can have a different ave regardless of trt type
#fixef(m1) #should give fixed effects of model


tmp = alldat_cont
tmp = tmp[!is.na(tmp$diff),]
cont_n_data = table(unique(tmp[,c("species_matched", "trt_type2", "trait")]))
cont_n_data = apply(cont_n_data, 2:3, sum)
write.csv(cont_n_data, "cont_n_data.csv")

##SLA
mSLA<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="SLA"))
summary(mSLA)

plot.sla<-as.data.frame(summary(mSLA)$coefficients)

plot.sla<-plot.sla%>%
  mutate(fixedef=row.names(plot.sla))

toplot.SLA<-plot.sla%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="SLA")

#not using SRL b/c data is no good on TRY.
# ###SRL
# 
# mSRL<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="SRL"))
# summary(mSRL)
# 
# plot.srl<-as.data.frame(summary(mSRL)$coefficients)%>%
#   mutate(fixedef=row.names(plot.srl))
# 
# toplot.SRL<-plot.srl%>%
#   separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
#   filter(!is.na(interaction))%>%
#   separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
#   select(-interaction, -drop)%>%
#   mutate(trait="SRL")

###seed mass
#this model failed to converge what does that mean?
mSM<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="seed_dry_mass"))
summary(mSM)

plot.sm<-as.data.frame(summary(mSM)$coefficients)

plot.sm<-plot.sm %>% mutate(fixedef=row.names(plot.sm))

toplot.SM<-plot.sm%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="Seed Mass")


# ###seed number - not doing this the data was weird.
# 
# mSN<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="seed_number"))
# summary(mSN)
# 
# plot.sn<-as.data.frame(summary(mSN)$coefficients)%>%
#   mutate(fixedef=row.names(plot.sn))
# 
# toplot.SN<-plot.sn%>%
#   separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
#   filter(!is.na(interaction))%>%
#   separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
#   select(-interaction, -drop)%>%
#   mutate(trait="Seed Number")

###rooting depth

mRD<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="rooting_depth"))
summary(mRD)

plot.rd<-as.data.frame(summary(mRD)$coefficients)

plot.rd<-plot.rd%>%
  mutate(fixedef=row.names(plot.rd))

toplot.RD<-plot.rd%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="Rooting Depth")

###LDMC

mLDMC<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="LDMC"))
summary(mLDMC)

plot.ldmc<-as.data.frame(summary(mLDMC)$coefficients)

plot.ldmc<-plot.ldmc%>%
  mutate(fixedef=row.names(plot.ldmc))

toplot.LDMC<-plot.ldmc%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="LDMC")

###plot height vegetative

mhght<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="plant_height_vegetative"))
summary(mhght)

plot.hght<-as.data.frame(summary(mhght)$coefficients)

plot.hght<-plot.hght%>%
  mutate(fixedef=row.names(plot.hght))

toplot.hgt<-plot.hght%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="PlantHeight")


##CtN
mCtN<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="leaf_C.N"))
summary(mCtN)

plot.ctn<-as.data.frame(summary(mCtN)$coefficients)

plot.ctn<-plot.ctn%>%
  mutate(fixedef=row.names(plot.ctn))

toplot.ctn<-plot.ctn%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="leaf_C.N")


toplot<-toplot.SLA%>%
  bind_rows(toplot.SM, toplot.RD, toplot.LDMC, toplot.hgt, toplot.ctn) %>% 
  rename(SE="Std. Error") %>% 
  mutate(trt_type2=factor(trt_type, levels=c("co2", "drought", "irrigation", "temp", "n", "p", "multnuts", "all mult"))) %>% 
  mutate(min=Estimate-SE, max=Estimate+SE) %>% 
  mutate(sig=ifelse(min>0&max>0, "*", ifelse(min<0&max<0, "*", "")))

#ggplot(data=toplot, aes(y=Estimate, x=1, label=sig))+
#  geom_point()+
#  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0.05)+
#  coord_flip()+
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
#  xlab("")+
#  scale_x_continuous(limits=c(0, 2))+
#  geom_hline(yintercept=0, linetype="dashed")+
#  geom_text(nudge_x = 0.2, nudge_y = 0.01, size=10, color="red")+
#  facet_grid(trt_type2~trait, labeller = labeller(trt_type2=trt.labels))

##############
# new plot
##############
make_boxplot = function(toplot_data = toplot, # data to plot
                                   trt.labels_data = trt.labels, # treatment labels
                                   trait.labels_data=trait.labels, # trait labels
                                   p_alpha = 0.025, # alpha value for significance tests
                                   groupbytrait = FALSE, # should individual panels include multiple traits?
                                   legend_line_length = 1.55, # length of lines in legend
                                   legend_cex = 0.9, # text size in legend
                                   lower_margin = 5.5, # margin for plotting category names
                                   sigadj = 0, # shift the asterisk up or down,
                                   legend_textwidth = 0.2, # text spacing for legend
                                   legend_yadj = 0.35, # margin spacing for legend
                                   legend_xadj = 0, # horizontal adjustment for the legend 
                                   flipaxes = TRUE, # plot with quantitative values on the x-axis
                                   traitorder = NULL, # optional order in which to plot traits
                                   group_colors = NULL,  # optional vector of colors for backgrounds of individual panels
                                   xtit = "Effect Size, DCi diff. vs. Standardized ln(Trait)", # axis title
                                   xlm = NULL,  # optional x-axis limits
                                   collst = NULL, # optional vector of colors for the bars
                                   autopar = TRUE, # should par settings be handled automatically?
                                   axisside = 1+flipaxes, # which side should the axis labels be plotted on?
                                   n_table = NULL, ncex = 0.6, yadj_text = 0) { # optional table of n values to be plotted; size for plotted text; vertical adjustment for text
  lwd_lines = 2.5
  lwd_diff_fact = 0.5
  
  # get significance based on estimate, se, and desired alpha
  tmp = pnorm(0, abs(toplot_data$Estimate), toplot_data$SE) <= p_alpha
  sig = tmp; sig[tmp] = "*"; sig[!tmp] = ""
  
  trt_type2_levels = as.character(sort(unique(toplot_data$trt_type2[!is.na(toplot_data$trt_type2)])))
  trt_type2_levels = trt_type2_levels[match(names(trt.labels_data), trt_type2_levels)]
  
  trait_levels = sort(unique(toplot_data$trait))
  trait_levels = trait_levels[match(names(trait.labels_data), trait_levels)]
  if(!is.null(traitorder) & !groupbytrait) {
    trait_levels = trait_levels[traitorder]
    trait.labels_data = trait.labels_data[traitorder]
  }
  
  if(groupbytrait) {
    v1 = trt_type2_levels
    v2 = trait_levels
  } else {
    v1 = trait_levels
    v2 = trt_type2_levels
  }
  if(is.null(collst)) {
    collst = RColorBrewer::brewer.pal(n = length(v2), name = "Dark2")
  }
  
  xps = seq(1, length(v1))
  dx = seq(0.35, -0.35, length=length(v2))
  whisker_length = min(diff(dx))/5
  
  yrng = c(min(toplot_data$Estimate-2*toplot_data$SE),
           max(toplot_data$Estimate+2*toplot_data$SE))
  
  if(autopar) {
    if(flipaxes) {
      par(mar=c(4,lower_margin,3.5,1))
    } else {
      par(mar=c(lower_margin,4,3,1))
    }
  }
  
  if(flipaxes) {
    if(!is.null(xlm)) {
      plot(yrng, c(min(xps)-0.5, max(xps)+0.5), xlab = "", ylab = "", type = "n", axes = FALSE, yaxs = "i", xlim = xlm)
    } else {
      plot(yrng, c(min(xps)-0.5, max(xps)+0.5), xlab = "", ylab = "", type = "n", axes = FALSE, yaxs = "i")
    }
    axis(1)
    if(groupbytrait) {
      axis(axisside, at = xps, labels = trt.labels_data, las = 2)
    } else {
      axis(axisside, at = xps, labels = trait.labels_data, las = 2)
    }
  } else {
    if(!is.null(xlm)) {
      plot(c(min(xps)-0.5, max(xps)+0.5), yrng, xlab = "", ylab = "", type = "n", axes = FALSE, xaxs = "i", xlim = xlm)
    } else {
      plot(c(min(xps)-0.5, max(xps)+0.5), yrng, xlab = "", ylab = "", type = "n", axes = FALSE, xaxs = "i")
    }
    axis(2)
    if(groupbytrait) {
      axis(axisside, at = xps, labels = trt.labels_data, las = 2)
    } else {
      axis(axisside, at = xps, labels = trait.labels_data, las = 2)
    }
  }
  box()
  
  if(!is.null(group_colors)) {
    for(i in 1:length(v1)) {
      if(flipaxes) {
        polygon(rep(c(yrng[1]-diff(range(yrng)), yrng[2]+diff(range(yrng))), each=2),
                xps[i]+c(-0.5, 0.5, 0.5, -0.5),
                col = group_colors[i], border = NA)
      }
    }
  }
  
  if(flipaxes) {
    abline(v=0, lty=2)
  } else {
    abline(h=0, lty=2)
  }
  
  for(i in 1:length(v1)) {
    for(j in 1:length(v2)) {
      if(groupbytrait) {
        ps = which(toplot_data$trait == trait_levels[j] & 
                   toplot_data$trt_type2 == trt_type2_levels[i])
      } else {
        ps = which(toplot_data$trait == trait_levels[i] & 
                     toplot_data$trt_type2 == trt_type2_levels[j])
      }
      
      xv = xps[i] + dx[j]
      yv = c(toplot_data$Estimate[ps]-toplot_data$SE[ps],
             toplot_data$Estimate[ps],
             toplot_data$Estimate[ps]+toplot_data$SE[ps],
             qnorm(p_alpha, toplot_data$Estimate[ps], toplot_data$SE[ps]),
             qnorm(1-p_alpha, toplot_data$Estimate[ps], toplot_data$SE[ps]))
      
      if(flipaxes) {
        segments(yv[1], xv, yv[3], xv,
                 col = collst[j], lwd = lwd_lines, lend = 2)
        segments(yv[4], xv, yv[5], xv,
                 col = collst[j], lwd = lwd_lines*lwd_diff_fact, lend = 2)
        segments(yv[c(4,5)], xv+whisker_length, yv[c(4,5)], xv-whisker_length,
                 col = collst[j], lwd = lwd_lines*lwd_diff_fact, lend = 2)
        points(yv[2], xv, col = collst[j], pch = 16, cex = 0.7)
      } else {
        segments(xv, yv[1], xv, yv[3],
                 col = collst[j], lwd = lwd_lines, lend = 2)
        segments(xv, yv[4], xv, yv[5],
                 col = collst[j], lwd = lwd_lines*lwd_diff_fact, lend = 2)
        segments(xv-whisker_length, yv[c(4,5)], xv+whisker_length, yv[c(4,5)],
                 col = collst[j], lwd = lwd_lines*lwd_diff_fact, lend = 2)
        points(xv, yv[2], col = collst[j], pch = 16, cex = 0.7)
      }
      ## add significance marker
      if(sum(ps)>0) {
        if(flipaxes) {
          if(all(yv[1:3]<0)) {
            text(yv[4], xv+sigadj, sig[ps], pos = 2,
                 cex = 0.8, offset = 0.1)
          } else {
            text(yv[5], xv+sigadj, sig[ps], pos = 4,
                 cex = 0.8, offset = 0.1)
          }
        } else {
          text(xv+sigadj, yv[3], sig[ps], pos = 3,
              cex = 0.8, offset = 0)
        }
      }
      
      #add n
      if(flipaxes) {
        if(!is.null(n_table)) {
          row_ps = which(row.names(n_table)==v2[j])
          col_ps = which(colnames(n_table)==v1[i])
          tmp_xrng = range(c(yrng, xlm))
          xps_text = tmp_xrng[2]-diff(tmp_xrng)*0.03
          yps_text = xv+yadj_text
            
          text(xps_text, yps_text, paste("(", n_table[row_ps, col_ps], ")", sep =""), cex=ncex)
        }
      }
      
    }
  }
  if(flipaxes) {
    abline(h = seq(1.5, length(v1)-0.5, by = 1), col = "darkgrey")
  } else {
    abline(v = seq(1.5, length(v1)-0.5, by = 1), col = "darkgrey")
  }
  
  if(flipaxes) {
    if(groupbytrait) {
      legend(mean(yrng)+legend_xadj, max(xps)+diff(range(xps))*legend_yadj, legend = trait.labels_data, lty = 1, lwd = lwd_lines*1.5,
             col = collst, ncol = ceiling(length(trait.labels_data)/2), xpd = NA, xjust = 0.5,
             seg.len = legend_line_length, cex = legend_cex, x.intersp = 0.5, text.width = legend_textwidth, bty = "n")
    } else {
      legend(mean(yrng)+legend_xadj, max(xps)+diff(range(xps))*legend_yadj, legend = trt.labels_data, lty = 1, lwd = lwd_lines*1.5,
             col = collst, ncol = ceiling(length(trt.labels_data)/2), xpd = NA, xjust = 0.5,
             seg.len = legend_line_length, cex = legend_cex, x.intersp = 0.5, text.width = legend_textwidth, bty = "n")
    }
    mtext(xtit, 1, line = 2.4)
  } else {
    if(groupbytrait) {
      legend(mean(xps)+legend_xadj, yrng[2]+diff(range(yrng))*legend_yadj, legend = trait.labels_data, lty = 1, lwd = lwd_lines*1.5,
             col = collst, ncol = ceiling(length(v2)/2), xpd = NA, xjust = 0.5,
             seg.len = legend_line_length, cex = legend_cex, x.intersp = 0.5, text.width = legend_textwidth, bty = "n")
    } else {
      legend(mean(xps)+legend_xadj, yrng[2]+diff(range(yrng))*legend_yadj, legend = trt.labels_data, lty = 1, lwd = lwd_lines*1.5,
             col = collst, ncol = ceiling(length(v2)/2), xpd = NA, xjust = 0.5,
             seg.len = legend_line_length, cex = legend_cex, x.intersp = 0.5, text.width = legend_textwidth, bty = "n")
    }
    mtext(xtit, 2, line = 2.4)
  }
}

trt.labels=c(co2="CO2", drought="Drt", irrigation="Irg.", temp="Temp.", n="N", p="P", multnuts="Mult. Nut.","all mult" ="Interact.")
trait.labels=c(LDMC="LDMC", leaf_C.N="C:N",
               PlantHeight="Plant Height", "Rooting Depth"="Rooting Depth",
               "Seed Mass"="Seed Mass", SLA="SLA")


pdf("traits_by_treat_cont.pdf", width = 5.2, height=10)
make_boxplot(toplot_data = toplot,
                        trt.labels_data = trt.labels,
                        trait.labels_data=trait.labels,
                        groupbytrait = FALSE,
                        legend_line_length = 1.4,
                        legend_textwidth = 0.012,
                        legend_yadj = 0.175,
                        legend_xadj = -0.002,
                        p_alpha = 0.05,
                        lower_margin = 6.8,
                        sigadj = -0.03,
                        traitorder = rev(c(2,1,6,4,3,5)),
                        group_colors = adjustcolor(rev(c("darkgreen",
                                                         "darkgreen",
                                                         "darkgreen",
                                                         "blue",
                                                         "blue",
                                                         "orange")), alpha.f = 0.08))
dev.off()


####Categorical data

table(catTraits$photo_path)#drop possible hybird, parasitic
table(catTraits$growth_form)#drop cactus
table(catTraits$clonal)#drop uncertain
table(catTraits$n_fixation)
table(catTraits$mycorrhizal)#drop unk
table(catTraits$lifespan)#drop uncertain

alldat_cat<-dcidiff_models%>%
  right_join(catTraits)%>%
  gather(growth_form:photo_path, key="trait", value="value")%>%
  na.omit()


tmp = alldat_cat
tmp = tmp[!is.na(tmp$diff),]
tmp$trait_value = paste(tmp$trait, tmp$value, sep = "_")
cat_n_data = table(unique(tmp[,c("species_matched", "trt_type2", "trait_value")]))
cat_n_data = apply(cat_n_data, 2:3, sum)

value2 = colnames(cat_n_data)
tmp=ifelse(value2=='clonal_no', "Non-Clonal",
  ifelse(value2=='clonal_yes', 'Clonal',
  ifelse(value2=='lifespan_ann.bien', 'Annual/Bienn.',
  ifelse(value2=='lifespan_perennial', 'Perennial',
  ifelse(value2=='mycorrhizal_yes', 'Mycorr.',
  ifelse(value2=='mycorrhizal_no', 'Non-Mycorr.',
  ifelse(value2=='n_fixation_yes', 'N-Fixer',
  ifelse(value2=='n_fixation_no', 'Non-N-Fixer',
  ifelse(value2=='photo_path_C3', "C3",
  ifelse(value2=='photo_path_C4', 'C4', value2))))))))))
tmp = gsub("growth_form_", "", tmp, fixed = TRUE)
#tmp = paste(toupper(substr(tmp,1,1)), substr(tmp,2,99), sep = "")
colnames(cat_n_data) = tmp
write.csv(cat_n_data, "cat_n_data.csv")



#photo path
mpp<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="photo_path"&value!="parasitic"&value!="possible hybrid"&value!="CAM"))
summary(mpp)

plot.mpp<-as.data.frame(emmeans(mpp, ~ value*trt_type2)) %>% 
  mutate(trait="photo_path")

#lifespan
ml<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="lifespan"&value!='uncertain'))
summary(ml)

plot.ml<-as.data.frame(emmeans(ml, ~ value*trt_type2)) %>% 
  mutate(trait="lifespan")

# #growth form
# mgf<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="growth_form"&value!="cactus"))
# summary(mgf)
# 
# plot.mgf<-as.data.frame(emmeans(mgf, ~ value*trt_type2)) %>% 
#   mutate(trait="growth_form")
# 
###clonaltiy
mc<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="clonal"&value!='uncertain'))

plot.mc<-as.data.frame(emmeans(mc, ~ value*trt_type2))%>%
  mutate(trait="Clonal")

##nfixer
mnf<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="n_fixation"))

plot.nf<-as.data.frame(emmeans(mnf, ~ value*trt_type2))%>%
  mutate(trait="Nfix")


##mycorhizal
mmyc<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="mycorrhizal"&value!="unk"))

plot.myc<-as.data.frame(emmeans(mmyc, ~ value*trt_type2))%>%
  mutate(trait="Myc")

##growthform
mgf<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="growth_form"&value!='cactus'&value!='fern'))
plot.mgf<-as.data.frame(emmeans(mgf, ~ value*trt_type2))%>%
  mutate(trait="GF")


toplot.cat<-plot.mpp%>%
  bind_rows(plot.ml, plot.mc, plot.nf, plot.myc, plot.mgf)


toplotesacat<-toplot.cat %>% 
  mutate(trt_type3=factor(trt_type2, levels=c("co2", "drought", "irrigation", "temp", "n", "p", "multnuts", "all mult"))) %>% 
  mutate(value2=paste(trait, value, sep="_" )) %>% 
  mutate(min=emmean-SE, max=emmean+SE) %>% 
  mutate(sig=ifelse(min>0&max>0, "*", ifelse(min<0&max<0, "*", ""))) %>% 
  mutate(Trait_name=ifelse(value2=='Clonal_no', "Non-Clonal",
                    ifelse(value2=='Clonal_yes', 'Clonal',
                    ifelse(value2=='lifespan_ann.bien', 'Annual/Bienn.',
                    ifelse(value2=='lifespan_perennial', 'Perennial',
                    ifelse(value2=='Myc_yes', 'Mycorr.',
                    ifelse(value2=='Myc_no', 'Non-Mycorr.',
                    ifelse(value2=='Nfix_yes', 'N-Fixer',
                    ifelse(value2=='Nfix_no', 'Non-N-Fixer',
                    ifelse(value2=='photo_path_C3', "C3",
                    ifelse(value2=='photo_path_C4', 'C4', "GF")))))))))))# %>% 
  #mutate(Traits2=factor(Trait_name,levels=c("Annual", 'Perennial', 'Clonal', 'Non-Clonal', 'C3', 'C4', 'Mycorr.', 'Non-Mycorr.', 'N-Fixer', 'LDMC', 'Plant Height', 'Rooting Depth', 'Seed Mass', 'SLA'))) #%>% 
  #na.omit()
head(toplotesacat)
toplotesacat$value = as.character(toplotesacat$value)
toplotesacat$Trait_name = as.character(toplotesacat$Trait_name)
toplotesacat$value[toplotesacat$Trait_name!="GF"] = toplotesacat$Trait_name[toplotesacat$Trait_name!="GF"]
toplotesacat$value = factor(toplotesacat$value)

unique(data.frame(toplotesacat$value, toplotesacat$Trait_name))


trt.labels=c(co2="CO2", drought="Drt", irrigation="Irg", temp="Temp", n="N", p="P", multnuts="Mult. Nut.","all mult" ="Interact.")
trait.labels = sort(unique(as.character(toplotesacat$value)))
tmp = trait.labels
trait.labels = paste(toupper(substr(trait.labels,1,1)),
      substr(trait.labels,2,99), sep = "")
names(trait.labels) = tmp
trait.categories = c("Lifespan", "Photo. Pathway", "Photo. Pathway",
                     "Clonality", "Growth Form",
                     "Growth Form", "Mycorr. Assoc.",
                     "N. Fixation", "Clonality", "Mycorr. Assoc.",
                     "N. Fixation", "Lifespan", "Growth Form", "Growth Form")
trait.super.categories = c("Reproduction", "Leaf Traits", "Leaf Traits",
                     "Reproduction", "Growth Form",
                     "Growth Form", "Symbiosis",
                     "Symbiosis", "Reproduction", "Symbiosis",
                     "Symbiosis", "Reproduction", "Growth Form", "Growth Form")
cbind(trait.labels, trait.categories, trait.super.categories)
toplotesacat$Estimate = toplotesacat$emmean
toplotesacat$trait = toplotesacat$value

##all cat in one fig

#ggplot(data=toplotesacat, aes(y=emmean, x=1, label=sig))+
#  geom_point()+
#  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.05)+
#  scale_y_continuous(breaks = c(-0.05, 0, 0.05))+
#  coord_flip()+
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
#  xlab("")+
#  ylab("Estimate")+
#  scale_x_continuous(limits=c(0, 2))+
#  geom_hline(yintercept=0, linetype="dashed")+
#  geom_text(nudge_x = 0.2, nudge_y = 0.01, size=10, color="red")+
#  facet_grid(trt_type3~Traits2, labeller = labeller(trt_type3=trt.labels))

##############
# new plot
##############
gcol = adjustcolor(rev(c(rep("darkgreen",2),
                         rep("blue",4),
                         rep("orange",4),
                         rep("red", 4))), alpha.f = 0.08)
tord = rev(c(2,3,
             7,10,
             8,11,
             1,12,
             4,9,
             5,6,
             13,14))

if(FALSE) {
pdf("traits_by_treat_cat.pdf", width = 5.2, height=10)
make_boxplot(toplot_data = toplotesacat,
                        trt.labels_data = trt.labels,
                        trait.labels_data=trait.labels,
                        groupbytrait = FALSE,
                        legend_line_length = 1.4,
                        legend_textwidth = 0.037,
                        legend_cex = 0.9,
                        legend_yadj = 0.1,
                        legend_xadj = -0.01,
                        p_alpha = 0.05,
                        lower_margin = 6.8,
                        sigadj = -0.09,
                        traitorder = tord,
                        group_colors = gcol,
             xtit = "Effect Size, Mean DCi diff. by. Trait",
             n_table = cat_n_data)
dev.off()
}


tord1 = rev(c(2,11,7,12,4,5,14))
tord2 = rev(c(3,8,10,1,9,6,13))
gcol_split = adjustcolor(rev(c(rep("darkgreen",1),
                         rep("blue",2),
                         rep("orange",2),
                         rep("red", 2))), alpha.f = 0.08)

pdf("traits_by_treat_cat_2col.pdf", width = 10.4, height=10)
par(mar=c(2,6.8,3.5,0.2), oma =c(3,1,0,0), mfrow=c(1,2))
make_boxplot(toplot_data = toplotesacat,
             trt.labels_data = trt.labels,
             trait.labels_data=trait.labels,
             groupbytrait = FALSE,
             legend_line_length = 1.4,
             legend_textwidth = 0.037,
             legend_cex = 1.2,
             legend_yadj = 0.185,
             legend_xadj = 0.14,
             p_alpha = 0.05,
             sigadj = -0.04,
             xlm = c(-0.12,0.12),
             traitorder = tord1,
             group_colors = gcol_split, 
             xtit = "", autopar = FALSE, n_table = cat_n_data)

par(mar=c(2,0.2,3.5,6.8))
make_boxplot(toplot_data = toplotesacat,
             trt.labels_data = trt.labels,
             trait.labels_data=trait.labels,
             groupbytrait = FALSE,
             legend_line_length = 1.4,
             legend_textwidth = 0.037,
             legend_cex = 0.9,
             legend_yadj = 10,
             legend_xadj = -0.01,
             p_alpha = 0.05,
             sigadj = -0.04,
             xlm = c(-0.12,0.12),
             traitorder = tord2,
             group_colors = gcol_split,
             xtit = "", axisside = 4, autopar = FALSE, n_table = cat_n_data)
mtext("Effect Size, Mean DCi diff. by. Trait", side = 1, outer = TRUE, line = 1.2, cex = 1.7)
dev.off()



