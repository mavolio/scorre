### Identify winners, losers, neutral, and colonizers using Dci
###
### Authors: Wilcox (wilcoxkr@gmail.com)
### Date created: Dec 13, 2022

### Set up workspace
library(tidyverse)

### Read in Dci data -- run DCi responses_treatmentsubset.R script from winners losers folder in github
### find CT_diff
rm(list=ls()[!ls() %in% c('CT_diff')])

CT_diff_working <- CT_diff %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type))

### Calculate quanitles of Dci diff
diff_quantile <- CT_diff_working %>%
  group_by(site_proj_comm, treatment) %>%
  summarize(p05 = quantile(diff, probs = 0.05, na.rm = FALSE),
            p10 = quantile(diff, probs = 0.10, na.rm = FALSE),
            p20 = quantile(diff, probs = 0.20, na.rm = FALSE),
            p80 = quantile(diff, probs = 0.80, na.rm = FALSE),
            p90 = quantile(diff, probs = 0.90, na.rm = FALSE),
            p95 = quantile(diff, probs = 0.95, na.rm = FALSE)
            )

CT_diff_full <- CT_diff_working %>%
  left_join(diff_quantile, by=c("site_proj_comm", "treatment")) %>%
  mutate(species_status_95 = ifelse(DCi==0 & treatDCi>0, "colonizer",
                                 ifelse(diff < CT_diff_full$p05, "loser",
                                        ifelse(diff > CT_diff_full$p95, "winner",
                                               "neutral")))) %>%
  mutate(species_status_90 = ifelse(DCi==0 & treatDCi>0, "colonizer",
                                    ifelse(diff < CT_diff_full$p10, "loser",
                                           ifelse(diff > CT_diff_full$p90, "winner",
                                                  "neutral")))) %>%
  mutate(species_status_80 = ifelse(DCi==0 & treatDCi>0, "colonizer",
                                    ifelse(diff < CT_diff_full$p20, "loser",
                                           ifelse(diff > CT_diff_full$p80, "winner",
                                                  "neutral")))) %>%
  mutate(site_proj_comm_trt = paste(site_proj_comm, treatment, sep="_"))


### Create tables for plotting
species_status_80_table <- as.data.frame(
  with(CT_diff_full, table(site_proj_comm_trt, species_status_80)))%>%
  spread(key=species_status_80, value=Freq)

species_status_80_table$richness <- with(species_status_80_table,
                                         colonizer+loser+neutral+winner)

species_status_90_table <- as.data.frame(
  with(CT_diff_full, table(site_proj_comm_trt, species_status_90)))%>%
  spread(key=species_status_90, value=Freq)

species_status_90_table$richness <- with(species_status_90_table,
                                         colonizer+loser+neutral+winner)


species_status_95_table <- as.data.frame(
  with(CT_diff_full, table(site_proj_comm_trt, species_status_95)))%>%
  spread(key=species_status_95, value=Freq)

species_status_95_table$richness <- with(species_status_95_table,
                                         colonizer+loser+neutral+winner)

### Plotting 
# winners versus richness
winner_richness_fig80 <- ggplot(species_status_80_table, aes(x=richness, y= winner)) +
  geom_point() +ylim(0,30) +ggtitle("p20, p80")

winner_richness_fig90 <- ggplot(species_status_90_table, aes(x=richness, y= winner)) +
  geom_point() +ylim(0,30)+ggtitle("p10, p90")

winner_richness_fig95 <- ggplot(species_status_95_table, aes(x=richness, y= winner)) +
  geom_point() +ylim(0,30)+ggtitle("p05, p95")


grid.arrange(winner_richness_fig80, winner_richness_fig90,winner_richness_fig95, nrow = 2)

## losers versus richness
loser_richness_fig80 <- ggplot(species_status_80_table, aes(x=richness, y= loser)) +
  geom_point() +ylim(0,30) +ggtitle("p20, p80")

loser_richness_fig90 <- ggplot(species_status_90_table, aes(x=richness, y= loser)) +
  geom_point() +ylim(0,30)+ggtitle("p10, p90")

loser_richness_fig95 <- ggplot(species_status_95_table, aes(x=richness, y= loser)) +
  geom_point() +ylim(0,30)+ggtitle("p05, p95")


grid.arrange(loser_richness_fig80, loser_richness_fig90,loser_richness_fig95, nrow = 2)


# winners versus richness
colonizer_richness_fig80 <- ggplot(species_status_80_table, aes(x=richness, y= colonizer)) +
  geom_point() +ylim(0,30) +ggtitle("p20, p80")

colonizer_richness_fig90 <- ggplot(species_status_90_table, aes(x=richness, y= colonizer)) +
  geom_point() +ylim(0,30)+ggtitle("p10, p90")

colonizer_richness_fig95 <- ggplot(species_status_95_table, aes(x=richness, y= colonizer)) +
  geom_point() +ylim(0,30)+ggtitle("p05, p95")


grid.arrange(colonizer_richness_fig80, colonizer_richness_fig90,colonizer_richness_fig95, nrow = 2)
### exploratory code below

print(zzz)
zz <- filter(CT_diff_working, site_proj_comm=="KNZ IRG l")


losers <- filter(CT_diff_temp, diff < diff_quantile_temp[1])$species_matched
winners <- filter(CT_diff_temp, diff > diff_quantile_temp[2])$species_matched

zz <- filter(alldat, site_code=="KNZ" & project_name=="IRG")
unique(zz$calendar_year)




CT_diff_temp <- filter(CT_diff_working, site_proj_comm == "KNZ IRG l")

ggplot(subset(CT_diff_full, site_proj_comm== 'KNZ IRG l'),
       aes(x=diff)) + geom_histogram(binwidth=.1) +
  geom_vline(xintercept=as.numeric(subset(diff_quantile, site_proj_comm=='KNZ IRG l')[c(3:8)]),col=c("red","blue","purple","purple","blue","red"))

ggplot(subset(CT_diff_full, site_proj_comm== 'KNZ IRG l'),
       aes(x=diff)) + geom_histogram(binwidth=.1) +
  geom_vline(xintercept=as.numeric(subset(diff_quantile, site_proj_comm=='KNZ IRG l')[c(3:8)]),col=c("red","blue","purple","purple","blue","red"))

