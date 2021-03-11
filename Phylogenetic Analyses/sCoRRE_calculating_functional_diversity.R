### 
### Calculate Functional dispersion for CoRRE data
### 

### Trait data
#read in data
contTraits <- read.csv('Trait Data\\TRY Data\\Gap_Filled\\TRY_new.csv')%>%
  rename(species_matched=Species)%>%
  select(-X.1, -X, -Family, -Genus, -ObservationID)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()

contTraitsSubset <- contTraits%>%
  rename(ssd=X4, rooting_depth=X6, SLA=X11, leaf_C_mass=X13, leaf_N_mass=X14, leaf_P_mass=X15, stem_diameter=X21, seed_mass=X26, seed_length=X27, leaf_thickness=X46, LDMC=X47, leaf_dry_mass=X55, germination_rate=X95, leaf_length=X144, leaf_width=X145, leaf_CN=X146, stem_conduit_density=X169, stem_conduit_diameter=X281, seed_number=X138, SRL=X1080)%>%
  select(-X18, -X50, -X78, -X163, -X223, -X224, -X237, -X282, -X289, -X3112, -X3113, -X3114, -X3120)

traits <- read.csv('CoRRE data\\CoRRE data\\trait data\\sCoRRE categorical trait data - traits_complete_pre spot check_03102021.csv')%>%
  full_join(contTraitsSubset) %>%
  drop_na()

traitsOutliersRemoved <- traits %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

traitsScaled <- traitsOutliersRemoved %>% ## only scales continuous traits
  mutate_at(vars(ssd:SRL), scale)


# Read in relative abund



C:\Users\wilco\Dropbox\shared working groups\sDiv_sCoRRE_shared\CoRRE data\CoRRE data\community composition\CoRRE_RelativeAbundanceMar2021.csv

#Experimental information
# C:\Users\wilco\Dropbox\shared working groups\sDiv_sCoRRE_shared\CoRRE data\CoRRE data\community composition\CoRRE_ExperimentInfoMar2021.csv

# for matched_names (key for matching CORRE and TRY species names)
#C:\Users\wilco\Dropbox\shared working groups\sDiv_sCoRRE_shared\CoRRE data\CoRRE data\trait data\corre2trykey.csv




