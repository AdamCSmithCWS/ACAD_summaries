library(tidyverse)
library(ebirdst)
library(terra)
library(sf)
library(bbsBayes2)

re_ebird <- TRUE

sp <- "Golden-winged Warbler"


## eBird seasons
## 



tmp <- ebirdst_runs

non_res <- tmp |> 
  filter(!is_resident) |> 
  mutate(length_breeding = breeding_end-breeding_start,
         length_nonbreeding = nonbreeding_end - nonbreeding_start,
         length_postbreeding_migration = postbreeding_migration_end - postbreeding_migration_start,
         length_prebreeding_migration = prebreeding_migration_end - prebreeding_migration_start)


# eBird relative abundance ------------------------------------------------

sp_ebird <- ebirdst::get_species(sp)

qual_sel <- ebirdst_runs[which(ebirdst_runs$species_code == sp_ebird),]
breed_qual <- unname(unlist(qual_sel[,"breeding_quality"]))
resident_qual <- unname(unlist(qual_sel[,"resident_quality"]))
resident <- unname(unlist(qual_sel[,"is_resident"]))

breeding_start <- unname(unlist(qual_sel[,"breeding_start"]))
breeding_end <- unname(unlist(qual_sel[,"breeding_end"]))

season <- ifelse(resident,"resident","breeding")
yr_ebird <- as.integer(unname(unlist(qual_sel[,"status_version_year"])))

if(re_ebird){
  
  
  down <- try(ebirdst::ebirdst_download_status(sp_ebird,
                                               download_ranges = FALSE,
                                               download_abundance = TRUE,
                                               download_occurrence = FALSE,
                                               force = FALSE,
                                               pattern = "abundance_(median|upper|lower)_3km"),
              silent = TRUE)
  
}

if(class(down) == "try-error"){
  
  warning(paste("Download fail", sp))
  
  next}




### insert code to pull weekly abundance for the relevant season
### 
df_full <- df_full %>% 
  mutate(doy = yday(as_date(paste(year,month,day,sep = "-"))))


survey_days <- df_full %>% 
  select(month,day) %>% 
  distinct() %>% 
  mutate(year = yr_ebird,
         doy = as_date(paste(year,month,day,sep = "-"))) %>% 
  drop_na() ## drops routes without day of survey information (or with incorrect information, e.g., April 31)

survey_day_range <- as_date(min(survey_days$doy):max(survey_days$doy))


abd_seasonal_abundance <- ebirdst::load_raster(species = sp_ebird,
                                               resolution = "3km",
                                               #period = "seasonal",
                                               product = "abundance",
                                               metric = "median")  #3km high resolution

weeks_keep <- names(abd_seasonal_abundance)
weeks_keep <- weeks_keep[weeks_keep %in% survey_day_range]

## keep only estimates for weeks with survey data
breed_abundance <- abd_seasonal_abundance[[weeks_keep]]

## translate the weeks into breaks for the doy values
doy_weeks_keep <- c((yday(weeks_keep)-3),(yday(weeks_keep)[length(weeks_keep)]+4))




## crop relative abundance to US-Canada
bbs_strata_buffer <- bbsBayes2::load_map("bbs") %>%
  st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata_buffer, st_crs(breed_abundance))

# bbs_strata_o <- bbs_strata

bbs_strata <- bbs_strata |>
  #   select(-area_sq_km) |>
  rename(bbs_strata = strata_name) 
#   st_intersection(strata_map) |> 
#   mutate(area_km2 = as.numeric(st_area(geom))/1e6) |> 
#   filter(area_km2 > 10) |> 
#   sf::st_make_valid()


## projected bbs strata
bbs_strata_proj <- st_transform(bbs_strata, st_crs(breed_abundance))%>% 
  st_make_valid()


strata_used_proj <-  st_transform(strata_used, st_crs(breed_abundance)) %>% 
  st_make_valid()

strata_map_proj <-  st_transform(strata_map, st_crs(breed_abundance)) %>% 
  st_make_valid()


replace_zeros <- function(x){
  min_x <- min(x[which(x>0)],na.rm = TRUE)
  y <- ifelse(x == 0,min_x*0.5, # currently forces a non-zero abundance value that is half of the lowest positive abundance
              x)
  y <- ifelse(is.na(y),
              min_x*0.1, # also forces a non-zero abundance value that is 1/10th of the lowest positive abundance if missing, but also has monitoring observations
              y)
  return(y)
} 


abundance_in_strata_used_wide <- terra::extract(breed_abundance,
                                                strata_used_proj,
                                                fun = sum,
                                                na.rm = TRUE,
                                                ID = FALSE,
                                                exact = TRUE) %>% 
  mutate(strata_name = strata_used_proj$strata_name,
         stratum = strata_used_proj$stratum,
         across(.cols = starts_with(as.character(yr_ebird)),
                .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data

if(any(!is.finite(colSums(abundance_in_strata_used_wide[,weeks_keep])))){
  weeks_drop <- which(!is.finite(colSums(abundance_in_strata_used_wide[,weeks_keep])))
  
  weeks_keep <- weeks_keep[-weeks_drop]
  
  # doy_keep <- min(yday(as_date(weeks_keep))):max(yday(as_date(weeks_keep)))
  # df_full <- df_full |> 
  #   filter(doy %in% doy_keep)
  # 
  
  ## keep only estimates for weeks with survey data
  breed_abundance <- abd_seasonal_abundance[[weeks_keep]]
  
  ## translate the weeks into breaks for the doy values
  doy_weeks_keep <- c((yday(weeks_keep)-3),(yday(weeks_keep)[length(weeks_keep)]+4))
  
  
  abundance_in_strata_used_wide <- terra::extract(breed_abundance,
                                                  strata_used_proj,
                                                  fun = sum,
                                                  na.rm = TRUE,
                                                  ID = FALSE,
                                                  exact = TRUE) %>% 
    mutate(strata_name = strata_used_proj$strata_name,
           stratum = strata_used_proj$stratum,
           across(.cols = starts_with(as.character(yr_ebird)),
                  .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data
  
  
}


saveRDS(breed_abundance,paste0("data/species_relative_abundance_",yr_ebird,"/",sp_ebird,"_derived_breeding_weekly_relative_abundance.rds"))



# abundance_in_strata_used_wide <- terra::extract(breed_abundance,
#                                            strata_used_proj,
#                                            fun = sum,
#                                            na.rm = TRUE,
#                                            ID = FALSE,
#                                            exact = TRUE) %>% 
#   mutate(strata_name = strata_used_proj$strata_name,
#          stratum = strata_used_proj$stratum,
#          across(.cols = starts_with(as.character(yr_ebird)),
#                 .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data
# 


prop_in_strata_used <- abundance_in_strata_used_wide %>% 
  mutate(across(.cols = starts_with(as.character(yr_ebird)), # calculates the proportion of the total abundance in each polygon
                .fns = ~.x/sum(.x,na.rm = TRUE))) %>% 
  arrange(stratum)

abundance_in_strata_used_summary <- abundance_in_strata_used_wide %>% 
  pivot_longer(cols = starts_with(as.character(yr_ebird)),
               values_to = "rel_abundance",
               names_to = "week") %>% 
  group_by(week) %>% 
  mutate(prop_population = rel_abundance/sum(rel_abundance,na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(strata_name,stratum) %>% 
  summarise(mean_rel_abundance = mean(rel_abundance),
            median_rel_abundance = median(rel_abundance),
            min_rel_abundance = min(rel_abundance),
            max_rel_abundance = max(rel_abundance),
            sd_rel_abundance = sd(rel_abundance),
            mean_prop_population = mean(prop_population),
            median_prop_population = median(prop_population),
            min_prop_population = min(prop_population),
            max_prop_population = max(prop_population),
            sd_prop_population = sd(prop_population)) %>% 
  mutate(cv_rel_abundance = sd_rel_abundance/median_rel_abundance,
         cv_prop_population = sd_prop_population/mean_prop_population)





# Abundance in BCR-prov ---------------------------------------------------


abundance_in_bcr_prov_wide <- terra::extract(breed_abundance,
                                             bbs_strata_proj,
                                             fun = sum,
                                             na.rm = TRUE,
                                             ID = FALSE,
                                             exact = TRUE) %>% 
  mutate(bbs_strata = bbs_strata$bbs_strata,
         country = bbs_strata_proj$country,
         prov_state = bbs_strata_proj$prov_state,
         bcr = bbs_strata_proj$bcr)#,
# across(.cols = starts_with(as.character(yr_ebird)),
#        .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data


prop_in_bcr_prov <- abundance_in_bcr_prov_wide %>% 
  mutate(across(.cols = starts_with(as.character(yr_ebird)),
                .fns = ~.x/sum(.x,na.rm = TRUE))) # calculates the proportion of the weekly abundance in each cell

abundance_in_bcr_prov_summary <- abundance_in_bcr_prov_wide %>% 
  pivot_longer(cols = starts_with(as.character(yr_ebird)),
               values_to = "rel_abundance",
               names_to = "week") %>% 
  group_by(week) %>% 
  mutate(prop_population = rel_abundance/sum(rel_abundance,na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(bbs_strata, 
           country,
           prov_state,
           bcr) %>% 
  summarise(mean_rel_abundance = mean(rel_abundance),
            median_rel_abundance = median(rel_abundance),
            min_rel_abundance = min(rel_abundance),
            max_rel_abundance = max(rel_abundance),
            sd_rel_abundance = sd(rel_abundance),
            mean_prop_population = mean(prop_population),
            median_prop_population = median(prop_population),
            min_prop_population = min(prop_population),
            max_prop_population = max(prop_population),
            sd_prop_population = sd(prop_population)) %>% 
  mutate(cv_rel_abundance = sd_rel_abundance/median_rel_abundance,
         cv_prop_population = sd_prop_population/mean_prop_population)




# Abundance in all strata -------------------------------------------------



abundance_in_strata_wide <- terra::extract(breed_abundance,
                                           strata_map_proj,
                                           fun = sum,
                                           na.rm = TRUE,
                                           ID = FALSE,
                                           exact = TRUE) %>% 
  mutate(strata_name = strata_map_proj$strata_name)#,
# across(.cols = starts_with(as.character(yr_ebird)),
#        .fns = ~replace_zeros(.x))) # forces a non-zero abundance value for regions with monitoring data


prop_in_strata <- abundance_in_strata_wide %>% 
  mutate(across(.cols = starts_with(as.character(yr_ebird)),
                .fns = ~.x/sum(.x,na.rm = TRUE))) 

abundance_in_strata_summary <- abundance_in_strata_wide %>% 
  pivot_longer(cols = starts_with(as.character(yr_ebird)),
               values_to = "rel_abundance",
               names_to = "week") %>% 
  group_by(week) %>% 
  mutate(prop_population = rel_abundance/sum(rel_abundance,na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(strata_name) %>% 
  summarise(mean_rel_abundance = mean(rel_abundance),
            median_rel_abundance = median(rel_abundance),
            min_rel_abundance = min(rel_abundance),
            max_rel_abundance = max(rel_abundance),
            sd_rel_abundance = sd(rel_abundance),
            mean_prop_population = mean(prop_population),
            median_prop_population = median(prop_population),
            min_prop_population = min(prop_population),
            max_prop_population = max(prop_population),
            sd_prop_population = sd(prop_population)) %>% 
  mutate(cv_rel_abundance = sd_rel_abundance/median_rel_abundance,
         cv_prop_population = sd_prop_population/mean_prop_population)



