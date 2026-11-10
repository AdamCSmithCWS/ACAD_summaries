## ebird trend summaries for ACAD regions
## 
## 



## spatial summary of eBird trends
## 
## following the methods described here to summarise trends with uncertainty
## https://ebird.github.io/ebirdst/articles/trends.html#applications-regional 
#


library(ebirdst)
library(tidyverse)
library(bbsBayes2)
library(terra)
library(sf)

# species-level summaries of available trends
tr <- ebirdst::ebirdst_runs |> 
  filter(has_trends, # selects only the species with trends
         !grepl("example",species_code)) # drops the example species data 

# regions for trend summary
regions1 <- read_sf("data/ACAD_non-marine_regions_2026.shp") |> 
  group_by(region) |> 
  summarise() 

#transform to equal area projection and re-calculate area
regions <- regions1 |> 
  st_transform(crs = "ESRI:102008") |> 
  mutate(area_km2 = as.numeric(st_area(regions1)/1e6))

region_names <- regions |> 
  st_drop_geometry() |> 
  as.data.frame()


# viewing a potential mismatch between some polygons ----------------------


tmp <- regions |> 
  filter(region %in% c("MXCL","Guatemala"))


bb <- st_bbox(tmp)

tst <- ggplot()+
  geom_sf(data = regions,
          aes(colour = region),
          fill = NA) +
  coord_sf(xlim = c(bb["xmin"]+500000,bb["xmax"]),
           ylim = c(bb["ymin"],bb["ymax"]-500000))

pdf("overlap.pdf",
    width = 17, height = 11)
print(tst)
dev.off()

## end mismatch viewing



plot_overlap <- TRUE # generates a pdf map showing the distribution of trend grid-cell centres to polygons

for(i in 1:nrow(tr)){
  
  species <- unname(unlist(tr[i,"species_code"]))
  
  
  path <- try(ebirdst_download_trends(species),
              silent = TRUE)
  
  if(class(path) == "try-error"){next}
  trends <- load_trends(species,
                        fold_estimates = TRUE) # fold estimates allow uncertainty
  


# selecting the point centres to simplify the spatial overlay -------------
  trend_points <- trends |> 
    group_by(srd_id,longitude,latitude) |> 
    summarise(.groups = "drop")
  
  trends_sf <- st_as_sf(trend_points, 
                        coords = c("longitude", "latitude"), 
                        crs = 4326) %>% 
    sf::st_transform(.,st_crs(regions))

  # Spatial join --------------------------------------------------------
  # joins the point centres
  trends_regional_join_tst <- trends_sf %>% 
    st_join(regions,join = st_within) |> 
    select(-area_km2)

  
  
# Linking coastal trends that are within 15km of polygon ------------------
  linked <- trends_regional_join_tst |> 
    filter(!is.na(region)) # trends that have been linked to a region
  
 missing <- trends_regional_join_tst |> 
    filter(is.na(region))  # trends that are not linked to a region
  
 nearby1 <- st_is_within_distance(missing,
                              regions,
                              dist = 15000,
                              sparse = TRUE) # identifying trend points within 15km of region polygon
 
 missing_in_15 <- missing[which(lengths(nearby1) > 0),] #dropping the points outside of 15km distance
 
 nearby <- st_nearest_feature(missing_in_15,
                       regions) # identifying the nearest polygon for each trend point
 
 linked_miss <- missing_in_15 |> 
   mutate(region = region_names[nearby,"region"]) #adding region information to trends that missed
 
 trends_regional_join <- linked |> # recombining the newly linked trends with the overlap trends
   bind_rows(linked_miss)
  
 # calculating the number of grid cells joined to each region
 # useful below for a coarse coverage estimate of the region (although not of the species range in the region)
 n_grids <- trends_regional_join |> 
   group_by(region) |> 
   summarise(n_grid_cells = n(),
             .groups = "drop")
 
 if(plot_overlap){
 pdf(paste0("check_plots/",species,"region_trend_overlap.pdf"),
     width = 8.5,
     height = 11)
 tst <- ggplot()+
    geom_sf(data = trends_regional_join,
            aes(colour = region))+
    geom_sf(data = regions,
            fill = NA)+
   scale_colour_viridis_d()
  print(tst)
  dev.off()
 } 
  
 trends_regional_join <- trends_regional_join |> #list of which region is associated with each grid cell
   st_drop_geometry()
  

# joining the spatial associations back to the full trend data ------------
 trends_region <- trends %>% 
   st_drop_geometry() |> 
    inner_join(trends_regional_join,
               by = "srd_id") |>  
    group_by(region,fold) %>%
    summarise(weighted_abd_ppy = sum(abd * abd_ppy) / sum(abd), # trend weighted by mean abundance in cell
              .groups = "drop") %>%                           # summarised for each of the 100 folds
    group_by(region) %>% # now summarising across all folds
    summarise(abd_ppy_median = median(weighted_abd_ppy, na.rm = TRUE),
              abd_ppy_lower = quantile(weighted_abd_ppy, 0.05, na.rm = TRUE), #90% CIs
              abd_ppy_upper = quantile(weighted_abd_ppy, 0.95, na.rm = TRUE),
              .groups = "drop") %>% 
    arrange(abd_ppy_median)
  
  
  trends_region <- trends_region |> 
    inner_join(n_grids,
               by = "region") |> 
    inner_join(region_names,
               by = "region") |> 
    mutate(approx_area_w_trend = n_grid_cells*(27^2), # area of 27x27km grid cell
           approx_p_w_trend = approx_area_w_trend/area_km2) # approx area w trends, but naive to (biased high) coastal regions
  
  
  saveRDS(trends_region,paste0("species_trend_summaries/",species,"_eBird_trends_by_region.rds"))
  print(round(i/nrow(tr),2))
}




# compile species values estiamted above --------------------------------------------------
# 
 tr_join <- tr |> 
  select(species_code,
           taxon_concept_id,
           scientific_name,
           common_name,
           taxonomic_order,
           contains("trend"),
         rsquared,
         beta0)

all_eBird_trends <- NULL

for(i in 1:nrow(tr)){
  
  species <- unname(unlist(tr[i,"species_code"]))
  
  tmp <- readRDS(paste0("species_trend_summaries/",species,"_eBird_trends_by_region.rds")) |> 
    mutate(species_code = species) |> 
    select(species_code,
           region,
           abd_ppy_median,
           abd_ppy_lower,
           abd_ppy_upper,
           n_grid_cells,
           area_km2,
           approx_area_w_trend,
           approx_p_w_trend) 
  
  all_eBird_trends <- bind_rows(all_eBird_trends,tmp)
  
}

all_eBird_trends <- all_eBird_trends |> 
  left_join(tr_join,
            by = "species_code") |> 
  rename(eBird_trend = abd_ppy_median,
         eBird_trend_lower_80p_ci = abd_ppy_lower,
         eBird_trend_upper_80p_ci = abd_ppy_upper,
         number_trend_grid_cells_in_region = n_grid_cells,
         area_of_region_km2 = area_km2,
         approximate_area_w_eBird_trends = approx_area_w_trend,
         approximate_proportion_region_area_w_eBird_trends = approx_p_w_trend) |> 
  rowwise() |> 
  mutate(approx_proportion_region_covered = min(1,approximate_proportion_region_area_w_eBird_trends)) |> 
  relocate(species_code,
           common_name,
           scientific_name,
           region,
           starts_with("eBird_trend"),
           approx_proportion_region_covered,
           taxon_concept_id,
           taxonomic_order,
           starts_with("trends"))



# exploring the area summaries as a coarse sanity check -------------------



reg_sum <- all_eBird_trends |> 
  group_by(region) |> 
  summarise(n_sp = n(),
            n_sp_over_100_percent_covered = length(which(approximate_proportion_region_area_w_eBird_trends > 1.1)),
            mean_p_covered = mean(approximate_proportion_region_area_w_eBird_trends),
            mean_area = mean(area_of_region_km2),
            max_area = max(area_of_region_km2),
            min_area = min(area_of_region_km2)) |> 
  mutate(p_sp_over_100_percent_covered = n_sp_over_100_percent_covered/n_sp)

tst <- ggplot(data = reg_sum)+
  geom_point(aes(x = mean_area,
                 y = p_sp_over_100_percent_covered,
                 colour = n_sp))
tst

tst2 <- ggplot(data = reg_sum)+
  geom_point(aes(x = mean_area,
                 y = mean_p_covered,
                 colour = n_sp))
tst2

plot_reg_sum <- regions |> 
  left_join(reg_sum,
            by = "region")

tst <- ggplot()+
  geom_sf(data = plot_reg_sum,
          aes(fill = p_sp_over_100_percent_covered))
tst


sp_sum <- all_eBird_trends |> 
  group_by(species_code) |> 
  summarise(n_regs = n(),
            n_regs_over_100_percent_covered = length(which(approximate_proportion_region_area_w_eBird_trends > 1.1)),
            mean_p_covered = mean(approximate_proportion_region_area_w_eBird_trends),
            min_p_covered_area = min(approximate_proportion_region_area_w_eBird_trends),
            max_p_covered_area = max(approximate_proportion_region_area_w_eBird_trends)) |> 
  mutate(p_regs_over_100_percent_covered = n_regs_over_100_percent_covered/n_regs)




# final export of rounded values ------------------------------------------


all_eBird_trends <- all_eBird_trends |> 
  select(-c(approximate_area_w_eBird_trends,approximate_proportion_region_area_w_eBird_trends)) |> 
  mutate(across(.cols = c(starts_with("eBird_trend"),
                          approx_proportion_region_covered),
                .fns = ~signif(.x,3)),
         across(.cols = c(area_of_region_km2),
                .fns = ~signif(.x,4)))

write_csv(all_eBird_trends,"All_ebird_trends_ACAD.csv")




