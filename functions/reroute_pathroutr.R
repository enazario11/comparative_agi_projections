#libraries
library(tidyverse)
library(pathroutr)
library(sf)


reroute_one <- function(d, land_region, vis_graph) {
  d <- arrange(d, date)
  d_trim <- prt_trim(d, land_region)
  rrt <- prt_reroute(d_trim, land_region, vis_graph, blend = TRUE)
  prt_update_points(rrt, d_trim)
}

safe_reroute <- function(d, land_region, vis_graph) {
  tryCatch(reroute_one(d, land_region, vis_graph),
           error = function(e) {
             message("Failed: id = ", d$id[1], ", rep = ", d$rep[1], ": ", conditionMessage(e))
             NULL
           })
}

run_reroute <- function(sp_filt) {
  sim_df <- sp_filt %>%
    unnest(cols = sims) %>%
    dplyr::select(id, rep, date, lon, lat)

  pts <- st_as_sf(sim_df, coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    st_transform(3857)

  bb <- st_bbox(c(xmin = -170, xmax = -100, ymin = 0, ymax = 55),
                crs = st_crs(4326)) %>%
    st_as_sfc() %>%
    st_segmentize(units::set_units(10, "km")) %>%
    st_transform(st_crs(pts))

  land <- ne_countries(scale = 10, returnclass = "sf") %>%
    st_transform(st_crs(pts))

  land_region <- st_intersection(st_union(land), bb) %>%
    st_cast("POLYGON") %>%
    st_sf()

  vis_graph <- prt_visgraph(land_region, centroids = TRUE)

  pts %>%
    group_by(id, rep) %>%
    group_split() %>%
    lapply(safe_reroute, land_region = land_region, vis_graph = vis_graph) %>%
    bind_rows()
}