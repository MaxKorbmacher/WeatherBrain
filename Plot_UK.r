# Load required libraries

pacman::p_load(ggplot2,sf, rnaturalearth, rnaturalearthdata,ggrepel)

# Get UK map data
uk <- ne_countries(scale = "large", country = "united kingdom", returnclass = "sf")

# Define the three cities with coordinates
cities <- data.frame(
  name = c("Cheadle", "Newcastle", "Reading"),
  lon = c(-2.2141, -1.6178, -0.9781),
  lat = c(53.3951, 54.9783, 51.4543)
)

# Convert to sf object
cities_sf <- st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)

# Create the map
UKp = ggplot() +
  # UK boundary
  geom_sf(data = uk, fill = "#f0f0f0", color = "#666666", linewidth = 0.5) +
  # Highlighted cities
  geom_sf(data = cities_sf, color = "#d62828", size = 4, alpha = 0.8) +
  geom_sf(data = cities_sf, color = "#d62828", size = 2, alpha = 1) +
  # City labels
  geom_text_repel(
    data = cities,
    aes(x = lon, y = lat, label = name),
    size = 4.5,
    fontface = "bold",
    color = "#2d3142",
    bg.color = "white",
    bg.r = 0.15,
    box.padding = 0.5,
    point.padding = 0.5
  ) +
  # Styling
  coord_sf(xlim = c(-8, 2), ylim = c(49.5, 59)) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#e8f4f8", color = NA),
    panel.grid = element_line(color = "#d0e8f0", linewidth = 0.3),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#555555", margin = margin(b = 15)),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) + theme_blank() + xlab("") + ylab("")
  # + labs(
  #   title = "Featured UK Cities",
  #   subtitle = "Cheadle, Newcastle, and Reading"
  # )


ggsave("UKplot.pdf",UKp, height = 6, width = 4)
