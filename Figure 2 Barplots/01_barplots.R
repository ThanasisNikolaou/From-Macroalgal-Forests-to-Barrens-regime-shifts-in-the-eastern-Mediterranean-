library(ggplot2)
library(dplyr)
library(patchwork)

dataset <- read.delim("clipboard")
dataset$id <- factor(dataset$id, levels = unique(dataset$id))

# shared theme and fill
my_theme <- theme(
  panel.background = element_rect(fill = "white", color = "black"),
  panel.grid.major = element_line(color = "grey", linewidth = 0.2),
  panel.grid.minor = element_line(color = "grey", linewidth = 0.2),
  axis.line = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.ticks.length = unit(0.3, "cm"),
  axis.ticks.y = element_line(color = "black", linewidth = 0.8),
  axis.text = element_text(size = 12, color = "black"),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  axis.title = element_text(size = 14, face = "bold")
)

region_fill <- scale_fill_manual(values = c("#7aa3cd", "#d17480"))

# function to build each bar plot
make_plot <- function(data, y_var, y_label, y_max) {
  ggplot(data, aes(x = id, y = .data[[y_var]], fill = region)) +
    geom_bar(stat = "identity", color = "black") +
    region_fill +
    my_theme +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(x = "Site", y = y_label, fill = "Region")
}

canopy <- make_plot(dataset, "canopy_avg", "Canopy cover (%)",           100)
pl     <- make_plot(dataset, "pl",         "Paracentrotus lividus (ind/m2)", 16)
al     <- make_plot(dataset, "al",         "Arbacia lixula (ind/m2)",     16)
ss     <- make_plot(dataset, "ss",         "Sarpa salpa (g/m2)",          10)
sc     <- make_plot(dataset, "sc",         "Sparisoma cretense (g/m2)",   10)
sr     <- make_plot(dataset, "sr",         "Siganus rivulatus (g/m2)",    10)
sl     <- make_plot(dataset, "sl",         "Siganus luridus (g/m2)",      10)
CHI    <- make_plot(dataset, "CHI",        "CHI",                       40)

# save individual plots
ggsave("canopy.png",            plot = canopy, dpi = 300, width = 8, height = 6, units = "in")
ggsave("paracentrotus_lividus.png", plot = pl, dpi = 300, width = 8, height = 6, units = "in")
ggsave("arbacia_lixula.png",    plot = al,  dpi = 300, width = 8, height = 6, units = "in")
ggsave("sarpa_salpa.png",       plot = ss,  dpi = 300, width = 8, height = 6, units = "in")
ggsave("sparisoma_cretense.png",plot = sc,  dpi = 300, width = 8, height = 6, units = "in")
ggsave("siganus_rivulatus.png", plot = sr,  dpi = 300, width = 8, height = 6, units = "in")
ggsave("siganus_luridus.png",   plot = sl,  dpi = 300, width = 8, height = 6, units = "in")
ggsave("CHI.png",               plot = CHI, dpi = 300, width = 8, height = 6, units = "in")

# combined figure
combined <- canopy / pl / al / ss / sc / sl / sr / CHI
ggsave("combined_species.png", plot = combined, dpi = 300, width = 15, height = 13, units = "in")

