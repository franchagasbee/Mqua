###############################################################
# MELIPONA TE PIZZA - GRÁFICO DE ELEMENTOS TRANSPONÍVEIS
# Script público/laboratorial limpo e funcional
###############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --------------------------
# Parâmetros de entrada/saída
# --------------------------
arquivo_bed <- "/home/f-chagas/Documents/Mquadrifasciata_genome/Plots/Mqua.filteredRepeats.withSize_ek.bed"
output_dir <- dirname(arquivo_bed)
base_name <- "Melipona_TE_Pizza_FINAL"

output_files <- list(
  pdf = file.path(output_dir, paste0(base_name, ".pdf")),
  jpg = file.path(output_dir, paste0(base_name, ".jpg"))
)

# --------------------------
# Cores por classe
# --------------------------
cores <- c(
  DNA="#65B0A2", LINE="#FFF59D", LTR="#D1C4E9", SINE="#FFC080",
  RC="#7FBDE3", Satellite="#BEE07F", Low_complexity="#FF8F70",
  Simple_repeat="#F4B6C8", Unknown="#A9A9A9"
)

# Classes cujas labels ficam externas
external_labels <- c("Low_complexity", "RC", "SINE", "Satellite", "Simple_repeat")

# --------------------------
# Ler e organizar dados
# --------------------------
data <- read.table(arquivo_bed, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(data) <- c("Chrom","Start","End","ClassType","Count","Strand","Extra")

data <- data %>%
  separate(ClassType, into=c("Class","Type"), sep="/", fill="right") %>%
  mutate(
    Class = ifelse(is.na(Class), Type, Class),
    Type  = ifelse(is.na(Type), Class, Type)
  )

# --------------------------
# Dados agregados por classe
# --------------------------
PD_class <- data %>%
  group_by(Class) %>%
  summarise(n=sum(Count), .groups="drop") %>%
  mutate(label_pizza = sprintf("%s %.1f%%", Class, n / sum(n) * 100)) %>%
  arrange(Class) %>%
  mutate(
    ymax = cumsum(n),
    ymin = lag(ymax, default=0),
    ymid = (ymin+ymax)/2,
    is_external = Class %in% external_labels
  )

total_n <- sum(PD_class$n)

# --------------------------
# Configurações da pizza
# --------------------------
pizza_radius <- 18.5
radius_offset <- 1.9

# Linhas e labels externas
PD_external <- PD_class %>%
  filter(is_external) %>%
  mutate(
    angle = 90 - 360*(ymid/total_n),
    angle_rad = angle * pi/180,
    x_start = pizza_radius,               # início da linha = borda da fatia
    y_start = ymid,
    x_end = pizza_radius + radius_offset, # fim da linha = label
    y_end = ymid
  )

# Aqui você adiciona o ajuste para deixar a linha da Low_complexity mais longa
PD_external <- PD_external %>%
  mutate(
    x_end = case_when(
      Class == "Low_complexity" ~ x_end + 0.5,
      Class == "RC" ~ x_end + 2.5,
      Class == "Simple_repeat" ~ x_end + 6.5, 
      Class == "Satellite" ~ x_end + 5.5,
      Class == "SINE" ~ x_end + 4.5, # linha mais longa
      TRUE                      ~ x_end
    )
  )

# --------------------------
# Construção da pizza
# --------------------------
gg_pie <- ggplot(PD_class, aes(ymax=ymax, ymin=ymin, xmax=pizza_radius, xmin=0, fill=Class)) +
  geom_rect(color="black", size =0.03) +
  coord_polar(theta="y") +
  scale_fill_manual(values=cores, guide="none") +
  theme_void() +
  ggtitle("Melipona quadrifasciata - Transposable Elements") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  # labels internas
  geom_text(
    data = PD_class %>% filter(!is_external),
    aes(x=pizza_radius/2, y=ymid, label=label_pizza),
    size=3.0
  ) +
  # linhas externas
  geom_segment(
    data = PD_external,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    color = "black", size = 0.1
  ) +
  # labels externas
  geom_text(
    data = PD_external,
    aes(x = x_end, y = y_end, label = label_pizza),
    size = 3,
    hjust = ifelse(PD_external$angle < 180, 0, 1),
    vjust = 0.5
  )

# --------------------------
# Legenda hierárquica
# --------------------------
PD <- data %>%
  group_by(Class, Type) %>%
  summarise(n = sum(Count), .groups="drop") %>%
  mutate(label_full = sprintf("%s %.1f%%", Type, n / sum(n) * 100))

line_height <- 1
y_start <- 0
legend_bars <- list()
legend_labels <- list()

for(cls in levels(factor(PD$Class))) {
  sub_df <- PD %>% filter(Class == cls) %>% arrange(Type)
  n_sub <- nrow(sub_df)
  
  legend_bars[[cls]] <- data.frame(
    Class = cls,
    y_min = y_start,
    y_max = y_start + n_sub * line_height
  )
  
  legend_labels[[cls]] <- data.frame(
    Label = cls,
    y = y_start + n_sub * line_height + line_height*0.3,
    Class = cls,
    Font = "bold"
  )
  
  legend_labels[[paste0(cls,"_sub")]] <- data.frame(
    Label = paste0("   ", sub_df$label_full),
    y = y_start + seq(line_height/2, n_sub*line_height, by=line_height),
    Class = cls,
    Font = "plain"
  )
  
  y_start <- y_start + n_sub * line_height + line_height
}

legend_bars_df <- bind_rows(legend_bars)
legend_labels_df <- bind_rows(legend_labels)

gg_legend <- ggplot() +
  geom_rect(data=legend_bars_df,
            aes(xmin=0.9, xmax=1, ymin=y_min, ymax=y_max, fill=Class),
            show.legend=FALSE) +
  geom_text(data=legend_labels_df,
            aes(x=1.05, y=y, label=Label, fontface=Font),
            hjust=0, size=3.8) +
  scale_fill_manual(values=cores, guide="none") +
  theme_void() +
  coord_cartesian(xlim=c(0.7,4))

# --------------------------
# Combinar pizza + legenda
# --------------------------
final_plot <- gg_pie + gg_legend + plot_layout(ncol=2, widths=c(3.0,1.5))

# --------------------------
# Salvar arquivos
# --------------------------
ggsave(output_files$pdf, final_plot, width=12, height=12, units="in")
ggsave(output_files$jpg, final_plot, width=12, height=12, units="in", dpi=300)
