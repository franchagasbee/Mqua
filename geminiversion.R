###############################################################
# MELIPONA TE PIZZA - GRÁFICO DE ELEMENTOS TRANSPONÍVEIS (VERSÃO FINAL CONSOLIDADA)
# Correção do posicionamento "natural" via start=0 no coord_polar
###############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --------------------------
# Parâmetros de Entrada/Saída
# --------------------------
arquivo_bed <- "/home/f-chagas/Documents/Mquadrifasciata_genome/Plots/Mqua.filteredRepeats.withSize_ek.bed"
output_dir <- dirname(arquivo_bed)
base_name <- "Melipona_TE_Pizza_FINAL_COM_GENOMA"

output_files <- list(
  pdf = file.path(output_dir, paste0(base_name, ".pdf")),
  jpg = file.path(output_dir, paste0(base_name, ".jpg"))
)

# Tamanho total do genoma de Melipona quadrifasciata (em bases)
TAMANHO_GENOMA_TOTAL <- 270207263

# --------------------------
# Paleta de Cores & Configurações Visuais
# --------------------------
cores <- c(
  DNA="#65B0A2", LINE="#FFF59D", LTR="#D1C4E9", SINE="#FFC080",
  RC="#7FBDE3", Satellite="#BEE07F", Low_complexity="#FF8F70",
  Simple_repeat="#F4B6C8", Unknown="#A9A9A9",
  TE_free="#1f78b4"
)

# Classes cujas labels são posicionadas externamente
external_labels <- c("Low_complexity", "RC", "SINE", "Satellite", "Simple_repeat")

# Parâmetros de deslocamento da linha externa (declarativo, escalável)
offsets_df <- tibble(
  Class = c("Low_complexity", "RC", "Simple_repeat", "Satellite", "SINE"),
  # Valores de ajuste adicionais à x_end (radius_offset)
  offset_x = c(0.5, 6.5, 12.5, 12.6, 12.5)
)

# --------------------------
# Leitura e Organização dos Dados
# --------------------------
data <- read.table(arquivo_bed, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(data) <- c("Chrom","Start","End","ClassType","Count","Strand","Extra")

data <- data %>%
  separate(ClassType, into=c("Class","Type"), sep="/", fill="right") %>%
  # Usa coalesce para garantir que Class e Type tenham valor mesmo que um seja NA
  mutate(
    Class = coalesce(Class, Type),
    Type = coalesce(Type, Class)
  )

# --------------------------
# Agregação e Cálculo de Fatias (Incluindo TE-free)
# --------------------------
PD_class_TE <- data %>%
  group_by(Class) %>%
  summarise(n=sum(Count), .groups="drop")

total_n_anotado <- sum(PD_class_TE$n)
# Garante que TE_free_n não é negativo
TE_free_n <- max(0, TAMANHO_GENOMA_TOTAL - total_n_anotado)

PD_class <- bind_rows(
  PD_class_TE,
  data.frame(Class="TE_free", n=TE_free_n)
) %>%
  arrange(Class) %>% # A ordenação alfabética é mantida ("natural")
  mutate(
    label_pizza = sprintf("%s %.2f%%", Class, n / sum(n) * 100),
    ymax = cumsum(n),
    # Cálculo limpo de ymin
    ymin = c(0, head(ymax, -1)),
    ymid = (ymin + ymax)/2,
    is_external = Class %in% external_labels
  )

total_n <- sum(PD_class$n)

# --------------------------
# Configuração dos Rótulos Externos
# --------------------------
pizza_radius <- 30
radius_offset <- 2

PD_external <- PD_class %>%
  filter(is_external) %>%
  mutate(
    # Cálculo do ângulo (necessário para o hjust)
    angle = 90 - 360*(ymid/total_n),
    x_start = pizza_radius,
    y_start = ymid,
    x_end = pizza_radius + radius_offset,
    y_end = ymid
  ) %>%
  # JUNÇÃO DECLARATIVA: Aplica offsets de linha usando left_join
  left_join(offsets_df, by="Class") %>%
  mutate(x_end = x_end + coalesce(offset_x, 0)) # Adiciona o offset_x, se existir

# --------------------------
# Construção do Gráfico de Pizza (gg_pie)
# --------------------------
gg_pie <- ggplot(PD_class, aes(ymax=ymax, ymin=ymin, xmax=pizza_radius, xmin=0, fill=Class)) +
  geom_rect(color="black", size=0.03) +
  # <<-- CORREÇÃO PARA GARANTIR PONTO DE PARTIDA CONSISTENTE (TOPO)
  coord_polar(theta="y", start=0) +
  scale_fill_manual(values=cores, guide="none") +
  theme_void() +
  ggtitle("Melipona quadrifasciata - Transposable Elements (Genome Coverage)") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  # Rótulos internos
  geom_text(
    data = PD_class %>% filter(!is_external),
    aes(x=pizza_radius/2, y=ymid, label=label_pizza),
    size=3.0
  ) +
  # Linhas de conexão
  geom_segment(
    data = PD_external,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
    color = "black", size = 0.1
  ) +
  # Rótulos externos
  geom_text(
    data = PD_external,
    aes(x = x_end, y = y_end, label = label_pizza),
    size = 3,
    hjust = ifelse(PD_external$angle < 180, 0, 1),
    vjust = 0.5
  )

# --------------------------
# Construção da Legenda Hierárquica (gg_legend)
# --------------------------
PD_legend_data <- data %>%
  group_by(Class, Type) %>%
  summarise(n = sum(Count), .groups="drop") %>%
  mutate(label_full = sprintf("%s %.2f%%", Type, n / sum(n) * 100))

line_height <- 1
y_start <- 0
legend_bars <- list()
legend_labels <- list()

# 1. Loop para construir a hierarquia dos TEs
for(cls in unique(PD_legend_data$Class)) {
  sub_df <- PD_legend_data %>% filter(Class == cls) %>% arrange(Type)
  n_sub <- nrow(sub_df)
  
  legend_bars[[cls]] <- data.frame(Class = cls, y_min = y_start, y_max = y_start + n_sub * line_height)
  
  # Rótulo da Classe Principal
  legend_labels[[cls]] <- data.frame(Label = cls, y = y_start + n_sub * line_height + line_height*0.3, Class = cls, Font = "bold")
  # Rótulos dos Sub-tipos
  legend_labels[[paste0(cls,"_sub")]] <- data.frame(
    Label = paste0("   ", sub_df$label_full), 
    y = y_start + seq(line_height/2, n_sub*line_height, by=line_height), 
    Class = cls, Font = "plain"
  )
  y_start <- y_start + n_sub * line_height + line_height
}

legend_bars_df <- bind_rows(legend_bars)
legend_labels_df <- bind_rows(legend_labels)

# 2. Adiciona o TE_free (Barra e Label)
TE_free_pct <- TE_free_n / TAMANHO_GENOMA_TOTAL * 100
legend_bars_df <- bind_rows(
  legend_bars_df,
  data.frame(Class = "TE_free", y_min = y_start, y_max = y_start + line_height)
)
legend_labels_df <- bind_rows(
  legend_labels_df,
  data.frame(
    Label = sprintf("TE_free %.2f%%", TE_free_pct),
    y = y_start + line_height*0.3,
    Class = "TE_free", Font = "bold"
  )
)

gg_legend <- ggplot() +
  geom_rect(
    data=legend_bars_df,
    aes(xmin=-0.5, xmax=0.0, ymin=y_min, ymax=y_max, fill=Class),
    show.legend=FALSE
  ) +
  geom_text(
    data=legend_labels_df,
    aes(x=0.05, y=y, label=Label, fontface=Font),
    hjust=0, size=3.8
  ) +
  scale_fill_manual(values=cores, guide="none") +
  theme_void() +
  coord_cartesian(xlim=c(-0.6, 7.0))

# --------------------------
# Combinação e Exportação Final
# --------------------------
final_plot <- gg_pie + gg_legend + plot_layout(ncol=2, widths=c(6.0, 1.0))

ggsave(output_files$pdf, final_plot, width=15, height=12, units="in")
ggsave(output_files$jpg, final_plot, width=15, height=12, units="in", dpi=300)