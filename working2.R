library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # mais moderno que cowplot

# -----------------------------
# 0. Arquivo e cores
# -----------------------------
arquivo_bed <- "/home/f-chagas/Documents/Mquadrifasciata_genome/Plots/Mqua.filteredRepeats.withSize_ek.bed"
output_dir <- dirname(arquivo_bed)
base_name <- "Melipona_TE_Pizza_Patchwork"

cores <- c(
  DNA="#65B0A2", LINE="#FFF59D", LTR="#D1C4E9", SINE="#FFC080",
  RC="#7FBDE3", Satellite="#BEE07F", Low_complexity="#FF8F70",
  Simple_repeat="#F4B6C8", Unknown="#A9A9A9"
)

pdf_out <- file.path(output_dir, paste0(base_name, ".pdf"))
jpg_out <- file.path(output_dir, paste0(base_name, ".jpg"))

# -----------------------------
# 1. Leitura e preparação
# -----------------------------
data <- read.table(arquivo_bed, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(data) <- c("Chrom","Start","End","ClassType","Count","Strand","Extra")

data <- data %>%
  separate(ClassType, into=c("Class","Type"), sep="/", fill="right") %>%
  mutate(
    Class = ifelse(is.na(Class), Type, Class),
    Type  = ifelse(is.na(Type), Class, Type)
  )

# -----------------------------
# 2. Contagem
# -----------------------------
PD <- data %>%
  group_by(Class, Type) %>%
  summarise(n = sum(Count), .groups="drop") %>%
  mutate(label_full = sprintf("%s %.1f%%", Type, n / sum(n) * 100))

PD_pie <- PD %>%
  group_by(Class) %>%
  summarise(n = sum(n), .groups="drop") %>%
  mutate(label_pizza = sprintf("%s %.1f%%", Class, n / sum(n) * 100))

PD$Class <- factor(PD$Class, levels=PD_pie$Class)
PD_pie$Class <- factor(PD_pie$Class, levels=PD_pie$Class)

# -----------------------------
# 3. Gráfico de pizza
# -----------------------------
gg_pie <- ggplot(PD_pie, aes(x="", y=n, fill=Class)) +
  geom_bar(stat="identity", color="black", linewidth=0, show.legend = FALSE) +
  coord_polar(theta="y") +
  scale_fill_manual(values=cores) +
  geom_text(aes(label=label_pizza), position=position_stack(vjust=0.5), size=3.5) +
  theme_void() +
  ggtitle("Melipona quadrifasciata - Transposable Elements") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))

# -----------------------------
# 4. Legenda hierárquica automática
# -----------------------------
line_height <- 1
legend_bars <- list()
legend_labels <- list()
y_start <- 0

for(cls in levels(PD$Class)) {
  sub_df <- PD %>% filter(Class == cls) %>% arrange(Type)
  n_sub <- nrow(sub_df)
  
  # barra contínua da classe
  legend_bars[[cls]] <- data.frame(
    Class = cls,
    y_min = y_start,
    y_max = y_start + n_sub * line_height
  )
  
  # label da classe
  legend_labels[[cls]] <- data.frame(
    Label = cls,
    y = y_start + n_sub * line_height + line_height*0.3,
    Class = cls,
    Font = "bold"
  )
  
  # labels dos subtipos
  legend_labels[[paste0(cls,"_sub")]] <- data.frame(
    Label = paste0("   ", sub_df$label_full),
    y = y_start + seq(line_height/2, n_sub*line_height, by = line_height),
    Class = cls,
    Font = "plain"
  )
  
  y_start <- y_start + n_sub * line_height + line_height
}

legend_bars_df <- bind_rows(legend_bars)
legend_labels_df <- bind_rows(legend_labels)

gg_legend <- ggplot() +
  geom_rect(data=legend_bars_df, aes(xmin=0.9, xmax=1, ymin=y_min, ymax=y_max, fill=Class),
            color=NA, show.legend=FALSE) +
  geom_text(data=legend_labels_df, aes(x=1.05, y=y, label=Label, fontface=Font),
            hjust=0, size=3.8) +
  scale_fill_manual(values=cores, guide="none") +
  theme_void() +
  coord_cartesian(xlim=c(0.7,4))

# -----------------------------
# 5. Combinar usando patchwork
# -----------------------------
final_plot <- gg_pie + gg_legend + plot_layout(ncol=2, widths = c(2.5, 1.5))

# -----------------------------
# 6. Salvar
# -----------------------------
ggsave(pdf_out, final_plot, width=12, height=12, units="in")
ggsave(jpg_out, final_plot, width=12, height=12, units="in", dpi=300)
