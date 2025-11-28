library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# -----------------------------
# 0. Arquivo e cores
# -----------------------------
arquivo_bed <- "/home/f-chagas/Documents/Mquadrifasciata_genome/Plots/Mqua.filteredRepeats.withSize_ek.bed"
output_dir <- dirname(arquivo_bed)
base_name <- "Melipona_TE_Pizza_Legenda_FINAL"

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
  geom_bar(stat="identity", color="black", linewidth=0.) +   # borda fina consistente
  coord_polar(theta="y") +
  scale_fill_manual(values=cores) +
  geom_text(aes(label=label_pizza), position=position_stack(vjust=0.5), size=3) +
  theme_void() +
  ggtitle("Melipona quadrifasciata - Transposable Elements") +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        legend.position="none")

# -----------------------------
# 4. Legenda hierárquica
# -----------------------------
legend_df <- lapply(levels(PD$Class), function(cls) {
  sub_df <- PD %>% filter(Class==cls) %>% arrange(Type)
  data.frame(
    Label = c(cls, sub_df$label_full),
    Font = c("bold", rep("plain", nrow(sub_df))),
    Color = cls,
    stringsAsFactors=FALSE
  )
}) %>% bind_rows()

legend_df$y <- rev(seq_len(nrow(legend_df)))

legend_grob <- ggplot(legend_df) +
  geom_tile(aes(x=1, y=y, fill=Color), width=0.15, height=1.2, show.legend=FALSE) +
  geom_text(aes(x=1.1, y=y, label=Label, fontface=Font), hjust=0, size=3.8) +
  scale_fill_manual(values=cores) +
  theme_void() +
  coord_cartesian(xlim=c(0.7,4))

# -----------------------------
# 5. Combinar gráficos
# -----------------------------
final_plot <- plot_grid(gg_pie, legend_grob, ncol=2, rel_widths=c(2.5,1.5))

# -----------------------------
# 6. Salvar com ggsave (PDF + JPG)
# -----------------------------
ggsave(pdf_out, final_plot, width=12, height=12, units="in")                 # PDF vetorial
ggsave(jpg_out, final_plot, width=12, height=12, units="in", dpi=300)      # JPG raster
