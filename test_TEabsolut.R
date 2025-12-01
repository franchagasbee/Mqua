# 1. Carregar as bibliotecas necessárias
library(tidyverse)
library(ggplot2)

# 2. Definir o caminho para o seu arquivo
arquivo_bed <- "/home/f-chagas/Documents/Mquadrifasciata_genome/Plots/Mqua.filteredRepeats.withSize_ek.bed"

# --- INÍCIO DA CORREÇÃO ---

# 3. Ler o arquivo BED (Corrigido)
# Lê o arquivo. Como o formato .bed do RepeatMasker (Earl Grey) tem 12 ou mais colunas,
# ler todas e depois selecionar (e renomear) as que interessam é mais seguro.
df_repeats <- read_delim(arquivo_bed,
                         delim = "\t",
                         col_names = FALSE,
                         # Definimos um número grande de colunas (15 c's) para garantir que todas sejam lidas como 'character'
                         col_types = paste(rep('c', 15), collapse = "")) 

# Renomear as colunas importantes e selecionar apenas elas (V1 e V4)
# V1 = Scaffold (Col. 1)
# V4 = Nome Completo da Repetição (Col. 4, ex: LTR/Gypsy)
df_repeats <- df_repeats %>%
  select(Scaffold = X1, RepeatInfo = X4) # X1 e X4 são os nomes padrão quando col_names=FALSE

# --- FIM DA CORREÇÃO ---


# 4. Preparação e Contagem dos Dados
# 4.1. Extrair o Tipo de Repetição Principal (ex: 'LTR' de 'LTR/Gypsy')
df_repeats_clean <- df_repeats %>%
  # A linha que estava dando erro de contexto (mas não a causa raiz)
  mutate(
    RepeatType = str_extract(RepeatInfo, "^[^/]+") # Extrai tudo antes do primeiro '/'
  )

# 4.2. Refinar Classes e Contar
df_summary <- df_repeats_clean %>%
  mutate(
    RepeatType = case_when(
      is.na(RepeatType) | RepeatType == "" ~ "Unknown", 
      RepeatType == "Simple" ~ "Simple_repeat", 
      RepeatType == "Low" ~ "Low_complexity", 
      TRUE ~ RepeatType
    )
  ) %>%
  group_by(Scaffold, RepeatType) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  ungroup()

# --- NOVO PASSO DE FILTRO (RECOMENDADO) ---
# 4.3. Filtrar para Manter Apenas os Scaffolds Desejados
scaffold_order_to_keep <- paste0("Mqua", str_pad(1:9, 3, pad = "0"))

df_summary <- df_summary %>%
  filter(Scaffold %in% scaffold_order_to_keep) 
# -------------------------------------------

# 5. Definir a ordem dos Scaffolds (como no seu gráfico) e Cores
scaffold_order <- paste0("Mqua", str_pad(1:9, 3, pad = "0"))
df_summary$Scaffold <- factor(df_summary$Scaffold, levels = scaffold_order)

# Cores aproximadas baseadas na sua imagem
custom_colors <- c(
  "DNA" = "#65B0A2",
  "LINE" = "#FFF59D",
  "LTR" = "#D1C4E9",
  "SINE" = "#FFC080",
  "RC" = "#7FBDE3",
  "Satellite" = "#BEE07F",
  "Low_complexity" = "#FF8F70",
  "Simple_repeat" = "#F4B6C8",
  "Unknown" = "#A9A9A9"
)

# 6. Gerar o Gráfico (ggplot2)
p <- ggplot(df_summary, aes(x = Scaffold, y = Count, fill = RepeatType)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  
  scale_fill_manual(values = custom_colors, name = "Repeat Type (Class)") +
  
  labs(
    title = "Absolute Abundance and Classification of Transposable Elements (TEs) per Scaffold",
    x = "Chromosome/Scaffold",
    y = "Amount of TEs"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_text(size = 10, face = "bold")
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# 7. Visualizar e Salvar
print(p)
ggsave("Class_Abund_TEs_bars_ekdata.pdf", units="in", width=10, height=12, dpi=300)