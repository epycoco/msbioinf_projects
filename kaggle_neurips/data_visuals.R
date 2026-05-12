train_desc_df <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/train_desc.csv")

test_df <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_clean.csv")
test_df$id <- as.character(train_clean_df$id)
summary(test_df)
head (test_df)


######################### Percentuale di valori nulli per target (grafico a barre)
plot_missing_percent <- function(df, targets = c("Tg", "FFV", "Tc", "Density", "Rg")) {
  # Calcolo percentuale NA per ciascun target
  na_df <- data.frame(
    Target = targets,
    MissingPercent = sapply(df[targets], function(x) mean(is.na(x)) * 100)
  )
  
  # Plot
  ggplot(na_df, aes(x = reorder(Target, -MissingPercent), y = MissingPercent)) +
    geom_bar(stat = "identity", fill = "#8893AA") +
    geom_text(aes(label = sprintf("%.1f%%", MissingPercent)), vjust = -0.5, size = 4) +
    labs(
      title = "Percentage of NA values",
      x = "Targets",
      y = "NA (%)"
    ) +
    theme_minimal(base_size = 14)
}

plot_missing_percent(train_clean_df)


######################### BOXPLOT
plot_boxplots <- function(dataset, targets = c("Tg", "FFV", "Tc", "Density", "Rg")) {
  library(ggplot2)
  
  # Palette di colori personalizzabile (una per target)
  palette <- c("#BABAD8", "#CCB2C4", "#BFC5D6", "#E9C2BC", "#D9D9D9")
  
  # Se ci sono più target che colori → riciclo
  if (length(palette) < length(targets)) {
    palette <- rep(palette, length.out = length(targets))
  }
  
  for (i in seq_along(targets)) {
    t <- targets[i]
    color <- palette[i]
    
    df_t <- data.frame(Target = t, Value = dataset[[t]])
    
    p <- ggplot(df_t, aes(x = Target, y = Value, fill = Target)) +
      geom_boxplot(
        alpha = 0.8,
        outlier.color = "red",
        outlier.size = 2,
        color = "black",
        linewidth = 0.8,
        width = 0.3,
        coef = 1.5
      ) +
      scale_fill_manual(values = setNames(color, t)) +
      labs(
        title = paste("Boxplot", t)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold")
      )
    
    print(p)
  }
}

calc_outlier_percent <- function(df, targets = c("Tg", "FFV", "Tc", "Density", "Rg")) {
  results <- data.frame(Target = character(), OutlierPercent = numeric(), stringsAsFactors = FALSE)
  
  for (target in targets) {
    values <- df[[target]]
    values <- values[!is.na(values)]  # rimuovi NA
    
    if (length(values) > 0) {
      Q1 <- quantile(values, 0.25, na.rm = TRUE)
      Q3 <- quantile(values, 0.75, na.rm = TRUE)
      IQR_val <- Q3 - Q1
      
      lower_bound <- Q1 - 1.5 * IQR_val
      upper_bound <- Q3 + 1.5 * IQR_val
      
      n_outliers <- sum(values < lower_bound | values > upper_bound)
      percent_outlier <- (n_outliers / length(values)) * 100
      
      results <- rbind(results, data.frame(Target = target, OutlierPercent = percent_outlier))
    } else {
      results <- rbind(results, data.frame(Target = target, OutlierPercent = NA))
    }
  }
  
  return(results)
}


plot_boxplots(train_clean_df)
calc_outlier_percent(train_clean_df)



############## hist
save_histograms <- function(df, targets = c("Tg", "FFV", "Tc", "Density", "Rg"),
                            output_dir = "C:/Users/andre/Desktop/machine_learning_prjct/FINAL/ppt/image") {
  library(ggplot2)
  library(tidyr)
  # Crea la cartella se non esiste
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Palette personalizzata
  palette <- c("#BABAD8", "#CCB2C4", "#BFC5D6", "#E9C2BC", "#D9D9D9")
  palette <- rep(palette, length.out = length(targets))
  
  # Loop su ogni target
  for (i in seq_along(targets)) {
    target <- targets[i]
    color <- palette[i]
    
    # Subset dei dati per il singolo target
    df_target <- data.frame(Value = df[[target]])
    
    # Plot singolo
    p <- ggplot(df_target, aes(x = Value)) +
      geom_histogram(
        bins = 30,
        color = "black",
        fill = color,
        alpha = 0.8,
        linewidth = 0.3
      ) +
      labs(
        title = paste(target),
        x = "Values",
        y = "Frequency"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text = element_text(color = "grey15"),
        axis.title = element_text(face = "bold"),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    # Salvataggio
    filename <- file.path(output_dir, paste0("hist_", target, ".png"))
    ggsave(filename, plot = p, width = 7, height = 5, dpi = 300)
    
    # Stampa anche a schermo (opzionale)
    print(p)
  }
}

save_histograms(train_clean_df)


##########################################################################################################
##########################################################################################################
##########################################################################################################

# Crea i singoli dataset per ciascun target
train_tg <- train_clean_df %>% select(id, SMILES, Tg) %>% drop_na(Tg) 
train_ffv <- train_clean_df %>% select(id, SMILES, FFV) %>% drop_na(FFV) 
train_tc <- train_clean_df %>% select(id, SMILES, Tc) %>% drop_na(Tc) 
train_dens <- train_clean_df %>% select(id, SMILES, Density) %>% drop_na(Density) 
train_rg <- train_clean_df %>% select(id, SMILES, Rg) %>% drop_na(Rg)

# Lista dei dataset per i diversi target
dfs_trains <- list(
  Tg = train_tg,
  FFV = train_ffv,
  Tc = train_tc,
  Density = train_dens,
  Rg = train_rg
)

qqnorm(train_tg$Tg)
qqline(train_tg$Tg, col = "red")

ks.test(train_tg$Tg, "pnorm", mean(train_tg$Tg), sd(train_tg$Tg))
ks.test(train_ffv$FFV, "pnorm", mean(train_ffv$FFV), sd(train_ffv$FFV))
ks.test(train_tc$Tc, "pnorm", mean(train_tc$Tc), sd(train_tc$Tc))
ks.test(train_dens$Density, "pnorm", mean(train_dens$Density), sd(train_dens$Density))
ks.test(train_rg$Rg, "pnorm", mean(train_rg$Rg), sd(train_rg$Rg))

#shapiro.test(train_ffv$FFV)
#shapiro.test(train_tc$Tc)
#shapiro.test(train_dens$Density)
#shapiro.test(train_rg$Rg)

train_tg <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/target_train/Tg_train.csv")
train_fvv <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/target_train/FFV_train.csv")
train_tc <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/target_train/Tc_train.csv")
train_dens <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/target_train/Density_train.csv")
train_Rg <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/target_train/Rg_train.csv")

trains_dfs <- list(train_tg, train_fvv, train_tc, train_dens, train_Rg)

for (df in trains_dfs){
  summary(df)
}



library(ggplot2)
library(patchwork)
library(dplyr)

plots <- lapply(descriptors, function(desc) {
  ggplot(df, aes_string(x = desc, y = target)) +
    geom_point(alpha = 0.5, color = "#2C7BB6") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    theme_minimal(base_size = 10) +
    labs(title = desc, x = desc, y = NULL)
})

plot_dotplots_by_batch <- function(data, 
                                   target_col = "Tg", 
                                   batch_size = 25, 
                                   ncol_grid = 5, 
                                   output_dir = "C:/Users/andre/Desktop/machine_learning_prjct/FINAL/ppt/dotplots") {
  # Crea cartella di output se non esiste
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Prendi tutti i descrittori (tutte le colonne tranne il target)
  descriptors <- setdiff(names(data), target_col)
  
  # Numero totale di batch
  n_batches <- ceiling(length(descriptors) / batch_size)
  
  message("Totale descrittori: ", length(descriptors))
  message("Numero di batch: ", n_batches)
  
  for (i in seq_len(n_batches)) {
    start <- (i - 1) * batch_size + 1
    end <- min(i * batch_size, length(descriptors))
    current_desc <- descriptors[start:end]
    
    # Genera i grafici per il batch corrente
    plots <- lapply(current_desc, function(desc) {
      ggplot(data, aes_string(x = desc, y = target_col)) +
        geom_point(alpha = 0.6, color = "#2C7BB6", size = 1.2) +
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", size = 0.5) +
        theme_minimal(base_size = 10) +
        labs(title = desc, x = desc, y = NULL) +
        theme(plot.title = element_text(size = 9, hjust = 0.5))
    })
    
    # Combina i 25 grafici in una griglia 5x5
    patch <- wrap_plots(plots, ncol = ncol_grid)
    
    # Salva il batch come immagine
    file_name <- sprintf("%s/dotplots_batch_%02d.png", output_dir, i)
    ggsave(file_name, patch, width = 16, height = 12, dpi = 300)
    
    message("saved: ", file_name)
  }
}

#mordred_df <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/train_mordred.csv")
#rdkit2_df <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/train_rdkit2.csv")
#rdkit3_df <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/train_rdkit3.csv")

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# Librerie
library(tidyverse)
library(reshape2)
library(ggplot2)

desc_df <- read.csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/train_desc.csv")
desc_df$id <- as.character(desc_df$id) 
train_Rg$id <- as.character(train_Rg$id)
num_desc <- ncol(desc_df)-1

rg_desc_df <- inner_join(train_Rg[, c("id", "Rg")], desc_df, by = "id")
rg_desc_df$id <- NULL
plot_dotplots_by_batch(rg_desc_df, target_col = "Rg")

# Calcolo delle correlazioni
cor_vec <- rg_desc_df %>%
  select(-Density) %>%
  summarise(across(everything(), ~ cor(.x, rg_desc_df$Density, use = "complete.obs"))) %>%
  pivot_longer(cols = everything(), names_to = "descriptor", values_to = "correlation")

# Ordina per correlazione (da più negativa a più positiva)
cor_vec <- cor_vec %>%
  arrange(correlation) %>%
  mutate(descriptor = factor(descriptor, levels = descriptor))

min_corr <- round(min(cor_vec$correlation, na.rm = TRUE), 3)
max_corr <- round(max(cor_vec$correlation, na.rm = TRUE), 3)

# Heatmap 1D con correlazione reale (-1 → blu, 0 → bianco, +1 → rosso)
ggplot(cor_vec, aes(x = descriptor, y = 1, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste0("Density - Correlation with descriptors (min = ", min_corr, ", max = ", max_corr, ")"),
    fill = "R"
  )


# ---- Variante per valore assoluto ----
cor_vec_abs <- cor_vec %>%
  mutate(abs_corr = abs(correlation)) %>%
  arrange(desc(abs_corr)) %>%
  mutate(descriptor = factor(descriptor, levels = descriptor))

min_corr <- round(min(cor_vec_abs$abs_corr, na.rm = TRUE), 3)
max_corr <- round(max(cor_vec_abs$abs_corr, na.rm = TRUE), 3)

# Heatmap 1D ordinata per valore assoluto della correlazione
ggplot(cor_vec_abs, aes(x = descriptor, y = 1, fill = abs_corr)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste0("Density - Absolute correlation with descriptors (min = ", min_corr, ", max = ", max_corr, ")"),
    fill = "|R|"
  )



library(FactoMineR)
library(factoextra)


# 2️⃣ Seleziona i descrittori con correlazione >= 0.60 (in valore assoluto)
selected_desc <- cor_vec %>%
  filter(abs(correlation) >= 0.60) %>%
  pull(descriptor)

cat("Descrittori selezionati (|r| >= 0.60):\n")
print(selected_desc)

# 3️⃣ Crea subset solo con i descrittori selezionati
subset_df <- tg_desc_df %>%
  select(all_of(selected_desc))

# 4️⃣ Standardizza (media=0, varianza=1)
subset_scaled <- scale(subset_df)

# 5️⃣ PCA
pca_res <- prcomp(subset_scaled, center = TRUE, scale. = TRUE)

# 6️⃣ Scree plot con soglia 90% varianza cumulata
var_explained <- summary(pca_res)$importance[2, ] * 100
cum_var <- summary(pca_res)$importance[3, ] * 100
pc_num <- seq_along(var_explained)

scree_data <- data.frame(
  PC = factor(pc_num),
  Var = var_explained,
  CumVar = cum_var
)

ggplot(scree_data, aes(x = PC, y = Var)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(y = CumVar), group = 1, color = "red") +
  geom_point(aes(y = CumVar), color = "red") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = which(cum_var >= 90)[1], linetype = "dotted", color = "darkgreen") +
  labs(
    title = "Density - Scree Plot",
    subtitle = "Cumulative Variance: 90%",
    x = "Principal Component",
    y = "Varianza (%)"
  ) +
  theme_minimal()

# 7️⃣ Matrice dei loadings (coeff. degli autovettori)
loadings_df <- as.data.frame(pca_res$rotation)
loadings_df$descriptor <- rownames(loadings_df)

# Ordina per contributo decrescente nella PC1
loadings_pc1 <- loadings_df %>%
  select(descriptor, PC1) %>%
  arrange(desc(abs(PC1)))

cat("\nLoadings ordinati per PC1:\n")
print(loadings_pc1)






















