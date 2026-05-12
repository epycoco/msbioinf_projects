# Librerie necessarie
library(tidyverse)
library(ggplot2)

plot_dotplots_by_batch <- function(data, 
                                   target, 
                                   batch_size = 25, 
                                   ncol_grid = 5) {
  
  output_dir <- paste0("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/ppt/dotplots/",target,"_dotplots")
  
  # Crea cartella di output se non esiste
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Prendi tutti i descrittori (tutte le colonne tranne il target)
  descriptors <- setdiff(names(data), target)
  
  # Numero totale di batch
  n_batches <- ceiling(length(descriptors) / batch_size)
  
  message("Total descriptors: ", length(descriptors))
  message("Batch Number: ", n_batches)
  
  for (i in seq_len(n_batches)) {
    start <- (i - 1) * batch_size + 1
    end <- min(i * batch_size, length(descriptors))
    current_desc <- descriptors[start:end]
    
    # Genera i grafici per il batch corrente
    plots <- lapply(current_desc, function(desc) {
      ggplot(data, aes_string(x = desc, y = target)) +
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

TARGETS <- c("FFV","Tc","Density","Rg")

for (t in TARGETS){
  # 0.1) Import descriptor data.frame
  target_desc_df_dir <- paste0("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/desc_corr_target/",t,"_desc_target.csv")
  target_desc_df <- read.csv(target_desc_df_dir)
  # 0.2) Draw dot plots
  #plot_dotplots_by_batch(target_desc_df, target = t)
  
  # 1) Import correlation value 
  cor_value_dir <- paste0("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/desc_corr_target/",t,"_desc_corr_with_target.csv")
  cor_values <- read.csv(cor_value_dir, row.names = 1, check.names = FALSE)
  
  # 2) Transpose the matrix
  cor_values <- t(cor_values)
  
  # 3) Convert in data.frame
  cor_val_df <- as.data.frame(cor_values)
  
  # 4) Convert in vector
  cor_vec <- tibble(
    descriptor = rownames(cor_val_df),
    correlation = as.numeric(cor_val_df[, 1])
  )
  cor_vec$descriptor <- factor(cor_vec$descriptor, levels = cor_vec$descriptor)
  
  # 5) Compute min and MAX of corr values
  min_corr <- min(cor_vec$correlation, na.rm = TRUE)
  max_corr <- max(cor_vec$correlation, na.rm = TRUE)
  
  # 6) Heatmap of correlation values
  print(ggplot(cor_vec, aes(x = descriptor, y = 1, fill = correlation)) +
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
      title = paste0(t," - Correlation with descriptors (min = ", min_corr, ", max = ", max_corr, ")"),
      fill = "R"
    ))
  
  # 7) Compute the absolute value of correlation values
  cor_vec_abs <- cor_vec %>%
    mutate(abs_corr = abs(correlation)) %>%
    arrange(desc(abs_corr)) %>%
    mutate(descriptor = factor(descriptor, levels = descriptor))
  
  # 8) Compute min and MAX of absolute corr values
  min_corr <- round(min(cor_vec_abs$abs_corr, na.rm = TRUE), 3)
  max_corr <- round(max(cor_vec_abs$abs_corr, na.rm = TRUE), 3)
  
  # 9) Heatmap of absolute correlation values
  print(ggplot(cor_vec_abs, aes(x = descriptor, y = 1, fill = abs_corr)) +
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
      title = paste0(t," - Absolute correlation with descriptors (min = ", min_corr, ", max = ", max_corr, ")"),
      fill = "|R|"
  ))
}
