# Load required libraries
library(jsonlite)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(tools)

# Get list of zip files
zip_files <- list.files(pattern = "\\.zip$", full.names = TRUE)

# Prepare to collect plots
all_plots <- list()

for (zip_file in zip_files) {
  
  # Create a fresh temporary folder for this zip file
  tmp_dir <- tempfile("unzipped_")
  dir.create(tmp_dir)
  
  # Unzip to temporary directory
  unzip(zip_file, exdir = tmp_dir)
  
  # Identify the full_data and job_request files
  full_data_file <- list.files(
    tmp_dir,
    pattern = "_full_data_0\\.json$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  job_request_file <- list.files(
    tmp_dir,
    pattern = "_job_request\\.json$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Proceed only if both required files are found
  if (length(full_data_file) >= 1 && length(job_request_file) >= 1) {
    
    # Read JSON data
    pae_json <- fromJSON(full_data_file[1], simplifyVector = FALSE)
    job_json <- fromJSON(job_request_file[1], simplifyVector = FALSE)
    
    # Get sequence lengths
    seq1_length <- nchar(job_json$sequences[[1]]$proteinChain$sequence)
    seq2_length <- nchar(job_json$sequences[[2]]$proteinChain$sequence)
    
    # Extract PAE matrix
    pae_matrix <- matrix(
      unlist(pae_json$pae),
      nrow = length(pae_json$pae),
      byrow = TRUE
    )
    
    n_res <- nrow(pae_matrix)
    rownames(pae_matrix) <- 1:n_res
    colnames(pae_matrix) <- 1:n_res
    
    # Melt matrix into long format
    pae_df <- melt(pae_matrix)
    colnames(pae_df) <- c("Residue1", "Residue2", "PAE")
    
    # Convert columns to numeric
    pae_df$Residue1 <- as.numeric(as.character(pae_df$Residue1))
    pae_df$Residue2 <- as.numeric(as.character(pae_df$Residue2))
    pae_df$PAE <- as.numeric(pae_df$PAE)
    
    # Generate plot title from zip file name
    plot_title <- paste("PAE:", file_path_sans_ext(basename(zip_file)))
    
    # Create plot
    p <- ggplot(pae_df, aes(x = Residue1, y = Residue2, fill = PAE)) +
      geom_tile() +
      scale_fill_viridis_c() +
      coord_fixed() +
      theme_minimal(base_size = 6) +
      labs(
        title = plot_title,
        x = "Residue",
        y = "Residue",
        fill = "PAE (Å)"
      ) +
      geom_hline(yintercept = seq1_length, color = "black", linewidth = 0.3) +
      geom_vline(xintercept = seq1_length, color = "black", linewidth = 0.3)
    
    # Store plot
    all_plots[[length(all_plots) + 1]] <- p
    
  } else {
    message("Skipping ", basename(zip_file), 
            ": required JSON files not found.")
  }
  
  # Clean up temporary directory for this zip file
  unlink(tmp_dir, recursive = TRUE)
}

# Save all plots into separate PNGs with 18 plots per image (6x3 layout)
plots_per_page <- 18
page_index <- 1

if (length(all_plots) > 0) {
  for (i in seq(1, length(all_plots), by = plots_per_page)) {
    
    page_plots <- all_plots[i:min(i + plots_per_page - 1, length(all_plots))]
    png_filename <- paste0("PAE_Page_", page_index, ".png")
    
    png(
      filename = png_filename,
      width = 11,
      height = 8.5,
      units = "in",
      res = 300
    )
    
    grid.arrange(grobs = page_plots, ncol = 6, nrow = 3)
    dev.off()
    
    page_index <- page_index + 1
  }
  
  message("Finished. Generated ", length(all_plots), " plots across ", page_index - 1, " PNG file(s).")
  
} else {
  message("No plots were generated. Check whether the zip files contain the expected JSON files.")
}