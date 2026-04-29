#Calculates "Interaction Score" from predictions in the same folder

library(jsonlite)
library(tools)

# Convert PAE JSON to CSV
convertJSONtoCSV <- function(json_file_path, csv_file_path) {
  json_data <- fromJSON(json_file_path, flatten = TRUE)
  write.csv(json_data$pae, file = csv_file_path, row.names = FALSE)
} 

# Calculate protein lengths from job_request.json
getProteinLengths <- function(job_request_path) {
  job <- fromJSON(job_request_path, flatten = TRUE)
  
  # Access the nested data frame
  if (!is.null(job$sequences) &&
      is.data.frame(job$sequences[[1]]) &&
      "proteinChain.sequence" %in% names(job$sequences[[1]])) {
    
    seq_df <- job$sequences[[1]]
    
    if (nrow(seq_df) < 2) {
      warning(paste("Only one sequence found in:", job_request_path))
      return(NULL)
    }
    
    prey_size <- nchar(seq_df$proteinChain.sequence[1])
    bait_size <- nchar(seq_df$proteinChain.sequence[2])
    
    return(list(prey = prey_size, bait = bait_size))
  }
  
  warning(paste("Unexpected sequence structure in:", job_request_path))
  return(NULL)
}


# Calculate the contacts metric
calculateContacts <- function(csv_file_path, prey_size, bait_size) {
  check <- read.csv(csv_file_path)
  size <- dim(check)
  check2 <- check[, 3:(size[2] - 2)]  # trims non-PAE columns if needed
  check3 <- check2[(prey_size + 1):size[1], 1:prey_size]
  background <- matrix(31.70, prey_size, bait_size)
  contacts <- sum(background) - sum(check3)
  return(contacts)
}

# Loop over zip files
zip_files <- list.files(pattern = "*.zip", full.names = TRUE)
results <- data.frame(Gene = character(), Contacts = numeric(), Prey = numeric(), Bait = numeric(), stringsAsFactors = FALSE)

for (zip_file in zip_files) {
  gene_name <- sub("fold_", "", file_path_sans_ext(basename(zip_file)))
  temp_dir <- file.path(tempdir(), gene_name)
  dir.create(temp_dir)
  
  # Unzip and locate files
  unzip(zip_file, exdir = temp_dir)
  json_file_path <- list.files(temp_dir, pattern = "_full_data_0\\.json$", full.names = TRUE)
  job_file_path <- list.files(temp_dir, pattern = "_job_request\\.json$", full.names = TRUE)
  csv_file_path <- file.path(temp_dir, "Rank1_Scores.csv")
  
  unzip(zip_file, exdir = temp_dir)
  
  all_files <- list.files(temp_dir, recursive = TRUE)
  
  print(zip_file)
  print(all_files)
  
  json_file_path <- list.files(
    temp_dir,
    pattern = "_full_data_0\\.json$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  job_file_path <- list.files(
    temp_dir,
    pattern = "(^|_)job_request\\.json$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  print(json_file_path)
  print(job_file_path)
  
  # Only process if both required files are present
  if (length(json_file_path) == 1 && length(job_file_path) == 1) {
    lengths <- getProteinLengths(job_file_path)
    convertJSONtoCSV(json_file_path, csv_file_path)
    contacts <- calculateContacts(csv_file_path, lengths$prey, lengths$bait)
    results <- rbind(results, data.frame(Gene = gene_name, Contacts = contacts, Prey = lengths$prey, Bait = lengths$bait))
    print(paste("Processed:", gene_name))
  } else {
    warning(paste("Missing files in zip:", zip_file))
  }
  
  unlink(temp_dir, recursive = TRUE)
}


# Export results
write.csv(results, "SPRYonly_Contacts.csv", row.names = FALSE)

results
