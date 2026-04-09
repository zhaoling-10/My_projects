## UnigeneFasta file for enrichment analysis
# We want to ensure that the header is on a separate line from the sequence and that the sequence lines are concatenated without spaces but still remain on a separate line from the header. 
# The script provided does exactly that, ensuring each header is followed by the corresponding sequence on a new line. 
# Here’s a detailed breakdown of the process:
# Each header line (starting with ">") will be followed by its full concatenated sequence on the next line.
# The sequence will be accumulated across multiple lines until the next header is encountered.
# Here’s the corrected R script to ensure this format:

# Step 1: Read the text file
file_path <- "Lepus_europaeus_shared_proteins.fa"
lines <- readLines(file_path)

# Step 2: Initialize variables for processing
output_lines <- c()
current_sequence <- ""

# Step 3: Process each line
for (line in lines) {
  if (startsWith(line, ">")) {
    # Add the previous sequence (if exists) to the output
    if (current_sequence != "") {
      output_lines <- c(output_lines, current_sequence)
      current_sequence <- ""
    }
    # Add the header to the output
    output_lines <- c(output_lines, line)
  } else {
    # Remove spaces from the sequence
    current_sequence <- paste0(current_sequence, gsub(" ", "", line))
  }
}

# Add the last sequence to the output
if (current_sequence != "") {
  output_lines <- c(output_lines, current_sequence)
}

# Step 4: Write the output to a new file
output_file_path <- "UnigeneLepus_europaeus_shared_proteins.fa"
writeLines(output_lines, output_file_path)

