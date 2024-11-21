#Done in R.
#Creates a file that contains the gene id "MSTRG" and the corresponding genome Id "Bany"
#Input is the transcriptome gtf file
#Output is a double column txt file
#Process: isolate lines with both terms, remove extra text, remove genome name extensions, remove the duplicate line

setwd("C:/Users/molen/Desktop/transcriptomics_analysis/Assemblies/Bany_MSTRGcorrespondence")
list.files("C:/Users/molen/Desktop/transcriptomics_analysis/Assemblies/Bany_MSTRGcorrespondence")

#1. Remove first column of the input file (all lines have "Bany" on the first column)
# Define the input and output file paths
input_file <- "stringtie_brain_assembly.gtf"
output_file <- "nofirstcolumn.txt"

#Open the input file and the output file
infile <- file(input_file, open = "r")
outfile <- file(output_file, open = "w")

# Process the file line by line
while (length(line <- readLines(infile, n = 1, warn = FALSE)) > 0) {
  # Split the line by tab
  fields <- unlist(strsplit(line, "\t"))
  
  # Remove the first element (first column)
  modified_line <- paste(fields[-1], collapse = "\t")
  
  # Write the modified line to the output file
  writeLines(modified_line, outfile)
}

# Close the input and output files
close(infile)
close(outfile)

#2.Keep lines with both Bany and MSTRG
# Define the input and output file paths
input_file <- "nofirstcolumn.txt"
output_file <- "MSTRG_Bany.txt"

# Read the file as a vector of lines
lines <- readLines(input_file)

# Filter lines that contain both "MSTRG" and "Bany"
filtered_lines <- grep("MSTRG", lines, ignore.case = TRUE, value = TRUE)
filtered_lines <- grep("Bany", filtered_lines, ignore.case = TRUE, value = TRUE)

# Write the filtered lines to the output file
writeLines(filtered_lines, output_file)

cat("Filtered lines saved to", output_file, "\n")

#3.Extract the "Bany" and "MSTRG" info
# Define the input and output file paths
input_file <- "MSTRG_Bany.txt"
output_file <- "MSTRG_Bany_extracted_info.txt"

# Open the input and output files
infile <- file(input_file, open = "r")
outfile <- file(output_file, open = "w")

# Process the file line by line
while (length(line <- readLines(infile, n = 1, warn = FALSE)) > 0) {
 
   # Check if the line contains both "MSTRG" and "Bany"
  if (grepl("MSTRG", line) & grepl("Bany", line)) {
    # Extract MSTRG and Bany related information using regex
    MSTRG_info <- regmatches(line, regexpr("MSTRG[^\\s]+", line))  # Extract "MSTRG" and any following characters until whitespace
    Bany_info <- regmatches(line, regexpr("Bany[^\\s]+", line))    # Extract "Bany" and any following characters until whitespace
    
    # Combine the extracted information and write to the output file
    output_line <- paste(MSTRG_info, Bany_info, sep = "\t")
    writeLines(output_line, outfile)
  }
}

# Close the input and output files
close(infile)
close(outfile)

cat("Extracted MSTRG and Bany information saved to", output_file, "\n")

#4. Extract only "Bany" and "mstrg"
# Define the input and output file paths
input_file <- "MSTRG_Bany_extracted_info.txt"  
output_file <- "MSTRG_Bany_only.txt"

# Open the input and output files
infile <- file(input_file, open = "r")
outfile <- file(output_file, open = "w")

# Process the file line by line
while (length(line <- readLines(infile, n = 1, warn = FALSE)) > 0) {
  # Check if the line contains both "MSTRG" and "Bany"
  if (grepl("MSTRG", line) & grepl("Bany", line)) {
    # Extract MSTRG information (up to the first occurrence of a semicolon or space)
    MSTRG_info <- regmatches(line, regexpr("MSTRG\\.\\d+", line))
    
    # Extract Bany information (up to the first semicolon or space)
    Bany_info <- regmatches(line, regexpr("Bany_[^\\s;\"]+", line))
    
    # Combine the extracted information and write to the output file
    output_line <- paste(MSTRG_info, Bany_info, sep = "\t")
    writeLines(output_line, outfile)
  }
}

# Close the input and output files
close(infile)
close(outfile)

cat("Extracted MSTRG and Bany information saved to", output_file, "\n")

#5.Remove the Bany name extensions (-RA to -RM)
# Define the input and output file paths
input_file <- "MSTRG_Bany_only.txt" 
output_file <- "no-Rs.txt"

# Open the input file
lines <- readLines(input_file)

# Remove the specified suffixes from all lines
cleaned_lines <- gsub("(-RA|-RB|-RC|-RD|-RE|-RF|-RG|-RH|-RI|-RJ|-RK|-RL|-RM)", "", lines)

# Write the cleaned lines to the output file
writeLines(cleaned_lines, output_file)

cat("Suffixes removed and result saved to", output_file, "\n")

#5.Remove all duplicate lines
# Define the input and output file paths
input_file <- "no-Rs.txt" 
output_file <- "MSTRG_Bany_unique.txt"

# Read all lines from the input file
lines <- readLines(input_file)

# Remove duplicate lines
unique_lines <- unique(lines)

# Write the unique lines to the output file
writeLines(unique_lines, output_file)

cat("Unique lines saved to", output_file, "\n")

