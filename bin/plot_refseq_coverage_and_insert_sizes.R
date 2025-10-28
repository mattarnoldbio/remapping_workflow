#!/usr/bin/env Rscript

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # lib_dir=args[1]
  depth_input           = args[1]
  insert_size_input     = args[2]
  R_lib_dir             = args[3]
  output_dir            = "./"
} else {
  # if running via RStudio
  depth_input           = "../results/save/collected_per_base_depth.tsv"
  insert_size_input     = "../results/save/collected_insert_sizes.txt"
  R_lib_dir             = NA
  output_dir            = "../results/plot/"
}

# this library will be available in the tidyverse singularity image we are using 
# (or analogous conda env)
library(tidyverse)

# these libraries are not part of the standard tidyverse, so may have to load it
# from a specified path
# either from pipeline's R lib dir or from R environment
if (!is.na(R_lib_dir)) {
  library(rstatix, lib.loc=R_lib_dir)
  library(ggpubr, lib.loc=R_lib_dir)
  library(patchwork, lib.loc=R_lib_dir)

} else {
  # in this case assuming these will be installed
  library(rstatix)
  library(ggpubr)
  library(patchwork)
}

# read in coverge depth input
depth_df <- read.delim(depth_input, header=F, sep="\t")
colnames(depth_df) <- c("dataset", "reference_sequence", "position", "depth")

# read in insert size input
# this is the IS data from samtools stats output
insert_sizes_in <- read.delim(insert_size_input, header=F, sep="\t")
colnames(insert_sizes_in) <- c("dataset", "reference_sequence", "insert_size", "pairs_total", "inward_oriented_pairs", "outward_oriented_pairs", "other_pairs")


# ----------------------------------
# calculate average depth in windows
# ----------------------------------
# TODO: make window size configurable via a command-line arg
window_size = 10
# %/% is the integer division operator
depth_df <- depth_df %>% mutate (window = position %/% window_size)

# calculate average coverage depth in each window
df_windowed <- depth_df %>% 
  group_by(dataset, reference_sequence, window)  %>% 
  summarize(depth = mean(depth), .groups = "drop") %>% 
  mutate(position = (window*window_size) + 1) %>% 
  ungroup()

##now plot coverage data on multiple pdf pages

# these are the datasets
datasets <- depth_df %>% group_by(dataset) %>% summarize(.groups="drop") %>% pull()

# these are the dataset/refseq combinations
dataset_refseqs <- depth_df %>% group_by(dataset, reference_sequence) %>% summarize(.groups="drop")

# a function to create coverage plots for all datasets
plot_all_refseqs <- function(dataset_refseqs){
  
  # TODO: make this configurable via a command-line arg
  plots_per_page <- 8
  
  page_number <- 1
  
  pdf_list <- c()
  
  # iterate through the datasets, doing up to datasets_per_page per page
  for (i in seq(1, nrow(dataset_refseqs), plots_per_page)) {
    
    # plots_per_page at a time
    subset_datasets <- dataset_refseqs %>% filter(row_number() >= i & row_number() <= (i+(plots_per_page-1)))
    
    pdf_name <- paste0(output_dir, "/coverage_plot_page_", page_number, ".pdf")
    
    pdf_list <- c(pdf_list, pdf_name)
    
    # generate & print plot
    plot_datasets(subset_datasets, pdf_name)
    
    page_number = page_number + 1
  }
  
  # return a list of PDF filenames
  pdf_list
}

# a function to create coverage plots for a certain number of datasets
plot_datasets <- function(dataset_refseq_names, pdf_name = "test.pdf"){
  
  plots <- apply(dataset_refseq_names, 1, plot_one_refseq)
  
  page_p <- NULL

  # this uses patchwork to add the plot to the page
  for (plot in plots) {
    page_p <- page_p / plot 
  }

  # add a single x axis label at the bottom
  # since all x axes are the same
  page_p <- page_p + xlab("genome position (nt)")
  
  # output plot to console
  # print(page_p)

  # save a 1-page PDF
  ggsave(pdf_name, page_p, height=10.5, width=7.5, units="in")
}

# a function to create on coverage plot
plot_one_refseq <- function (dataset_refseq_names, number_plot_cols = 1) {
  
  dataset_to_plot <- dataset_refseq_names[1]
  refseq_to_plot  <- dataset_refseq_names[2] 
  
  debug <- 0
  if (debug) {
    dataset_to_plot <- "NY-1901"
    refseq_to_plot  <- "partitivirus_NY_1901_RNA_2"
  }

  # subset the main dataframes to get the data just for this dataset/reference_sequence
  subset_df <- df_windowed %>% 
    filter(dataset  == dataset_to_plot & reference_sequence == refseq_to_plot)
  
  # convert any depth of 0 to 1, since plotting on a log10 y scale...
  subset_df <- subset_df %>% mutate(depth = if_else(depth == 0, 1, depth))
  
  # select only necessary columns
  subset_df <- subset_df %>% select(reference_sequence, depth, position)
  
  cov_p <- ggplot(subset_df) + 
    geom_line(aes(x=position, y=depth), linewidth=0.5) +
    geom_area(aes(x=position, y=depth), fill="darkslateblue", alpha=0.5) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_y_log10() +
    xlab("Genome position (nt)") +
    ylab ("Coverage depth (x)") +
    ggtitle(NULL, subtitle = paste0(dataset_to_plot, "-", refseq_to_plot)) +
    theme(strip.text.y = element_text(angle = 0)) 
  
  cov_p
  
  # ----------------
  # insert size plot
  # ----------------
  is_to_plot <- filter(insert_sizes_in, 
    dataset  == dataset_to_plot & reference_sequence == refseq_to_plot)
  
  is_plot <- ggplot(filter(is_to_plot, pairs_total > 0)) +
    # geom_point(aes(x=insert_size, y=pairs_total)) +
    geom_col(aes(x=insert_size, y=pairs_total), width=1, color="darkslateblue") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    # scale_y_log10() +
    xlab("Insert size (nt)") +
    ylab ("Read pairs") 
    
  is_plot
   
  # return the combined plot
  cov_p + is_plot + plot_layout(widths = c(3, 1))
}

# plot all reference sequences
plot_all_refseqs(dataset_refseqs)

