library(tidyverse)

# plot insert sizes from mapped reads samtools stats output
#
# Mark Stenglein  6/2025

# ------------------------
# import tsv
# ------------------------
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  tsv_input=args[1]
  output_directory="./"
} else {
  # if running via RStudio
  tsv_input="../results/all_insert_sizes.txt"
  output_directory="../results/"
}

df <- read.delim(tsv_input, sep="\t", header=F)
# The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs      
colnames(df) <- c("sample_id", "insert_size", "pairs_total", "inward_oriented_pairs", "outward_oriented_pairs", "other_pairs")


df_to_plot <- df 
df_to_plot <- df %>% filter(str_detect(sample_id, "PA-196|NY-"))

ggplot(df_to_plot) +
  geom_point(aes(x=insert_size, y=pairs_total), 
             shape=21, fill="darkslateblue", color="black", stroke=0.25, size=1.5) +
  theme_bw(base_size = 12) +
  scale_y_log10() +
  facet_wrap(~sample_id, ncol=1, strip.position = "right") +
  theme(strip.text.y = element_text(angle = 0)) 

ggplot(df_to_plot) +
  geom_violin(aes(x=insert_size, y=pairs_total), 
             shape=21, fill="darkslateblue", color="black", stroke=0.25, size=1.5) +
  theme_bw(base_size = 12) +
  scale_y_log10() +
  facet_wrap(~sample_id, ncol=1)

df_avg <- df %>% 
  group_by(sample_id) %>%
  summarize(mean_insert_size = mean(insert_size),
            median_insert_size = median(insert_size),
            sd_insert_size = sd(insert_size))

ggplot(df) +
  # geom_boxplot(aes(x=sample_id, y=insert_size),
               # width=0.5, outlier.shape = NA, fill=NA) +
  geom_violin(aes(x=sample_id, y=insert_size),
               fill=NA) +
  geom_jitter(aes(x=sample_id, y=insert_size),
               width=0.25, height=0, 
              shape=21, fill="darkslateblue", color="black", stroke=0.25, size=1, alpha=0.5) +
  theme_bw()




