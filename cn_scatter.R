# Generate a linear scatterplot of copy number data for each gene
# 9/25/2014 - D Sims
library(ggplot2)
library(stringr)

args <- commandArgs( trailingOnly = TRUE )
if ( length(args) != 1 ) stop( "ERROR: You must load 1 table into this script" )
input_file <- args[1]

# Read table in so that we can capture the data from the table title
in_fh <- file( input_file, 'r' )
title <- readLines( in_fh, n=1 )

# Capture elements from the first line to be used later in the graph
regex <- '(Gender: \\w+), (Cellularity: \\d\\.\\d+), (MAPD: \\d\\.\\d+)'
elems <- str_extract( title, perl(regex) )
gtitle <- str_extract( title, perl( '(CNVs Found in [_-\\w]+)') )
gtitle <- gsub( 'Found in', 'Report for', gtitle )

# Now read the rest of the table into a df
cn_data <- data.frame( read.table( in_fh, header=TRUE ) ) 

# Create the plot
plot <- ggplot( cn_data, aes( x=Gene, y=CN )) +
        geom_point() +
        geom_hline( aes( yintercept = 2 ), linetype = "dashed", width = 0.3 ) +
        geom_hline( aes( yintercept = c(1,4)), width = 0.3 ) +
        geom_errorbar( data = cn_data, aes( ymin = CI_05, ymax = CI_95 )) +
        ggtitle( paste(gtitle, "\n", elems) ) + 
        theme_bw() +
        theme( axis.text.x = element_text( size = 6, face = "bold", angle = 90, hjust = 1 ))

filename <- gsub( 'txt', 'pdf', input_file )
ggsave( plot = plot, file = filename )
#plot
