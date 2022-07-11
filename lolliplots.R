library(plotly)

get_exon_data <- function(){
  XML <- '<?xml version="1.0" encoding="UTF-8"?>
	<!DOCTYPE Query>
	<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
				
		<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
			<Attribute name = "ensembl_gene_id" />
			<Attribute name = "ensembl_transcript_id" />
			<Attribute name = "exon_chrom_start" />
			<Attribute name = "exon_chrom_end" />
			<Attribute name = "start_position" />
			<Attribute name = "end_position" />
			<Attribute name = "external_gene_name" />
			<Attribute name = "ensembl_exon_id" />
		</Dataset>
	</Query>'
  XML <- gsub('`n','', XML)
  XML <- gsub('\t','', XML)
  
  download_url <- paste0('http://www.ensembl.org/biomart/martservice?query=', XML)
  file <- 'ensembl_exon_positions/b38_downloaded'
  exists <- file.exists(file)
  
  if(!exists){
    print('Ensembl Exon positions file does not exist')
    rept <- TRUE
    while(rept){
      print(paste0('Current folder: ', getwd()))
      download_yes <- readline(prompt="Would you like to download them to folder:
                           ensembl_exon_positions/b38_downloaded
                            (if folder doesn't exist, it will be created [y/n]")
      if(tolower(download_yes) == 'y'){
        if(!file.exists('ensembl_exon_positions')){
          dir.create('ensembl_exon_positions')
        }
        download.file(download_url, destfile=file, method='wget')
        rept <- FALSE
      }
      else if(tolower(download_yes) == 'n'){
        rept <- FALSE
      }
      else{
        print('Input not understood please input "y" or "n" only')
      }
    }
  }
  return(NULL)
}

read_exon_locs <- function(file){
  HEADERS <- c('Gene stable ID',
               'Transcript stable ID',
               'Exon start (bp)',
               'Exon stop (bp)',
               'Gene start (bp)',
               'Gene end(bp)',
               'Gene name',
               'Exon Name')
  HEADERS <- gsub(' \\(', '_', HEADERS)
  HEADERS <- gsub(' ', '\\.', HEADERS)
  HEADERS <- gsub(')', '', HEADERS)
  exons <- read.table(file, sep = "\t", col.names = HEADERS)
  return(exons)
}

lolliplot_raw <- function(results, exon_info,
                          title,
                          fig_height = 500,
                          fig_width = 1400,
                          ex_start_col = 'Exon.start_bp',
                          ex_stop_col = 'Exon.stop_bp',
                          lolli_x = 'GENPOS',
                          lolli_y = 'LOG10P',
                          lolli_size = 'BETA',
                          lolli_col = 'MASK',
                          lolli_direction = 'BETA',
                          lollipop_max_size = 5,
                          lollipop_stem_width = 0.1){
  # Make the Exome rectangle points
  exon_x <- c()
  exon_y <- c()
  for(row in row.names(exon_info)){
    exon_x <- exon_x %>% c(exon_info[row, ex_start_col], exon_info[row, ex_start_col])
    exon_x <- exon_x %>% c(exon_info[row, ex_stop_col], exon_info[row, ex_stop_col])
    exon_x <- exon_x %>% c(exon_info[row, ex_start_col])
    exon_x <- exon_x %>% c(NA)
    exon_y <- exon_y %>% c(0.1, -0.1, -0.1, 0.1, 0.1, NA)
  }
  
  # Draw the figure and rectangles
  fig <- plot_ly(width=fig_width, height=fig_height,
                 type = 'scatter',
                 x = exon_x,
                 y = exon_y,
                 mode = 'lines',
                 fill = 'toself',
                 showlegend=FALSE,
                 name='Exon regions',
                 line=list(color='black', width=0.5),
                 fillcolor='rgba(100,100,255,1)')
  
  # Draw lollipops onto the exomes 
  for(color in unique(results[,lolli_col])){
    df_filt = results[which(results[,lolli_col] == color),]
    # Create data to draw the lines to the exome
    lines_x = c()
    lines_y = c()
    
    for(row in row.names(df_filt)){
      lines_y <- lines_y %>% c(df_filt[row, lolli_y] * sign(df_filt[row, lolli_direction]), 0, NA)
      lines_x <- lines_x %>% c(df_filt[row, lolli_x], df_filt[row, lolli_x], NA)
    }
    
    # Add bubble scatter and lines
    fig <- fig %>% add_trace(type = 'scatter', inherit=FALSE,
                             x = df_filt[, lolli_x],
                             y = df_filt[, lolli_y]*sign(df_filt[, lolli_direction]),
                             mode='markers',
                             marker=list(size=abs(df_filt[, lolli_size]),
                                         sizeref=2.*max(results[, lolli_size])/(lollipop_max_size^2)),
                             name=color,
                             legendgroup=color
                             )
    
    fig <- fig %>% add_trace(type = 'scatter', inherit=FALSE,
                           x = lines_x,
                           y = lines_y,
                           mode='lines',
                           showlegend=FALSE,
                           line=list(color='black', width=lollipop_stem_width),
                           hoverinfo='skip',
                           name = color,
                           legendgroup = color
                           )
  }
  # Update the figure so legend colour is visible
  #fig <- fig %>% layout(legend= list(itemsizing='constant'))
  
  # Reverse the order so that exome locations and the bubbles are on top
  fig$x$attrs <- rev(fig$x$attrs)
  
  # Add titles
  fig <- fig %>% layout(
    title=title,
    xaxis=list(title='Exon Positions'),
    yaxis=list(title=lolli_y),
    legend=list(title=lolli_col)
  )
  
  return(fig)
}

reduce_gaps <- function(gene_results, exon_info, new_gap = 10,
                        ex_start_col = 'Exon.start_bp',
                        ex_stop_col = 'Exon.stop_bp',
                        gene_pos_col = 'GENPOS'){
  # Calculate gap distances
  exon_sorted <- exon_info[order(exon_info[,ex_start_col]),]
  row.names(exon_sorted) <- NULL
  gap <- c(tail(exon_sorted[, ex_start_col], -1), head(exon_sorted[, ex_start_col], 1)) - exon_sorted[, ex_stop_col]
  gap <- gap[1:length(gap)-1]
  results_out <- gene_results
  both_col <- c(ex_start_col, ex_stop_col)
  for(i in 1:length(gap)){
    exon_sorted[(i+1):nrow(exon_sorted), both_col] <- exon_sorted[(i+1):nrow(exon_sorted), both_col] %>% apply(2, function(x){x - gap[i] + new_gap})
    my_min <- exon_sorted[i+1, ex_start_col]
    results_out[which(results_out[, gene_pos_col] > my_min), gene_pos_col] <- results_out[which(results_out[, gene_pos_col] > my_min), gene_pos_col] %>% lapply(function(x){x - gap[i] + new_gap}) %>% unlist()
  }
  out <- list()
  out$results <- results_out
  out$exons <- exon_sorted
  return(out)
}

lp.example <- function(){
  ## Important necessary files:
  exon_file <-  "ensembl_exon_positions/b38_downloaded"
  exon_file <- '/slade/home/beml201/programs/lolliplots/ensembl_exon_positions/b38_downloaded'
  associations <-  "/slade/projects/UKBB/DNA_Nexus/bw_raw_raw_regenie_burden_dnanexus_2022-02-17/Single_Variant_bw_raw_Step2_Chr15_bw_raw.regenie"
  masks <-  "/slade/projects/UKBB/DNA_Nexus/set_lists_450k_v2/annotations_chr15.txt"
  gene_name <- 'IGF1R'
  transcript_id <- 'ENST00000649865'
  
  exons <-  read_exon_locs(exon_file)
  
  selected_exons_df <- exons[which(exons[, 'Transcript.stable.ID'] == transcript_id),]
  assoc_df <- read.table(associations, sep="", comment='#', header=TRUE)
  masks_df <- read.table(masks, sep="", col.names=c('ID','TID','MASK'))
  
  # Merge gene results with masks (so that only selected exomes are used)
  assoc_results <- merge(masks_df[grep(transcript_id, masks_df[, 'TID']),], assoc_df, how='inner', on='ID' )
  
  fig1 <- lolliplot_raw(assoc_results, selected_exons_df, gene_name)
  reduced_gap <-  reduce_gaps(assoc_results, selected_exons_df)
  fig2 <- lolliplot_raw(reduced_gap$results, reduced_gap$exons, gene_name)
  
  print('Example figure of exomes in "real" size')
  print(fig1)
  print('Example figure of exomes only (introns removed)')
  print(fig2)
  #htmlwidgets::saveWidget(as_widget(fig2), "lolliplots_example.html")
}

if(sys.nframe() == 0){
  lp.example()
}