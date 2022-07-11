# r-lolliplots
This is the lolliplots tool, converted to be able to be run in R.

A tool to create lollipop plots (lolliplots) for coding region association results

## Setup
The plotly package will need to be installed to run this script.
To install plotly, use the command ```install.packages("plotly")```

## Loading it into your own script
Lolliplots is not currently a CRAN-available package, and so needs to be imported by sourcing the *lolliplots.R* file to your current system. This can be done with the following:
```
source("lolliplots.R")
```

An absolute file path to lolliplots.R may need to be given to source the file correctly.

You should then be able to use the functions as normal.

Exon location data will need to be downloaded only once. This can be run using ```get_exon_data()```

This will open a prompt to download the biomart data.

```read_exon_locs()``` should be given the file that was downloaded using ```get_exon_data()```.

This data should then be filtered to the desired transcript before plotting:
``` exons[which(exons[, 'Transcript.stable.ID'] == transcript_id),] ```
Where **transcript_id** is the name of your desired transcript

Use ``` reduced_gap <- reduce_gaps(assoc, selected_exons_df)``` to remove introns from the association data and exon location data before passing the new data to ``` lolliplot_raw() ```.
Data should be accessed from the ```reduce_gaps()``` output using ```reduced_gap$results``` and ```reduced_gap$exons```

## Functions

### get_exon_data():
Arguments: **None**

Run this command (No need to put into a variable) to download the exon location data.
Downloads the data to a new folder in the current working directory: *ensembl_exon_positions/b38_downloaded*
This only need to be run once. In subsequent runs, this should be loaded in using ```read_exon_locs()``` through your own hard coded location.
User inputs will be requested.

### read_exon_locs(file):
|Arguments| About|
|---|---|
|**file** | string: the file location of the downloaded exon data (eg "ensembl_exon_positions/b38_downloaded") |

Function reads in the downloaded file data using pandas and assigns the appropriate headers to work with the default options of ``` lolliplot_raw() ```
Column names are as follows:
- Gene.stable.ID
- Transcript.stable.ID
- Exon.start_bp
- Exon.stop_bp
- Gene.start_bp
- Gene.end_bp
- Gene.name
- Exon.Name


### reduce_gaps(gene_results, exon_info)
|Arguments|About|Default|
|---|---|---|
|**gene_results** | data.frame: A dataframe of the association results to plot| - |
|**exon_info** | data.frame: A dataframe of the filtered exon locations| - |
|new_gap| numeric: New gap size between the exons | 10|
|ex_start_col | string: Name of column of exon_info of exome start locations | 'Exon.start_bp'|
|ex_stop_col | string: Name of column of exon_info of exome stop locations | 'Exon.stop_bp'|
|gene_pos_col | string: Name of column of gene_results of the variant positions | 'GENPOS'|

Reduce the spaces between the coding regions to a smaller gap (designated by **new_gap**). Returns one lsit with two named attributes, **results** and **exons**, the modified **gene_results** and the modified **exon_info** respectively.
The resulting dataframes can be used to make the lolliplots with the intronic regions removed.

> Note: the exon_info should be filtered to your desired transcript before inputting into ```reduce_gaps()``` or ```lolliplot_raw()```.

### lolliplot_raw(results, exon_info, title)
|Arguments|About|Default|
|---|---|---|
|**results**| data.frame: : A dataframe of the association results to plot, or the results output of ```reduce_gaps()```| - |
|**exon_info**| data.frame: : A pandas dataframe of the filtered exon locations, or the exons output of ```reduce_gaps()```| - |
|**title**| string: Title of the plot| - |
|fig_height| int: Figure height (can be adjusted using update_layout()) | 500|
|fig_width | int: Figure height (can be adjusted using update_layout()) | 1400|
|ex_start_col | string: Name of column of exon_info of exome start locations | 'Exon.start_bp'|
|ex_stop_col | string: Name of column of exon_info of exome stop locations| 'Exon.stop_bp'|
|lolli_x | string: Name of column of results df for x axis | 'GENPOS'|
|lolli_y | string: Name of column of results df for y axis | 'LOG10P'|
|lolli_size | string: Name of column of results df for size of bubbles | 'BETA'|
|lolli_col | string:  Name of column of results df for colour of bubbles | 'MASK'|
|lolli_direction | string: Name of column of reults df for direction of bubbles | 'BETA'|
|lollipop_max_size | float: Maximum size of bubbles (Note this need to be a decimal, not an integer) | 5.|
|lollipop_stem_width | float: Width of the lollipop stems | 0.1|

Returns a lolliplot as a plotly figure.

> Note: the exon_info should be filtered to your desired transcript before inputting into ```reduce_gaps()``` or ```lolliplot_raw()```.


### lp.example():
Arguments: **None**

This is an example lolliplot script, using available associations on the Exeter server. This function can be modified and run to get lolliplots out quickly and easily. To save the figures, the ```print()``` parts should be commented out, and the ```htmlwidgets::saveWidget()```should be uncommented. The name of the output file should be changed from *lolliplots_example.html*.