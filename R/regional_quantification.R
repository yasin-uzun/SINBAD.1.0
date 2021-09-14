library(dplyr)
library(data.table)
library(GenomicRanges)

read_region_annot <- function(region_annot_file, format_file, type_filtering =  'none')
{
    df_region_annot = read.table(region_annot_file, sep ='\t')
    df_format = read.table(format_file, sep ='\t')
    colnames(df_region_annot) = as.character(df_format$V2)
    head(df_region_annot)
    dim(df_region_annot)
    #df_region_annot %>% distinct(region_name, .keep_all = TRUE)

    df_region_annot <- dplyr::distinct(df_region_annot, region_name, .keep_all = TRUE)
    dim(df_region_annot)
    head(df_region_annot)
    length(unique(df_region_annot$region_name))


    if(type_filtering != 'none')
    {
      df_region_annot  = df_region_annot[df_region_annot$region_type == type_filtering, ]
    }

    rownames(df_region_annot) =  df_region_annot$region_name

    head(df_region_annot)

    return(df_region_annot)
}

get_promoters <- function(df_gene_annot)
{

    head(df_gene_annot)
    df_promoters = df_gene_annot
    head(df_promoters)
    attach(df_gene_annot)
    df_promoters$start[strand == '+'] = df_gene_annot$start[strand == '+'] - 2000
    df_promoters$end[strand == '+'] = df_gene_annot$start[strand == '+'] + 500

    df_promoters$end[strand == '-'] = df_gene_annot$end[strand == '-'] + 2000
    df_promoters$start[strand == '-'] = df_gene_annot$end[strand == '-'] - 500
    detach(df_gene_annot)
    head(df_promoters)

    df_promoters$start[df_promoters$start < 0] = 0

    return(df_promoters)
}

convert_to_granges <- function(df_region)
{

    gr_region <- with(df_region, GRanges(chrom,
                                           IRanges(start+1, end),
                                           strand = '*',
                                           region_name,
                                           region_type))

    return(gr_region)
}


compute_region_met_matrix <- function(  list_call_count_matrices,
                                       min_call_count_threshold = 10)
{

  full_met_call_count_matrix = list_call_count_matrices$full_met_call_count_matrix
  full_total_call_count_matrix = list_call_count_matrices$full_total_call_count_matrix

  met_rate_matrix = full_met_call_count_matrix / full_total_call_count_matrix

  met_rate_matrix[full_total_call_count_matrix < min_call_count_threshold] = NA

  return(met_rate_matrix)

}


compute_aggr_met_rate <- function( list_call_count_matrices )
{

  aggr_met_counts = rowSums(list_call_count_matrices$full_met_call_count_matrix, na.rm = T)
  aggr_total_counts = rowSums(list_call_count_matrices$full_total_call_count_matrix, na.rm = T)

  head(aggr_met_counts)
  head(aggr_total_counts)

  aggr_rate = aggr_met_counts / aggr_total_counts

  head(aggr_rate)

  return(aggr_rate)

}



compute_call_count_matrices <- function(  df_region,
                               methylation_calls_dir,
                               methylation_type = 'CpG',
                               exclude_cells = c()
                              )
{


  gr_region = convert_to_granges(df_region)
  setwd(methylation_calls_dir)
  pattern  = paste0(methylation_type, '.*organism.cov.gz')
  cov_files = list.files(methylation_calls_dir, pattern)



  result_list = list()

  cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
  registerDoSNOW(cl)

  #for(i in  1:length(cov_files))
  #result_list <- foreach(i=1:length(cov_files)) %dopar%
  met_hits_list <- foreach(i=1:length(cov_files)) %dopar%
  {
    library(data.table)
    library(GenomicRanges)

    print('************************')
    print(i)
    cov_file = cov_files[i]

    print(cov_file)
    print('************************')
    cell_id = gsub('.organism.cov.gz', '', cov_file)
    cell_id = gsub(paste0(methylation_type,'_calls.'), '', cell_id)


    #dt_cov = fread(paste0(methylation_calls_dir, cov_file) )
    dt_cov.dummy = data.table(data.frame(chrom = 'chrZz',
                              start = 1,
                              end = 2,
                              met_rate = -.5,
                              met = 1,
                              demet = 1   ))


    dt_cov = tryCatch({
      fread(paste0(methylation_calls_dir, cov_file) )
    }, warning = function(w) {
      dt_cov.dummy
    }, error = function(e) {
      dt_cov.dummy
    }, finally = {

    }
    )

    colnames(dt_cov) = c('chrom', 'start', 'end', 'met_rate', 'met', 'demet')
    #dt_cov$chr = paste0('chr', dt_cov$chr)
    head(dt_cov)
    unique(dt_cov$chr)

    gr_cov <- with(dt_cov, GRanges(chrom, IRanges(start+1, end), strand = '*', met_rate, met, demet)  )
    gr_cov

    #dt_inter = data.table(intersect_bed(gr_region, gr_cov))
    #dim(dt_inter)
    #head(dt_inter)

    hits_obj <- findOverlaps(gr_region, gr_cov)
    class(hits_obj)
    da = as.data.frame(gr_region[queryHits(hits_obj)])
    db = as.data.frame(gr_cov[subjectHits(hits_obj)])
    dt_inter <- data.table(cbind(da, db))


    quant_cols = c('met', 'demet')

    dt_aggr <- dt_inter[, lapply(.SD, sum), by = .(region_name), .SDcols = quant_cols  ]
    length(unique(dt_aggr$region_name))
    head(dt_aggr)
    dim(dt_aggr)

    dt_aggr$call_count = dt_aggr$met + dt_aggr$demet



    #result_list[[cell_id]] = met_rate_vector

    df_aggr = data.frame(dt_aggr)

    head(df_aggr)

    rownames(df_aggr) = df_aggr$region_name

    df_aggr_x = base::merge(df_region, df_aggr, by.x = 'region_name', by.y = 'region_name', all.x = T)

    rownames(df_aggr_x) = df_aggr_x$region_name

    df_aggr_x
  }

  stopCluster(cl)  #not reached

  temp = met_hits_list
  cell_ids = cov_files
  cell_ids = gsub('.organism.cov.gz', '', cell_ids)
  cell_ids = gsub(paste0(methylation_type,'_calls.'), '', cell_ids)
  names(met_hits_list) = cell_ids


  cell_ids = names(met_hits_list)

  temp = met_hits_list
  for(cell_id in cell_ids)
  {
    met_hits_list[[cell_id]] = temp[[cell_id]][df_region$region_name, ]
  }

  include_cells = setdiff(cell_ids, exclude_cells)

  met_hits_list = met_hits_list[include_cells]


  #met_hits_list = lapply(met_call_count_list, '[', df_region$region_name)


  met_call_count_list = lapply(met_hits_list, '[[', 'met')
  demet_call_count_list = lapply(met_hits_list, '[[', 'demet')
  total_call_count_list = lapply(met_hits_list, '[[', 'call_count')

  head(lapply(met_call_count_list, length))
  head(met_call_count_list[[1]])
  length(met_call_count_list[[1]])

  class(met_call_count_list[[1]])

  full_met_call_count_list = met_call_count_list
  full_total_call_count_list = total_call_count_list


  full_met_call_count_matrix = do.call('cbind', full_met_call_count_list)
  full_total_call_count_matrix = do.call('cbind', full_total_call_count_list)

  rownames(full_met_call_count_matrix) = df_region$region_name
  rownames(full_total_call_count_matrix) = df_region$region_name

  full_met_call_count_matrix[1:5, 1:5]
  full_total_call_count_matrix[1:5, 1:5]

  list_call_count_matrices = list('full_met_call_count_matrix' = full_met_call_count_matrix,
                     'full_total_call_count_matrix' = full_total_call_count_matrix
                     )


  return(list_call_count_matrices)


}



impute_nas <- function(met_mat, max_ratio_of_na_cells = 0.25)
{
  count_na_by_region = rowSums(is.na(met_mat))
  head(count_na_by_region)
  ratio_na_by_region = count_na_by_region / ncol(met_mat)
  head(ratio_na_by_region)
  sum(ratio_na_by_region > max_ratio_of_na_cells)
  use_regions = ratio_na_by_region <= max_ratio_of_na_cells
  met_mat_filtered = met_mat[use_regions, ]
  print(dim(met_mat_filtered))
  #met_mat_filtered[1:5, 1:5]

  met_mat_imputed = met_mat_filtered
  na_index <- which(is.na(met_mat_imputed), arr.ind=TRUE)
  met_mat_imputed[na_index] <- rowMeans(met_mat_imputed, na.rm=TRUE)[na_index[,1]]

  return(met_mat_imputed)
}

