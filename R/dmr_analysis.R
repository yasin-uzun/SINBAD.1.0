require(openxlsx)

dm_stat_test <- function(met_mat, cells_of_interest, bg_cells = NULL, test_type = 't', minimum_cell_count = 20)
{

  if(is.null(bg_cells))
  {
    all_cells = colnames(met_mat)
    bg_cells = setdiff(all_cells, cells_of_interest)
  }

  non_na_cell_counts_for_cells_of_interest = rowSums(!is.na(met_mat[, cells_of_interest]))
  head(non_na_cell_counts_for_cells_of_interest)
  non_na_cell_counts_for_bg_cells = rowSums(!is.na(met_mat[, bg_cells]))
  head(non_na_cell_counts_for_bg_cells)

  var.g = apply(met_mat, MARGIN = 1, var, na.rm = T)
  head(var.g)
  length(var.g)
  var.g[is.na(var.g)] = 0
  sum(var.g > 0)

  filtering1 = non_na_cell_counts_for_cells_of_interest > minimum_cell_count
  filtering2 = non_na_cell_counts_for_bg_cells > minimum_cell_count
  filtering3 = var.g > 0

  #filtering3['TBR1']

  filtering = filtering1 & filtering2 & filtering3

  head(filtering)
  sum(filtering)

  met_mat_filtered = met_mat[filtering, ]
  dim(met_mat_filtered)

  if(test_type == 't'| test_type == 'student')
  {
      p_values <- sapply(1:nrow(met_mat_filtered),
                     FUN = function(i)
                       t.test(met_mat_filtered[i, cells_of_interest],
                              met_mat_filtered[i, bg_cells])$p.value)
  }
  else if(test_type == 'wilcox')
  {
    p_values <- sapply(1:nrow(met_mat_filtered),
                       FUN = function(i)
                         wilcox.test(met_mat_filtered[i, cells_of_interest],
                                met_mat_filtered[i, bg_cells])$p.value)
  }else
  {
    msg = paste0('Test type provided (', test_type,  ') is not valid. Please specify t or wilcox.')
    stop(msg)
  }

  length(p_values)
  adjusted_p_values = p.adjust(p_values, method = 'BH')
  head(adjusted_p_values )
  sum(adjusted_p_values < 0.05, na.rm = T)


  pseudo_number = 0.05
  fold_changes <- sapply(1:nrow(met_mat_filtered),
                         FUN = function(i)
                           mean(met_mat_filtered[i, cells_of_interest] + pseudo_number, na.rm =  T)
                         / mean(met_mat_filtered[i, bg_cells] + pseudo_number, na.rm =  T) )

  mean_interest = rowMeans(met_mat_filtered[, cells_of_interest], na.rm = T)
  mean_bg = rowMeans(met_mat_filtered[, bg_cells], na.rm = T)
  log2_FC = round(log2(fold_changes), 2)

  df_result = data.frame(region = rownames(met_mat_filtered),
                         mean_bg = mean_bg,
                         mean_group = mean_interest,
                         log2_FC = log2_FC,
                         p.value = p_values,
                        adjusted.p.value = adjusted_p_values
                         )

  rownames(df_result) = rownames(met_mat_filtered)

  head(df_result)

  df_result_sorted = df_result[order(df_result$p.value), ]

  return(df_result_sorted)

}


dm_stat_test_for_clusters <- function(cluster_vec, feature_mat, minimum_cell_count = 20, dmr_adj_p_value_cutoff = 0.01, dmr_log2_fc_cutoff = 0.25)
{

  all_cells = colnames(feature_mat)
  cluster_ids = sort(unique(cluster_vec))
  dm_result_list = list()
  sig_region_list = list()
  sig_count_list = list()

  for(cluster_id in cluster_ids)
  {
    print(cluster_id)
    cluster_cells = all_cells[cluster_vec == cluster_id]

    dm_test = dm_stat_test(met_mat = feature_mat,
                           cells_of_interest = cluster_cells,
                           minimum_cell_count = minimum_cell_count)

    dm_result = data.frame(cluster = cluster_id, dm_test)
    sig_ids = as.character(dm_result$region[dm_result$adjusted.p.value < dmr_adj_p_value_cutoff &
                                              abs(dm_result$log2_FC) > dmr_log2_fc_cutoff])
    sig_region_list[[cluster_id]] = sig_ids

    dm_result_up = dm_result[dm_result$log2_FC > 0, ]
    dm_result_down = dm_result[dm_result$log2_FC < 0, ]

    dm_result_list[[paste0('C_', cluster_id, '.Incr.')]] = dm_result_up
    dm_result_list[[paste0('C_', cluster_id, '.Dec.')]] = dm_result_down

    sig_up_count = sum(dm_result_up$adjusted.p.value < dmr_adj_p_value_cutoff)
    sig_down_count = sum(dm_result_down$adjusted.p.value < dmr_adj_p_value_cutoff)

    sig_count_list[[paste0('C_', cluster_id)]] = c(sig_up_count, sig_down_count)
  }

  sig_ids =  as.character(unlist(sig_region_list))

  DM_summary = do.call(rbind.data.frame, sig_count_list)
  colnames(DM_summary) = c('Incr.', 'Decr.')
  DM_summary = data.frame(Cluster = cluster_ids, DM_summary)


  dm_result_list_with_summary = list()
  dm_result_list_with_summary[[1]] = DM_summary
  dm_result_list_with_summary[2:(length(dm_result_list) + 1)] = dm_result_list
  names(dm_result_list_with_summary) = c('Summary',  names(dm_result_list))
  names(dm_result_list_with_summary)
  length(dm_result_list_with_summary)

  analyis_result = list('sig_ids' = sig_ids, 'dm_result_list_with_summary' = dm_result_list_with_summary)

  return(analyis_result)
}



