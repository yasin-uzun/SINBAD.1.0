library(doSNOW)

#source('/mnt/isilon/tan_lab/uzuny/scripts/R_utilities/biology.v07.R')
source('/mnt/isilon/tan_lab/uzuny/scripts/R_utilities/seurat_plotting.v03.R')

test <- function()
{
  cat('SINBAD installation is ok.\n')
}
#

read_configs <- function(config_dir)
{
  print('Reading config.general.R for program paths...')
  source(paste0(config_dir, '/config.general.R') )
  print('Reading config.genome.R for genome paths...')
  source(paste0(config_dir, '/config.genome.R') )
  print('Reading config.project.R for project settings...')
  source(paste0(config_dir, '/config.project.R') )
  #system('echo $PATH')

}



#sample_name = 'Test'
#readRenviron(path = 'variables.env')
#system('source ~/.bashrc')

#image_file = paste0(working_dir, '/Image_2020_06_03.img')
#image_file = paste0(working_dir, '/Image_2020_06_11.img')
#save.image(image_file)

construct_sinbad_object <- function(raw_fastq_dir  = NA,
                                    demux_index_file = NA,
                                    working_dir ,
                                    sample_name)
{

  sample_working_dir = paste0(working_dir, '/', sample_name, '/')
  dir.create(sample_working_dir, recursive = T)

  main_log_dir <<- paste0(sample_working_dir, '/logs/')
  demux_fastq_dir <<- paste0(sample_working_dir, '/demux_fastq/')
  trimmed_fastq_dir <<- paste0(sample_working_dir, '/trimmed_fastq/')
  merged_alignment_dir <<- paste0(sample_working_dir, '/alignments_merged/')



  methylation_calls_dir <<- paste0(sample_working_dir, '/methylation_calls/')
  summary_dir <<- paste0(sample_working_dir, '/summary/')
  plot_dir <<- paste0(sample_working_dir, '/plots/')
  matrix_dir <<- paste0(sample_working_dir, '/matrices/')
  object_dir <<- paste0(sample_working_dir, '/objects/')
  alignment_summary_file <<- paste0(sample_working_dir, '/Alignment_Summary.tsv')

  if(trimmer == "NoTrimming")
  {
    trimmed_fastq_dir = demux_fastq_dir
  }


  dir.create(main_log_dir, showWarnings = F, recursive = T)
  dir.create(trimmed_fastq_dir, showWarnings = F, recursive = T)
  dir.create(merged_alignment_dir, showWarnings = F, recursive = T)
  dir.create(methylation_calls_dir, showWarnings = F, recursive = T)
  dir.create(summary_dir, showWarnings = F, recursive = T)
  dir.create(plot_dir, showWarnings = F, recursive = T)
  dir.create(demux_fastq_dir, showWarnings = F, recursive = T)
  dir.create(matrix_dir, showWarnings = F, recursive = T)
  dir.create(object_dir, showWarnings = F, recursive = T)

  sinbad_object = list('sample_working_dir' = sample_working_dir,
                       'main_log_dir' = main_log_dir,
                       'raw_fastq_dir' = raw_fastq_dir,
                       'demux_index_file' = demux_index_file,
                       'demux_fastq_dir' = demux_fastq_dir,
                       'trimmed_fastq_dir' =  trimmed_fastq_dir,
                       'merged_alignment_dir' = merged_alignment_dir,
                       'methylation_calls_dir'  = methylation_calls_dir,
                       'summary_dir' = summary_dir,
                       'plot_dir' = plot_dir,
                       'sample_name' = sample_name,
                       'matrix_dir' = matrix_dir,
                       'object_dir' = object_dir,
                       'annot_list' = list(),
                       'met_matrices' = list(),
                       'list_of_list_call_count_matrices' = list(),
                       'list_aggr_rates' = list()
                       )


  if(sequencing_type == 'paired')
  {
    r2_meta_fastq_dir <<- paste0(sample_working_dir, '/r2_meta_fastq/')
    r1_and_r2_alignments_dir <<- paste0(sample_working_dir, '/alignments_r1_and_r2/')

    dir.create(r2_meta_fastq_dir, showWarnings = F, recursive = T)
    dir.create(r1_and_r2_alignments_dir, showWarnings = F, recursive = T)

    sinbad_object$r2_meta_fastq_dir = r2_meta_fastq_dir
    sinbad_object$r1_and_r2_alignments_dir = r1_and_r2_alignments_dir

  }


  class(sinbad_object) = 'Sinbad'

  return(sinbad_object)


}


wrap_demux_fastq_files <- function(sinbad_object, flag_r2_index_embedded_in_r1_reads = FALSE)
{
  #Demultiplex fastq files
  #TODO: Check the input fastq dir and demux file exists, give error and stop otherwise


  if(sequencing_type == 'paired' )
  { #Pair ended
    demux_fastq_files(raw_fastq_dir = sinbad_object$raw_fastq_dir,
                      demux_index_file = sinbad_object$demux_index_file,
                      demux_index_length = demux_index_length,
                      demux_fastq_dir = sinbad_object$demux_fastq_dir,
                      main_log_dir = sinbad_object$main_log_dir,
                      sample_name = sinbad_object$sample_name,
                      read_type = 'R1')

    if(flag_r2_index_embedded_in_r1_reads)
    {

      demux_fastq_files(sinbad_object$r2_meta_fastq_dir,
                      sinbad_object$demux_index_file,
                      demux_index_length,
                      sinbad_object$demux_fastq_dir,
                      sinbad_object$main_log_dir,
                      sinbad_object$sample_name,
                      read_type = 'R2')
    }else{
      demux_fastq_files(sinbad_object$raw_fastq_dir,
                        sinbad_object$demux_index_file,
                        demux_index_length,
                        sinbad_object$demux_fastq_dir,
                        sinbad_object$main_log_dir,
                        sinbad_object$sample_name,
                        read_type = 'R2')
    }# if(flag_r2_index_embedded_in_r1_reads)
  }else
  { #single ended
    demux_fastq_files(sinbad_object$raw_fastq_dir,
                      sinbad_object$demux_index_file,
                      demux_index_length,
                      sinbad_object$demux_fastq_dir,
                      sinbad_object$main_log_dir,
                      sinbad_object$sample_name)
  }

  return(sinbad_object)

}

wrap_demux_stats <- function(sinbad_object)
{

  sinbad_object$df_demux_reports = read_demux_logs(sinbad_object$main_log_dir)

  demux_summary_file = paste0(summary_dir, '/Demux_statistics.tsv')
  write.table(sinbad_object$df_demux_reports,
              file = demux_summary_file,
              sep = '\t', quote = F,
              row.names = F, col.names = T)

  #Count demuxd reads
  sinbad_object$demux_read_counts =  count_fastq_reads(sinbad_object$demux_fastq_dir)

  return(sinbad_object)

}

wrap_trim_fastq_files <- function(sinbad_object)#Trim adapters
{
  trim_fastq_files(demux_fastq_dir = sinbad_object$demux_fastq_dir,
                   trimmed_fastq_dir = sinbad_object$trimmed_fastq_dir,
                   main_log_dir = sinbad_object$main_log_dir)

  return(sinbad_object)
}

wrap_trim_stats <- function(sinbad_object)#Trim adapters
{
  sinbad_object$trimmed_read_counts = count_fastq_reads(sinbad_object$trimmed_fastq_dir)



  return(sinbad_object)
}

wrap_plot_preprocessing_stats <- function(sinbad_object)
{
  dir.create(paste0(sinbad_object$plot_dir, '/QC/'), recursive = T)

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Preprocessing_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = T, title = sample_name)
  #postscript(plot_file, paper = 'special', horizontal = F, width = 8, height = 8, title = sample_name)

  stat1 = plot_preprocessing_results(sample_name = sinbad_object$sample_name,
                             demux_reports = sinbad_object$df_demux_reports,
                             demux_read_counts = sinbad_object$demux_read_counts,
                             trimmed_read_counts = sinbad_object$trimmed_read_counts)
  dev.off()

  df.stat = data.frame(names(stat1), stat1)
  tsv_file = paste0(sinbad_object$summary_dir, '/Overall_preprocessing_statistics.tsv')
  write.table(df.stat, file = tsv_file, quote = F, row.names = F, col.names = F, sep = '\t')

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Preprocessing_statistics.png')
  png(plot_file, width = 800, height = 600)
  plot_preprocessing_results(sample_name = sinbad_object$sample_name,
                             demux_reports = sinbad_object$df_demux_reports,
                             demux_read_counts = sinbad_object$demux_read_counts,
                             trimmed_read_counts = sinbad_object$trimmed_read_counts)
  dev.off()

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Preprocessing_statistics.pdf')
  pdf(plot_file, width = 10, height = 7)
  plot_preprocessing_results(sample_name = sinbad_object$sample_name,
                             demux_reports = sinbad_object$df_demux_reports,
                             demux_read_counts = sinbad_object$demux_read_counts,
                             trimmed_read_counts = sinbad_object$trimmed_read_counts)
  dev.off()

}

wrap_align_sample <- function(sinbad_object, pattern = '')
{
  pbat_flag = 0
  #if(sequencing_type == 'paired' | protocol == 'snmc'){

  #}
  alignment_dir = sinbad_object$r1_and_r2_alignments_dir


  #Run aligner
  #if(sequencing_type == 'paired')
  #{
    if(grepl('--pbat', bismark_aligner_param_settings) ) {pbat_flag = 1}
  #}

  #pattern = '*_001_ACTTGA.bam'

  if(pbat_flag == 0) #If no pbat, align as usual
  {
    align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
               genomic_sequence_path = genomic_sequence_path,
               alignment_dir = alignment_dir,
               aligner = aligner,
               num_cores= num_cores,
               mapq_threshold =mapq_threshold,
               main_log_dir = sinbad_object$main_log_dir,
               aligner_param_settings = bismark_aligner_param_settings,
               pattern = pattern)
  }else
  {
    print('Bismark-paired-pbat')
    bismark_aligner_param_settings.r1 = bismark_aligner_param_settings
    bismark_aligner_param_settings.r2 = gsub('--pbat', '', bismark_aligner_param_settings)

    #Left reads
    aln_result = align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
                 genomic_sequence_path = genomic_sequence_path,
                 alignment_dir = alignment_dir,
                 aligner = aligner,
                 num_cores= num_cores,
                 mapq_threshold =mapq_threshold,
                 main_log_dir = sinbad_object$main_log_dir,
                 aligner_param_settings = bismark_aligner_param_settings.r1,
                 pattern = 'R1')

    count_attempts = 0
    while(aln_result > 0)
    {
      count_attempts = count_attempts + 1
      if(count_attempts > 10)
      {
        msg1 = '!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!\n'
        msg2 = 'Alignment failed for some cells \n'
        msg3 = paste(aln_result, 'cells yet to be aligned.')
        msg4 = 'Please consider increasing memory.\n'
        msg5 = '!!!!!!!!!!!!!!!!!Exiting!!!!!!!!!!!!!!!!!!!!!!\n'
        msg = paste0(msg1, msg2, msg3, msg4, msg5)
        stop(paste(msg1, msg2, msg3))

      }

      aln_result = align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
                   genomic_sequence_path = genomic_sequence_path,
                   alignment_dir = alignment_dir,
                   aligner = aligner,
                   num_cores= num_cores,
                   mapq_threshold =mapq_threshold,
                   main_log_dir = sinbad_object$main_log_dir,
                   aligner_param_settings = bismark_aligner_param_settings.r1,
                   pattern = 'R1')
    }



    #Right reads

    aln_result = align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
                 genomic_sequence_path = genomic_sequence_path,
                 alignment_dir = alignment_dir,
                 aligner = aligner,
                 num_cores= num_cores,
                 mapq_threshold =mapq_threshold,
                 main_log_dir = sinbad_object$main_log_dir,
                 aligner_param_settings = bismark_aligner_param_settings.r2,
                 pattern = 'R2')


    count_attempts = 0
    while(aln_result > 0)
    {
      count_attempts = count_attempts + 1
      if(count_attempts > 10)
      {
        msg1 = '!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!\n'
        msg2 = 'Alignment failed for some cells \n'
        msg3 = paste(aln_result, 'cells yet to be aligned.')
        msg4 = 'Please consider increasing memory.\n'
        msg5 = '!!!!!!!!!!!!!!!!!Exiting!!!!!!!!!!!!!!!!!!!!!!\n'
        msg = paste0(msg1, msg2, msg3, msg4, msg5)
        stop(paste(msg1, msg2, msg3))

      }

      aln_result = align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
                   genomic_sequence_path = genomic_sequence_path,
                   alignment_dir = alignment_dir,
                   aligner = aligner,
                   num_cores= num_cores,
                   mapq_threshold =mapq_threshold,
                   main_log_dir = sinbad_object$main_log_dir,
                   aligner_param_settings = bismark_aligner_param_settings.r2,
                   pattern = 'R2')
    }





  }


}



wrap_generate_alignment_stats <- function(sinbad_object)
{


  sinbad_object$df_alignment_reports = process_bismark_alignment_reports(alignment_dir = sinbad_object$r1_and_r2_alignments_dir)
  sinbad_object$df_bam_read_counts = count_bam_files(alignment_dir = sinbad_object$r1_and_r2_alignments_dir)
  dim(sinbad_object$df_alignment_reports)
  dim(sinbad_object$df_bam_read_counts)


  sinbad_object$df_alignment_stats_r1_and_r2 = base::merge(sinbad_object$df_alignment_reports,
                                                 sinbad_object$df_bam_read_counts, by = 0)
  dim(sinbad_object$df_alignment_stats_r1_and_r2)
  head(sinbad_object$df_alignment_stats_r1_and_r2)

  sinbad_object$df_alignment_stats = merge_r1_and_r2_alignment_stats(sinbad_object$df_alignment_stats_r1_and_r2)

  return(sinbad_object)

}

#temp1 = sinbad_object$df_alignment_stats
#sinbad_object$df_alignment_stats_r1_and_r2 = sinbad_object$df_alignment_stats

wrap_compute_coverage_rates <- function(sinbad_object, parallel = T)
{
  log_file = paste0(sinbad_object$main_log_dir, '/alignment/coverage.log')
  df_coverage_rates = compute_coverage_rates(alignment_dir = sinbad_object$merged_alignment_dir,
                                             log_file = log_file, parallel = parallel)
  sinbad_object$df_coverage_rates =df_coverage_rates

  df_coverage_rates_ordered = sinbad_object$df_coverage_rates[as.character(sinbad_object$df_alignment_stats$Cell_ID), ]

  sinbad_object$df_alignment_stats$base_count = df_coverage_rates_ordered[,1]
  sinbad_object$df_alignment_stats$coverage_rate = df_coverage_rates_ordered[,2]



  return(sinbad_object)
}



wrap_plot_alignment_stats <- function(sinbad_object)
{
  ###Scatter plot filtering
  filter_aln_rate =  sinbad_object$df_alignment_stats$Alignment_rate  > alignment_rate_threshold
  #filter_read_count = df_alignment_stats$organism_read_counts > organism_minimum_filtered_read_count
  filter_read_count = sinbad_object$df_alignment_stats$nc_filtered_read_counts > minimum_filtered_read_count

  passing = filter_aln_rate &  filter_read_count

  sinbad_object$df_alignment_stats$pass = 0
  sinbad_object$df_alignment_stats$pass[passing] = 1

  sinbad_object$df_alignment_stats$Row.names = NULL
  rownames(sinbad_object$df_alignment_stats) = sinbad_object$df_alignment_stats$Cell_ID
  sinbad_object$alignment_summary_file = paste0(sinbad_object$summary_dir, '/Alignment_statistics.tsv')
  write.table(sinbad_object$df_alignment_stats,
              file = sinbad_object$alignment_summary_file,
              sep = '\t', quote = F, row.names = F, col.names = T)

  dir.create(paste0(sinbad_object$plot_dir, '/QC/') )

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Alignment_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = T, title = sinbad_object$sample_name)
  stat1 = plot_alignment_stats(sample_name = sinbad_object$sample_name,
                       df_alignment_stats = sinbad_object$df_alignment_stats)
  dev.off()

  df.stat = data.frame(names(stat1), stat1)
  print(df.stat)
  tsv_file = paste0(sinbad_object$summary_dir, '/Overall_alignment_statistics.tsv')
  write.table(df.stat, file = tsv_file, quote = F, row.names = F, col.names = F, sep = '\t')


  plot_file = paste0(sinbad_object$plot_dir, '/QC/Alignment_statistics.png')
  png(plot_file, width = 800, height = 600)
  plot_alignment_stats(sinbad_object$sample_name, sinbad_object$df_alignment_stats)
  dev.off()

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Alignment_statistics.pdf')
  pdf(plot_file, width = 10, height = 7)
  plot_alignment_stats(sinbad_object$sample_name, sinbad_object$df_alignment_stats)
  dev.off()

  sinbad_object$df_overall_alignment_stat = df.stat
  return(sinbad_object)

}



wrap_merge_r1_and_r2_bam <- function(sinbad_object)
{
  print('Merging lambda')
  merge_r1_and_r2_bam_for_sample(alignment_dir_in = sinbad_object$r1_and_r2_alignments_dir,
                                 alignment_dir_out = sinbad_object$merged_alignment_dir,
                                 bam_type = 'lambda')

  print('Merging organism')

  merge_r1_and_r2_bam_for_sample(alignment_dir_in = sinbad_object$r1_and_r2_alignments_dir,
                                 alignment_dir_out = sinbad_object$merged_alignment_dir,
                                  bam_type = 'organism')
}

wrap_call_methylation_sites <- function(sinbad_object)
{

  call_methylation_sites_for_sample(alignment_dir = sinbad_object$merged_alignment_dir,
                                    methylation_calls_dir = sinbad_object$methylation_calls_dir,
                                    log_dir = sinbad_object$main_log_dir,
                                    bme_param_settings = bme_param_settings )

  return(sinbad_object)


}

wrap_generate_methylation_stats <- function(sinbad_object)
{

  sinbad_object$df_org_split_reports = process_bismark_split_reports(methylation_calls_dir =  sinbad_object$methylation_calls_dir,
                                                                     genome_type = 'organism')


  sinbad_object$list_org_bias_reports = process_bismark_bias_reports(methylation_calls_dir =  sinbad_object$methylation_calls_dir,
                                                                     genome_type = 'organism')

  sinbad_object$df_lambda_split_reports = NULL
  sinbad_object$list_lambda_bias_reports = NULL

  if(exists('lambda_control') )
  {
    if(lambda_control)
    {
      sinbad_object$df_lambda_split_reports = process_bismark_split_reports(methylation_calls_dir =  sinbad_object$methylation_calls_dir,
                                                                            genome_type = 'lambda')

      sinbad_object$list_lambda_bias_reports = process_bismark_bias_reports(methylation_calls_dir =  sinbad_object$methylation_calls_dir,
                                                                            genome_type = 'lambda')
    }
  }

  list_met_call_counts = list()
  met_types = c('CpG', 'CHG', 'CHH')

  for(met_type in met_types)
  {
    list_met_call_counts[[met_type]] = get_met_call_counts(methylation_calls_dir = sinbad_object$methylation_calls_dir, met_type)
  }

  sinbad_object$list_met_call_counts = list_met_call_counts

  return(sinbad_object)
}

wrap_plot_met_stats <- function(sinbad_object)
{
  dir.create(paste0(sinbad_object$plot_dir, '/QC/') )

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Met_call_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = T, title = sinbad_object$sample_name)
  plot_split_reports(df_org_split_reports = sinbad_object$df_org_split_reports,
                     df_lambda_split_reports = sinbad_object$df_lambda_split_reports,
                     list_org_bias_reports = sinbad_object$list_org_bias_reports,
                     list_met_call_counts = sinbad_object$list_met_call_counts,
                     lambda_flag = lambda_control)
  dev.off()

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Met_call_statistics.png')
  png(plot_file, width = 800, height = 600)
  plot_split_reports(sinbad_object$df_org_split_reports,
                     sinbad_object$df_lambda_split_reports,
                     sinbad_object$list_org_bias_reports,
                     list_met_call_counts = sinbad_object$list_met_call_counts,
                     lambda_flag = lambda_control)
  dev.off()

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Met_call_statistics.pdf')
  pdf(plot_file, width = 10, height = 7)
  plot_split_reports(sinbad_object$df_org_split_reports,
                     sinbad_object$df_lambda_split_reports,
                     sinbad_object$list_org_bias_reports,
                     list_met_call_counts = sinbad_object$list_met_call_counts,
                     lambda_flag = lambda_control)
  dev.off()

  return(sinbad_object)

}

wrap_read_annot <- function(sinbad_object, annot_file, annot_format_file, annot_name = NULL)
{

  if(! 'annot_list' %in% names(sinbad_object)  )
  {
    sinbad_object$annot_list = list()
  }

  #Read regions
  print('Reading genomic coordinates')
  df_annot = read_region_annot(annot_file, annot_format_file)
  head(df_annot)
  df_annot$start[df_annot$start < 0] = 0
  num_regions = nrow(df_annot)
  print(paste('Found ',num_regions,' regions'))
  sizes = df_annot$end - df_annot$start
  mean_size = round(mean(sizes))
  print(paste('Mean region size is ',mean_size,' bp'))

  if(is.null(annot_name))
  {
    #No annot name is given, use file name
    annot_name = basename(annot_file)
    annot_name = sub('regions.', '', annot_name)
    annot_name = sub('.bed', '', annot_name)

  }

  sinbad_object$annot_list[[annot_name]] = df_annot

  print(paste('Saved the genomic coordinates to of the sinbad_object$annot_list with the entry name: ',annot_name))


  return(sinbad_object)

}


wrap_quantify_regions <- function(sinbad_object, annot_name = 'Bins_100Kb', methylation_type = 'CpG', min_call_count_threshold_for_region = 10)
{


  if(  !('annot_list' %in% names(sinbad_object) )  | length(sinbad_object$annot_list) == 0 )
  {
    print('Annotation regions are not computed. Call wrap_read_annots() function first. Exiting.')
    return(sinbad_object)
  }

  if(  !(annot_name %in% names(sinbad_object$annot_list) )   )
  {
    print(paste0('Annotation regions are not computed for: ', annot_name) )
    return(sinbad_object)
  }


  if(! 'met_matrices' %in% names(sinbad_object)  )
  {
    sinbad_object$met_matrices = list()
  }

  exclude_cells = c('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')

  if( 'df_alignment_stats' %in% names(sinbad_object)  )
  {
    attach(sinbad_object$df_alignment_stats)
    exclude_cells = rownames(sinbad_object$df_alignment_stats[pass == 0, ])
    detach(sinbad_object$df_alignment_stats)
  }


  if(! 'list_of_list_call_count_matrices' %in% names(sinbad_object)  )
  {
    sinbad_object$list_of_list_call_count_matrices = list()
  }

  if(! 'list_aggr_rates' %in% names(sinbad_object)  )
  {
    sinbad_object$list_aggr_rates = list()
  }




  print(annot_name)

  list_call_count_matrices = compute_call_count_matrices(
              df_region = sinbad_object$annot_list[[annot_name]],
              methylation_calls_dir = sinbad_object$methylation_calls_dir,
              methylation_type = methylation_type,
              exclude_cells = exclude_cells)

  names(list_call_count_matrices)
  list_call_count_matrices$full_met_call_count_matrix[1:5,1:5]
  sum(!is.na(list_call_count_matrices$full_met_call_count_matrix))


  aggr_rate = compute_aggr_met_rate(list_call_count_matrices)

  sinbad_object$list_of_list_call_count_matrices[[methylation_type]][[annot_name]] = list_call_count_matrices
  sinbad_object$list_aggr_rates[[methylation_type]][[annot_name]] = aggr_rate

  met_mat = compute_region_met_matrix(list_call_count_matrices = list_call_count_matrices,
                                     min_call_count_threshold = min_call_count_threshold_for_region)
  sinbad_object$met_matrices[[annot_name]] = met_mat


  met_rate_file = paste0(sinbad_object$matrix_dir, 'met_rate_mat.' , annot_name, '.',methylation_type ,'.rds')
  saveRDS(met_mat, file = met_rate_file)

  aggr_rate_file = paste0(sinbad_object$matrix_dir, 'aggr_rate.' , annot_name, '.',methylation_type ,'.rds')
  saveRDS(aggr_rate, file = aggr_rate_file)

  met_count_mat_file = paste0(sinbad_object$matrix_dir, 'met_count_mat.' , annot_name, '.',methylation_type ,'.rds')
  saveRDS(list_call_count_matrices$full_met_call_count_matrix, file = met_count_mat_file)

  total_count_mat_file = paste0(sinbad_object$matrix_dir, 'total_count_mat.' , annot_name, '.',methylation_type ,'.rds')
  saveRDS(list_call_count_matrices$full_total_call_count_matrix, file = total_count_mat_file)




  return(sinbad_object)

}



wrap_impute_nas <- function(sinbad_object, annot_name = 'Bins_100kb',  max_ratio_of_na_cells = 0.25)
{

    met_mat = sinbad_object$met_matrices[[annot_name]]
    names(sinbad_object$met_matrices)
    #sinbad_object$met_matrices$Bins_10Kb[1:5, 1:5]
    imputed_matrix = impute_nas(met_mat = met_mat
                                , max_ratio_of_na_cells = max_ratio_of_na_cells)

    imputed_file = paste0(sinbad_object$matrix_dir, 'imputed.' , annot_name, '.rds')
    saveRDS(imputed_matrix, file = imputed_file)

    sinbad_object$imputed_matrices[[annot_name]] = imputed_matrix


  return(sinbad_object)

}



wrap_dim_red <- function(sinbad_object,
                         annot_type = 'Bins_100Kb',
                         legend_title = 'Methylation',
                         methylation_type = 'CpG',
                         rho_threshold = 1, delta_threshold = 3.5)
{
  #Reduce dimensionality
  #marker_genes = get_marker_genes('leuk')

  if( ! ('dim_red_objects' %in% names(sinbad_object)  )   )
  {
    sinbad_object$dim_red_objects = list()
  }

  name_for_dim_red = annot_type

  dim_red_object = reduce_dims_for_sample(
                         met_mat_for_dim_red = sinbad_object$met_matrices[[annot_type]]  ,
                         name_for_dim_red = name_for_dim_red  ,
                         plot_dir = sinbad_object$plot_dir  ,
                         methylation_type = methylation_type
                         )


  plot_dir = sinbad_object$plot_dir

  title = paste0(sample_name, ' - ' , methylation_type ,
                 '\nDR region: ', annot_type
                 #, '\nFeature region: ', name_for_features
                 )

  dir.create(paste0(plot_dir, '/DimRed/') )

  gg1 = plot_dim_red(dim_red = dim_red_object$umap, title = title )
  print(gg1)

  plot_file = paste0(plot_dir, '/DimRed/UMAP.',name_for_dim_red,'.eps')
  ggsave(gg1, filename = plot_file, device = 'eps', width = 16, height = 16, units = 'cm')

  plot_file = paste0(plot_dir, '/DimRed/UMAP.',name_for_dim_red,'.png')
  ggsave(gg1, filename = plot_file, device = 'png', width = 16, height = 16, units = 'cm')



  clusters = compute_clusters(umap = dim_red_object$umap, rho_threshold = rho_threshold, delta_threshold = delta_threshold)
  dim_red_object$clusters = clusters
  sinbad_object$dim_red_objects[[annot_type]] = dim_red_object


  gg1 = plot_dim_red(dim_red = dim_red_object$umap, title = title, cell_groups = clusters )
  print(gg1)

  plot_file = paste0(plot_dir, '/DimRed/UMAP.',name_for_dim_red,'.clusters.eps')
  ggsave(gg1, filename = plot_file, device = 'eps', width = 16, height = 16, units = 'cm')

  plot_file = paste0(plot_dir, '/DimRed/UMAP.',name_for_dim_red,'.clusters.png')
  ggsave(gg1, filename = plot_file, device = 'png', width = 16, height = 16, units = 'cm')


  return(sinbad_object)

}


wrap_plot_features <- function(sinbad_object, features,
                               dim_red_annot_type = 'Bins_100Kb',
                               features_annot_type = 'promoters',
                               legend_title = 'Met. Rate')
{

  dim_red_object = sinbad_object$dim_red_objects[[dim_red_annot_type]]

  feature_matrix = sinbad_object$met_matrices[[features_annot_type]]
  feature_matrix[1:5, 1:5]


  #if(features_annot_type == 'promoters')
  #{
  #  feature_matrix = 1 - feature_matrix
  #  legend_title = 'De-methylation'
  #}
  #if(features_annot_type == 'MAPLE')
  #{

    #legend_title = 'MAPLE'
  #}

  plot_features(umap = dim_red_object$umap,
                feature_matrix = feature_matrix,
                features = features,
                name_for_dim_red = dim_red_annot_type,
                name_for_features = features_annot_type,
                legend_title = legend_title,
                plot_dir = sinbad_object$plot_dir,
                clusters = dim_red_object$clusters)


}


wrap_dmr_analysis <- function(sinbad_object,
                              dim_red_annot_type = 'Bins_100Kb',
                              feature_annot_type = 'Gene_Body',
                              remove_noncoding_genes = T,
                              heatmap_gene_count_per_cluster = 10
                              )
{
  #DM Analysis
  dim_red_object = sinbad_object$dim_red_objects[[dim_red_annot_type]]
  cluster_vec = dim_red_object$clusters
  feature_mat = sinbad_object$met_matrices[[feature_annot_type]]

  head(cluster_vec)
  feature_mat[1:5, 1:5]
  sum(!is.na(feature_mat))

  sinbad_object$dm_stat_list_for_clusters = dm_stat_test_for_clusters(cluster_vec = cluster_vec ,
                                                                      feature_mat = feature_mat,
                                                                      minimum_cell_count = 10,
                                                                      dmr_adj_p_value_cutoff = 0.01,
                                                                      dmr_log2_fc_cutoff = 0.25)
  sinbad_object$dm_result_file = paste0(sinbad_object$summary_dir, '/DM_Analysis.',
                                        '.clusters.', dim_red_annot_type,
                                        '.feature_annot_type.',feature_annot_type,'.xlsx')
  write.xlsx(sinbad_object$dm_stat_list_for_clusters$dm_result_list_with_summary,
             file = sinbad_object$dm_result_file)


  de_list = sinbad_object$dm_stat_list_for_clusters$dm_result_list_with_summary[2:length(sinbad_object$dm_stat_list_for_clusters$dm_result_list_with_summary)]

  short_list = list()
  counter = 1
  for(dm_stat in de_list)
  {
    dm_stat.short = dm_stat[abs(dm_stat$log2_FC) >= 0.25 & dm_stat$adjusted.p.value < 0.01 ,]
    if(remove_noncoding_genes)
    {
      dm_stat.short = dm_stat.short[!grepl('^AC0', dm_stat.short$region),]
      dm_stat.short = dm_stat.short[!grepl('^AC1', dm_stat.short$region),]
      dm_stat.short = dm_stat.short[!grepl('^AL0', dm_stat.short$region),]
      dm_stat.short = dm_stat.short[!grepl('-AS1$', dm_stat.short$region),]
      dm_stat.short = dm_stat.short[!grepl('^LINC', dm_stat.short$region),]
      dm_stat.short = dm_stat.short[!grepl('\\.', dm_stat.short$region),]

    }

    dm_stat.short = dm_stat.short[1:min(heatmap_gene_count_per_cluster, nrow(dm_stat.short)),]
    short_list[[counter]] =dm_stat.short
    counter = counter  + 1
  }

  all_dms = do.call("rbind", short_list )
  all_dms.uniq = all_dms[!duplicated(all_dms$region), ]
  print(all_dms.uniq)
  rownames(all_dms.uniq) = all_dms.uniq$region

  #all_dms.uniq[all_dms.uniq$region == 'TYRO3', ]
  #all_dms.uniq['TBR1', ]
  #all_dms.uniq[all_dms.uniq$region == 'SLC6A1', ]
  #all_dms.uniq[all_dms.uniq$region == 'SLC17A7', ]

  head(all_dms.uniq)

  #all_dms.uniq = all_dms.uniq[abs(all_dms.uniq$log2_FC) >= 0.5, ]
  head(all_dms.uniq)

  dim(all_dms)
  dim(all_dms.uniq)

  #all_dms.uniq['BAIAP2',]
  #all_dms.uniq['CELF2',]

  heatmap_dm_regions = rownames(all_dms.uniq)


  #sig_regions = dm_stat_list_for_clusters$sig_ids
  sorted_clusters = sort(cluster_vec)

  feature_mat = replace_nas_by_column_mean(feature_mat)
  sig_mat = as.matrix(feature_mat[heatmap_dm_regions, names(sorted_clusters) ])
  df_annot = data.frame(cluster =  sorted_clusters)

  head(df_annot)
  sig_mat[1:5, 1:5]
  num_clus = length(unique(sorted_clusters))
  colors1 = hue_pal()(num_clus)
  names(colors1) = as.character(1:num_clus)

  ann_colors = list(
    cluster = colors1
  )

  ph <- pheatmap::pheatmap( sig_mat, main = 'Differentially Methylated Regions',
                           cluster_cols = F,
                           cluster_rows = T,
                           show_colnames = F,
                           fontsize = 13,
                           annotation_col = df_annot,
                           #color = viridis(100),
                           color = colorRampPalette(c( "lightblue", "navy"))(100),
                           annotation_colors = ann_colors[1],
                           fontsize_row = 9)

  dir.create(paste0(sinbad_object$plot_dir, '/DMR/'), showWarnings = F)
  plot_file = paste0(sinbad_object$plot_dir, '/DMR/Heatmap.DMR',
                     '.clusters.', dim_red_annot_type,
                     '.feature_annot_type.',feature_annot_type,'.eps')
  ggsave(ph, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')

  return(sinbad_object)


}


generate_Seurat_object <- function(sinbad_object,
                                   dim_red_annot_type = 'Bins_100Kb',
                                   feature_annot_type = 'Gene_Body',
                                   features_to_plot = c(),
                                   vmr_count = 2000, max_ratio_of_na_cells = 0.75)
{

  dim_red_object = sinbad_object$dim_red_objects[[dim_red_annot_type]]
  cluster_vec = dim_red_object$clusters
  dim_red_mat = sinbad_object$met_matrices[[dim_red_annot_type]]

  feature_mat = sinbad_object$met_matrices[[feature_annot_type]]
  met_mat_for_features = feature_mat

  library(Seurat)

  if(feature_annot_type == 'Promoters')
  {
    met_mat_for_features = 1 - feature_mat
    legend_title = 'De-methylation'
  }
  if(feature_annot_type == 'MAPLE')
  {
    met_mat_for_features = feature_mat
    legend_title = 'MAPLE'
  }

  met_mat_for_dim_red = impute_nas(dim_red_mat)
  met_mat_for_features = impute_nas(met_mat_for_features)

  seurat_object = CreateSeuratObject(met_mat_for_dim_red, assay = dim_red_annot_type)

  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = vmr_count)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), npcs = 10)

  DimPlot(seurat_object, reduction = "pca")

  seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
  seurat_object <- FindClusters(seurat_object, resolution = 0.1)

  seurat_object <- RunUMAP(seurat_object, dims = 1:10, min.dist = 0.01)

  head(seurat_object@reductions$umap@cell.embeddings)
  head(dim_red_object$umap)

  seurat_object@reductions$umap@cell.embeddings = dim_red_object$umap
  colnames(seurat_object@reductions$umap@cell.embeddings) = c('UMAP_1', 'UMAP_2')
  seurat_object@meta.data$sinbad_clusters = dim_red_object$clusters


  #gg1 = DimPlot(seurat_object, reduction = "umap", label = T)
  #print(gg1)



  gg1 = DimPlot(seurat_object, reduction = "umap", label = T, group.by = 'sinbad_clusters')
  print(gg1)


  dir.create(paste0(sinbad_object$plot_dir, '/Seurat/'))

  plot_file = paste0(sinbad_object$plot_dir, '/Seurat/Seurat.UMAP.DR_',dim_red_annot_type,'.clusters.eps')
  ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')

  plot_file = paste0(sinbad_object$plot_dir, '/Seurat/Seurat.UMAP.DR_',dim_red_annot_type,'.clusters.png')
  ggsave(gg1, filename = plot_file, device = 'png', width = 20, height = 20, units = 'cm')


  feature_assay = CreateAssayObject(met_mat_for_features)

  feature_annot_type = gsub('-', 'u', feature_annot_type)
  feature_annot_type = gsub('\\+', 'd', feature_annot_type)

  seurat_object[[feature_annot_type]] <- feature_assay



  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000, assay = feature_annot_type)
  seurat_object <- ScaleData(seurat_object, assay = feature_annot_type)

  DefaultAssay(seurat_object) = feature_annot_type
  #FeaturePlot(seurat_object, features = 'GNB1')


  feature_plot_dir = paste0(sinbad_object$plot_dir, '/Seurat/gene_expression.dim_red_annot_type.',dim_red_annot_type,
                            '.feature_annot_type.',feature_annot_type,'/')
  dir.create(feature_plot_dir, recursive = T, showWarnings = F)
  dir.create(paste0(feature_plot_dir, '/eps/'))
  dir.create(paste0(feature_plot_dir, '/png/'))


  counts_genes = rownames(seurat_object@assays[[feature_annot_type]]@data)
  head(counts_genes)
  features = intersect(features_to_plot, counts_genes)

  n = length(features)


  if( n > 0 )
  {
    for(feature in features)
    {
      #print(feature)
      #cairo_ps(plot_file, fallback_resolution = 2400)
      #postscript(plot_file, onefile = F, width = 7, height = 6)
      gg1 = FeaturePlot(seurat_object, features = feature)
      print(gg1)
      plot_file = paste0(feature_plot_dir, '/eps/', feature, '.eps')
      ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')
      plot_file = paste0(feature_plot_dir, '/png/', feature, '.png')
      ggsave(gg1, filename = plot_file, device = 'png', width = 20, height = 20, units = 'cm')
      #dev.off()
    }#for(gene in marker_gene s)


    if(n <= 10)
    {
      gg1 = StackedVlnPlot(obj = seurat_object, features = features, group.by = 'seurat_clusters') #, group.by = 'cell_type', assay = 'SCT')
      print(gg1)
      plot_file = paste0(feature_plot_dir, '/eps/', 'Violin', '.eps')
      ggsave(gg1, filename = plot_file, device = 'eps', width = 10, height = 28, units = 'cm')
      plot_file = paste0(feature_plot_dir, '/png/', 'Violin', '.png')
      ggsave(gg1, filename = plot_file, device = 'png', width = 10, height = 28, units = 'cm')
    }else
    {
      for(i in 1:ceiling(n/10))
      {
        start = (i-1)*10 + 1
        end = min(i*10, n)
        gg1 = StackedVlnPlot(obj = seurat_object,
                             features = features[start:end],
                             group.by = 'seurat_clusters', y.max = 0.5) #, group.by = 'cell_type', assay = 'SCT')
        print(gg1)
        plot_file = paste0(feature_plot_dir, '/eps/', 'Violin.',i ,'.eps')
        ggsave(gg1, filename = plot_file, device = 'eps', width = 10, height = 28, units = 'cm')
        plot_file = paste0(feature_plot_dir, '/png/', 'Violin.',i, '.png')
        ggsave(gg1, filename = plot_file, device = 'png', width = 10, height = 28, units = 'cm')
      }

    }#else



  }#if



  return(seurat_object)

}




add_feature_to_Seurat_object <- function(seurat_object, sinbad_object,
                                         feature_annot_type = 'Gene_Body',
                                         features_to_plot = c()
)
{

  feature_mat = sinbad_object$met_matrices[[feature_annot_type]]
  met_mat_for_features = feature_mat

  met_mat_for_features = impute_nas(met_mat_for_features)


  feature_assay = CreateAssayObject(met_mat_for_features)

  feature_annot_type = gsub('-', 'u', feature_annot_type)
  feature_annot_type = gsub('\\+', 'd', feature_annot_type)

  seurat_object[[feature_annot_type]] <- feature_assay



  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000,
                                 assay = feature_annot_type)
  seurat_object <- ScaleData(seurat_object, assay = feature_annot_type)

  DefaultAssay(seurat_object) = feature_annot_type
  #FeaturePlot(seurat_object, features = 'GNB1')


  feature_plot_dir = paste0(sinbad_object$plot_dir, '/Seurat/gene_expression.dim_red_annot_type.',dim_red_annot_type,
                            '.feature_annot_type.',feature_annot_type,'/')
  dir.create(feature_plot_dir, recursive = T, showWarnings = F)
  dir.create(paste0(feature_plot_dir, '/eps/'))
  dir.create(paste0(feature_plot_dir, '/png/'))


  counts_genes = rownames(seurat_object@assays[[feature_annot_type]]@data)
  head(counts_genes)
  features = intersect(features_to_plot, counts_genes)

  n = length(features)


  if( n > 0 )
  {
    for(feature in features)
    {
      #print(feature)
      #cairo_ps(plot_file, fallback_resolution = 2400)
      #postscript(plot_file, onefile = F, width = 7, height = 6)
      gg1 = FeaturePlot(seurat_object, features = feature)
      print(gg1)
      plot_file = paste0(feature_plot_dir, '/eps/', feature, '.eps')
      ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')
      plot_file = paste0(feature_plot_dir, '/png/', feature, '.png')
      ggsave(gg1, filename = plot_file, device = 'png', width = 20, height = 20, units = 'cm')
      #dev.off()
    }#for(gene in marker_gene s)


    if(n <= 10)
    {
      gg1 = StackedVlnPlot(obj = seurat_object, features = features, group.by = 'seurat_clusters') #, group.by = 'cell_type', assay = 'SCT')
      print(gg1)
      plot_file = paste0(feature_plot_dir, '/eps/', 'Violin', '.eps')
      ggsave(gg1, filename = plot_file, device = 'eps', width = 10, height = 28, units = 'cm')
      plot_file = paste0(feature_plot_dir, '/png/', 'Violin', '.png')
      ggsave(gg1, filename = plot_file, device = 'png', width = 10, height = 28, units = 'cm')
    }else
    {
      for(i in 1:ceiling(n/10))
      {
        start = (i-1)*10 + 1
        end = min(i*10, n)
        gg1 = StackedVlnPlot(obj = seurat_object,
                             features = features[start:end],
                             group.by = 'seurat_clusters', y.max = 0.5) #, group.by = 'cell_type', assay = 'SCT')
        print(gg1)
        plot_file = paste0(feature_plot_dir, '/eps/', 'Violin.',i ,'.eps')
        ggsave(gg1, filename = plot_file, device = 'eps', width = 10, height = 28, units = 'cm')
        plot_file = paste0(feature_plot_dir, '/png/', 'Violin.',i, '.png')
        ggsave(gg1, filename = plot_file, device = 'png', width = 10, height = 28, units = 'cm')
      }

    }#else



  }#if



  return(seurat_object)

}




wrap_remove_temporary_files <- function(sinbad_object)
{

  print('Deleting temporary files .... ')
  sys_command = paste0('rm ', sinbad_object$r2_meta_fastq_dir, '/*')
  print(sys_command)
  system(sys_command)

  sys_command = paste0('rm ', sinbad_object$demux_fastq_dir, '/*')
  print(sys_command)
  system(sys_command)

  sys_command = paste0('rm ', sinbad_object$r1_and_r2_alignments_dir, '/*bam')
  print(sys_command)
  system(sys_command)

  sys_command = paste0('rm ', sinbad_object$methylation_calls_dir, '/*txt.gz')
  print(sys_command)
  system(sys_command)



  print('Deleted temporary files .... ')




}

process_sample_wrapper <- function(sinbad_object)
{


  par(mfrow = c(2,2))

  sinbad_object = wrap_demux_fastq_files(sinbad_object)

  sinbad_object = wrap_demux_stats(sinbad_object)

  sinbad_object = wrap_trim_fastq_files(sinbad_object)

  sinbad_object = wrap_trim_stats(sinbad_object)

  sinbad_object = wrap_align_sample(sinbad_object)

  sinbad_object = wrap_generate_alignment_stats(sinbad_object)

  sinbad_object = wrap_call_methylation_sites(sinbad_object)

  sinbad_object = wrap_generate_methylation_stats(sinbad_object)

  sinbad_object = wrap_read_regions(sinbad_object)

  sinbad_object = wrap_quantify_regions(sinbad_object)

  sinbad_object = wrap_dim_red(sinbad_object)

  sinbad_object = wrap_dmr_analysis(sinbad_object)


  return(sinbad_object)


}


