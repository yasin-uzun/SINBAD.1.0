library(doSNOW)
#library(vioplot)
library(readr)
library(Rsamtools)
library(ShortRead)


align_sample <- function(read_dir,
                         genomic_sequence_path,
                         alignment_dir,
                         aligner,
                         num_cores,
                         mapq_threshold,
                         main_log_dir,
                         aligner_param_settings = '',
                         pattern = '')
{
  alignment_log_dir= paste0(main_log_dir, '/alignment/')
  dir.create(alignment_log_dir, recursive = T)
  setwd(alignment_dir)

  fastq_files_1 = list.files(read_dir, pattern = "*.fastq.gz")
  fastq_files_2 = list.files(read_dir, pattern = "*.fastq")

  fastq_files = union(fastq_files_1, fastq_files_2)
  fastq_files = fastq_files[grepl(pattern, fastq_files)]

  num_fastq = length(fastq_files)

  print('***********************************************************************')
  print(paste('There are', num_fastq, 'fastq files.'))
  print('***********************************************************************')

  if(num_cores > 1)
  {
    cl <- makeCluster(num_cores, outfile="", type = 'SOCK')

    clusterExport(cl, ls(.GlobalEnv))
    registerDoSNOW(cl)
    clusterExport(cl, ls(.GlobalEnv))

    foreach(i=1:length(fastq_files)) %dopar%
      {
        fastq_file = fastq_files[i]
        print( fastq_file)
        bam_file = paste0(alignment_dir, fastq_file)
        bam_file = gsub('.fastq.gz', '.organism.bam', bam_file)
        print(bam_file)

        if(!file.exists(bam_file))
        {
          align_cell(read_dir = read_dir,
                     fastq_file = fastq_file,
                     aligner = aligner,
                     genomic_sequence_path = genomic_sequence_path,
                     alignment_dir = alignment_dir,
                     log_dir = alignment_log_dir,
                     aligner_param_settings = aligner_param_settings)
        }

      }#foreach
  }else#if(num_cores > 1)
  {
    for(i in 1:length(fastq_files))
    {
      fastq_file = fastq_files[i]
      print('*****************************************')

      #bam_file = paste0(alignment_dir, fastq_file)
      #bam_file = gsub('\\.fastq\\.gz', '\\.organism\\.bam', bam_file)
      print(fastq_files)
      #print(bam_file)

      print('*****************************************')
      if(!file.exists(bam_file))
      {
        align_cell(read_dir = read_dir,
                   fastq_file = fastq_file,
                   aligner = aligner,
                   genomic_sequence_path = genomic_sequence_path,
                   alignment_dir = alignment_dir,
                   aligner_param_settings = aligner_param_settings,
                   log_dir = alignment_log_dir)
      }

    }#foreach



  }#else


  aligner_log_dir = paste0(alignment_log_dir, '/bismark_aligner/')
  failed_alignments = find_failed_alignments(aligner_log_dir, read_dir, pattern = pattern)

  count_failed_alignments = length(failed_alignments)

  return(count_failed_alignments)

}#run_alignment

#aligner_log_dir = paste0(sinbad_object$main_log_dir, '/alignment/bismark_aligner/')
#trimmed_fastq_dir = sinbad_object$trimmed_fastq_dir

find_failed_alignments <- function(aligner_log_dir, read_dir, pattern = '')
{
  log_files = list.files(aligner_log_dir, "*.log")
  log_files = log_files[grepl(pattern, log_files)]
  length(log_files)
  setwd(aligner_log_dir)

  log_matches <- sapply(log_files, FUN=function(x){
    grep("Alignment is successful", readLines(x))
  })
  length(log_matches)
  matches = unlist(log_matches)
  head(log_matches)
  sum(log_matches == 41)
  completed =  names(log_matches)[log_matches>0]
  successfull_fastq = gsub('.log', '', completed)

  all_fastq_files = list.files(read_dir, full.names = FALSE)
  all_fastq_files = all_fastq_files[grepl(pattern, all_fastq_files)]

  failed_fastq = setdiff(all_fastq_files, successfull_fastq)

  return(failed_fastq)

}#find_failed_alignments


align_cell <- function(read_dir, fastq_file, aligner, genomic_sequence_path,
                       alignment_dir, log_dir, aligner_param_settings = ''){

  print(aligner_param_settings)
  print('align_cell')
  cat(' -- fastq_file: ', fastq_file, '\n')

  cell_id = gsub('.fastq.gz', '', fastq_file)

  #Align
  if(aligner == 'bismark')
  {

    run_bismark_aligner(read_dir = read_dir,
                        fastq_file_left = fastq_file,
                        genomic_sequence_path = genomic_sequence_path,
                        alignment_dir = alignment_dir,
                        cell_id = cell_id, log_dir = log_dir,
                        aligner_param_settings = aligner_param_settings)

  }else if(aligner == 'bsmap')
  {

    run_bsmap_aligner(read_dir = read_dir,
                      fastq_file = fastq_file,
                      genomic_sequence_path = genomic_sequence_path,
                      alignment_dir = alignment_dir,
                      cell_id = cell_id, log_dir = log_dir)
  }else if(aligner == 'bs_seeker')
  {

    run_bs_seeker_aligner(read_dir = read_dir,
                          fastq_file = fastq_file,
                          genomic_sequence_path = genomic_sequence_path,
                          alignment_dir = alignment_dir,
                          cell_id = cell_id, log_dir = log_dir)
  }else
  {
    msg = paste0('Unrecognized aligner name: ', aligner, ' . Should be one of: bismark, bs_map, bs_seeker')
    stop(msg)
  }

  gc()
  #
  #
  #Mapq filtering
  if(mapq_threshold > 0)
  {
    filter_mapq(alignment_dir, cell_id, mapq_threshold = mapq_threshold, log_dir = log_dir)
  }#if(mapq_threshold > 0)
  #
  use_mapq_filtered_reads = F
  if(mapq_threshold > 0)
  {
    use_mapq_filtered_reads = T
  }

  gc()


  remove_duplicate_reads(alignment_dir = alignment_dir,
                         cell_id = cell_id,
                         use_mapq_filtered_reads = use_mapq_filtered_reads,
                         duplicate_remover = duplicate_remover, log_dir = log_dir)
  gc()


  filter_non_conversion(alignment_dir = alignment_dir,
                        cell_id = cell_id, log_dir = log_dir)

  gc()


  split_lambda(alignment_dir = alignment_dir, cell_id = cell_id, log_dir = log_dir)

  gc()




}



filter_cell <- function(read_dir, fastq_file, aligner, genomic_sequence_path,
                        alignment_dir, log_dir){

  use_mapq_filtered_reads = T
  print('Filtering cell')
  cat(' -- fastq_file: ', fastq_file, '\n')

  cell_id = gsub('.fastq.gz', '', fastq_file)


  print(paste('*Duplicate remover: ', duplicate_remover))


  remove_duplicate_reads(alignment_dir = alignment_dir,
                         cell_id = cell_id,
                         use_mapq_filtered_reads = use_mapq_filtered_reads,
                         duplicate_remover = duplicate_remover, log_dir = log_dir)
  gc()


  filter_non_conversion(alignment_dir = alignment_dir,
                        cell_id = cell_id, log_dir = log_dir)

  gc()


  split_lambda(alignment_dir = alignment_dir, cell_id = cell_id, log_dir = log_dir)

  gc()




}


run_bismark_aligner <- function(read_dir, fastq_file_left, fastq_file_right = NULL,
                                genomic_sequence_path, alignment_dir,
                                aligner_param_settings = '',
                                cell_id, log_dir)
{

  setwd(alignment_dir)

  log_sub_dir = paste0(log_dir, '/bismark_aligner/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, fastq_file_left, '.log')

  sys_command = paste0('bismark --bowtie2 ', aligner_param_settings
                       ,' --fastq  '
                       ,' --basename ', cell_id
                       ,' --bam ', genomic_sequence_path
                       ,'  ', read_dir, fastq_file_left
                       ,' > ', log_file
  )

  cat('Aligning ', fastq_file_left, '\n')
  print(sys_command)
  command_result = system(sys_command)

  sink(log_file, append = T)
  if(command_result == 0)
  {
    cat('Alignment is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR:Bismark aligner failed to run for cell ', cell_id ,'. Exiting the pipeline. Please see the output log.')
  }
  sink()

  return(command_result)
}#run_alignment



run_bs_seeker_aligner <- function(read_dir, fastq_file_left, fastq_file_right = NULL,
                                  genomic_sequence_path, alignment_dir,
                                  cell_id, log_dir, is_paired = F)
{

  #setwd(alignment_dir)
  setwd(bs_seeker_path)

  log_sub_dir = paste0(log_dir, '/align/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, fastq_file_left, '.log')

  if(is.null(fastq_file_right)) #single ended
  {
    sys_command = paste0('./bs3-align  '
                         ,' -i ', read_dir, fastq_file_left
                         ,' -o ', alignment_dir, cell_id, '.bam'
                         ,' ', bs_seeker_aligner_param_settings
                         ,' -g  ', bs_seeker_genome_dir, bs_seeker_genome_fasta
                         ,' --db=', bs_seeker_genome_dir
                         #,' > ', log_file
    )



  }else #pair ended
  {
    sys_command = paste0('bs3-align  '
                         ,' -1', read_dir, fastq_file_left
                         ,' -2', read_dir, fastq_file_right
                         ,' -o ', alignment_dir, cell_id, '.bam'
                         ,' ', bs_seeker_aligner_param_settings
                         ,' -g  ', bs_seeker_genome_dir, bs_seeker_genome_fasta
                         ,' --db=', bs_seeker_genome_dir
                         #,' > ', log_file
    )


  }




  cat('Aligning ', fastq_file, '\n')
  print(sys_command)
  command_result = system(sys_command)

  sink(log_file, append = T)
  if(command_result == 0)
  {
    cat('Alignment is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR:BS_Seeker aligner failed to run for cell ', cell_id ,'. Exiting the pipeline. Please see the output log.')
  }
  sink()

  return(command_result)
}#run_alignment

run_bsmap_aligner <- function(read_dir, fastq_file_left, fastq_file_right = NULL,
                              genomic_sequence_path, alignment_dir,
                              cell_id, log_dir, is_paired = F)
{

  #setwd(alignment_dir)

  log_sub_dir = paste0(log_dir, '/align/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, fastq_file_left, '.log')

  if(is.null(fastq_file_right)) #single ended
  {
    sys_command = paste0(bsmap_path, '/bsmap  '
                         ,' -a ', read_dir, fastq_file_left
                         ,' -d  ', reference_genome_fasta
                         ,' -o ', alignment_dir, cell_id, '.bam'
                         ,' ', bsmap_aligner_param_settings

                         #,' > ', log_file
    )





  }else #pair ended
  {
    sys_command = paste0(bsmap_path, '/bsmap  '
                         ,' -a ', read_dir, fastq_file_left
                         ,' -b ', read_dir, fastq_file_right
                         ,' -d  ', reference_genome_fasta
                         ,' -o ', alignment_dir, cell_id, '.bam'
                         ,' ', bsmap_aligner_param_settings

                         #,' > ', log_file
    )



  }




  cat('Aligning ', fastq_file, '\n')
  print(sys_command)
  command_result = system(sys_command)

  sink(log_file, append = T)
  if(command_result == 0)
  {
    cat('Alignment is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR:Bismark aligner failed to run for cell ', cell_id ,'. Exiting the pipeline. Please see the output log.')
  }
  sink()

  return(command_result)
}#run_alignment



filter_mapq <- function(alignment_dir, cell_id, mapq_threshold = 10, log_dir)
{
  setwd(alignment_dir)
  command_result = 0
  log_sub_dir = paste0(log_dir, '/filter_mapq/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, cell_id, '.log')

  sink(log_file)
  if(mapq_threshold > 0)
  {
    cat('Filtering alignments mapq > ', mapq_threshold, ' :' , cell_id, '\n')
    sys_command = paste0('samtools view -bhq ',mapq_threshold,' ',cell_id,'.bam > ',cell_id,'.mapq_filtered.bam')
    command_result = system(sys_command)
  }else
  {
    print(paste0('Mapq threshold is 0. No mapq filtering done for ', cell_id))
  }
  sink()

  if(command_result == 0)
  {
    cat('Read filtering is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR: Read filtering failed to run ', cell_id, '. Exiting the pipeline. Please see the output log.')
  }


  return(command_result)
}





remove_duplicate_reads <- function(alignment_dir, cell_id,
                                   use_mapq_filtered_reads = T,
                                   duplicate_remover = 'samtools'
                                   , log_dir)
{
  bam_file = paste0(cell_id,'.bam')
  if(use_mapq_filtered_reads)
  {
    bam_file = paste0(cell_id,'.mapq_filtered.bam')
  }
  rmdup_file = paste0(cell_id,'.rmdup.bam')

  log_sub_dir = paste0(log_dir, '/rmdup/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, cell_id, '.log')

  cat('Using ', duplicate_remover, ' to remove duplicates...')

  if(duplicate_remover == 'picard')
  {
    sorted_file = paste0(cell_id,'.sorted.bam')
    metric_file = paste0(cell_id,'.picard_metrics.txt')

    sys_command = paste0('java -jar ',picard_path,' SortSam ',
                         ' INPUT=', bam_file,
                         ' OUTPUT=', sorted_file,
                         ' SORT_ORDER=coordinate ',
                         ' >> ', log_file)

    command_result = system(sys_command)

    if(command_result == 0)
    {
      cat('Sorted alignments for cell ', cell_id, '\n')
    }else
    {
      stop(paste('ERROR: Picard sort failed. Could not sort alignments (bam file) for cell ', cell_id, '\n'))
    }

    sys_command = paste0('java -jar ',picard_path,' MarkDuplicates ',
                         ' INPUT=', sorted_file,
                         ' OUTPUT=', rmdup_file,
                         ' REMOVE_DUPLICATES=true ',
                         ' METRICS_FILE=',metric_file)

    command_result = system(sys_command)

    if(command_result == 0)
    {
      cat('Removed duplicate reads for cell ', cell_id, '\n')
    }else
    {
      stop(paste('ERROR: Picard remove duplicates failed. Could not sort alignments (bam file) for cell ', cell_id, '\n'))
    }

  }#if(duplicate_remover == 'picard')

  if(duplicate_remover == 'samtools')
  {
    sys_command = paste0('samtools rmdup -S ', bam_file, '  ', rmdup_file)
    command_result = system(sys_command)

    if(command_result == 0)
    {
      cat('Removed duplicate reads for cell ', cell_id, '\n')
    }else
    {
      stop(paste('ERROR: Samtools remove duplicates failed. Could not sort alignments (bam file) for cell ', cell_id, '\n'))
    }

  }

}#remove_duplicate_reads



filter_non_conversion <- function(alignment_dir, cell_id, log_dir)
{
  setwd(alignment_dir)
  command_result = 0
  rmdup_file = paste0(cell_id,'.rmdup.bam')

  log_sub_dir = paste0(log_dir, '/filter_non_conversion/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)
  log_file = paste0(log_sub_dir, cell_id, '.log')

  cat('Filtering nonconversion reads  :' , cell_id, '\n')
  sys_command = paste0(bismark_path, '/filter_non_conversion  '
                       , ' --samtools_path ', samtools_path
                       , ' ', filter_non_conversion_param_settings
                       , ' ', rmdup_file
                       , ' > ', log_file)
  command_result = system(sys_command)



  if(command_result == 0)
  {
    cat('Non-conversion filtering is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR: Non-conversion filtering failed to run for ', cell_id, '.\n Exiting the pipeline. Please see the output log.')
  }

  return(command_result)
}


split_lambda <- function(alignment_dir, cell_id, log_dir)
{
  setwd(alignment_dir)
  command_result = 0
  noncon_file = paste0(cell_id,'.rmdup.nonCG_filtered.bam')
  sorted_file = paste0(cell_id,'.rmdup.nonCG_filtered.sorted.bam')

  lambda_file = paste0(cell_id,'.lambda.bam')
  organism_file = paste0(cell_id,'.organism.bam')

  log_sub_dir = paste0(log_dir, '/split_lambda/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)
  log_file = paste0(log_sub_dir, cell_id, '.log')

  sink(log_file)
  cat('Splitting lambda control  :' , cell_id, '\n')
  sys_command = paste0('samtools sort -o ', sorted_file, ' -T tmp_sort_', cell_id, ' ' , noncon_file)
  print(sys_command)
  command_result_0 = system(sys_command)
  sys_command = paste0('samtools index -b ', sorted_file)
  command_result_1 = system(sys_command)
  sys_command = paste0('samtools view -hb  ', sorted_file, ' ' , lambda_chrom_numbers , ' > ', lambda_file)
  command_result_2 = system(sys_command)
  sys_command = paste0('samtools view -hb  ', sorted_file, ' ' , organism_chrom_numbers , ' > ', organism_file)
  command_result_3 = system(sys_command)
  sink()

  sorted_organism_file = gsub('organism.bam', 'organism.sorted.bam', organism_file)
  sys_command = paste0('samtools sort -o ', sorted_organism_file, ' -T tmp_sort_', cell_id, ' ' , organism_file)
  print(sys_command)
  command_result_4 = system(sys_command)

  sys_command = paste0('samtools index -b ', sorted_organism_file)
  print(sys_command)
  command_result_5 = system(sys_command)


  command_result = command_result_1 + command_result_2 + command_result_3

  if(command_result == 0)
  {
    cat('Lambda control split is successful for cell ', cell_id, '\n')
    cat('Output (organism): ', organism_file, '\n')
  }else
  {
    stop('Lambda control splitting failed to run for ', cell_id, '. Exiting the pipeline. Please see the output log.')
  }

  return(command_result)

}



split_lambda_old <- function(alignment_dir, cell_id, log_dir)
{
  setwd(alignment_dir)
  command_result = 0
  noncon_file = paste0(cell_id,'.rmdup.nonCG_filtered.bam')
  sorted_file = paste0(cell_id,'.rmdup.nonCG_filtered.sorted.bam')

  lambda_file = paste0(cell_id,'.lambda.bam')
  organism_file = paste0(cell_id,'.organism.bam')

  log_sub_dir = paste0(log_dir, '/split_lambda/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)
  log_file = paste0(log_sub_dir, cell_id, '.log')

  sink(log_file)
  cat('Splitting lambda control  :' , cell_id, '\n')
  sys_command = paste0('samtools sort ', sorted_file, ' -T tmp_sort_', cell_id, ' ' , noncon_file)
  command_result_0 = system(sys_command)
  sys_command = paste0('samtools index -b ', noncon_file)
  command_result_1 = system(sys_command)
  sys_command = paste0('samtools view -hb  ', noncon_file, ' ' , lambda_chrom_numbers , ' > ', lambda_file)
  command_result_2 = system(sys_command)
  sys_command = paste0('samtools view -hb  ', noncon_file, ' ' , organism_chrom_numbers , ' > ', organism_file)
  command_result_3 = system(sys_command)
  sink()

  command_result = command_result_1 + command_result_2 + command_result_3

  if(command_result == 0)
  {
    cat('Lambda control split is successful for cell ', cell_id, '\n')
  }else
  {
    stop('Lambda control splitting failed to run for ', cell_id, '. Exiting the pipeline. Please see the output log.')
  }

  return(command_result)

}


process_bismark_alignment_reports <- function(alignment_dir)
{
  setwd(alignment_dir)
  report_files = list.files(alignment_dir, pattern = "*SE_report.txt")
  row_names = c()
  result_list = list()
  for(report_file in report_files)
  {
    #report_file = "GSM3444126_R2_SE_report.txt"
    print(paste('Reading ', report_file))

    cell_id = gsub('_SE_report.txt', '', report_file)

    lines = readLines(report_file)
    informative_lines = lines[grep('\t', lines)]
    informative_lines = informative_lines[!grepl('strand', informative_lines)]
    informative_lines = informative_lines[!grepl('C methylated in Unknown context', informative_lines)]

    if( length(informative_lines) > 0 )
    {
      df_report = read.delim(text = informative_lines, sep = '\t',  header = F)
      class(df_report)
      row_names = gsub(':', '', df_report$V1)
      result_list[[cell_id]] = as.character(df_report$V2)
    }#if
  }#for

  row_names[1] = "Total_reads"
  row_names[2] = "Alignments"
  row_names[3] = "Alignment_rate"

  row_names[4] = "Sequences_with_no_alignments"
  row_names[6] = "Discarded_sequences"

  #unlist(lapply(result_list, length))
  table(unlist(lapply(result_list, length)))


  df_result = do.call(cbind, result_list)
  rownames(df_result) = row_names
  t_df_result =  t(df_result)
  final_table = data.frame(Cell_ID = rownames(t_df_result),  t_df_result, stringsAsFactors = F)

  return(final_table)

}#process_bismark_alignment_reports

#df_alignment_stats = final_table

plot_alignment_stats <- function(sample_name, df_alignment_stats)
{

  par(mfrow = c(2,3))

  font.face = 'Helvetica'
  font.size = 1.5
  color_vec = hue_pal()(3)
  sep.lwd = 1
  #sep.color = "#33a02c"

  sep.color = "#ccebc5"

  overall_stat = c()
  ###Stats##
  total_reads = as.numeric(df_alignment_stats$Total_reads)
  aligned_reads = as.numeric(df_alignment_stats$Alignments)
  mapq_read_counts = df_alignment_stats$mapq_read_counts
  rmdup_read_counts = df_alignment_stats$rmdup_read_counts
  nc_filtered_read_counts = df_alignment_stats$nc_filtered_read_counts
  organism_read_counts = df_alignment_stats$organism_read_counts

  total_reads.sum = sum(total_reads)
  aligned_reads.sum = sum(aligned_reads)
  mapq_reads.sum = sum(mapq_read_counts)
  overall.alignment.rate = round(aligned_reads.sum / total_reads.sum * 100, 1)

  overall_stat = c( 'Sample_Name' = sample_name
                   ,  'Total_reads' = total_reads.sum
                   , 'Aligned_reads' = aligned_reads.sum
                   , 'Mapq_reads' = mapq_reads.sum
                   , 'Overall_alignment_rate' = overall.alignment.rate

                   )

  par(mar = c(2,4,6,2))
  plot.new()

  mtext(side =3 ,line=1, cex.main=font.size, family =font.face, font=2, adj=0, "Alignment Statistics:")

  lbl = paste0("All reads: ",format(total_reads.sum, scientific = F, big.mark = ","))
  mtext(side =3 , line=-2, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  lbl = paste0("Aligned reads: ",format(aligned_reads.sum, scientific = F, big.mark = ","))
  mtext(side =3 , line=-4, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  lbl = paste0("Overall alignment rate: ",format(overall.alignment.rate, scientific = F, big.mark = ","), '%')
  mtext(side =3 , line=-6, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  lbl = paste0("High quality reads: ",format(mapq_reads.sum, scientific = F, big.mark = ","))
  mtext(side =3 , line=-8, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)

  box("figure", col=sep.color,  lwd = sep.lwd)


  ##################Plot alignment rates##################

  par(mar = c(6,6,6,3))

  if(is.character( df_alignment_stats$Alignment_rate) )
  {
    df_alignment_stats$Alignment_rate =   parse_number(df_alignment_stats$Alignment_rate)
  }

  alignment_rates = df_alignment_stats$Alignment_rate
  roof = max(alignment_rates) * 1.2

  par(xpd = F)
  dens_aln = density(alignment_rates)
  plot(dens_aln, main="", xlab = "% Alignment rate", ylab = "", las = 1,
       cex.lab = font.size, cex.axis = font.size)
  polygon(dens_aln, col="#a6cee3", border="#a6cee3")

  title(ylab="Density", line=4, cex.lab=font.size, family=font.face)

  median_alignment_rate = median(alignment_rates)
  abline(v = median_alignment_rate, col = 'red', lty = 2, lwd = 2)

  #legend('topleft', legend = c('Median'),
  #       col = c('red'), cex = font.size,
  #       lty = 2, bty = "n", inset = c(-0.15, 0))

  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Alignment Rate Distr.")

  box("figure", col=sep.color,  lwd = sep.lwd)


  #############Plot number of reads#######

  bb = boxplot(
    total_reads/1000000,
    #aligned_reads/1000000,
    plot = F)

  roof = max(bb$stats) * 1.2

  x_labels = c('Total', 'Aligned', 'Mapq.', 'Rm-dup.', 'NC filt.', 'Org.')

  par(mar = c(8,6,6,6))

  boxplot(
    total_reads/1000000,
    aligned_reads/1000000,
    mapq_read_counts/1000000,
    rmdup_read_counts/1000000,
    nc_filtered_read_counts/1000000,
    organism_read_counts/1000000,
    ylim = c(0, roof),
    col = "#ffff99", names = x_labels,
    ylab = 'million reads',

    #main = 'Read statistics',
    las = 2,
    cex.axis = font.size, cex.names = font.size,
    cex.main = font.size, cex.lab = font.size,
    outline =F,  medlwd = 1
  )    #Plot number of reads


  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Read counts")

  box("figure", col=sep.color,  lwd = sep.lwd)




  ###Coverage rates#################

  par(mar = c(8,9,3,2))



  coverage_rates =  df_alignment_stats$coverage_rate
  coverage_rates = coverage_rates[!is.na(coverage_rates)]
  median_coverage_rate = median(coverage_rates)
  par(xpd = F)

  dens_cov = density(coverage_rates)
  plot(dens_cov, main="", xlab = "% Coverage", ylab = "", las = 1,
       cex.lab = font.size, cex.axis = font.size)
  polygon(dens_cov, col="#a6cee3", border="#a6cee3")

  title(ylab="Density", line=4, cex.lab=font.size, family=font.face)

  abline(v = median_coverage_rate, col = 'red', lty = 2, lwd = 2)

  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Genomic Coverage")

  box("figure", col=sep.color,  lwd = sep.lwd)

  #df_alignment_stats[, c('Cell_ID', 'coverage_rate')]


  ###Scatter plot filtering#################
  par(mar = c(8,6,3,2))

  filter_aln_rate =  df_alignment_stats$Alignment_rate  > alignment_rate_threshold
  #filter_read_count = df_alignment_stats$organism_read_counts > organism_minimum_filtered_read_count
  filter_read_count = df_alignment_stats$nc_filtered_read_counts > minimum_filtered_read_count

  passing = filter_aln_rate &  filter_read_count

  df_alignment_stats$pass = 0
  df_alignment_stats$pass[passing] = 1

  Alignment.rate = df_alignment_stats$Alignment_rate[passing]
  Filtered.reads = df_alignment_stats$nc_filtered_read_counts[passing] / 1000000
  plot(Alignment.rate,
       Filtered.reads,
       xlim = c(0, max(Alignment.rate)),
       ylim = c(0, max(Filtered.reads)),
       pch = 20, col = 'darkgreen',
       xlab = '% Alignment', ylab = 'Filtered reads (million)'
       , cex.lab = font.size, cex.axis = font.size
       , main = ""
       , las = 1
  )

  Alignment.rate = df_alignment_stats$Alignment_rate[!passing]
  Filtered.reads = df_alignment_stats$nc_filtered_read_counts[!passing] / 1000000
  points(Alignment.rate,
         Filtered.reads,
         pch = 20, col = 'red')

  abline(v = alignment_rate_threshold, lty = 2)
  abline(h = minimum_filtered_read_count/ 1000000 + 0.01, lty = 2)

  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Cell QC")

  box("figure", col=sep.color,  lwd = sep.lwd)

  #####################Filtering pie chart#############

  par(mar = c(2,2,2,2))

  cnt_failed_low_aln = sum(!filter_aln_rate & filter_read_count)
  cnt_failed_low_reads = sum(filter_aln_rate & !filter_read_count)
  cnt_failed_low_aln_and_low_reads = sum(!filter_aln_rate & !filter_read_count)

  cnt_passed = sum(passing)
  cnt_failed = sum(!passing)

  #slices = c(cnt_passed, cnt_failed_low_aln, cnt_failed_low_reads, cnt_failed_low_aln_and_low_reads)
  slices = c(cnt_passed, cnt_failed)
  #lbls = c('Passed', 'Low Aln', 'Low Reads', 'Low Aln and Reads')
  lbls = c('Passed', 'Discard')

  lbls <- paste0(lbls, "\n(", slices, ")") # add percents to labels
  lbls <- paste(lbls,"",sep="") # ad % to labels

  nice.pie(slices, labels = lbls, col=c('#b2df8a', '#fb9a99'), cex = font.size,
           text_col ='black', main = ''  )

  mtext(line=-1, cex.main=font.size, family =font.face, font=2, "Cell filtering")
  box("figure", col=sep.color,  lwd = sep.lwd)

  overall_stat['Passed'] = cnt_passed
  overall_stat['Discard'] = cnt_failed

  #############OUTER BOX############

  box("outer", col="#ccebc5",  lwd = 75)
  box("outer", col="black",  lwd = 3)

  title(paste(sample_name, ': Alignment results'), line = -1.5, outer = TRUE,
        cex.main=font.size * 1.2, family=font.face)


  return(overall_stat)
}



merge_r1_and_r2_alignment_stats <- function(df_alignment_stats)
{
  df_alignment_stats$Row.names = NULL
  stat_columns = colnames(df_alignment_stats)
  rate_columns = stat_columns[grepl('rate', stat_columns) | grepl('^C\\.', stat_columns)]
  count_columns = setdiff(stat_columns, rate_columns)
  df_count_stats = df_alignment_stats[, count_columns]

  cell_id_wo_r   = df_count_stats$Cell_ID
  cell_id_wo_r = gsub('_R1', '', cell_id_wo_r)
  cell_id_wo_r = gsub('_R2', '', cell_id_wo_r)

  df_count_stats$Cell_ID = cell_id_wo_r

  for(count_column in setdiff(count_columns, 'Cell_ID' ) )
  {
    df_count_stats[, count_column] = as.integer(df_count_stats[, count_column])
  }

  attach(df_count_stats)
  df_count_stats_sum = aggregate(as.matrix(df_count_stats[, 2:ncol(df_count_stats)]),
                                 by=list(Cell_ID), FUN = sum, na.rm = TRUE)
  detach(df_count_stats)
  colnames(df_count_stats_sum)[1] = 'Cell_ID'

  df_count_stats_sum$Alignment_rate  = round(df_count_stats_sum$Alignments / df_count_stats_sum$Total_reads * 100)
  df_chrom_sizes = read.table(chrom_sizes_file, header = F)
  genome_size = sum(df_chrom_sizes$V2)

  df_count_stats_sum$coverage_rate  = -1

  attach(df_count_stats_sum)
  df_count_stats_sum$C.methylated.in.CpG.context = round(100 * Total.methylated.C.s.in.CpG.context / (Total.methylated.C.s.in.CpG.context + Total.unmethylated.C.s.in.CpG.context), 1)
  df_count_stats_sum$C.methylated.in.CHG.context = round(100 * Total.methylated.C.s.in.CHG.context / (Total.methylated.C.s.in.CHG.context + Total.unmethylated.C.s.in.CHG.context), 1)
  df_count_stats_sum$C.methylated.in.CHH.context = round(100 * Total.methylated.C.s.in.CHH.context / (Total.methylated.C.s.in.CHH.context + Total.unmethylated.C.s.in.CHH.context), 1)

  detach(df_count_stats_sum)

  rownames(df_count_stats_sum) = df_count_stats_sum$Cell_ID
  head(df_count_stats_sum)



  ###Scatter plot filtering
  filter_aln_rate =  df_count_stats_sum$Alignment_rate  > alignment_rate_threshold
  #filter_read_count = df_alignment_stats$organism_read_counts > organism_minimum_filtered_read_count
  filter_read_count = df_count_stats_sum$nc_filtered_read_counts > minimum_filtered_read_count

  passing = filter_aln_rate &  filter_read_count

  df_count_stats_sum$pass = 0
  df_count_stats_sum$pass[passing] = 1

  return(df_count_stats_sum)
}


count_bam_files <- function(alignment_dir)
{
  setwd(alignment_dir)
  mapq_files = list.files(pattern = '*.mapq_filtered.bam')
  rmdup_files = list.files(pattern = '*.rmdup.bam')
  nc_filtered_files = list.files(pattern = '.*.nonCG_filtered.bam$')
  organism_bam_files = list.files(pattern = '.*.organism.bam$')

  #cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
  #registerDoSNOW(cl)

  #temp = mapq_files[1:10]

  #t1 = Sys.time()
  #temp1 = mcsapply(temp, FUN = countBam, mc.cores = 10)
  #t2 = Sys.time()
  #print(t2-t1)

  library(ShortRead)

  print('Counting mapq filtered bam files')
  df_mapq_read_counts = mcsapply(mapq_files, FUN = countBam, mc.cores = num_cores)
  mapq_read_counts = unlist(df_mapq_read_counts['records', ] )
  names(mapq_read_counts) = gsub('.mapq_filtered.bam', '', names(mapq_read_counts))


  print('Counting rmdup filtered bam files')
  df_rmdup_read_counts = mcsapply(rmdup_files, FUN = countBam, mc.cores = num_cores)
  rmdup_read_counts = unlist(df_rmdup_read_counts['records', ] )
  names(rmdup_read_counts) = gsub('.rmdup.bam', '', names(rmdup_read_counts))


  print('Counting nonconversion filtered bam files')
  df_nc_filtered_read_counts = mcsapply(nc_filtered_files, FUN = countBam, mc.cores = num_cores)
  nc_filtered_read_counts = unlist(df_nc_filtered_read_counts['records', ] )
  names(nc_filtered_read_counts) = gsub('.rmdup.nonCG_filtered.bam', '', names(nc_filtered_read_counts))


  print('Counting organism bam files')
  df_organism_read_counts = mcsapply(organism_bam_files, FUN = countBam, mc.cores = num_cores)
  organism_read_counts = unlist(df_organism_read_counts['records', ] )
  names(organism_read_counts) = gsub('.organism.bam', '', names(organism_read_counts))
  print('Finished counting bam files')

  all_names = Reduce(union,
                     list(
                       names(mapq_read_counts),
                       names(rmdup_read_counts),
                       names(nc_filtered_read_counts),
                       names(organism_read_counts)
                     )
  )

  df_bam_read_counts <- data.frame(mapq_read_counts = mapq_read_counts[all_names],
                                   rmdup_read_counts = rmdup_read_counts[all_names],
                                   nc_filtered_read_counts = nc_filtered_read_counts[all_names],
                                   organism_read_counts = organism_read_counts[all_names]
  )
  head(df_bam_read_counts)

  rownames(df_bam_read_counts) = gsub('.mapq_filtered.bam', '', rownames(df_bam_read_counts))

  #df_bam_read_counts = df_bam_read_counts[!grepl('sorted', rownames(df_bam_read_counts)), ]

  return(df_bam_read_counts)
}#count_bam_files <- function(alignment_dir)


#bam_file = 'Lane1_ACTTGA.rmdup.nonCG_filtered.bam'
compute_coverage_rates <- function(alignment_dir, parallel = T, log_file)
{
  setwd(alignment_dir)
  print('***********************')
  print('Computing coverage rates')

  df_chrom_sizes = read.table(chrom_sizes_file, header = F)

  df_chrom_ranges = data.frame(chrom = df_chrom_sizes$V1, start = 1, end = df_chrom_sizes$V2)
  total_genome_size = sum(df_chrom_ranges$end)

  chrom_ranges = makeGRangesFromDataFrame(df_chrom_ranges,
                                          keep.extra.columns=FALSE,
                                          ignore.strand=T,
                                          seqinfo=NULL,
                                          starts.in.df.are.0based=FALSE)

  sbp <- ScanBamParam(which=chrom_ranges )

  p_param <- PileupParam(distinguish_nucleotides=FALSE,distinguish_strands=FALSE,
                         min_base_quality=10, min_nucleotide_depth=1)


  bam_files = list.files(alignment_dir,  pattern = '*organism.sorted.bam$')
  #   bam_files = bam_files[1:250]

  base_counts = c()

  count_bam = length(bam_files)

  print(paste('There are ', count_bam, ' files in the alignment directory.')  )

  if(parallel)
  {
    print('Running coverage in parallel')
    cl <- makeCluster(num_cores, outfile=log_file, type = 'SOCK')
    #clusterExport(cl, varlist = ls(), envir = environment())
    #clusterExport(cl, list = ls(), envir = environment())
    #clusterExport(cl, ls(.GlobalEnv))

    registerDoSNOW(cl)
    #clusterExport(cl, varlist = ls(), envir = environment())
    clusterExport(cl, ls(.GlobalEnv))

    print('Starting')

    #base_counts = foreach(i=1:10, .export= ls(globalenv()) ) %dopar%
    base_counts = foreach(i=1:length(bam_files), .export= ls(globalenv()) ) %dopar%
      {

        bam_file = bam_files[i]
        print(paste(i, 'Computing coverage rate for', bam_file))
        #res <- pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)

        library(Rsamtools)

        #res = NA

        res = tryCatch({
          pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)
        }, warning = function(w) {

        }, error = function(e) {
          NA
        }, finally = {

        }
        )

        print(paste(i, 'Computed coverage rate for', bam_file))
        if(is.null(res))
        {
          base_count = NA

        }else
        {
          base_count = nrow(res)

        }

        print(base_count)
        base_count
      }#foreach

  }else{

    base_counts = c()

    for(i in 1:length(bam_files) )
    {
      library(Rsamtools)

      bam_file = bam_files[i]
      print(paste(i, 'Computing coverage rate for', bam_file))
      #res <- pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)

      #res = NA

      res = tryCatch({
        pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)
      }, warning = function(w) {

      }, error = function(e) {
        NA
      }, finally = {

      }
      )

      print(paste(i, 'Computed coverage rate for', bam_file))
      if(is.null(res))
      {
        base_count = NA

      }else
      {
        base_count = nrow(res)

      }

      print(base_count)
      base_counts[i] = base_count

    }#for

  }#else





  cell_ids = bam_files
  cell_ids = sub('.organism.sorted.bam', '', cell_ids)
  cell_ids = sub('.organism.bam', '', cell_ids)
  cell_ids = sub('.bam', '', cell_ids)

  names(base_counts) = cell_ids
  base_counts = unlist(base_counts)

  # for(bam_file in bam_files)
  # {
  #   print(bam_file)
  #   res <- pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)
  #   class(res)
  #   dim(res)
  #   head(res)
  #   tail(res)
  #   cell_id = sub('.rmdup.nonCG_filtered.bam', '', bam_file)
  #   base_counts[cell_id] = nrow(res)
  # }

  coverage_rates = round(base_counts / total_genome_size * 100, 2)
  #names(coverage_rates) = gsub('.rmdup.nonCG_filtered.bam', '', names(coverage_rates))

  #hist(coverage_rates, breaks = 20)
  df_coverage_rates = data.frame(base_count = base_counts, coverage_rate = coverage_rates)
  print('Computed coverage rates.')

  return(df_coverage_rates)
}


merge_r1_and_r2_bam_for_cell <- function(r1_bam, r2_bam, merged_bam)
{
  sys_command = paste0('samtools merge '
                       , ' -h  ', r1_bam
                       #, ' --samtools_path ', samtools_path
                       , ' ', merged_bam
                       , ' ', r1_bam
                       , ' ', r2_bam
  )
  print(sys_command)
  system(sys_command)


  temp_dir = gsub('.bam', '.tmp_sort/', merged_bam)
  sorted_bam = gsub('.bam', '.sorted.bam', merged_bam)
  sys_command = paste0('mkdir -p ', temp_dir)
  print(sys_command)
  system(sys_command)

  sys_command = paste0('samtools sort -o ', sorted_bam, ' -T ', temp_dir, ' ' , merged_bam)
  print(sys_command)
  system(sys_command)

  sys_command = paste0('samtools index ', sorted_bam)
  print(sys_command)
  system(sys_command)

}



merge_r1_and_r2_bam_for_sample <- function(alignment_dir_in,  alignment_dir_out, bam_type = 'organism' )
{

  r1_bam_files = list.files(alignment_dir_in, pattern = 'R1')
  bam_type_str = paste0(bam_type, '.bam$')
  r1_bam_files = r1_bam_files[grepl(bam_type_str, r1_bam_files)]

  if(num_cores > 1)
  {
    cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
    clusterExport(cl, ls(.GlobalEnv))
    registerDoSNOW(cl)
    clusterExport(cl, ls(.GlobalEnv))


    foreach(i=1:length(r1_bam_files)) %dopar%
      {
        r1_bam = r1_bam_files[i]
        r2_bam = paste0(alignment_dir_in,  gsub('R1', 'R2', r1_bam) )
        merged_bam = paste0(alignment_dir_out,  gsub('_R1', '', r1_bam) )
        r1_bam = paste0(alignment_dir_in, r1_bam)
        msg = paste('Merging', r1_bam, 'and', r2_bam)
        print(msg)
        merge_r1_and_r2_bam_for_cell(r1_bam, r2_bam, merged_bam)
      }#foreach

    stopCluster(cl)

  }else
  {

    for(r1_bam in r1_bam_files)
    {

      r2_bam = paste0(alignment_dir_in,  gsub('R1', 'R2', r1_bam) )
      merged_bam = paste0(alignment_dir_out,  gsub('_R1', '', r1_bam) )
      r1_bam = paste0(alignment_dir_in, r1_bam)

      msg = paste('Merging', r1_bam, 'and', r2_bam)
      print(msg)
      merge_r1_and_r2_bam_for_cell(r1_bam, r2_bam, merged_bam)
    }#for

  }####else

}#merge_r1_and_r2_bam_for_sample






