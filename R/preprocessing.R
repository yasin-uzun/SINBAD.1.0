library(ShortRead)
library(scales)


get_fast_files <- function(fastq_dir, pattern = '')
{
  fastq_files_1 = list.files(fastq_dir, pattern = "\\.fastq.gz$")
  fastq_files_2 = list.files(fastq_dir, pattern = "\\.fastq$")

  fastq_files = union(fastq_files_1, fastq_files_2)
  fastq_files = fastq_files[grepl(pattern = pattern, fastq_files)]

  return(fastq_files)
}

get_r2_indeces_from_r1 <- function(r1_fastq_dir, r2_input_fastq_dir, r2_output_fastq_dir, sample_name = '')
{
  setwd(r1_fastq_dir)
  dir.create(r2_output_fastq_dir, recursive = T)

  r1_fastq_files = get_fast_files(fastq_dir = r1_fastq_dir, pattern = '_R1')
  r1_fastq_files = r1_fastq_files[!grepl(pattern = 'Undetermined', r1_fastq_files)]

  r2_fastq_files = get_fast_files(fastq_dir = r2_input_fastq_dir, pattern = '_R2')
  r2_fastq_files = r2_fastq_files[!grepl(pattern = 'Undetermined', r2_fastq_files)]

  r1_fastq_files = r1_fastq_files[grepl(sample_name, r1_fastq_files)]
  r2_fastq_files = r2_fastq_files[grepl(sample_name, r2_fastq_files)]


  r1_fastq_file_count = length(r1_fastq_files)




  if(num_cores > 1)
  {
    thread_count = min(r1_fastq_file_count, num_cores)
    cl <- makeCluster(thread_count, outfile="", type = 'SOCK')
    registerDoSNOW(cl)

    foreach(i=1:r1_fastq_file_count, .export = ls(globalenv()) ) %dopar%
    {
      r1_fastq_file = r1_fastq_files[i]
      print(r1_fastq_file)

      r2_fastq_file = gsub('_R1', '_R2', r1_fastq_file)

      command = paste0('perl ', perl_index_transfer_path, ' ',
                       ' --input_r1_fastq_file ',  r1_fastq_dir, '/', r1_fastq_file ,
                       ' --input_r2_fastq_file ', r2_input_fastq_dir, '/', r2_fastq_file ,
                       ' --demux_index_length ', demux_index_length,
                       ' --output_r2_fastq_file ', r2_output_fastq_dir, '/', r2_fastq_file )

      print(command)
      system(command)

    }#foreach(i=1:length(raw_fastq_files))

    stopCluster(cl)



  }else
  {

    for(i in 1:r1_fastq_file_count)
    {
      r1_fastq_file = r1_fastq_files[i]
      print(r1_fastq_file)

      r2_fastq_file = gsub('_R1', '_R2', r1_fastq_file)

      command = paste0('perl ', perl_index_transfer_path, ' ',
                       ' --input_r1_fastq_file ',  r1_fastq_dir, '/', r1_fastq_file ,
                       ' --input_r2_fastq_file ', r2_input_fastq_dir, '/', r2_fastq_file ,
                       ' --demux_index_length ', demux_index_length,
                       ' --output_r2_fastq_file ', r2_output_fastq_dir, '/', r2_fastq_file )

      print(command)
      system(command)

    }#foreach(i=1:length(raw_fastq_files))


  }



}#demux_fastq_files

demux_fastq_files <- function(raw_fastq_dir, demux_index_file, demux_index_length, demux_fastq_dir,
                              main_log_dir, read_type = '', sample_name = '')
{
  setwd(raw_fastq_dir)

  fastq_files_1 = list.files(raw_fastq_dir, pattern = "\\.fastq.gz$")
  fastq_files_2 = list.files(raw_fastq_dir, pattern = "\\.fastq$")

  raw_fastq_files = union(fastq_files_1, fastq_files_2)
  raw_fastq_files = raw_fastq_files[!grepl('Undetermined', raw_fastq_files)]

  raw_fastq_files = raw_fastq_files[grepl(sample_name, raw_fastq_files)]
  raw_fastq_files = raw_fastq_files[grepl(read_type, raw_fastq_files)]

  demux_log_dir = paste0(main_log_dir, '/demux/')
  dir.create(demux_log_dir, showWarnings = F, recursive = T)

  if(num_cores > 1)
  {
    raw_fastq_file_count = length(raw_fastq_files)
    thread_count = min(raw_fastq_file_count, num_cores)
    cl <- makeCluster(thread_count, outfile="", type = 'SOCK')
    registerDoSNOW(cl)
    foreach(i=1:raw_fastq_file_count, .export = ls(globalenv()) ) %dopar%
    {
      raw_fastq_file = raw_fastq_files[i]
      print(raw_fastq_file)
      output_prefix = gsub('.fastq.gz', '', raw_fastq_file)
      #ssdemux_log_file = paste0(demux_log_dir, '/', raw_fastq_file, '.log')

      demux_command = paste0('perl ', perl_demux_path, ' ',
                             ' --demux_index_file ',  demux_index_file ,
                             ' --raw_fastq_file ', raw_fastq_dir, '/', raw_fastq_file ,
                             ' --demux_index_length ', demux_index_length,
                             ' --output_dir ', demux_fastq_dir ,
                             ' --output_prefix ', output_prefix,
                             ' --log_dir ', demux_log_dir)

      print(demux_command)
      system(demux_command)

    }#foreach(i=1:length(raw_fastq_files))

    stopCluster(cl)

  }#if(num_cores > 1)



}#demux_fastq_files

read_demux_logs <- function(main_log_dir)
{

  demux_log_dir = paste0(main_log_dir, '/demux/')
  demux_log_files = list.files(demux_log_dir, pattern = "\\.log$")
  demux_log_files = demux_log_files[!grepl('Undetermined', demux_log_files)]
  setwd(demux_log_dir)
  list_df_demux_combined = list()

  for(demux_log_file in demux_log_files)
  {
    print(demux_log_file)
    lane_id = gsub('.log', '', demux_log_file)
    list_df_demux_combined[[lane_id]] = read.table(demux_log_file, header = T)
  }#for

  df_demux_combined = do.call('rbind', list_df_demux_combined)
  df_demux_combined = data.frame(Group = rownames(df_demux_combined), df_demux_combined)
  head(df_demux_combined)

  return(df_demux_combined)
}#read_demux_logs

#source("/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/preprocessing.R")

trim_fastq_files <- function(demux_fastq_dir, trimmed_fastq_dir, main_log_dir)
{

  print('Trimming fastq files')
  setwd(demux_fastq_dir)

  fastq_files_1 = list.files(demux_fastq_dir, pattern = "*.fastq.gz")
  fastq_files_2 = list.files(demux_fastq_dir, pattern = "*.fastq")

  demux_fastq_files = union(fastq_files_1, fastq_files_2)
  demux_fastq_files = demux_fastq_files[!grepl('No_matching_index', demux_fastq_files)]
  print('demux_fastq_files :')
  #print( demux_fastq_files)
  trimming_log_dir = paste0(main_log_dir, '/trimming/')
  dir.create(trimming_log_dir, showWarnings = F, recursive = T)

  if(num_cores > 1)
  {
    cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
    registerDoSNOW(cl)
    clusterExport(cl, ls(), envir = environment())


    if(sequencing_type  == 'paired')
    {
      demux_fastq_files = demux_fastq_files[grepl('R1', demux_fastq_files)]
    }

    print('*******************Trimming....')
    foreach(i=1:length(demux_fastq_files), .export = ls(globalenv())) %dopar%
    {
      fastq_file = demux_fastq_files[i]
      print(paste('Input fastq: ', fastq_file) )
      log_file = paste0(trimming_log_dir, fastq_file, '.log')

      if(trimmer == 'cutadapt')
      {
        if(sequencing_type  == 'paired')
        {
          fastq_file.left = fastq_file
          fastq_file.right = gsub('R1', 'R2', fastq_file)
          command = paste0(cutadapt_path, '/cutadapt ',cutadapt_param_settings,
                           ' -o ', trimmed_fastq_dir, '/', fastq_file.left,
                           ' -p ', trimmed_fastq_dir, '/', fastq_file.right,
                           '  ',  demux_fastq_dir, '/',fastq_file.left ,
                           '  ',  demux_fastq_dir, '/',fastq_file.right ,
                           ' >  ',  log_file
          )

        }else
        {
           command = paste0(cutadapt_path, '/cutadapt ',cutadapt_param_settings,
                             ' -o ', trimmed_fastq_dir, '/', fastq_file,
                             '  ',  demux_fastq_dir, '/',fastq_file ,
                             ' >  ',  log_file
                             )
      }
          print(command)
          system(command)


      }else if(trimmer == 'trim_galore')
      {

        command = paste0(trim_galore_path, '/trim_galore ',trim_galore_param_settings,
                           ' -o ', trimmed_fastq_dir,
                           '  ',  demux_fastq_dir, '/',fastq_file ,
                           ' >  ',  log_file
                          )
        print(command)
        system(command)
        trimmed_fq_name = sub('.fastq.gz', '_trimmed.fq.gz',  fastq_file)
        mv_command = paste0('mv ', trimmed_fastq_dir, trimmed_fq_name, ' ',  trimmed_fastq_dir, fastq_file)
        system(mv_command)

      }else if(trimmer == 'Trimmomatic')
      {

        command = paste0('java -jar ', Trimmomatic_jar_path,
                         ' ', Trimmomatic_param_settings,
                         ' -trimlog ', trimming_log_dir, '/',log_file,
                         ' ', demux_fastq_dir, '/',fastq_file ,
                         ' ',  fastq_file
                         )

        system(command)


      }else
      {
        stop('Invalid  trimmer name, exiting. Trimmer must be one of these: cutadapt trim_galore Trimmomatic')
      }


      #ILLUMINACLIP:${ADAPTER_SEQ}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:25



    }#foreach(i=1:length(raw_fastq_files))

  }#if(num_cores > 1)



}


process_cutadapt_logs <- function(main_log_dir)
{
  cutadapt_log_dir = paste0(log_dir, '/cutadapt/')
  cutadapt_log_files = list.files(cutadapt_log_dir, pattern = "*.log")
  setwd(cutadapt_log_dir)

  row_names = c()
  result_list = list()
  for(report_file in cutadapt_log_files)
  {
    print(report_file)
    cell_id = gsub('.fastq.gz.log', '', report_file)
    lines = readLines(report_file)
    informative_lines = lines[grep(':', lines)]
    informative_lines = informative_lines[!grepl('strand', informative_lines)]
    informative_lines = gsub('\\(', ':', informative_lines)
    informative_lines = gsub('\\)', ':', informative_lines)

    df_report = read.delim(text = informative_lines, sep = '\t',  header = F)
    class(df_report)
    row_names = gsub(':', '', df_report$V1)
    result_list[[cell_id]] = as.character(df_report$V2)
  }

  row_names[1] = "Total_reads"
  row_names[2] = "Alignments"
  row_names[3] = "Alignment_rate"

  row_names[4] = "Sequences_with_no_alignments"
  row_names[6] = "Discarded_sequences"

  df_result = do.call(cbind, result_list)
  rownames(df_result) = row_names
  t_df_result =  t(df_result)
  final_table = data.frame(Cell_ID = rownames(t_df_result),  t_df_result, stringsAsFactors = F)

}

count_fastq_reads <- function(fastq_dir)
{
  print(paste('Reading fastq files from the directory:', fastq_dir))
  fastq_files_1 = list.files(fastq_dir, pattern = "*.fastq.gz")
  fastq_files_2 = list.files(fastq_dir, pattern = "*.fastq")

  fastq_files = union(fastq_files_1, fastq_files_2)
  fastq_files = fastq_files[!grepl('No_matching_index', fastq_files)]

  num_fastq = length(fastq_files)
  print(paste('I found:', num_fastq, ' fastq files: '))
  print(fastq_files)
  print(paste('Now counting reads inside fastqs '))


  cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
  registerDoSNOW(cl)
  clusterExport(cl, ls(), envir = environment())
  #clusterExport(cl, varlist = ls(), envir = environment())

  #read_counts = c()
  line_counts_list = foreach(i=1:length(fastq_files), .export = ls(globalenv())) %dopar%
  {
    library(ShortRead)
    fastq_file = fastq_files[i]
    print(fastq_file)
    countLines(fastq_dir, fastq_file)
  }

  line_counts = unlist(line_counts_list)
  read_counts = line_counts / 4
  names(read_counts) = gsub('.fastq.gz', '', names(read_counts))
  names(read_counts) = gsub('.fastq', '', names(read_counts))

  return(read_counts)

}

plot_preprocessing_results <- function(sample_name, demux_reports, demux_read_counts, trimmed_read_counts)
{
  font.face = 'Helvetica'
  font.size = 1.5
  color_vec = hue_pal()(3)
  sep.lwd = 1

  #sep.color = "#33a02c"

  sep.color = "#fdbf6f"


  #layout(mat = matrix(c(1, 2, 1, 3),
  #                    2, 2, byrow = TRUE))

  #par(oma = c(2,2,2,2))

  par(mfrow = c(2,3))


  ###########DEMUX STATS####################

  par(mar = c(7,7,8,5))

  count_multiplexed_fastq = nrow(demux_reports)
  count_demultiplexed_fastq = length(demux_read_counts)

  Total_reads.sum = sum(demux_reports$Total_reads)
  Reads_with_matching_index.sum = sum(demux_reports$Reads_with_matching_index)
  Reads_without_matching_index.sum = sum(demux_reports$Reads_without_matching_index)

  par(mar = c(2,4,6,2))

  plot.new()
  overall_stat = c()

  overall_stat["Sample_Name"]  = sample_name

  #title("Demultiplexing Statistics:", line=1, cex = font.size, family =font.face, adj=0)
  mtext(side =3 ,line=1, cex.main=font.size, family =font.face, font=2, adj=0, "Demultiplexing Statistics:")

  lbl = paste0("Input count: ",format(count_multiplexed_fastq, scientific = F, big.mark = ","))
  mtext(side =3 , line=-2, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  overall_stat["Input_count"]  = count_multiplexed_fastq

  lbl = paste0("Output count: ",format(count_demultiplexed_fastq, scientific = F, big.mark = ","))
  mtext(side =3 , line=-4, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  overall_stat["Output_count"]  = count_demultiplexed_fastq

  lbl = paste0("Total reads: ",format(Total_reads.sum, scientific = F, big.mark = ","))
  mtext(side =3 , line=-6, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  overall_stat["Total_reads"]  = Total_reads.sum


  lbl = paste0("Matching index: ",format(Reads_with_matching_index.sum, scientific = F, big.mark = ","))
  mtext(side =3 , line=-8, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  overall_stat["Matching_index"]  = Reads_with_matching_index.sum


  lbl = paste0("No matching index: ",format(Reads_without_matching_index.sum, scientific = F, big.mark = ","))
  mtext(side =3 , line=-10, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  overall_stat["No_matching_index"]  = Reads_without_matching_index.sum



  box("figure", col=sep.color,  lwd = sep.lwd)


  ###########DEMUX BOXPLOT####################
  par(mar = c(6,7,6,5))


  list_demux_counts = demux_reports[c('Total_reads', 'Reads_with_matching_index')] /1000000
  roof = max(list_demux_counts) * 1.2

  boxplot(list_demux_counts, names = c('Total', 'With index'), #main = 'Demultiplexing statistics',
           col = color_vec[2], outline = F, las = 1, ylab = 'million reads',
          cex.lab = font.size, cex.axis = font.size, ylim = c(0, roof)
          ,  medlwd = 1)
  #title("Demultiplexing", line=1, cex.lab=font.size, cex.main = font.size,  family=font.face)
  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Demultiplexing")

  box("figure", col=sep.color,  lwd = sep.lwd)



  ###########DEMUX LOSS RATES####################

  par(mar = c(6,6,6,4))


  index_loss_rates = (1 - list_demux_counts$Reads_with_matching_index / list_demux_counts$Total_reads ) * 100

  #barplot(sort(index_loss_rates))
  roof = max(10, max(index_loss_rates * 1.2))
  barplot(sort(index_loss_rates, decreasing = T),
       #xlab = 'Multiplexed input',
       ylab = '% of total reads',  cex.lab = font.size, cex.axis = font.size,
       col = 'brown1',

       las = 2,
       #main = 'Pct of reads without index',
       ylim = c(0, roof)

       #, cex.main = 0.8
  )
  title(xlab="Multiplexed input", line=1, cex.lab=font.size, family=font.face)
  #title("Pct. of reads without index", line=1, cex.lab=font.size, family=font.face)
  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Pct. of reads without index")

  box("figure", col=sep.color,  lwd = sep.lwd)



  ###########TRIM STATS####################



  count_multiplexed_fastq = nrow(demux_reports)
  count_demultiplexed_fastq = length(demux_read_counts)

  Total_reads.trim = sum(trimmed_read_counts)
  Mean_reads.trim = mean(trimmed_read_counts)
  Pct.read.trimmed.out = (demux_read_counts - trimmed_read_counts) / demux_read_counts * 100
  Pct.read.trimmed.out.mean = round(mean(Pct.read.trimmed.out), 1)

  par(mar = c(2,4,6,2))

  plot.new()
  #title("Demultiplexing Statistics:", line=1, cex = font.size, family =font.face, adj=0)
  mtext(side =3 ,line=1, cex.main=font.size, family =font.face, font=2, adj=0, "Trimming Statistics:")

  lbl = paste0("Total trimmed reads: ",format(Total_reads.trim, scientific = F, big.mark = ","))
  mtext(side =3 , line=-2, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  overall_stat["Total_trimmed_reads"]  = Total_reads.trim


  lbl = paste0("Mean trimmed reads: ",format(Mean_reads.trim, scientific = F, big.mark = ","))
  mtext(side =3 , line=-4, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  overall_stat["Mean_trimmed_reads"]  = round(Mean_reads.trim)

  lbl = paste0("Mean pct. of reads filtered out: ",format(Pct.read.trimmed.out.mean, scientific = F, big.mark = ","), "%")
  mtext(side =3 , line=-6, cex.main=font.size,adj=0,
        family=font.face,  font.main =1, lbl)
  box("figure", col=sep.color,  lwd = sep.lwd)
  overall_stat["Pct_removed_reads"]  = Pct.read.trimmed.out.mean



  ###########TRIM LINE PLOT####################

  par(mar = c(8,7,4,4))

  #trimmed_read_counts = trimmed_read_counts[!grepl('R2', names(trimmed_read_counts))]
  #demux_read_counts = demux_read_counts[!grepl('R2', names(demux_read_counts))]


  trim_loss_rates = (1 - trimmed_read_counts / demux_read_counts[names(trimmed_read_counts)] ) * 100

  #barplot(sort(loss_rates))
  plot(sort(trim_loss_rates), type = 'l',
       xlab = 'Demultiplexed input', ylab = '% reads filtered out', col = 'blue',
       #main = 'Adapter trimming'
       #, cex.main = 0.8
       , cex.lab=font.size
       , las =1
       , cex.axis = font.size
  )
  #title("Adapter trimming", line=1, cex.lab=font.size, family=font.face)
  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Adapter trimming")

  box("figure", col=sep.color,  lwd = sep.lwd)




  ###########TRIM BOXPLOT####################

  par(mar = c(8,7,4,5))


  x_labels = c('Demux', 'Trimmed')
  bb = boxplot(
    demux_read_counts/1000000,
    trimmed_read_counts/1000000, plot = F)

  roof = max(bb$stats) * 1.2

  boxplot(
    demux_read_counts/1000000,
    trimmed_read_counts/1000000,
    col = color_vec[3], names = x_labels,
    ylab = 'million reads'
    ,cex.lab = font.size, cex.axis = font.size
    #,main = 'Read statistics'
    , outpch = 19, outcex = 0.2
    , las = 1
    , ylim = c(0, roof)
    #,cex.axis = 1.25, cex.names = 1.25,
    #cex.main = 1.25, cex.lab = 1.25
    , outline =F
    ,  medlwd = 1
  )    #Plot number of reads
  #title("Trimming", line=1, cex.lab=font.size, family=font.face)
  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Trimming")

  box("figure", col=sep.color,  lwd = sep.lwd)



   ##############OUTER BOX########################


  #box("inner", col="black",  lwd = 1)


  box("outer", col=sep.color,  lwd = 75)
  box("outer", col="black",  lwd = 3)

  title(paste(sample_name, ': Preprocessing results'), line = -1.5,
        outer = TRUE, cex.main=font.size * 1.2, family=font.face)

  return(overall_stat)


}


