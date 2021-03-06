
call_methylation_sites_for_sample <- function(alignment_dir, methylation_calls_dir, log_dir, bme_param_settings)
{
  dir.create(methylation_calls_dir, showWarnings = F, recursive = T)
  setwd(alignment_dir)

  bam_files = list.files(alignment_dir, pattern = "*.organism.bam$")

  if(num_cores > 1)
  {
    cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
    registerDoSNOW(cl)
    foreach(i=1:length(bam_files), .export= ls(globalenv())       ) %dopar%
    {
      bam_file = bam_files[i]
      cell_id = gsub('.organism.bam', '', bam_file)

      call_methylation_sites_for_cell(alignment_dir = alignment_dir,
                                      methylation_calls_dir = methylation_calls_dir,
                                      cell_id = cell_id,
                                      bme_param_settings = bme_param_settings,
                                      log_dir = log_dir)

    }#foreach
  }else#if(num_cores > 1)
  {
    for(i in 1:length(bam_files))
    {
      bam_file = bam_files[i]
      cell_id = gsub('.organism.bam', '', bam_file)

    }#foreach

  }#else



}

call_methylation_sites_for_cell <- function(alignment_dir, methylation_calls_dir, cell_id, bme_param_settings, log_dir)
{
  setwd(alignment_dir)

  lambda_file = paste0(cell_id,'.lambda.bam')
  organism_file = paste0(cell_id,'.organism.bam')

  log_sub_dir = paste0(log_dir, '/call_methylation_sites/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)
  log_file = paste0(log_sub_dir, cell_id, '.log')


  #Extract methylation sites
  sys_command = paste0('bismark_methylation_extractor ',bme_param_settings, ' ',
                       lambda_file,' --output  ', methylation_calls_dir, '/'
                       , ' > ', log_file)
  command_result_1 = system(sys_command)

  sys_command = paste0('bismark_methylation_extractor ',bme_param_settings, ' ',
                       organism_file,' --output  ', methylation_calls_dir, '/'
                       , ' >> ', log_file)
  command_result_2 = system(sys_command)

  #Convert to bismark cov format (Only for organism, not for lambda control)
  setwd(methylation_calls_dir)

  sink(log_file, append = T)

  sys_command = paste0('bismark2bedGraph -o CpG_calls.',cell_id,'.organism ', 'CpG_context_',cell_id,'.organism.txt.gz')
  command_result_3 = system(sys_command)


  if(file.exists(paste0('Non_CpG_context_',cell_id,'.organism.txt.gz')))
  {
    sys_command = paste0('bismark2bedGraph --CX -o CH_calls.',cell_id,'.organism ', 'Non_CpG_context_',cell_id,'.organism.txt.gz')
    command_result_4 = system(sys_command)
  }

  if(file.exists(paste0('CHG_context_',cell_id,'.organism.txt.gz')))
  {
    sys_command = paste0('bismark2bedGraph --CX -o CHG_calls.',cell_id,'.organism ', 'CHG_context_',cell_id,'.organism.txt.gz')
    command_result_5 = system(sys_command)

    sys_command = paste0('bismark2bedGraph --CX -o CHH_calls.',cell_id,'.organism ', 'CHH_context_',cell_id,'.organism.txt.gz')
    command_result_6 = system(sys_command)
  }

  file.rename(paste0('CpG_calls.',cell_id,'.organism.gz.bismark.cov.gz'),  paste0('CpG_calls.',cell_id,'.organism.cov.gz'))
  file.rename(paste0('CH_calls.',cell_id,'.organism.gz.bismark.cov.gz'),  paste0('CH_calls.',cell_id,'.organism.cov.gz'))
  file.rename(paste0('CHG_calls.',cell_id,'.organism.gz.bismark.cov.gz'),  paste0('CHG_calls.',cell_id,'.organism.cov.gz'))
  file.rename(paste0('CHH_calls.',cell_id,'.organism.gz.bismark.cov.gz'),  paste0('CHH_calls.',cell_id,'.organism.cov.gz'))

  delete_command = paste0('rm *CpG_context_',cell_id, '.organism.txt.gz')
  #command_result_9 = system(delete_command)

  sink()

  command_result = command_result_1 + command_result_2 + command_result_3


  if(command_result == 0)
  {
    cat('Methylation calling is successful for cell ', cell_id, '\n')
  }else
  {
    stop('Methylation calling failed to run for ', cell_id, '.\n Exiting the pipeline. Please see the output log.')
  }

}



process_bismark_split_reports <- function(methylation_calls_dir, genome_type = 'organism')
{
  setwd(methylation_calls_dir)
  pattern = paste0('.', genome_type, "_splitting_report.txt")
  report_files = list.files(methylation_calls_dir, pattern = paste0("*", pattern) )
  row_names = c()
  result_list = list()
  for(report_file in report_files)
  {
    print(report_file)
    cell_id = gsub(pattern = pattern, '', report_file)
    lines = readLines(report_file)
    informative_lines = lines[grep('\t', lines)]
    informative_lines = informative_lines[!grepl('strand', informative_lines)]

    df_report = read.delim(text = informative_lines, sep = '\t',  header = F)
    class(df_report)
    row_names = gsub(':', '', df_report$V1)
    result_list[[cell_id]] = as.character(df_report$V2)
  }

  df_result = do.call(cbind, result_list)
  rownames(df_result) = row_names
  t_df_result =  t(df_result)
  final_table = data.frame(Cell_ID = rownames(t_df_result),  t_df_result, stringsAsFactors = F)

  return(final_table)
}



process_bismark_bias_reports <- function(methylation_calls_dir, genome_type = 'organism')
{
  setwd(methylation_calls_dir)
  pattern = paste0('.', genome_type, ".M-bias.txt")
  report_files = list.files(methylation_calls_dir, pattern = paste0("*", pattern) )
  row_names = c()
  result_list = list()
  CpG_met_rates_list = list()
  CH_met_rates_list = list()
  for(report_file in report_files)
  {
    print(report_file)
    cell_id = gsub(pattern = pattern, '', report_file)
    result_list[[cell_id]] = process_bismark_bias_report(cell_id, report_file)
    CpG_met_rates_list[[cell_id]]  = result_list[[cell_id]]$CpG_met_rate
    CH_met_rates_list[[cell_id]]  = result_list[[cell_id]]$CH_met_rate
  }

  CpG_met_rates  = do.call("cbind", CpG_met_rates_list)
  CH_met_rates  = do.call("cbind", CH_met_rates_list)

  return(list(CpG_met_rates = CpG_met_rates, CH_met_rates = CH_met_rates))
}

#report_file = 'Test___Lane1_CAGATC.lambda.M-bias.txt'
#cell_id = 'Test___Lane1_CAGATC'
process_bismark_bias_report <- function(cell_id, report_file)
{

  x <- readLines(report_file)

  start <- grep("C", x)
  mark <- vector('integer', length(x))
  mark[start] <- 1
  # determine limits of each table
  mark <- cumsum(mark)
  # split the data for reading
  df_list <- lapply(split(x, mark), function(.data){
    .input <- read.table(textConnection(.data), skip=2, header=TRUE, sep = '\t')
    attr(.input, 'name') <- .data[1]  # save the name
    .input
  })
  # rename the list
  names(df_list) <- sapply(df_list, attr, 'name')
  df_list
  df_list$`CH context` = df_list$`CHG context` + df_list$`CHH context`
  df_list$`CH context`$X..methylation = 100 * df_list$`CH context`$count.methylated / (df_list$`CH context`$count.methylated + df_list$`CH context`$count.unmethylated)
  df_list$`CH context`$X..methylation = round(df_list$`CH context`$X..methylation, 2)

  result_df = data.frame(position = df_list$`CpG context`$position,
                         CpG_met_rate = df_list$`CpG context`$X..methylation,
                         CH_met_rate = df_list$`CH context`$X..methylation
                         )

  result_df$position = NULL
  return(result_df)
}




plot_split_reports <- function(df_org_split_reports, df_lambda_split_reports, list_org_bias_reports)
{
  par(mfrow = c(3,2))
  #Conversion rate
  conv_CpG = 100 - parse_number(df_lambda_split_reports$C.methylated.in.CpG.context)
  met_CpG = parse_number(df_org_split_reports$C.methylated.in.CpG.context)
  x_labels = c('CpG')
  pcts_conv = list(conv_CpG)
  pcts_met = list(met_CpG)

  if(!is.null(df_lambda_split_reports$C.methylated.in.non.CpG.context))
  {
    conv_CH = 100 - parse_number(df_lambda_split_reports$C.methylated.in.non.CpG.context)
    met_CH = 100 - parse_number(df_org_split_reports$C.methylated.in.non.CpG.context)
    x_labels = c(x_labels, 'CH')
    other_conv  = list(conv_CH)
    pcts_conv[[length(pcts_conv)+1]] = conv_CH
    pcts_met[[length(pcts_met)+1]] = met_CH


  }

  if(!is.null(df_lambda_split_reports$C.methylated.in.CHG.context))
  {
    conv_CHG = 100 - parse_number(df_lambda_split_reports$C.methylated.in.CHG.context)
    conv_CHH = 100 - parse_number(df_lambda_split_reports$C.methylated.in.CHH.context)

    met_CHG = parse_number(df_org_split_reports$C.methylated.in.CHG.context)
    met_CHH = parse_number(df_org_split_reports$C.methylated.in.CHH.context)

    x_labels = c(x_labels, 'CHG', 'CHH')

    pcts_conv[[length(pcts_conv)+1]] = conv_CHG
    pcts_conv[[length(pcts_conv)+1]] = conv_CHH

    pcts_met[[length(pcts_met)+1]] = met_CHG
    pcts_met[[length(pcts_met)+1]] = met_CHH

  }



  roof = max(unlist(pcts_conv))
  floor = min(unlist(pcts_conv))

  if(floor > 99) {floor = 99}

  vioplot(pcts_conv, ylim = c(floor, roof),
          col = "#fb8072", names = x_labels, cex.main = 1.25, cex.sub = 1.5,
          main = '% Conversion rate', cex.lab = 1.50, cex.names = 1.5, cex.axis = 1.5, cex = 1.5,
          ylab = ''
  )

  #Bias

  org_bias_reports_CpG = list_org_bias_reports$CpG_met_rates
  org_bias_reports_CH = list_org_bias_reports$CH_met_rates

  x = 1:nrow(org_bias_reports_CpG)
  y = rowMeans(org_bias_reports_CpG)
  roof = max(y) * 1.05
  floot = min(y) * 1.05

  plot(x, y, ylim = c(0,100),
       type = 'l', xlab = 'position', ylab = 'Met. rate (%)', col = 'red',
       main = 'Position bias', legend = c('CpG', 'CH') , cex.lab = 1.50, cex.axis = 1.5 )

  x = 1:nrow(org_bias_reports_CH)
  y = rowMeans(org_bias_reports_CH)
  roof = max(y) * 1.05
  floot = min(y) * 1.05
  lines(x, y, ylim = c(0,10),
        type = 'l', xlab = 'position', ylab = 'Met. rate (%)', col = 'blue')

  legend('right', legend = c('CpG', 'CH'), col = c('red', 'blue'), lty = 1, bty = "n")


  #Methylation rate
  met_CpG = parse_number(df_org_split_reports$C.methylated.in.CpG.context)
  #met_CH = parse_number(df_org_split_reports$C.methylated.in.non.CpG.context)


  roof = max(unlist(pcts_met))

  vioplot(pcts_met,
          col = "lightgreen", names = x_labels, cex.main = 1.25,
          main = '% Methylation rate',
          ylab = '', cex.lab = 1.50, cex.axis = 1.5, cex.names = 1.5
  )

  #CpG Call Counts
  call_CpG = parse_number(df_org_split_reports$Total.methylated.C.s.in.CpG.context)
             + parse_number(df_org_split_reports$Total.C.to.T.conversions.in.CpG.context)

  hist(call_CpG/1000000, col = '#bebada',
       main = 'CpG sites called',
       xlab = 'million CpG sites',
       ylab = 'Number of cells', cex.lab = 1.50, cex.axis = 1.5)




  title(paste(sample_name, ': Methylation calling results'), line = -1, outer = TRUE)


}



