library(doSNOW)

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

      call_methylation_sites_for_cell(alignment_dir = alignment_dir,
                                      methylation_calls_dir = methylation_calls_dir,
                                      cell_id = cell_id,
                                      bme_param_settings = bme_param_settings,
                                      log_dir = log_dir)

    }#foreach

  }#else

  failed_cells = find_failed_met_calls(methylation_calls_dir)
  for(cell_id in failed_cells)
  {
    call_methylation_sites_for_cell(alignment_dir, methylation_calls_dir, cell_id, bme_param_settings, log_dir)
  }


}

find_failed_met_calls <- function(methylation_calls_dir)
{
  cpg_met_files = list.files(methylation_calls_dir, pattern = "CpG_calls.*cov.gz")
  cpg_met_paths = paste0(methylation_calls_dir, cpg_met_files)
  cell_ids = cpg_met_files
  cell_ids = gsub("CpG_calls.", "", cell_ids)
  cell_ids = gsub(".organism.cov.gz", "", cell_ids)
  names(cpg_met_paths) = cell_ids
  file.sizes = file.info(cpg_met_paths)$size
  names(file.sizes) = cell_ids
  zero.size.files = cell_ids[file.sizes == 0]


  setwd(methylation_calls_dir)
  pattern = paste0('.', 'organism', "_splitting_report.txt")
  report_files = list.files(methylation_calls_dir, pattern = paste0("*", pattern) )
  if(length(report_files) == 0 ) {return(NULL) }
  row_names = c()
  result_list = list()
  inf_line_counts = c()
  for(report_file in report_files)
  {
    print(report_file)
    cell_id = gsub(pattern = pattern, '', report_file)
    lines = readLines(report_file)
    informative_lines = lines[grep('\t', lines)]
    informative_lines = informative_lines[!grepl('strand', informative_lines)]
    inf_line_counts[cell_id] = length(informative_lines)
  }

  zero.inf_lines = names(inf_line_counts[inf_line_counts == 0])

  failed_cells = union(zero.size.files, zero.inf_lines)
  return(failed_cells)
}

call_methylation_sites_for_cell <- function(alignment_dir, methylation_calls_dir, cell_id, bme_param_settings, log_dir)
{
  setwd(alignment_dir)

  lambda_file = paste0(cell_id,'.lambda.bam')
  organism_file = paste0(cell_id,'.organism.bam')

  log_sub_dir = paste0(log_dir, '/call_methylation_sites/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)
  log_file = paste0(log_sub_dir, cell_id, '.log')


  sys_command = paste0('rm ', methylation_calls_dir, '/*', cell_id, '*')
  command_result_0 = system(sys_command)

  #Extract methylation sites
  sys_command = paste0(bismark_path, '/bismark_methylation_extractor ',bme_param_settings, ' ',
                       lambda_file,' --output  ', methylation_calls_dir, '/'
                       , ' > ', log_file)
  print(sys_command)
  command_result_1 = system(sys_command)

  sys_command = paste0(bismark_path, '/bismark_methylation_extractor ',bme_param_settings, ' ',
                       organism_file,' --output  ', methylation_calls_dir, '/'
                       , ' >> ', log_file)
  print(sys_command)
  command_result_2 = system(sys_command)

  #Convert to bismark cov format (Only for organism, not for lambda control)
  setwd(methylation_calls_dir)

  sink(log_file, append = T)

  sys_command = paste0(bismark_path, '/bismark2bedGraph -o CpG_calls.',cell_id,'.organism ', 'CpG_context_',cell_id,'.organism.txt.gz')
  command_result_3 = system(sys_command)

  file.rename(paste0('CpG_calls.',cell_id,'.organism.gz.bismark.cov.gz'),  paste0('CpG_calls.',cell_id,'.organism.cov.gz'))

  if(file.exists(paste0('Non_CpG_context_',cell_id,'.organism.txt.gz')))
  {
    sys_command = paste0(bismark_path, '/bismark2bedGraph --CX -o CH_calls.',cell_id,'.organism ', 'Non_CpG_context_',cell_id,'.organism.txt.gz')
    command_result_4 = system(sys_command)
    file.rename(paste0('CH_calls.',cell_id,'.organism.gz.bismark.cov.gz'),  paste0('CH_calls.',cell_id,'.organism.cov.gz'))

  }

  if(file.exists(paste0('CHG_context_',cell_id,'.organism.txt.gz')))
  {
    sys_command = paste0(bismark_path, '/bismark2bedGraph --CX -o CHG_calls.',cell_id,'.organism ', 'CHG_context_',cell_id,'.organism.txt.gz')
    command_result_5 = system(sys_command)

    sys_command = paste0(bismark_path, '/bismark2bedGraph --CX -o CHH_calls.',cell_id,'.organism ', 'CHH_context_',cell_id,'.organism.txt.gz')
    command_result_6 = system(sys_command)

    file.rename(paste0('CHG_calls.',cell_id,'.organism.gz.bismark.cov.gz'),  paste0('CHG_calls.',cell_id,'.organism.cov.gz'))
    file.rename(paste0('CHH_calls.',cell_id,'.organism.gz.bismark.cov.gz'),  paste0('CHH_calls.',cell_id,'.organism.cov.gz'))

  }


  #delete_command = paste0('rm *CpG_context_',cell_id, '.organism.txt.gz')
  #command_result_9 = system(delete_command)

  sink()

  command_result = command_result_1 + command_result_2 + command_result_3


  if(command_result == 0)
  {
    cat('Methylation calling is successful for cell ', cell_id, '\n')
  }else
  {
    stop('Methylation calling failed to run for ', cell_id, '.\n Exiting the pipeline. Please see the output log file :', log_file)
  }

}



process_bismark_split_reports <- function(methylation_calls_dir, genome_type = 'organism')
{
  setwd(methylation_calls_dir)
  pattern = paste0('.', genome_type, "_splitting_report.txt")
  report_files = list.files(methylation_calls_dir, pattern = paste0("*", pattern) )
  if(length(report_files) == 0 ) {return(NULL) }
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


get_met_call_counts <- function(methylation_calls_dir, met_type = 'CpG')
{

  print('--------------------------------------------------------------------')

  cov_files= list.files(methylation_calls_dir, pattern = paste0( met_type, '.*', 'cov.gz'))
  met_call_counts = c()
  counter = 0

  print('***************************************')
  print(methylation_calls_dir)
  print('***************************************')

  #cl <- makeCluster(num_cores, type = 'SOCK')
  #clusterExport(cl, ls(.GlobalEnv))
  #registerDoSNOW(cl)
  #clusterExport(cl, ls(.GlobalEnv))


  #met_call_counts = foreach(i=1:length(cov_files), .export= ls(globalenv()) ) %dopar%
  for(cov_file in cov_files)
  {
    counter = counter + 1
    print(counter)
    print(cov_file)
    cell_id = cov_file
    cell_id = sub('.organism.cov.gz', '', cell_id)
    cell_id = sub(paste0(met_type, '_calls.'), '', cell_id)
    print(cell_id)

    df_cov.dummy = data.frame(chrom = 'chrZz',
                        start = 1,
                        end = 2,
                        met_rate = -.5,
                        met = 1,
                        demet = 1   )

    num_calls = 0

    num_calls = tryCatch({
     nrow(fread(paste0(methylation_calls_dir,  cov_file), tmpdir = methylation_calls_dir ))

    }, warning = function(w) {
      0
    }, error = function(e) {
      0
    }, finally = {

    }
    )



    met_call_counts[cell_id] = num_calls
  }

  return(met_call_counts)
}


plot_site_counts <- function(call_counts, log_scale = F)
{
  #boxplot(list_met_call_counts, outline = F)

  call_count = list_met_call_counts$CpG

  if(log_scale)
  {
    hist(log10(call_count + 1), col = '#bebada',
         main = 'Sites called',
         xlab = 'Sites (log10)',
         ylab = 'Number of cells',
         breaks = 25,
         cex.lab = 1.50, cex.axis = 1.5)

    abline(v = median(log10(call_count + 1)), col = 'red', lty = 2)

  }else
  {
    hist(call_count/1000000, col = '#bebada',
         main = 'Sites called',
         xlab = 'million sites',
         ylab = 'Number of cells',
         breaks = 25,
         cex.lab = 1.50, cex.axis = 1.5)

    abline(v = median(call_count/1000000), col = 'red', lty = 2)

  }




}


plot_split_reports <- function(df_org_split_reports, df_lambda_split_reports, list_org_bias_reports,
                               list_met_call_counts, log_scale = F, lambda_flag = F)
{

  numeric_column_index = 2:ncol(df_org_split_reports)

  lambda_flag = !is.null(df_lambda_split_reports)

  for(i in numeric_column_index)
  {
    if(is.character(df_org_split_reports[,i]) )
    {
      df_org_split_reports[,i] = parse_number(df_org_split_reports[,i])
    }
    if(is.character(df_lambda_split_reports[,i]) & lambda_flag )
    {
      df_lambda_split_reports[,i] = parse_number(df_lambda_split_reports[,i])
    }
  }

  #Conversion rate

  if(lambda_flag)
  {
    conv_CpG = 100 - df_lambda_split_reports$C.methylated.in.CpG.context
    pcts_conv = list(conv_CpG)

  }

  met_CpG = df_org_split_reports$C.methylated.in.CpG.context
  x_labels = c('CpG')
  pcts_met = list(met_CpG)

  if(!is.null(df_org_split_reports$C.methylated.in.non.CpG.context))
  {
    met_CH = 100 - df_org_split_reports$C.methylated.in.non.CpG.context
    x_labels = c(x_labels, 'CH')
    pcts_met[[length(pcts_met)+1]] = met_CH

    if(lambda_flag)
    {
      conv_CH = 100 - df_lambda_split_reports$C.methylated.in.non.CpG.context
      other_conv  = list(conv_CH)
      pcts_conv[[length(pcts_conv)+1]] = conv_CH
    }

  }


  if(!is.null(df_org_split_reports$C.methylated.in.CHG.context))
  {


    met_CHG = df_org_split_reports$C.methylated.in.CHG.context
    met_CHH = df_org_split_reports$C.methylated.in.CHH.context

    pcts_met[[length(pcts_met)+1]] = met_CHG
    pcts_met[[length(pcts_met)+1]] = met_CHH

    x_labels = c(x_labels, 'CHG', 'CHH')


    if(lambda_flag)
    {

      conv_CHG = 100 - df_lambda_split_reports$C.methylated.in.CHG.context
      conv_CHH = 100 - df_lambda_split_reports$C.methylated.in.CHH.context

      pcts_conv[[length(pcts_conv)+1]] = conv_CHG
      pcts_conv[[length(pcts_conv)+1]] = conv_CHH

    }


  }

  #######################Potting#############
  par(mfrow = c(2,3))

  font.face = 'Helvetica'
  font.size = 1.5
  color_vec = hue_pal()(3)
  sep.lwd = 1
  #sep.color = "#33a02c"

  sep.color = "#bebada"

  par(mar = c(4,8,6,4))


  #############Conversion rates################

  if(lambda_flag)
  {

    roof = max(unlist(pcts_conv))
    floor = min(unlist(pcts_conv))

    if(floor > 95) {floor = 95}


    #vioplot(
    boxplot(
      pcts_conv, ylim = c(floor, roof),
      col = "#fb8072", names = x_labels, cex.main = font.size, cex.sub = font.size,
      #main = '% Conversion rate',
      cex.lab = font.size, cex.names = font.size, cex.axis = font.size, cex = font.size,
      ylab = '% conversion', outline = F
      , las = 1
    )

    mtext(line=1, cex.main=font.size, family =font.face, font=2, "Conversion rate")

    box("figure", col=sep.color,  lwd = sep.lwd)


  }


  ##################Bias#####################

  org_bias_reports_CpG = list_org_bias_reports$CpG_met_rates
  org_bias_reports_CH = list_org_bias_reports$CH_met_rates

  class(org_bias_reports_CpG)
  dim(org_bias_reports_CpG)

  par(family = 'Helvetica')

  x = 1:nrow(org_bias_reports_CpG)
  y = rowMeans(org_bias_reports_CpG, na.rm = T)
  roof = max(y) * 1.05
  floot = min(y) * 1.05

  plot(x, y, ylim = c(0,100),
       type = 'l', xlab = 'Position in read', ylab = 'Met. rate (%)', col = 'red',
       main = '',
       #legend = c('CpG', 'CH') ,
       cex.lab = font.size,
       cex.axis = font.size,
       las = 1,


  )

  x = 1:nrow(org_bias_reports_CH)
  y = rowMeans(org_bias_reports_CH)
  roof = max(y) * 1.05
  floot = min(y) * 1.05
  lines(x, y, ylim = c(0,10),
        type = 'l', xlab = 'position', ylab = 'Met. rate (%)', col = 'blue')

  legend('right', legend = c('CpG', 'CH'), col = c('red', 'blue'), lty = 1, bty = "n", cex = font.size)

  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Position bias")

  box("figure", col=sep.color,  lwd = sep.lwd)


  ####Methylation rate##########################3


  roof = max(unlist(pcts_met))

  #vioplot(
  boxplot(
    pcts_met,
    col = "lightgreen", names = x_labels,
    cex.main = font.size,
    #main = '% Methylation rate',
    ylab = '% Methylated C', cex.lab = font.size,
    cex.axis = font.size, cex.names = font.size
    , outline = F
    , las = 1
  )

  mtext(line=1, cex.main=font.size, family =font.face, font=2, "Methylation rate")


  ##########Call Counts###################

  par(mar = c(7,8.5,3,4))

  #c_type = 'CpG'

  c_types = c('CpG', 'CHG', 'CHH')

  for(c_type in c_types)
  {

    call_counts =list_met_call_counts[[c_type]] / 1000000
    dens = density(call_counts)

    plot(dens, main="", xlab = "million sites called", ylab = "", las = 1,
         cex.lab = font.size, cex.axis = font.size)
    polygon(dens, col="#a6cee3", border="#a6cee3")

    title(ylab="Density", line=4, cex.lab=font.size, family=font.face)

    median_call_count = median(call_counts)
    abline(v = median_call_count, col = 'red', lty = 2, lwd = 2)

    legend('topright', legend = c('Med.'),
           col = c('red'), cex = font.size * 0.9,
           lty = 2, bty = "n", inset = c(0, 0))

    box("figure", col=sep.color,  lwd = sep.lwd)


    mtext(line=1, cex.main=font.size, family =font.face, font=2, c_type)


  }



  box("outer", col=sep.color,  lwd = 75)
  box("outer", col="black",  lwd = 3)

  title(paste(sample_name, ': Methylation calling results'), line = -1.5, outer = TRUE,
        cex.main=font.size * 1.2, family=font.face)

}



