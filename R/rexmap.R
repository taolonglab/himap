#' Import these functions for everything
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table setkey
#' @importFrom data.table setcolorder
#' @importFrom data.table setorder
#' @importFrom data.table key
#' @importFrom data.table copy
#' @importFrom dada2 filterAndTrim
#' @importFrom stringr str_sub
#' @useDynLib rexmap, .registration = TRUE
NULL

pname = 'RExMap'
pname_l = 'rexmap'
# If you rename the pipeline (again), don't forget to re-run
# Rcpp::compileAttributes(). This will update package name in the compiled
# C++ files as well.


#' Detect operating system
detect_os = function () {
  switch(
    Sys.info()[['sysname']],
      Windows= {'win'},
      Linux  = {'linux'},
      Darwin = {'macos'}
  )
}

#' Auto-select precompiled binaries for the detected OS
exec_file = function (filename) {
  platform = detect_os()
  if (platform == 'macos') return(paste0(filename, '_macos'))
  else if (platform == 'linux') return(paste0(filename, '_linux'))
  else if (platform == 'win') return(paste0(sub('.exe', '', filename),
                                            '_win.exe'))
}


#-- Default options -----------------------------------------------------------
himap_opts = new.env()
allowed_options = c('aln_params', 'blast_max_seqs', 'blast_word_size',
                    'ncpu', 'print_strains_nmax', 'string_maxwidth',
                    'timing', 'verbose')
assign('himap_path', '', env=himap_opts)
# FUll current path is obtained by:
# assign()
# BLAST blast() output format
assign('blast_out_fmt',
       'qseqid sseqid qlen length qstart qend sstart send slen qseq sseq',
       env=himap_opts)
# BLAST collapse() format. Changing this might break collapse function
# used only for debugging.
assign('blast_coll_fmt',
       'qseqid sseqid qlen slen length qstart qend sstart send pident',
       env=himap_opts)
# BLAST executables paths
assign('path_makeblastdb',
       system.file('exec', exec_file('makeblastdb'), package=pname_l),
       env=himap_opts)
assign('path_blastn',
       system.file('exec', exec_file('blastn'), package=pname_l),
       env=himap_opts)
# BLAST alignment parameters
assign('aln_params', c(5L, -4L, -8L, -6L), env=himap_opts)
# Autodetect number of available threads for multithreading parts
# Set to 1 if you always  want to use only 1 thread.
# These (5,-5,-8,-6) are good ones for collapsing sequences using overlap
# alignments, if the lengths of query and subject aren't too different.
assign('osu_offset', 1000000L, env=himap_opts)
assign('ncpu', parallel::detectCores(), env=himap_opts)

# BLAST databases
assign('blast_dbs',
       data.table::fread(
         system.file('extdata', 'pcr_primers_table.txt', package=pname_l)
       ),
       env=himap_opts)
assign('blast_max_seqs', 500, env=himap_opts)
assign('blast_word_size', 50, env=himap_opts)
assign('database_build_version', '2020-01-29', env=himap_opts)

# Read merging
assign('mergepairs_matchqs',
       system.file('merge_tables', 'himap_mergepairs_match_qs.txt',
                   package=pname_l),
       env=himap_opts)
assign('mergepairs_mismatchqs',
       system.file('merge_tables', 'himap_mergepairs_mismatch_qs.txt',
                   package=pname_l),
       env=himap_opts)
# Taxonomy
assign('taxonomy_file',
       system.file('database', 'lineages-2020-02-14_bacteria_archaea.Rdata',
                   package=pname_l),
       env=himap_opts)

# Printing and Data.table adjustments
assign('string_maxwidth', 50, env=himap_opts)
assign('print_strains_nmax', 10, env=himap_opts)

# Display progress and timing of each function that supports these arguments
assign('verbose', FALSE, env=himap_opts)
assign('timing', FALSE, env=himap_opts)




# assign('maxrows', 12, env=himap_opts)

#' RExMap options
#'
#' @param option_names A string or a vector of strings with available options. If
#' not given, then the function lists available options.
#'
#' @export
rexmap_option = function (option_names=NULL) {
  if (is.null(option_names)) {
    cat(paste0(pname, ': available options'), fill=T)
    return(sapply(ls(himap_opts), rexmap_option))
  }
  if(!all(option_names %in% ls(himap_opts))) {
    warning("Invalid  option: ", option_names[!(option_names %in% ls(himap_opts))])
    option_names = option_names[option_names %in% ls(himap_opts)]
  }
  if (length(option_names) == 0) stop("Invalid options.")
  # get(option_names, env=himap_opts)
  if (length(option_names) == 1) get(option_names, env=himap_opts)
  else sapply(option_names, get, env=himap_opts)
}

reload_blast_dbs = function () {
  assign('blast_dbs',
         data.table::fread(
           system.file('extdata', 'pcr_primers_table.txt', package=pname_l)
         ),
         env=himap_opts)
}


# rexmap_option = rexmap_option


#' RExMap set default options
#'
#' @param option_name Name of the option to be changed.
#' @param value New value.
#'
#' Options that can be changed:
#'
#'
#' @export
rexmap_setoption = function (option_name, value) {
  # Simply set option_name to value. Used to change RExMap defaults. Not finished yet.
  # Allowed options
  if (option_name == 'ncpu') {
    # Check that it is an integer
    if (!(class(value) == 'integer')) stop('ncpu must be an integer.')
  } else if (option_name == 'string_maxwidth') {
    # Check that value is integer
    options('datatable.prettyprint.char'=value)
    assign('string_maxwidth', value, env=himap_opts)
  } else if (option_name %in% c('verbose', 'timing')) {
    if (!(value %in% c(TRUE, FALSE))) stop(option_name, ' can only be TRUE or FALSE.')
    else assign(option_name, value, env=himap_opts)
  } else if (option_name == 'aln_params') {
    # Check each possible combination for megablast
  } else {
    # Else just assign stuff
    assign(option_name, value, env=himap_opts)
  }
}

himap_setoption = rexmap_setoption


check_database_files = function () {
  # Check if all files are present for each hypervariable region in the PCR
  # primers table.
  pcr.dt = rexmap_option('blast_dbs')
  if (nrow(pcr.dt) == 0) return(NA)

  local_db_path = file.path(find.package(pname_l), 'inst', 'database')
  local_db_files = dir(local_db_path, 'V[0-9].*\\.(?:txt|nin|nhr|nsq)',
                       full.names=TRUE)
  local_db_filenames = basename(local_db_files)
  local_db_filenames_tables = local_db_filenames[local_db_filenames %like% '\\.txt$']
  local_db_filenames_dbnhr = local_db_filenames[local_db_filenames %like% '\\.nhr$']
  local_db_filenames_dbnsq = local_db_filenames[local_db_filenames %like% '\\.nsq$']
  local_db_filenames_dbnin = local_db_filenames[local_db_filenames %like% '\\.nin$']

  # For each PCR primer pair, check 1) copy number table, 2) blast db files
  pcr_valid.dt = pcr.dt[
    , .(valid = if (.BY[[2]] %in% local_db_filenames_tables &
                    paste0(.BY[[1]], '.nhr') %in% local_db_filenames_dbnhr &
                    paste0(.BY[[1]], '.nsq') %in% local_db_filenames_dbnsq &
                    paste0(.BY[[1]], '.nin') %in% local_db_filenames_dbnin) TRUE
                else FALSE)
    , by=.(DB, table)]

  if (!all(pcr_valid.dt[, valid])) {
    return(FALSE)
  } else {
    return(TRUE)
  }

}

#' Run this when the package is loaded and attached in R
.onAttach = function (libname, pkgname) {

  start_time = Sys.time()
  options('datatable.prettyprint.char'=rexmap_option('string_maxwidth'))

  # Check if the blastn and makeblastdb executable files are missing for the running
  # system platform. If they are missing, download them.
  blastn_file = system.file('exec', exec_file('blastn'), package=pname_l)
  makedb_file = system.file('exec', exec_file('makeblastdb'), package=pname_l)
  if (blastn_file=='' | makedb_file=='') {
    download_blast()
  }

  # Check if database files are missing. Need at least one set of blastdb
  # files and 1 table with matching primers...
  db_check = check_database_files()
  if (is.na(db_check) | !db_check) {
    cat('RExMap: errors in database.')
    download_database()
  }

  # Generate a unique list of all available hypervariable regions
  hregions = sub('_$', '', sub('[0-9]{4}-[0-9]{2}-[0-9]{2}$', '',
                 rexmap_option('blast_dbs')[, Hypervariable_region]))
  hregions = sub('_mux', '*', hregions, fixed=T)
  hregions = hregions[hregions != '']

  hregions_str = paste(unique(sort(hregions)), collapse=', ')
  last_date = rexmap_option('blast_dbs')[, max(date)]

  # startup_message = paste0(
  #    '| RExMap',
  #   # ' | Updated ', rexmap_option('database_build_version'),
  #   ' | ', last_date,
  #   ' | ',
  #   hregions_str,
  #   # paste(unique(sort(rexmap_option('blast_dbs')[
  #   #   , sub('_$', '', sub('^(V[0-9]+(\\-|_)(V[0-9]+)?).*$', '\\1', Hypervariable_region))]
  #   # )), collapse=', '),
  #   ' | Threads: ', rexmap_option('ncpu'), ' |'
  # )
  # # Positions of separators
  # startup_msg_separators = as.integer(gregexpr('|', startup_message, fixed=T)[[1]])
  # startup_msg_topbotframe = paste(
  #   mapply(
  #     function (x0, x1) paste0('+', paste(rep('-', x1-x0-1), collapse='')),
  #     startup_msg_separators[1:(length(startup_msg_separators)-1)],
  #     startup_msg_separators[2:length(startup_msg_separators)]
  #   ),
  # collapse='')
  # # Put together final startup string with line endings
  # startup_message_full = paste0(
  #   startup_msg_topbotframe, '+\n',
  #   startup_message, '\n',
  #   startup_msg_topbotframe, '+'
  # )

  startup_message_full = paste0(pname, ' | ', length(hregions),
                                ' regions | DB: ', last_date)

  end_time = Sys.time()
  dt = end_time - start_time

  packageStartupMessage(startup_message_full)
  # print(dt)
}

#' Show RExMap database version based on the last modification date
#'
#' Version == last modified date of the *_unique_variants_R.txt table.
#' @export
rexmap_db_version = function (region=NULL) {
  if (is.null(region)) {
    return(as.Date(file.info(system.file(
      'database', rexmap_option('blast_dbs')[, table], package=pname_l))$mtime))
  } else {
    return(as.Date(file.info(system.file(
      'database', rexmap_option('blast_dbs')[Hypervariable_region==region, table],
      package=pname_l))$mtime))
  }
}
himap_db_version = rexmap_db_version


#' Frequency table of sequence lengths
#'
#' Generate a frequency table of sequence lengths from FASTQ files.
#' The resulting table object can be visuaized with base plot function.
#'
#' @param fastq_files A character vector of FASTQ filenames.
#'
#' @export
sequence_length_table = function (fastq_files, ncpu=rexmap_option('ncpu')) {
  return(
    table(unlist(parallel::mclapply(
      fastq_files,
      function (f) nchar(sfastq_reader(f)$seqs),
      mc.cores=ncpu
    )))
  )
}

#' Quantile from frequency table
#'
#' Outputs a quantile corresponding to a given probability \code{prob}
#' from a frequency table \code{ft}. We use this to find the minimum sequence length
#' above which we have certain fraction of reads.
#'
#' Similar to the base R code \code{\link{quantile}}, but works with
#' frequency table instead of the raw list.
#'
#' @param ft Frequency table of sequence lengths. Output from \code{link{sequence_length_table}}.
#' @param prob Probability.
#'
#' @examples
#' # Find minimum sequence length, above which we have 99% reads
#' # (this example requires valid FASTQ files in fq_tri variable)
#' seqlen.ft = sequence_length_table(fq_tri)
#' ftquantile(seqlen.ft, 0.01)
#'
#' @export
ftquantile = function (ft, prob) {
  # Return a quantile from a frequency table ft at given probability prob.
  ft_relsum = cumsum(ft)/sum(ft)
  return(as.integer(names(ft_relsum[ft_relsum >= prob])[1]))
}


#' DADA2 denoising
#'
#' Uses DADA2 algorithm to denoise trimmed FASTQ files \code{fastq_trimmed}.
#' DADA2 partitions sequences based on their quality score profiles, using a
#' \code{pvalue_adjusted} parameter (\code{OMEGA_A} in the paper) as a cutoff
#' for forming new partitions. (See DADA2 paper
#' [Callahan et al. Nature Methods 2016] for details.)
#' If \code{pvalue}
#' is given (default), then the \code{pvalue_threshold} is calculated by dividing
#' the \code{pvalue} by the maximum number of unique sequences from all samples.
#' Otherwise \code{pvalue_threshold} is used directly.
#'
#' After denoising, the untrimmed parts of sequences are retrieved from
#' \code{fastq_untrimmed} files and for each partition a consensus sequence
#' of that part is concatenaded to its end, omitting everything after
#' (and including) any indel and N.
#'
#' @param pvalue P-value threshold before Bonferroni correction for DADA2.
#' @param pvalue_adjusted P-value threshold after Bonferroni correction for DADA2.
#' If left as NULL, adjusted pvalue is calculated by taking \code{pvalue} and
#' dividing it by the maximum number of unique sequences across all samples.
#' @param multithread TRUE/FALSE or the number of CPU threads to use for multithreading.
#' Does not work on Windows due to \code{parallel} package implementation.
#' @param verbose TRUE/FALSE: display of status messages.
#' @param dada_errors_path Path to the Rdata object containing previously calcualted
#' DADA2 error matrix
#'
#' @export
# @param use_intermediate Use saved intermediate files. Used for resuming very
# long runs.
dada_denoise = function (fastq_trimmed, fastq_untrimmed,
                         pvalue=1e-4,
                         pvalue_adjusted=NULL,
                         multithread=rexmap_option('ncpu'),
                         verbose=rexmap_option('verbose'),
                         timing=rexmap_option('timing'),
                         # error_estimation_nsamples=3,
                         dada_errors_path=NULL) {
  # Dereplicate reads into a derep object
  # Check whether intermediate folder exists
  # Load 1 sample at a time, process it, then combine

  # Time execution
  start_time = Sys.time()

  # Learn errors
  #if (verbose) cat('* learn errors')
  m('  Learn errors...', fill=F, time_stamp=T, verbose=verbose)
  #dada_derep = dada2::derepFastq(
  #  fastq_trimmed[1:min(length(fastq_trimmed), error_estimation_nsamples)])
  #if (class(dada_derep) != 'list') dada_derep = list(dada_derep)
  #if (verbose) cat('...')

  # Get or calculate the Bonferroni adjusted p-value
  if (is.null(pvalue_adjusted)) {
    # Find the maximum number of unique sequences across all samples
    max_num_uniques = max(unlist(
      parallel::mclapply(
        fastq_trimmed,
        function (f) {
          l = length(suppressMessages(suppressWarnings(dada2::derepFastq(f)$uniques)))
          return(l)
        }, mc.cores=multithread)
    ))
    pvalue_adjusted_calc = pvalue/max_num_uniques^2
  } else {
    pvalue_adjusted_calc = pvalue_adjusted
  }

  # Learn errors first, then use those downstream
  if (is.null(dada_errors_path)) {
    dada_errors = suppressWarnings(dada2::learnErrors(
      fastq_trimmed, multithread=multithread, OMEGA_A=pvalue_adjusted_calc
    ))
  } else {
    # cat(' loading pre-existing error file... ')
    m(' (loaded exiting)', fill=F, time_stamp=F, verbose=verbose)
    dada_errors = readRDS(dada_errors_path)
  }
  # if (verbose) cat(' OK.', fill=T)
  m(' OK.', fill=T, time_stamp=F, verbose=verbose)
  dada_results = list()
  dada_empty_result = NA

  # Keep track of

  # Iterate over each sample sequentially
  for (i in 1:length(fastq_trimmed)) {
    start_time_i = Sys.time()
    fq = fastq_trimmed[i]
    fqt = fastq_untrimmed[i]
    # if (verbose) cat('* processing ', basename(fq), fill=T)
    m('  -', basename(fq), fill=F, time_stamp=T, verbose=verbose)

    # Processed i/length(fastq_trimmed) files,
    # pct_processed = round(100*i/length(fastq_trimmed), 1)

    # Try loading the file
    dada_derep = tryCatch(
      list(dada2::derepFastq(fq)),
      error = function (e) NA)
    if (class(dada_derep) == 'logical') {
      if (is.na(dada_derep)) {
        m(' Error in derepFastq. Skipping.', fill=T, time_stamp=F,
          verbose=verbose)
        next
      }
    }
    pvalue_adjusted_calc = ie(is.null(pvalue_adjusted),
                              pvalue/length(dada_derep[[1]]$uniques)^2,
                              pvalue_adjusted)
    # if (is.null(pvalue_adjusted)) {
    #   # Adjust the p-value by an estimate for the total number of comparisons in dada2
    #   pvalue_adjusted_calc = pvalue/length(dada_derep[[1]]$uniques)^2
    # } else {
    #   pvalue_adjusted_calc = pvalue_adjusted
    # }
    dada_res = tryCatch(
      list(suppressMessages(suppressWarnings(
        dada2::dada(dada_derep, err=dada_errors, OMEGA_A=pvalue_adjusted_calc,
                    multithread=multithread)))),
      error = function (e) NA
    )
    if (class(dada_res) == 'logical') {
      if (is.na(dada_res)) {
        m(' Error in dada. Skipping.', fill=T, time_stamp=F, verbose=verbose)
        next
      }
    }


    # Extract length to which reads have been truncated
    trunclen = nchar(dada_res[[1]]$sequence[1])
    #dada_res = add_consensus(dada_res, dada_derep, fqt, fq, truncLen=trunclen,
    #                          ncpu=max(1, as.integer(multithread)), verbose=verbose)
    dada_res = untrim(dada_res, fq, fqt, ncpu=multithread, verbose=verbose)
    # Add the dada2 denoised object of this sample to the list
    dada_results[[length(dada_results)+1]] = dada_res[[1]]
    names(dada_results)[length(dada_results)] = basename(fq)
    end_time_i = Sys.time()
    dt_i = end_time_i - start_time_i
    m(' OK. [', round(dt_i, 1), ' ', attr(dt_i, 'units'), ']', fill=T,
      time_stamp=F, verbose=verbose)
  }

  end_time = Sys.time()
  dt = end_time - start_time
  m('  Done. [', round(dt, 1), ' ', attr(dt, 'units'), ']', fill=T,
    time_stamp=T, verbose=verbose)
  # diff_time = difftime(end_time, start_time, units='secs')
  #   cat('Finished in ', round(as.numeric(diff_time)/60), ' m ',
  #     round(as.numeric(diff_time)%%60, 1), ' s.\n')
  # }

  return(dada_results)
}


#' Untrim trimmed sequences in a DADA2 object
#'
#' @param dada_res DADA2 object with sequences to un-trim.
#' @param derep DADA2 derep object that was used to obtain dada_res.
#' @param fq_tri A character vector of FASTQ filenames pre-global-trimming
#' (but after PCR primer trim).
#' @param fq_fil A character vector of FASTQ filenames post-global-trimming and filtering.
#' @param truncLen Length used to trim.
#' @param verbose Boolean specifying whether to display progress bar.
#' @param ncpu Integer specifying number of CPU threads to use. This uses R package "parallel"
#' so works only on macOS and Linux.
#'
add_consensus = function (dada_res, derep, fq_tri, fq_fil, truncLen,
                          verbose=T, ncpu=rexmap_option('ncpu')) {
  if (verbose) cat('Retrieving full-length sequences...\n')
  for (s_id in 1:length(dada_res)) { # For each sample s_id
    if (verbose) cat('Sample ', s_id, '. Load...')
    x = partid_to_fastqid(dada_res[[s_id]]$map-1, derep[[s_id]]$map-1)
    # Load both merged (for retrieval) and filtered sequences for referencing
    # Filtered sequences are enumerated in dada2, so this file is used to
    # extract meta-data information for each read that we need to look up
    # in the merged file.
    fq_fil_data = sfastq_reader(fq_fil[s_id])
    fq_mer_data = sfastq_reader(fq_tri[s_id])
    if (verbose) cat('OK. Consensus...')
    pid_to_seq = unlist(parallel::mclapply(1:length(x), function (i) {
      metas = fq_fil_data[['meta']][x[i][[1]]+1]
      mer_seqs = fq_mer_data[['seqs']][which(fq_mer_data[['meta']] %in% metas)]
      # Select only the right hanging part
      mer_seqs_right = str_sub(mer_seqs, start=truncLen+1)
      mer_seqs_cons = gsub('[-|N]{1,}.*', '', consensus_sequence(mer_seqs_right))
      names(mer_seqs_cons) = as.character(i)
      return(mer_seqs_cons)
    }, mc.cores=ncpu))
    pid_to_seq = pid_to_seq[order(as.integer(names(pid_to_seq)))]
    if (verbose) cat('OK. Update...')
    # Concatenate dada2 middle partition sequence and right consensus
    dada_res[[s_id]]$sequence = paste(dada_res[[s_id]]$sequence,
                                     pid_to_seq,
                                     sep='')
    names(dada_res[[s_id]]$denoised) = dada_res[[s_id]]$sequence
    if (verbose) cat('OK.\n')
  }
  return(dada_res)
}


#' Un-trim trimmed sequences
#'
#' Use all pre-trimmed sequences to find the consensus of the trimmed part.
#'
#' This function modifies dada result \code{dada_res}.
#'
untrim = function (dada_res, fq_trimmed, fq_pretrimmed,
                   verbose=rexmap_option('verbose'),
                   ncpu=rexmap_option('ncpu')) {
  trim_len = nchar(dada_res[[1]]$sequence[1])
  if (verbose) cat('Trimmed length: ', trim_len, ' nt.', fill=T)
  if (verbose) cat('Retrieving full-length sequences...\n')
  for (s_id in 1:length(dada_res)) { # For each sample s_id
    seqs = dada_res[[s_id]]$sequence
    num_partitions = length(seqs)
    if (verbose) cat('Sample ', s_id, '. Load...')
    # x = partid_to_fastqid(dada_res[[s_id]]$map-1, derep[[s_id]]$map-1)
    # Load both merged (for retrieval) and filtered sequences for referencing
    # Filtered sequences are enumerated in dada2, so this file is used to
    # extract meta-data information for each read that we need to look up
    # in the merged file.
    fq_trimmed_data = sfastq_reader(fq_trimmed[s_id])
    fq_pretrimmed_data = sfastq_reader(fq_pretrimmed[s_id])
    if (verbose) cat('OK. Consensus...')
    part_to_seq = unlist(parallel::mclapply(1:length(seqs), function (i) {
      # Retrieve meta-data for each exact trimmed seq that maps to the same unique
      metas = fq_trimmed_data[['meta']][fq_trimmed_data[['seqs']]==seqs[i]]
      # Now get the untrimmed sequences with these IDs, they will differ only untrim
      # part
      untrim_seqs = fq_pretrimmed_data[['seqs']][fq_pretrimmed_data[['meta']] %in% metas]
      untrim_seqs_right = str_sub(untrim_seqs, start=trim_len+1)
      untrim_seqs_right_cons = gsub('[-|N]{1,}.*', '',
                                    consensus_sequence(untrim_seqs_right))
      # Retrieve the untrimmed extesions
      names(untrim_seqs_right_cons) = as.character(i)
      return(untrim_seqs_right_cons)
    }, mc.cores=ncpu))

    part_to_seq = part_to_seq[order(as.integer(names(part_to_seq)))]
    if (verbose) cat('OK. Update...')
    # Concatenate dada2 middle partition sequence and right consensus
    dada_res[[s_id]]$sequence = paste(dada_res[[s_id]]$sequence,
                                     part_to_seq,
                                     sep='')
    names(dada_res[[s_id]]$denoised) = dada_res[[s_id]]$sequence
    if (verbose) cat('OK.\n')
  }
  return(dada_res)
}


#' Convert abundance matrix to an abundance data table
#'
#' Each row in abundance matrix is a different sample, each column is a different
#' sequence.
#'
ab_mat_to_dt = function (ab_tab_nochim, fq_prefix_split='_',
                         fq_prefix_split_n=1) {
   ab_tab_nochim.dt = as.data.table(unname(ab_tab_nochim))
   ab_tab_nochim.dt[, sample_id := sapply(strsplit(
     dimnames(ab_tab_nochim)[[1]], fq_prefix_split, fixed=T
     ), `[`, fq_prefix_split_n)]
   ab_tab_nochim_m.dt = data.table::melt(ab_tab_nochim.dt, id.vars = 'sample_id',
                               variable.name = 'dada2_seqid', value.name = 'raw_count')
   ab_tab_nochim_m.dt[, qseqid := as.integer(gsub('V', '', dada2_seqid))]
   ab_tab_nochim_m.dt[, dada2_seqid := NULL]
   ab_tab_nochim_m.dt = merge(
     ab_tab_nochim_m.dt,
     data.table(qseqid=1:ncol(ab_tab_nochim), sequence=dimnames(ab_tab_nochim)[[2]])
   )
   setcolorder(ab_tab_nochim_m.dt, c('sample_id', 'qseqid', 'raw_count', 'sequence'))
   return(ab_tab_nochim_m.dt[])
}


pctsim_range = function (p) return(max(p, na.rm=T))


#' OSU abundance table
#'
#' Combine BLAST object and a sequence abundance table into a final OSU
#' abundance table.
#'
#' @param abundance_table Sequence abundance table. Data table with 3 columns
#' (in this exact order): sample_id, qseqid, raw_count.
#' @param blast_object BLAST output class, output from \code{\link{blast}} function.
#' @param verbose TRUE/FALSE: display status messages
#' @param raw_strains TRUE/FALSE: whether to report full strain information for each OSU.
#' If this is true "species" column will display a full list of best matching strains
#' for each OSU. Since this last may be very long, use \code{\link{print_strains}}
#' to make it more readable.
#' @param ncpu Integer specifying number of CPU threads to use. This uses R package "parallel"
#' (TRUE) or to simplify the output when multiple strains of the same species are in the
#' same OSU (FALSE).
#' @param custom_sampleids process only specified sampleids (separated by comma)
#' If NULL (default) process all sample ids.
#' @param pso_n Integer specifying number of times to run Particle Swarm
#' Optimizer for an OSU abundance estimation algorithm. (Default: 1000)
#'
#' @importFrom igraph make_empty_graph
#' @importFrom igraph add_vertices
#' @importFrom igraph add_edges
#' @importFrom igraph as.undirected
#' @importFrom igraph groups
#' @importFrom igraph clusters
#' @importFrom limSolve lsei
#' @importFrom pso psoptim
#' @importFrom data.table dcast
#'
#' @export
abundance = function (abundance_table, blast_object,
                      ncpu=rexmap_option('ncpu'),
                      ncpu_sample=1,
                      raw_strains=TRUE,
                      pso_n=1000,
                      custom_sampleids=NULL,
                      verbose=rexmap_option('verbose'),
                      debug=FALSE,
                      timing=FALSE,
                      lm_solver='lsei') {

  # Check that abundance table has only 1 qseqid per sample_id; this issue
  # may occur in case of incorrectly extracting sample_ids from dada2
  # object (check the separator).
  # if (nrow(abundance_table[, if (.N > 1) .SD, by=.(qseqid, sample_id)]) > 0) {
  #   stop('Error: Multiple sequence (qseqid) counts for identical sample_ids.
  #         Check that each file map to a unique sample_id in the sequence
  #         abundance table.')
  # }

  debug_print = function (...) {
    if (debug) {
      cat(...)
    }
  }

  # Initial timestamp
  t0 = Sys.time()
  # Generate OSU data table first
  debug_print('* Generating strain spectra (osu_data_m.dt)...')
  osu_data_m.dt = blast_cp_to_osu_dt(
    blast_best.dt=blast_object$alignments,
    cp.dt=blast_object$cp,
    ab_tab_nochim_m.dt=abundance_table[, 1:3],
    # ab_tab_nochim_m.dt=abundance_table,
    ncpu=max(ncpu, ncpu_sample),
    verbose=verbose
  )
  t1 = Sys.time()
  t1mt0 = t1-t0
  if (debug) {
    cat(' OK. [', round(t1mt0, 1), ' ', attr(t1mt0, 'units'), ']\n', sep='')
    print(head(osu_data_m.dt))
  }
  # Now generate osu abundance table
  osu_ab.dt = osu_cp_to_all_abs(abundance_table[, 1:3],
                                # abundance_table,
                                blast_object$alignments,
                                blast_object$cp,
                                osu_data_m.dt,
                                ncpu=ncpu,
                                ncpu_sample=ncpu_sample,
                                verbose=verbose,
                                pso_n=pso_n,
                                raw=raw_strains,
                                custom_sampleids=custom_sampleids,
                                debug=debug,
                                lm_solver=lm_solver)

  setcolorder(osu_ab.dt, c('sample_id', 'osu_id', 'osu_count', 'pctsim', 'species'))
  setorder(osu_ab.dt, sample_id, -osu_count)
  end_time = Sys.time()
  dt = end_time - t0
  m(' [', round(dt, 1), ' ', attr(dt, 'units'), ']', fill=T, time_stamp=F, verbose=verbose)
  if (verbose) cat('Abundance estimation complete.\n')
  return(osu_ab.dt)
}

#'
#'
#' @importFrom igraph make_empty_graph
#' @importFrom igraph add_vertices
#' @importFrom igraph add_edges
#' @importFrom igraph as.undirected
#' @importFrom igraph groups
#' @importFrom igraph clusters
#' @importFrom limSolve lsei
#' @importFrom pso psoptim
#' @importFrom data.table dcast
#'
osu_cp_to_all_abs = function (ab_tab_nochim_m.dt,
                              blast_best.dt,
                              cp.dt,
                              osu_data_m.dt,
                              ncpu=rexmap_option('ncpu'),
                              ncpu_sample=1,
                              verbose=T,
                              timing=F,
                              pso_n=1000,
                              raw=TRUE, debug=FALSE,
                              custom_sampleids=NULL,
                              lm_solver='lsei') {

  # lm_solver either 'lsei' or 'quadprogpp'

  debug_print = function (...) {
    if (debug) {
      cat(...)
    }
  }
  debug_print_head = function (dt) {
    if (debug) {
      print(head(dt))
    }
  }

  t0 = Sys.time()
  pctsim_min = 100
  osu_offset = rexmap_option('osu_offset')


  m('  Preparing blast tables...', fill=F, time_stamp=T, verbose=verbose)
  var_count.dt = unique(osu_data_m.dt[, .(variant_id, raw_count, sample_id)])
  common_variant_ids = var_count.dt[
    , if (all(raw_count > 0)) .SD, by=variant_id][, unique(variant_id)]
  blast_best2.dt = blast_best.dt[
    pctsim < pctsim_min, .(pctsim=pctsim_range(pctsim),
                           species=print_strains(strain, raw=raw)), by=qseqid]
  t1 = Sys.time()
  t1mt0 = t1-t0
  m(paste0('OK. [', round(t1mt0, 1), ' ', attr(t1mt0, 'units'), ']'),
    fill=T, time_stamp=F, verbose=verbose)

  if (debug) {
    # cat('DEBUG: var_count.dt \n')
    # print(head(var_count.dt))

    # cat('DEBUG: blast_best2.dt \n')
    # print(head(blast_best2.dt))
  }

  # Define optimization function
  H = function (x) as.numeric(x>0)
  f = function (x, A, B, x0) {
    # Return an optimization function
    eps = A %*% x - B
    (t(eps) %*% eps) + ((x0^2) %*% H(x))
  }
  if (debug) {
    m('  osu_data_m.dt:', fill=T, time_stamp=T, verbose=verbose)
    print(head(osu_data_m.dt))
    m('  osu_data_m.dt[raw_count > 0]:', fill=T, time_stamp=T, verbose=verbose)
    print(head(osu_data_m.dt[raw_count > 0]))
  }
  m('  Preparing initial OSU tables...', fill=F, time_stamp=T, verbose=verbose)
  if (nrow(osu_data_m.dt[raw_count > 0]) > 0) {
    osu_data_single.dt = osu_data_m.dt[
      , if (length(unique(variant_id)) == 1) .SD, by=osu_id]
    if (nrow(osu_data_single.dt) > 0) {
      osu_data_m_single.dt = unique(
        osu_data_single.dt[raw_count > 0, .(osu_id, variant_id, copy_number)]
      )
    }
  }
  t2 = Sys.time()

  if (debug) {
    # cat('DEBUG: osu_data_m_single.dt \n')
    # print(head(osu_data_m_single.dt))
  }
  osu_sp.dt = cp.dt[, .(species = print_strains(strain, raw=raw)), by=osu_id]
  if (debug) {
    # cat('DEBUG: osu_sp.dt \n')
    # print(head(osu_sp.dt))
  }
  if (!is.null(custom_sampleids)) {
    sample_ids = strsplit(custom_sampleids, split=',', fixed=T)[[1]]
  } else {
    sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]
  }
  t3 = Sys.time()
  t3mt1 = t3 - t1
  m(paste0('OK. [', round(t3mt1, 1), ' ', attr(t3mt1, 'units'), ']'),
    fill=T, time_stamp=F, verbose=verbose)
  m(paste0('  Calculating abundances (', ncpu, ' threads)'),
    fill=T, time_stamp=T, verbose=verbose)

  all_abs.dt = data.table::rbindlist(parallel::mclapply(sample_ids, function (s) {

    all_ab.dt = tryCatch({
      # if (debug) print()
      # Numbers of rows for optimized and non-optimized OSUs
      n_opt = 0
      n_non = 0
      t00 = Sys.time()
      debug_print('sample_id: ', s, '\n')

      # First prepare < 100% matches
      debug_print('  Preparing matrices for 100% matches...')
      if (nrow(osu_data_m.dt[sample_id==s & raw_count > 0]) > 0) {

        Ab.dt = dcast(
           osu_data_m.dt[sample_id==s], variant_id ~ osu_id,
           value.var='copy_number', fill=0
        )
        Ab.dt = merge(
           Ab.dt,
           unique(osu_data_m.dt[sample_id==s, .(variant_id, raw_count)]),
           by='variant_id'
        )

        A = as.matrix(Ab.dt[, 2:(ncol(Ab.dt)-1)])
        dimnames(A)[[1]] = Ab.dt[, variant_id]

        # Generate a matrix with B coefficients
        B = as.matrix(Ab.dt[, raw_count])
        dimnames(B)[[1]] = Ab.dt[, variant_id]

        t01 = Sys.time()
        t01mt00 = t01 - t00
        debug_print(' OK. [', round(t01mt00, 1), ' ', attr(t01mt00, 'units'), ']\n',
                    sep='')

        # Solve
        debug_print('  Solving linear model:\n')
        if (debug) {
          saveRDS(A, paste0(s, '_matrix_A.Rdata'))
          saveRDS(B, paste0(s, '_matrix_B.Rdata'))
          saveRDS(Ab.dt, paste0(s, '_Ab.dt.Rdata'))
          saveRDS(osu_data_m.dt, paste0(s, '_osu_data_m.dt.Rdata'))
        }
        if (lm_solver == 'lsei') {
          debug_print('  - Using limSolve::lsei() as least square solver.\n')
          sol = tryCatch(
            lsei(A, B, fulloutput=T, G=diag(ncol(A)), H=matrix(c(0), nrow=ncol(A),
                                                               ncol=1), type=2),
            error = function (x) NA
          )
        } else if (lm_solver == 'quadprogpp') {
          debug_print('  - Using quadprogpp as least square solver.\n')
          t01a = Sys.time()
          debug_print('  - Transforming matrices... ')
          dvec = crossprod(A, B)
          Dmat = crossprod(A, A)
          diag(Dmat) = diag(Dmat) + 1e-8
          Amat = t(diag(ncol(A)))
          bvec = matrix(c(0), nrow=ncol(A), ncol=1)
          t01b = Sys.time()
          t01bmt01a = t01b - t01a
          debug_print(' OK. [', round(t01bmt01a), ' ', attr(t01bmt01a, 'units'), ']\n')
          sol_vector = tryCatch(
            quadprogpp::quadprog.solve.QP(Dmat, dvec, Amat, bvec),
            error = function (x) NA
          )
          sol = list(X=sol_vector$solution)
          names(sol$X) = colnames(A)
        }
        t02 = Sys.time()
        t02mt01 = t02 - t01
        debug_print('  Done. [', round(t02mt01, 1), ' ', attr(t02mt01, 'units'), ']\n',
                    sep='')


        if (class(sol) == 'logical') {
          debug_print('  Linear model failed. [class(sol) == logical]')
          if (is.na(sol)) {
            # Least square model failed for some reason; in which case just
            # return an empty data table. Likely problematic inpjut sequene counts.
            return(data.table())
          }
        }

        debug_print('  Generating graph from the Ar matrix...')

        osu_th = 1e-1
        osu_ab = sol$X
        osu_count = as.integer(osu_ab[which(osu_ab > osu_th)])
        osu_ids   = as.integer(names(osu_ab[which(osu_ab > osu_th)]))
        osu_ab.dt = data.table(osu_id=osu_ids, osu_count=osu_count)
        data.table::setorder(osu_ab.dt, osu_id)

        # Try optimizing the full table first?
        x0 = osu_ab.dt[, osu_count]
        A_columns = as.character(osu_ab.dt[, osu_id])
        Ar = A[, A_columns, drop=F]
        Ar = Ar[setdiff(names(rowSums(Ar) > 0), names(B[B==0,])), , drop=F]
        Br = B[rownames(Ar), , drop=F]

        # Generate a graph from Ar matrix, then identify all
        # connected clusters.
        g = make_empty_graph()
        A_rows = rownames(Ar)
        # Add vertices like this
        #    1  2  3
        # 4  .  .  .
        # 5  .  .  .
        # 6  .  .  .
        #  nrow=3, ncol=4
        g = add_vertices(g, ncol(Ar))
        g = add_vertices(g, nrow(Ar))
        vxs = c(colnames(Ar), rownames(Ar))
        for (j in 1:ncol(Ar)) {
           g = add_vertices(g, 1)
        }
        for (i in seq(ncol(Ar)+1, ncol(Ar)+1+nrow(Ar))) {
           g = add_vertices(g, 1)
        }
        for (i in 1:nrow(Ar)) {
           for (j in 1:ncol(Ar)) {
              if (Ar[i,j] > 0) {
                 # cat('add: ', i+ncol(Ar), '-', j, '\n')
                 g = add_edges(g, c(i+ncol(Ar), j))
                 g = add_edges(g, c(j, i+ncol(Ar)))
              }
           }
        }
        g = as.undirected(g)
        t03 = Sys.time()
        t03mt02 = t03 - t02
        debug_print(' OK. [', round(t03mt02, 1), ' ', attr(t03mt02, 'units'), ']\n',
                    sep='')


        debug_print('  Finding connected clusters...')
        # Find all connected clusters
        cls = parallel::mclapply(
          groups(clusters(g)),
          function (x) if (length(x[x<=ncol(Ar)]) > 1) x else NA,
          mc.cores=ncpu_sample)
        cls = cls[!is.na(cls)]
        osu_ab2.dt = copy(osu_ab.dt)
        t04 = Sys.time()
        t04mt03 = t04 - t03
        debug_print(' OK. ', length(cls), ' clusters found. [',
                    round(t04mt03, 1), ' ', attr(t04mt03, 'units'), ']\n',
                    sep='')

        # For each cluster i
        if (length(cls) > 0) {
          debug_print('  Solving individual clusters\n')
          # for (i in 1:length(cls)) {
          osu_ab3_list = parallel::mclapply(1:length(cls), function (i) {
            # Add back any osu with a single variant_id
            # Sometimes, optimization will omit osu_ids with mapping to single
            # variant_ids so we bring those back here manually.
            debug_print('  cluster #', i, '\n')
            cl = cls[i][[1]]
            osu_ids = vxs[cl[cl<=ncol(Ar)]]
            variant_ids = vxs[cl[cl>ncol(Ar)]]
            debug_print('  - osu_ids:', paste(osu_ids, collapse=','), '\n')
            debug_print('  - variant_ids: ', paste(variant_ids, collapse=','), '\n')
            # osu_ids = c(osu_ids, osu_data_m_single.dt[variant_id %in% variant_ids][, osu_id])
            # variant_ids = c(variant_ids, osu_data_m_single.dt[variant_id %in% variant_ids][, osu_id])
            Ar2 = Ar[variant_ids, osu_ids, drop=F]
            Br2 = Br[variant_ids, drop=F]

            if (nrow(osu_data_m_single.dt[variant_id %in% variant_ids]) > 0) {
              # If we have any that need to be added
              x = as.matrix(dcast(osu_data_m_single.dt[variant_id %in% variant_ids],
                                  variant_id ~ as.character(osu_id), value.var='copy_number')[, -1])
              rownames(x) = osu_data_m_single.dt[variant_id %in% variant_ids, variant_id]
              x[is.na(x)] = 0
              Ar3.dt = merge(as.data.table(Ar2, keep.rownames=T),
                             as.data.table(x, keep.rownames=T), all.x=T)
              Ar3.dt[is.na(Ar3.dt)] = 0
              Ar3 = as.matrix(Ar3.dt[, -1])
              Ar3 = Ar3[, order(as.integer(colnames(Ar3)))]
              rownames(Ar3) = Ar3.dt[, rn]
            } else {
              # If not just use the Ar2 matrix
              Ar3 = Ar2
            }

            Br3 = as.matrix(B[rownames(Ar3),])
            debug_print('  - osu_ab3.dt creation:\n')
            debug_print_head(osu_ab.dt)

            osu_ab3.dt = merge(osu_ab.dt,
                               data.table(osu_id=as.integer(colnames(Ar3))),
                               by='osu_id', all.y=T)
            osu_ab3.dt[is.na(osu_count), osu_count := 0L]
            debug_print_head(osu_ab3.dt)
            # setorder(osu_ab3.dt, osu_id)

            x0 = osu_ab3.dt[, osu_count]

            x0_w = rep(0.3, length(x0))
            # x0_w = x0

            ar3 = copy(Ar3)
            for (r in 1:nrow(Ar3)) ar3[r,] = ar3[r,] / Br3[r]
            debug_print('  - pso running...')
            tmp = lapply(1:pso_n, function (it) {
              pso = psoptim(
                rep(0, length(x0)),
                f, lower=rep(0, length(x0)), upper=rep(max(Br3), length(x0)),
                A=ar3, B=rep(1, nrow(ar3)),
                x0 = as.numeric(x0_w),
                control=list('maxit'=10,
                             'vectorize'=T
                             # 'hybrid'=F, 'rand.order'=F
                ))
              list(round(pso$par, 2),
                   pso$value)
            })
            res = tmp[which.min(sapply(tmp, function (x) x[[2]]))]
            debug_print(' OK.\n')
            osu_ab3.dt[, osu_count := as.integer(res[[1]][[1]])]
            # osu_ab2.dt = merge(
            #   osu_ab2.dt[!(osu_id %in% colnames(Ar3))],
            #   osu_ab3.dt,
            #   all=T,
            #   by=names(osu_ab3.dt)
            # )

            # rm(Ar2, ar3, tmp, osu_ab3.dt, Ar3.dt)
            if (exists('Ar2')) rm(Ar2)
            if (exists('ar3')) rm(ar3)
            if (exists('tmp')) rm(tmp)
            # if (exists('osu_ab3.dt')) rm(osu_ab3.dt)
            if (exists('Ar3.dt')) rm(Ar3.dt)

            # New for mclapply
            return(list(osu_ab3.dt, colnames(Ar3)))

          }, mc.cores=ncpu_sample)
          t05 = Sys.time()
          t05mt04 = t05 - t04
          debug_print(' OK. [', round(t05mt04, 1), ' ', attr(t05mt04, 'units'), ']\n',
                      sep='')

          debug_print('  Merging tables...')
          for (osu_ab3_item in osu_ab3_list) {
            osu_ab2.dt = merge(
              osu_ab2.dt[!(osu_id %in% osu_ab3_item[[2]])],
              osu_ab3_item[[1]],
              all=T,
              by=names(osu_ab3_item[[1]])
            )
          }
          t06 = Sys.time()
          t06mt05 = t06 - t05
          debug_print(' OK. [', round(t06mt05, 1), ' ', attr(t06mt05, 'units'), ']\n',
                      sep='')


        } # end if length(cls) > 0

        if (!exists('t05')) {
          t05 = Sys.time()
        }

        # rm(g, cls, Ar, Br, sol, A, B)
        if (exists('g')) rm(g)
        if (exists('cls')) rm(cls)
        if (exists('Ar')) rm(Ar)
        if (exists('Br')) rm(Br)
        if (exists('sol')) rm(sol)
        if (exists('A')) rm(A)
        if (exists('B')) rm(B)

        debug_print('  Finalizing < 100% OSU table...')
        osu_ab2.dt = osu_ab2.dt[osu_count > 0]

        # Add OSU labels (species_label from osu_sp.dt)
        osu_ab4.dt = merge(osu_ab2.dt, osu_sp.dt, by='osu_id',
                          all.x=T)
        data.table::setorder(osu_ab4.dt, -osu_count)
        osu_ab5.dt = merge(osu_ab4.dt,
                          unique(osu_data_m.dt[, .(osu_id, pctsim=pctsim_range(pctsim))]),
                          by='osu_id')

        # Cleanup other variables
        # rm(osu_ab.dt, osu_ab2.dt, osu_ab4.dt, Ab.dt)
        if (exists('osu_ab.dt')) rm(osu_ab.dt)
        if (exists('osu_ab2.dt')) rm(osu_ab2.dt)
        if (exists('osu_ab4.dt')) rm(osu_ab4.dt)
        if (exists('Ab.dt')) rm(Ab.dt)
        n_opt = nrow(osu_ab5.dt)

        if (exists('t06')) {
          t07mt06 = t06 - t04
          debug_print(' OK. [', round(t05mt04, 1), ' ', attr(t05mt04, 'units'), ']\n',
                      sep='')

        } else {
          t07mt06 = t05 - t04
          debug_print(' OK. [', round(t05mt04, 1), ' ', attr(t05mt04, 'units'), ']\n',
                      sep='')
        }



      } else {
        debug_print(' No < 100% matches found.\n')
      }
      t08 = Sys.time()
      debug_print('  Processing < 100% matches...')
      # Check < 100% matches
      if (nrow(blast_best2.dt) > 0) {
        ab_tab2.dt = merge(
          ab_tab_nochim_m.dt[sample_id==s],
          blast_best2.dt,
          by='qseqid'
        )
        ab_tab2.dt[, osu_id := osu_offset+qseqid]
        ab_tab2.dt[, osu_count := raw_count]
        ab_tab2.dt[, c('raw_count', 'qseqid') := NULL]
        n_non = nrow(ab_tab2.dt)
      }

      if (n_opt > 0 & n_non > 0) {
        # Merge OSU analysis with low sim sequences
        # all_ab.dt = merge(osu_ab5.dt,
        #                   ab_tab2.dt,
        #                   by=intersect(names(osu_ab5.dt), names(ab_tab2.dt)),
        #                   all=T)
        all_ab.dt = data.table::rbindlist(list(
          osu_ab5.dt[, .(osu_id, osu_count, species, pctsim)],
          ab_tab2.dt[, .(osu_id, osu_count, species, pctsim)]
        ))
      } else if (n_opt == 0 & n_non > 0) {
        all_ab.dt = ab_tab2.dt
      } else if (n_opt > 0 & n_non == 0) {
        all_ab.dt = osu_ab5.dt
      } else {
        # Both optimized and non-optimized sequences are zero
        # Return empty data.table
        return(data.table())
      }
      all_ab.dt[, sample_id := s]
      if (debug) {
        print(head(all_ab.dt))
      }

      data.table::setcolorder(all_ab.dt, c('sample_id', 'osu_id', 'osu_count', 'species',
                               'pctsim'))
      # if (verbose) cat('- sample_id:', s, ' completed.\n')
      t09 = Sys.time()
      t09mt08 = t09 - t08
      debug_print('OK. [', round(t09mt08, 1), ' ', attr(t09mt08, 'units'), ']\n',
                  sep='')
      dt_total = t09 - t00
      m('  Sample', s, 'OK. [', round(dt_total, 1), ' ', attr(dt_total, 'units'),
        ']', fill=T, time_stamp=T, verbose=verbose)

      # Return value
      all_ab.dt[osu_count>0]
    }, error = function (e) {
      m('  Sample', s, 'Failed at OSU estimation.', fill=T, time_stamp=T, verbose=verbose)
      data.table(sample_id=NULL, osu_id=NULL, osu_count=NULL, species=NULL, pctsim=NULL)
    }
  )

  return(all_ab.dt)

  }, mc.cores=ncpu))

  m('  Done.', fill=T, time_stamp=T, verbose=verbose)
  return(unique(all_abs.dt))
}

#' Sequence abundance table
#'
#' Return abundances of each sequence from a DADA2 result object.
#'
#' @param dada_result dada() result
#' @param remove_bimeras Check and remove bimeric reads? (default: TRUE)
#' @param collapse_sequences Should we check for sequences that differ up to
#' shifts and add up the counts? (default: TRUE)
#' @param remove_bimeras_method Method to remove bimeras. See ?dada2::removeBimeraDenovo
#' for more options.
#'
#' @export
sequence_abundance = function (dada_result, remove_bimeras=T, collapse_sequences=T,
                               verbose=T, remove_bimeras_method='consensus',
                               remove_bimeras_oneoff=T, fq_prefix_split='.',
                               fq_prefix_split_n=1,
                               ncpu=TRUE, debug=F) {
  # Extract sequences and their counts from the dada class result
  # if (verbose) cat('* generating sequence table...')
  start_time = Sys.time()
  m('* DadaToTab.', fill=F, time_stamp=T, verbose=verbose)
  ab_tab = dada2::makeSequenceTable(dada_result)
  # if (verbose) cat(' OK.', fill=T)
  if (remove_bimeras) {
    # if (verbose) cat('* removing bimeras...')
    ab_tab_nochim = dada2::removeBimeraDenovo(
      ab_tab, method=remove_bimeras_method, allowOneOff=remove_bimeras_oneoff,
      # multithread=ie(rexmap_option('ncpu') > 1, T, F),
      multithread=ncpu,
      verbose=verbose)
    m(' DeBimered.', fill=F, time_stamp=F, verbose=verbose)
    # if (verbose) cat(' OK.', fill=T)
  } else {
    ab_tab_nochim = ab_tab
  }
  if (collapse_sequences) {
    # if (verbose) cat('* adding together sequences that differ in shifts on lengths...')
    ab_tab_nochim_coll = collapse(ab_tab_nochim, verbose=F)
    # if (verbose) cat(' OK.', fill=T)
    m(' Collapsed.', fill=F, time_stamp=F, verbose=debug)
  }
  else {
    ab_tab_nochim_coll = ab_tab_nochim
  }
  seq.dt = ab_mat_to_dt(ab_tab_nochim_coll, fq_prefix_split=fq_prefix_split,
                        fq_prefix_split_n=fq_prefix_split_n)
  m(' Tabulated.', fill=F, time_stamp=F, verbose=verbose)
  end_time = Sys.time()
  dt = end_time - start_time
  m(' [', round(dt, 1), ' ', attr(dt, 'units'), ']', fill=T, time_stamp=F,
    verbose=verbose)
  return(seq.dt)
}

#' Generate a table with sequences assigned to each OSU
#'
#' @param osu_abundances Data table with OSU abundances; output from the
#' \code{\link{abundance}} function.
#' @param blast_output Blast output class. Output from \code{\link{blast}}
#' function.
#'
#' @export
osu_sequences = function (osu_abundances, blast_output) {
  # First do 100% OSUs. For these combine blast output tables
  osu_var.dt = merge(
    unique(osu_abundances[, .(osu_id, species, pctsim)]),
    blast_output$sequences,
    by='osu_id'
  )
  osu_var.dt[, c('spectrum', 'variant_id') := NULL]
  osu_var.dt = merge(
    osu_var.dt,
    unique(blast_output$alignments[, .(qseqid, pctsim)]),
    all.x=T,
    by=c('qseqid', 'pctsim')
  )
  return(osu_var.dt[order(osu_id), .(osu_id, species, qseqid, copy_number, pctsim, sequence)])
}

#' Generate a table with a full list of all best-matched strains
#'
#' @param osu_abundances Data table with OSU abundances; output from the
#' \code{\link{abundance}} function.
#' @param blast_output Blast output class. Output from \code{\link{blast}}
#' function.
#'
#' @export
osu_matches = function (blast_output,
                        osu_offset=rexmap_option('osu_offset')) {

  # 100% hits are in the cp table
  osu_strain_exact.dt = unique(blast_output$cp[, .(osu_id, strain)])
  osu_strain_exact.dt[, pctsim := 100.00]

  # OSU offset is by default 1000000, i.e. the first OSU which is non-exact
  # hit.
  nonexact_osus_ref.dt = unique(
    blast_output$sequences[osu_id >= osu_offset, .(qseqid, osu_id)]
  )
  nonexact_seqs_ref.dt = unique(
    blast_output$alignments[, .(qseqid, strain, pctsim)]
  )
  osu_strain_nonexact.dt = unique(merge(
    nonexact_osus_ref.dt,
    nonexact_seqs_ref.dt,
    by='qseqid'
  )[, .(osu_id, strain, pctsim)])

  # Combine
  osu_strain.dt = data.table::rbindlist(list(
    osu_strain_exact.dt,
    osu_strain_nonexact.dt
  ))
  # Non-100% hits are in the
  return(osu_strain.dt[order(osu_id, strain)])
}




