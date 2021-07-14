#' Auto-detect PCR primers in a merged FASTQ file
#'
#' Samples `nseqs` from an input FASTQ file with merged forward and reverse
#' reads `fq`, then aligns all combinations of primer pairs from RExMapDB.
#' These can be seen using `rexmap_option('blast_dbs')`. Based on these alignments
#' and the parameter `min_seqs_pct`, it will try to guess the possible hypervariable
#' region by matching a forward and reverse primer.
#'
#' Returns a most likely hypervariable region as a character string.
#'
#' @param fq Path to the FASTQ file with merged reads to analyze.
#' @param pr_fwd_maxoff Maximum offset for the FORWARD primer from the start of
#' the read.
#' @param pr_rev_maxoff Maximum offset for the REVERSE primer from the start of
#' the read.
#' @param max_mismatch  Maximum number of nucleotide mismatch for the primer
#' alignment.
#' @param nseqs Number of sequences to randomly sample from the FASTQ file.
#' @param nseqs_seed Random seed for the sequence sampling.
#' @param min_seqs_pct Minimum percent of reads that need to align to the primer
#' sequence, so that we would consider this primer as a potential PCR primer.
#' @param ncpu Number of primer sequences to run in parallel.
#' @param ncpu_seqs Number read data sequences to run in parallel. Total number
#' of CPU threads willbe ncpu * ncpu_seqs.
#' @param verbose Verbose output.
#' @param debug Print debugging information to stdout.
#'
#' @export
detect_pcr_primers = function (fq,
                               # pr_fwd=NULL, pr_rev=NULL,
                               pr_fwd_maxoff=35, pr_rev_maxoff=35, max_mismatch=3,
                               nseqs=1000, nseqs_seed=42, min_seqs_pct=50,
                               ncpu=4, ncpu_seqs=1,
                               verbose=T, debug=F
                               ) {
   # Detect PCR primers in a single or multiple FASTQ files, given as an input
   # character vector 'fq'
   # - nseqs: Randomly sample nseqs from the input fastq files
   if (debug) {
      cat('------- Detect PCR primers -------', fill=T)
   }

   #---------------------- Check FASTQ files ---------------------------------
   # Check if all files in fq exist
   fq_exist = file.exists(fq)
   fq_good = fq[fq_exist]
   if (sum(!fq_exist) > 0) {
      if (debug) {
         cat('* Warning: Files: ', paste(fq[!fq_exist], collapse=','),
             ' do not exist.')
      }
      if (length(fq_good) == 0) {
         stop('No files can be processed.')
      }
   }

   #------------------- Prepare primer table ---------------------------------
   # if (is.null(pr_fwd) & is.null(pr_rev)) {
     primers.dt = rexmap_option('blast_dbs')[
       , .(Primer1, Primer2, `Primer1_sequence_5to3`, `Primer2_sequence_3to5`,
           Primer1_reg=sub('^(V[0-9])[_]?(V[0-9])?.*$', '\\1', Hypervariable_region),
           Primer2_reg=sub('^(V[0-9])[-]?(V[0-9])?.*$', '\\2', Hypervariable_region))]
     primers.dt[Primer2_reg=='', Primer2_reg := Primer1_reg]
     primers.dt = unique(primers.dt[Primer1 != ''])
     primers.dt[, Primer1_sequence_3to5 := sapply(Primer1_sequence_5to3, reverse_complement)]
     primers.dt[, Primer2_sequence_5to3 := sapply(Primer2_sequence_3to5, reverse_complement)]

     pr_m.dt = data.table::melt(
       primers.dt,
       measure.vars=list(
         c('Primer1', 'Primer2'),
         c('Primer1_sequence_5to3', 'Primer2_sequence_5to3'),
         c('Primer1_sequence_3to5', 'Primer2_sequence_3to5'),
         c('Primer1_reg', 'Primer2_reg'))
     )
     pr_m.dt[, variable := NULL]
     names(pr_m.dt) = c('pr_name', 'pr_5to3', 'pr_3to5', 'pr_lab')
     pr_m.dt[, pr_id := paste(pr_lab, sub('^[^RF]+([RF])$', '\\1', pr_name), sep='')]
     pr.dt = unique(data.table::melt(pr_m.dt, id.vars='pr_id', measure.vars=c('pr_5to3', 'pr_3to5'),
                         variable.name='pr_dir', value.name='seq_ext'))
     pr.dt[, seq_id := paste(pr_id, sub('pr_', '', pr_dir), sep='_')]
   # } else {
   #   pr.dt = data.table(
   #      pr_id=c('fwd_primer', 'rev_primer'),
   #      pr_dir=c('pr_5to3', 'pr_5to3'),
   #      seq_ext=c(pr_fwd, pr_rev),
   #      seq_id=c('fwd_primer_5to3', 'rev_primer_5to3')
   #   )
   # }


   # Load the files 'fq_good'
   # if (verbose) {
   #   t0 = Sys.time()
   # }
   fq_list = lapply(fq_good, fastq_reader)
   # if (verbose) {
   #   t1 = Sys.time
   #   t1mt0 = t1-t0
   # }
   set.seed(nseqs_seed)
   seqs_all = unlist(lapply(fq_list, function (x) x$seqs))
   nseqs_max = min(nseqs, length(seqs_all))
   seqs_sampled = sample(seqs_all, nseqs_max)
   if (debug) {
     cat('* Sampled: ', nseqs_max, 'sequences.', fill=T)
   }

   # Search of any primer, 5'->3' or 3'->5' directions
   primer_alignments.dt = data.table::rbindlist(parallel::mcmapply(
     function (seq_ext, seq_id, pr_id) {
       # Convert seq_ext -> seq_ext_n (replace all extd nt codes w/ N )
       # Align seq_ext_n vs each of the sequences seqs_sampled
       seq_ext_n = gsub('[^ACGT]', 'N', seq_ext)

       # Run alignments in C
       alns = parallel::mcmapply(function (read_seq) {

         # Run local sequence alignment between primer sequence (seq_ext_n)
         # and a read sequence read_seq
         aln = C_nwalign(seq_ext_n, read_seq, match=1, mismatch=-1, indel=-1)
         # This will result in a character vector of length 2 showing the
         # alignment, e.g.
         # ----------------ACACAACA----------------------   (primer)
         # ACTGACTAGCTACGTAACACAACATGCTACGACTGATCGACTACGT   (read)

         # aln_left_offset measures how much primer hangs over the beginning
         # of the read, e.g.
         # ACACAACA----------------------
         # ---CAACATGCTACGACTGATCGACTACGT
         # aln_left_offset = 4 (read begins at position 4 of the alignment)
         aln_left_offset = regexpr('[^-]', aln[2])[1]

         # pr_left and pr_right are the coordinates of the aligned part of the
         # sequences (without trailing -----s) in the **alignment** coordinate
         # system.
         # Primer left alignment offset. In the original case this is position
         # 17 (how far primer is into the read)
         pr_left = max(regexpr('[^-]', aln[1])[1], aln_left_offset)
         pr_right = regexpr('[-]{1,}$', aln[1])[1]-1

         # Number of Ns?
         pr_n = lengths(regmatches(seq_ext_n, gregexpr('N', seq_ext_n)))

         # Compare only the aligned part
         aln_stat = compare_alignment(
           str_sub(aln[1], start=pr_left, end=pr_right),
           str_sub(aln[2], start=pr_left, end=pr_right)
         )
         # Forward primer alignment is acceptable if it's near beginning (within first 5 nts)
         # and if it doesn't have more than 2 mismatches for non-N symbols, which includes indels.
         # pr_fwd_found = FALSE
         # pr_rev_found = FALSE
         pr_pot_dir = 'NA'
         if (pr_left < pr_fwd_maxoff &
             aln_stat[2]+aln_stat[3]+aln_stat[4]-pr_n <= max_mismatch) {
           # aln_stat[2]+aln_stat[3]+aln_stat[4]-pr_n = number of mismatches +
           # gap openings + gap extensions in the alignment
           # pr_fwd_found = TRUE
           # Potential direction for this primer is FORWARD (its near the begin)
           pr_pot_dir = 'F'
         }
         else if (pr_left > nchar(read_seq)-nchar(seq_ext_n)-pr_rev_maxoff &
             aln_stat[2]+aln_stat[3]+aln_stat[4]-pr_n <= max_mismatch) {
           # Reverse primer found
           pr_pot_dir = 'R'
         }

         # Trim if needed
         # start = 1L
         # end = -1L
         # if (pr_pot_dir != 'NA') {
         #   start = pr_right - aln_left_offset + 2
         #   end = pr_left - 1
         # } else {
         #
         # }
         return(pr_pot_dir)
       }, seqs_sampled, mc.cores=ncpu_seqs, USE.NAMES=F)
       # Generate statistics
       alns_stat_pct = round(
         100*table(factor(alns, levels=c('F', 'R', 'NA')))/nseqs_max,
         1)

       seq.dt = data.table(
         pr_id=pr_id, seq_id=seq_id, seq_ext=seq_ext,
         fwd_pct=alns_stat_pct['F'],
         rev_pct=alns_stat_pct['R'],
         na_pct=alns_stat_pct['NA']
       )
       return(seq.dt)

     }, pr.dt[, seq_ext], pr.dt[, seq_id], pr.dt[, pr_id],
     mc.cores=ncpu, USE.NAMES=F, SIMPLIFY=F
   ))

   # Use this statistics to find primer pair where either fwd or rev pct
   # matching sequence is greater than e.g. 5%
   primer_alignments.dt[, fwd_pct_bin := cut(fwd_pct, seq(0, 100, 10), right=T,
                                             include.lowest=T)]
   primer_alignments.dt[, rev_pct_bin := cut(rev_pct, seq(0, 100, 10), right=T,
                                             include.lowest=T)]
   good_alignments.dt = primer_alignments.dt[
     fwd_pct > min_seqs_pct | rev_pct > min_seqs_pct]


   if (nrow(good_alignments.dt) == 0) {
     if (debug) {
       cat('* No primer alignments found.', fill=T)
     }
     return(NA)
   } else {

     if (debug) {
       cat('* Found the following primer alignments:', fill=T)
       print(good_alignments.dt)
     }

     fwd_primer_best.dt = primer_alignments.dt[
       fwd_pct_bin==sort(fwd_pct_bin, decreasing=T)[1] & grepl('F$', pr_id)]

     if (verbose & debug) {
       cat('* fwd_primer_best.dt:', fill=T)
       print(fwd_primer_best.dt)
     }

     if (nrow(fwd_primer_best.dt) > 0) {
       # if (nrow(fwd_primer_best.dt[grepl('F$', pr_id)]) > 0) {
       #   fwd_primer_best.dt = fwd_primer_best.dt[grepl('F$', pr_id)][1]
       # } else {
       #   # Multiple of the same primer in the primer table
         fwd_primer_best.dt = fwd_primer_best.dt[1]
       # }
     } else {
       fwd_primer_best.dt = primer_alignments.dt[
         fwd_pct_bin==sort(fwd_pct_bin, decreasing=T)[1] & grepl('R$', pr_id)][1]
     }
     if (verbose & debug) {
       cat('* fwd_primer_best.dt:', fill=T)
       print(fwd_primer_best.dt)
     }

     # Extract F or R for the primer detected to be the forward
     fwd_primer_best_dir = fwd_primer_best.dt[
       , sub('^.*(R|F)$', '\\1', pr_id)]
     if (debug) {
        cat('* Forward primer best direction: ', fwd_primer_best_dir, fill=T)
     }
     if (is.na(fwd_primer_best_dir)) return(NA)

     # Is this primer marked as "forward" ? This would not happen in the
     # case where forward and reverse reads are swapped (or the primers)
     # so we want to catch both cases.
     if (fwd_primer_best_dir == 'F') {
       # OK first one is already F, so find R among the others
       rev_primer_best.dt = good_alignments.dt[
         rev_pct_bin==sort(rev_pct_bin, decreasing=T)[1] & grepl('R$', pr_id)]

       if (nrow(rev_primer_best.dt[grepl('R$', pr_id)]) > 0) {
         rev_primer_best.dt = rev_primer_best.dt[grepl('R$', pr_id)][1]
       } else {
         rev_primer_best.dt = rev_primer_best.dt[1]
       }
     } else if (fwd_primer_best_dir == 'R') {
       # Forward primer is one of the reverse primers
       # This means the sequences are stored in 3' -> 5' direction
       rev_primer_best.dt = primer_alignments.dt[
         rev_pct_bin==sort(rev_pct_bin, decreasing=T)[1] & grepl('F$', pr_id)]

       if (nrow(rev_primer_best.dt[grepl('F$', pr_id)]) > 0) {
         rev_primer_best.dt = rev_primer_best.dt[grepl('F$', pr_id)][1]
       } else {
         rev_primer_best.dt = rev_primer_best.dt[1]
       }
     } else {
       # fwd_primer_best_dir == 'NA' ?
     }
     rev_primer_best_dir = rev_primer_best.dt[
       , sub('^[^RF]+([RF])$', '\\1', pr_id)]
     if (debug) {
       cat('* Reverse primer best direction: ', rev_primer_best_dir, fill=T)
     }
     if (is.na(rev_primer_best_dir)) return(NA)

     if (debug) {
       cat('* Primers with highest % of hits:', fill=T)
       cat('  Forward: ', fwd_primer_best.dt[, pr_id],
           'in', sub('^.*_(3|5)to(3|5)$', '\\1\\\' -> \\2\\\'',
                     fwd_primer_best.dt[, seq_id]), fill=T)
       cat('  Reverse: ', rev_primer_best.dt[, pr_id],
           'in', sub('^.*_(3|5)to(3|5)$', '\\1\\\' -> \\2\\\'',
                     rev_primer_best.dt[, seq_id]), fill=T)

       cat('', fill=T)
     }
     if (
       (fwd_primer_best_dir=='F' & rev_primer_best_dir=='R') |
       (fwd_primer_best_dir=='R' & rev_primer_best_dir=='F')
     ) {
       if (fwd_primer_best.dt[, sub('(F|R)', '', pr_id)] ==
           rev_primer_best.dt[, sub('(F|R)', '', pr_id)]) {
         most_likely_region = fwd_primer_best.dt[, sub('(F|R)', '', pr_id)]
       } else {
         most_likely_region = paste0(
           fwd_primer_best.dt[, sub('(F|R)', '', pr_id)],
           '-',
           rev_primer_best.dt[, sub('(F|R)', '', pr_id)]
         )
       }
       if (debug) {
         cat('* Detected: ', most_likely_region, 'region.', fill=T)
       }
     } else {
       most_likely_region = NA
       if (debug) {
         cat('* No known region in the database with these primers.', fill=T)
       }
     }
     return(most_likely_region)
   }

}



#' PCR primer trimmer
#'
#' Trims PCR primers from FASTQ sequences.
#'
#' @param fq_in A character vector of input file names.
#' @param fq_out A character vector of output file names.
#' @param region Hypervariable region used. Used to automatically retrieve PCR
#' primer sequences for trimming. If this parameter is used (not NULL), then
#' \code{pr_fwd} and \code{pr_rev} are ignored.
#' @param pr_fwd Sequence for the \strong{forward} PCR primer (5' -> 3') in the extended
#' nucleotide code.
#' @param pr_rev Sequence for the \strong{reverse} PCR primer (5' -> 3') in the extended
#' nucleotide code. If you already have 3' -> 5' sequence, use \code{\link{reverse_complement}}
#' to obtain the 5' -> 3' sequence.
#' @param pr_fwd_maxoff Maximum allowed offset from the 5' sequence end for the
#' forward PCR primer. Primers aligned further into the read ignored.
#' @param pr_rev_maxoff Same but for reverse PCR primer, from the 3' end.
#' @param ncpu Number of CPU cores (threads) to use for multithreading. By default
#' all available CPU cores are used (\code{parallel::detectCores()}).
#' @param max_mismatch Maximum allowed number of mismatches for each PCR primer,
#' ignoring any ambiguous nucleotides (not A,C,G or T).
#' @param check_rc Check reverse complements of primer pairs? (default: FALSE)
#' @param overwrite Overwrite target files if they exist (default: TRUE)
#'
#' @export
remove_pcr_primers = function (
  fq_in, fq_out, region=NULL,
  pr_fwd=NULL, pr_rev=NULL, pr_fwd_maxoff=35, pr_rev_maxoff=35,
  max_mismatch=2, check_rc=FALSE,
  return_noprimer=T, ncpu=1, ncpu_sample=rexmap_option('ncpu'),
  verbose=rexmap_option('verbose'), overwrite=TRUE,
  extended_code_alignment=FALSE) {
  # fq_in = input fastq file
  # fq_out = output fastq file (without primers)
  # pr_fwd = forward primer 5'->3'
  # pr_rev = reverse primer 5'->3'
  #
  empty_result = list('fwd_trim'=0, 'rev_trim'=0)
  nts = ifelse(extended_code_alignment, 'ACGTYRWSKMDVHBX', 'ACGT')

  m('------- PCR primer removal ------', time_stamp=T, fill=T, verbose=verbose)

  if (is.null(region) & (is.null(pr_fwd) | is.null(pr_rev))) {
    stop('PCR primer remover: either region of pr_fwd and pr_rev need to be specified.')
  }
  # if region is specified, check that primers exist in the reference table
  if (!is.null(region)) {
    if (nrow(rexmap_option('blast_dbs')[Hypervariable_region==region]) == 0) {
      # Missing hypervariable region
      stop('PCR primer remover: hypervariable region \"', region, '\" not found.')
    }
    pr_fwd = rexmap_option('blast_dbs')[Hypervariable_region==region, Primer1_sequence_5to3]
    pr_rev = reverse_complement(
      rexmap_option('blast_dbs')[Hypervariable_region==region, Primer2_sequence_3to5]
    )
    m('* Using REGION', region,
      paste0('(fwd: ', pr_fwd, ', rev: ', pr_rev, ')'), fill=T, time_stamp=T,
      verbose=verbose)
  } else {
    m('* Using PRIMERS',
      paste0('(fwd: ', pr_fwd, ', rev: ', pr_rev, ')'), fill=T, time_stamp=T,
      verbose=verbose)

  }
  # Count extended DNA symbols
  # Ignore any extended DNA symbols. Shouldn't have too many of them anyway.
  pr_fwd = gsub(paste0('[^', nts, '-]'), 'N', pr_fwd)
  pr_rev = gsub(paste0('[^', nts, '-]'), 'N', pr_rev)
    # pr_fwd = gsub('[^ACGT-]', 'N', pr_fwd)
    # pr_rev = gsub('[^ACGT-]', 'N', pr_rev)

  pr_fwd_rc = reverse_complement(pr_fwd)
  pr_rev_rc = reverse_complement(pr_rev)


  fastq_trimmer = function (
    meta, seq, qual, pr_fwd, pr_rev,
    return_noprimer=return_noprimer
  ) {

    # Search for forward primer.
    aln = C_nwalign(pr_fwd, seq, match=1, mismatch=-1, indel=-1)

    aln_left_offset = regexpr('[^-]', aln[2])[1]
    pr_fwd_left = max(regexpr('[^-]', aln[1])[1],
                      aln_left_offset)
    # pr_fwd_left = regexpr('[^-]', aln[1])[1]
    pr_fwd_right = regexpr('[-]{1,}$', aln[1])[1]-1
    pr_fwd_n = lengths(regmatches(pr_fwd, gregexpr('N', pr_fwd)))
    pr_rev_n = lengths(regmatches(pr_rev, gregexpr('N', pr_rev)))
    # Alignment statistics: match, mismatch, gapopen, gapextend
    if (extended_code_alignment) {
      aln_stat = compare_alignment_ext(
        str_sub(aln[1], start=pr_fwd_left, end=pr_fwd_right),
        str_sub(aln[2], start=pr_fwd_left, end=pr_fwd_right)
      )
    } else {
      aln_stat = compare_alignment(
        str_sub(aln[1], start=pr_fwd_left, end=pr_fwd_right),
        str_sub(aln[2], start=pr_fwd_left, end=pr_fwd_right)
      )
    }
    # Forward primer alignment is acceptable if it's near beginning (within first 5 nts)
    # and if it doesn't have more than 2 mismatches for non-N symbols, which includes indels.
    pr_fwd_found = FALSE
    pr_rev_found = FALSE
    if (pr_fwd_left < pr_fwd_maxoff & aln_stat[2]+aln_stat[3]+aln_stat[4]-pr_fwd_n <= max_mismatch) {
      pr_fwd_found = TRUE
    }
    # Search for reverse primer
    aln2 = C_nwalign(pr_rev, seq, match=1, mismatch=-1, indel=-1)
    pr_rev_left = regexpr('[^-]', aln2[1])[1]
    # pr_rev_right = min(regexpr('[ACGTN][^ACGTN]*$', aln2[1])[1],
    #                    regexpr('[ACGTN][^ACGTN]*$', aln2[2])[1])
    pr_rev_right = min(regexpr(paste0('[', nts, 'N][^', nts, 'N]*$'), aln2[1])[1],
                       regexpr(paste0('[', nts, 'N][^', nts, 'N]*$'), aln2[2])[1])
    if (extended_code_alignment) {
      aln2_stat = compare_alignment_ext(
        str_sub(aln2[1], start=pr_rev_left, end=pr_rev_right),
        str_sub(aln2[2], start=pr_rev_left, end=pr_rev_right)
      )
    } else {
      aln2_stat = compare_alignment(
        str_sub(aln2[1], start=pr_rev_left, end=pr_rev_right),
        str_sub(aln2[2], start=pr_rev_left, end=pr_rev_right)
      )
    }
    if (pr_rev_left > nchar(seq)-nchar(pr_rev)-pr_rev_maxoff & aln2_stat[2]+aln2_stat[3]+aln2_stat[4]-pr_rev_n <= max_mismatch) {
      # Reverse primer found
      pr_rev_found = TRUE
    }

    # Trim if needed
    start = 1L
    end = -1L
    if (pr_fwd_found) start = pr_fwd_right - aln_left_offset + 2
    if (pr_rev_found) end = pr_rev_left - 1

    return(list('meta' = meta,
                'seqs' = str_sub(seq, start=start, end=end),
                'qual' = str_sub(qual, start=start, end=end),
                'trim_fwd' = pr_fwd_found,
                'trim_rev' = pr_rev_found))
  }

  m('  Processing samples:', fill=T, time_stamp=T, verbose=verbose)
  out_per_sample = parallel::mcmapply(function (fq_in_i, fq_out_i) {

    # Various checks
    start_time = Sys.time()
    m_buffer = ''
    if (ncpu == 1) {
      m('*', basename(fq_in_i), ':', time_stamp=T, fill=F, verbose=verbose)
    } else {
      m_buffer = paste0(m_buffer, '* ', basename(fq_in_i), ':')
    }
    if (!file.exists(fq_in_i)) {
      return(empty_result)
    }
    if (file.size(fq_in_i) == 0) {
      # Input file has zero size, so no sequences. This can happen if all the reads
      # are too noisy and have been filtered out in the merging step, if we were
      # unable to merge anything.
      # file.create(fq_out)
      if (ncpu > 1) {
        m(m_buffer, time_stamp=T, fill=F, verbose=verbose)
      }
      m(' (0 file size, skipping).', fill=T, verbose=verbose)
      return(empty_result)
    }
    if (!overwrite & file.exists(fq_out_i)) {
      # We are not overwriting files and output already exists
      if (ncpu > 1) {
        m(m_buffer, time_stamp=T, fill=F, verbose=verbose)
      }
      m(' (Output file exist, skipping).', fill=T, verbose=verbose)
      return(empty_result)
    }
    # Check if output folder exist, and if not create them
    output_folder = dirname(fq_out_i)
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive=T)
      if (ncpu > 1) {
        m_buffer = paste0(m_buffer, ' Create dir.')
      } else {
        m(' Create dir.', fill=F, verbose=verbose, time_stamp=F)
      }
    }

    # Load file
    in_fq = sfastq_reader(fq_in_i)
    if (ncpu > 1) {
      m_buffer = paste0(m_buffer, ' Load.')
    } else {
      m(' Load.', fill=F, verbose=verbose, time_stamp=F)
    }

    # seq = in_fq[['seqs']][1]
    # qual = in_fq[['qual']][1]
    ncpus = min(ncpu_sample, length(in_fq[['seqs']]))

    out_trimmed = parallel::mcmapply(
      function (meta_i, seqs_i, qual_i) {
        result = fastq_trimmer(
          meta=meta_i, seq=seqs_i, qual=qual_i,
          pr_fwd=pr_fwd, pr_rev=pr_rev
        )
        if (check_rc) {
          #if (!result$trim_fwd & !result$trim_rev) {
            # Check RC of primer pair if we couldnt find any
            result2 = fastq_trimmer(
              meta=meta_i, seq=seqs_i, qual=qual_i,
              pr_fwd=pr_rev_rc, pr_rev=pr_fwd_rc)
          #}
            # Take the alignment that managed to get better result
            if (sum(c(result$trim_fwd, result$trim_rev)) <=
                sum(c(result2$trim_fwd, result2$trim_rev))) {
              result = result2
            }
        }
        return(result)
      }, in_fq[['meta']], in_fq[['seqs']], in_fq[['qual']],
      mc.cores=ncpus, SIMPLIFY=F, USE.NAMES=F
    )
    if (ncpu > 1) {
      m_buffer = paste0(m_buffer, ' Trim.')
    } else {
      m(' Trim.', fill=F, verbose=verbose, time_stamp=F)
    }
    # if (verbose) cat('OK.', fill=T)

    fwd_trimmed_list = sapply(out_trimmed, function (x) x$trim_fwd)
    rev_trimmed_list = sapply(out_trimmed, function (x) x$trim_rev)
    fwd_trimmed = sum(fwd_trimmed_list)
    rev_trimmed = sum(rev_trimmed_list)
    any_trimmed = sum(fwd_trimmed_list | rev_trimmed_list)
    both_trimmed = sum(fwd_trimmed_list & rev_trimmed_list)

    # Save results in a new file
    fastq_list_writer(out_trimmed, fq_out_i, ncpu=ncpu_sample)
    pct_trimmed = 100*both_trimmed/length(fwd_trimmed_list)
    if (ncpu > 1) {
      m_buffer = paste0(
        m_buffer, ' OK. Trimmed: ',
        sprintf('% 3d', round(pct_trimmed)),
        '%|',
        sprintf('% 3d', round(100*fwd_trimmed/length(fwd_trimmed_list))),
        '%|',
        sprintf('% 3d', round(100*rev_trimmed/length(fwd_trimmed_list))),
        '% both|fwd|rev ')
    } else {
      m(' OK. Trimmed: ',
        paste0(
          sprintf('% 3d', round(pct_trimmed)), '%|',
          sprintf('% 3d', round(100*fwd_trimmed/length(fwd_trimmed_list))), '%|',
          sprintf('% 3d', round(100*rev_trimmed/length(fwd_trimmed_list))), '% both|fwd|rev '
        ),
        fill=F, verbose=verbose,
        time_stamp=F)
    }

    end_time = Sys.time()
    dt = end_time - start_time
    if (ncpu > 1) {
      m_buffer = paste0(
        m_buffer, ' [', round(dt), ' ', attr(dt, 'units'), ']'
      )
      m(m_buffer, time_stamp=T, fill=T, verbose=verbose)
    } else {
      m(' [', round(dt), ' ', attr(dt, 'units'), ']', fill=T, verbose=verbose,
        time_stamp=F)
    }

    # if (verbose) cat('OK.', fill=T)
    # if (timing) {
    #   end_time = Sys.time()
    #   diff_time = difftime(end_time, start_time, units='secs')
    #   cat('Finished in ', round(as.numeric(diff_time)/60), ' m ',
    #       round(as.numeric(diff_time)%%60, 1), ' s.\n')
    # }
    return(list(
      'fq_in'=basename(fq_in_i),
      'total'=length(out_trimmed),
      'both_trim'=both_trimmed,
      'any_trim'=any_trimmed,
      'fwd_trim'=fwd_trimmed,
      'rev_trim'=rev_trimmed))

  }, fq_in, fq_out, SIMPLIFY=F, USE.NAMES=F, mc.cores=ncpu)

  out.dt = data.table::rbindlist(lapply(out_per_sample, function (x) {
    return(data.table(
      fq_in=x$fq_in,
      total=x$total, both_trim=x$both_trim, any_trim=x$any_trim,
      fwd_trim=x$fwd_trim, rev_trim=x$rev_trim
    ))
  }))
  # out.dt[, fq_in := fq_in_i]
  # setcolorder(out.dt, x('fq_in', 'total', 'both_trim', 'any_trim',
  #                       'fwd_trim', 'rev_trim'))
  return(out.dt[])
}






