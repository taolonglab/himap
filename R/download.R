#' Download blastn and makeblastdb executables
#'
#' Downloads NCBI blastn and makeblastdb executables for the Operating System
#' running R and copies it to the RExMap package folder. This is executed automatically
#' if blastn and makeblastdb executables are missing upon attaching RExMap package.
#'
#' @param verbose (default: TRUE) Whether or not to print progress during download.
#'
#' @export
download_blast = function (verbose=T, override_os=NULL) {
   # First find the platform
   if (verbose) cat('RExMap: Download blastn+makeblastdb command line utilities\n')

   if (is.null(override_os)) {
      sys_os = detect_os()
   } else {
      sys_os = override_os
   }
   if (verbose) cat('  Operating system:', sys_os, '\n')

   # Check internet connection
   if (!internet_connection()) {
      stop('Error: Internet connection unavailable.')
   }

   # Check newest version
   if (verbose) cat('  Checking latest version...')
   curlpath = curl_path()
   ncbi_blast_url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/'
   latest_files = system2(curlpath, c(ncbi_blast_url), stdout=T, stderr=F)

   if (sys_os == 'linux') {
      pattern = 'ncbi-blast-[0-9]+\\.[0-9]+\\.[0-9]+\\+-x64-linux\\.tar\\.gz$'
   } else if (sys_os == 'macos') {
      pattern = 'ncbi-blast-[0-9]+\\.[0-9]+\\.[0-9]+\\+-x64-macosx\\.tar\\.gz$'
   } else {
      stop('Error: Unsupported OS.')
   }
   latest_file = latest_files[grep(pattern, latest_files)]
   latest_file2 = sub(paste0('^.*(', pattern, ')'), '\\1', latest_file)
   blast_version = sub('.*ncbi-blast-([0-9]+\\.[0-9]+\\.[0-9]+)\\+-.*', '\\1',
                       latest_file2)
   if (verbose) cat('', blast_version, 'found.\n')
   blast_url = paste0(ncbi_blast_url, latest_file2)

   # Download into temp folder
   if (verbose) cat('  Downloading executables...')
   blast_file_gz = file.path(tempdir(), basename(blast_url))
   download.file(blast_url, blast_file_gz, quiet=T)
   if (verbose) cat('OK.\n')

   # Unpack
   if (verbose) cat('  Unpacking...')
   blast_files = untar(blast_file_gz, verbose=F, list=T, exdir = tempdir())
   untar(blast_file_gz, verbose=F, exdir = tempdir())
   if (verbose) cat(' OK.\n')

   # Copy to RExMap folder
   if (verbose) cat('  Copying...')
   blast_exec = file.path(dirname(blast_file_gz),
                          grep('/bin/blastn$', blast_files, value=T))
   makedb_exec = file.path(dirname(blast_file_gz),
                          grep('/bin/makeblastdb$', blast_files, value=T))
   himap_path = find.package('himap')
   himap_blastn_path = file.path(himap_path, 'exec', paste0('blastn_', sys_os))
   himap_makedb_path = file.path(himap_path, 'exec', paste0('makeblastdb_',
                                                            sys_os))
   copy_success = file.copy(blast_exec, himap_blastn_path, overwrite=T)
   if (!copy_success) stop('Error. Failed: ', blast_exec, ' --[copy]--> ',
                           himap_blastn_path)
   copy_success = file.copy(makedb_exec, himap_makedb_path, overwrite=T)
   if (!copy_success) stop('Error. Failed: ', makedb_exec, ' --[copy]--> ',
                           himap_makedb_path)
   if (verbose) cat(' OK.\n')

   # Clean up files in the temp folder
   if (verbose) cat('  Cleaning temporary files...')
   temp_dirs = grep('/$', blast_files, value=T)
   temp_files = blast_files[!(blast_files %in% temp_dirs)]
   file.remove(file.path(tempdir(), temp_files))
   unlink(file.path(tempdir(), temp_dirs), recursive=T)
   file.remove(blast_file_gz)
   if (verbose) cat(' OK.\nDone.\n')
}

#' Check if there's an active internet connection
#'
#' Return TRUE or FALSE
internet_connection = function(url='https://cran.r-project.org') {
    out = tryCatch(
       readChar(con=url, nchars=1),
       error = function(e) return(FALSE),
       warning = function(w) return(FALSE)
    )
    if (class(out) == 'logical') return(FALSE)
    else return(TRUE)
}

#' curl path for system2 call
curl_path = function () {
   # Try system path first
   curl_sys_path = suppressWarnings(system2(c('which', 'curl'), stdout=T,
                                            stderr=F))
   if (length(curl_sys_path) > 0) return(curl_sys_path)
   # Curl not found in system path use the one bundled in RExMap
   return(file.path(find.package('himap'), 'exec', paste0('curl_', detect_os())))
}


#' Get latest BLAST executable path
#' (currently only for MacOS)
latest_blast_url = function () {
   curlpath = curl_path()
   latest_files = system2(
      curlpath,
      c('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/'),
      stdout=T, stderr=F
   )
   pattern = 'ncbi-blast-[0-9]+\\.[0-9]+\\.[0-9]+\\+-x64-macosx\\.tar\\.gz$'
   latest_file = latest_files[grep(pattern, latest_files)]
   latest_file2 = sub(paste0('^.*(', pattern, ')'), '\\1', latest_file)
}

#' Construct the pcr_primers_table after database update
#' (this is not finished yet)
#'
#' Automatically generates and verified inst/extdata/pcr_primers_table.txt
#' file based on the files in database/ folder.
#'
#' @param backup_old Whether to backup the current pcr_primers_table before
#' overwritting. Default: TRUE
#' @param verbose Verbose output, print status.
#'
# #' @export
update_pcr_primers_table = function (backup_old=T, verbose=T) {
   himap_database_path = file.path(find.package(pname_l), 'inst', 'database')
   if (!dir.exists(himap_database_path)) {
      stop('\nError: RExMap database folder ', himap_database_path, ' not found!')
   }
   pcr_sequences_filename = system.file('extdata', 'pcr_primer_sequences.txt',
                                        package=pname_l)
   if (pcr_sequences_filename == '') {
      stop('\nError: RExMap file with PCR sequences pcr_primer_sequences.txt not found!')
   }
   pcr_sequences.dt = data.table::fread(pcr_sequences_filename)
   if (verbose) {
      cat('* scanning ', himap_database_path, ' folder...')
   }
   himap_database_path_files = dir(himap_database_path, full.names=T)
   himap_database_path_files_base = basename(himap_database_path_files)
   if (verbose) cat(' OK.\n')
   # Find pairs of files matching regex:
   # ^(V[^_]+)_([0-9]+F\\-[0-9]+R)_hang([0-9]+).*\\.([^\\.]+)$
   # Regex capture groups:
   # \1 : hypervariable region(s) sequences
   # \2 : hypervariable region(s) universal PCR primer names
   # \3 : sequence overhang +- XX nts included for each side to ensure query
   #      sequences do not overextend ends of reference sequences
   # \4 : database version/date
   if (verbose) cat('* inspecting files.')
   regex_prefix = '^(V[^_]+)_([0-9]+F\\-[0-9]+R)_hang([0-9]+)[_]?([0-9]{4}\\-[0-9]{2}\\-[0-9]{2})?.*'
   # Find potential database entries
   potential_db_entries = unique(sub(
      regex_prefix, '\\1_\\2_hang\\3_\\4', basename(himap_database_path_files)
   ))
   if (verbose) cat('.')
   potential_db_entries = sub('_$', '', potential_db_entries)
   potential_db_entries = grep('^V[0-9]+', potential_db_entries, value=T)
   if (verbose) cat('.')

   # Now find 4 files matching potential_db_entries[i].* that end with
   # _R.txt, .nhr, .nsq, .nin
   # out of those 4, the n* extensions need to have identical full name besides
   # file extension. Also each primer needs to already be present in the
   # pcr_sequences.dt table.
   final_db_entries = list()

   for (db_entry in potential_db_entries) {
      if (verbose) cat('.')
      nhr_file = grep(paste0(db_entry, '.*\\.nhr$'), himap_database_path_files_base, value=T)
      nin_file = grep(paste0(db_entry, '.*\\.nin$'), himap_database_path_files_base, value=T)
      nsq_file = grep(paste0(db_entry, '.*\\.nsq$'), himap_database_path_files_base, value=T)
      # Do all three files exist?
      if (length(nhr_file) == 0 | length(nin_file) == 0 | length(nsq_file) == 0) {
         next
      }
      # Do all three files match apart from the extension?
      if (
         (sub('\\.nhr$', '', nhr_file) != sub('\\.nin$', '', nin_file)) |
         (sub('\\.nin$', '', nin_file) != sub('\\.nsq$', '', nsq_file))
      ) {
         next
      }
      # Do we find _R.txt file with the same database_prefix??
      r_file = grep(paste0(db_entry, '.*_R\\.txt$'), himap_database_path_files_base, value=T)
      if (length(r_file) == 0) {
         next
      }

      # OK we have all 4 files, generate a table entry if we have sequences
      hypervar_region = sub(regex_prefix, '\\1', db_entry)
      fr_primers = strsplit(sub(regex_prefix, '\\2', db_entry), '-', fixed=T)[[1]]

      # Forward and reverse primer sequences
      f_seq = pcr_sequences.dt[Primer==fr_primers[1], Primer_sequence_5to3]
      r_seq = pcr_sequences.dt[Primer==fr_primers[2], Primer_sequence_3to5]
      if (length(f_seq) == 0 | length(r_seq) == 0) {
         # Fwd or rev primer sequence not found, skip
         next
      }

      db_date = sub(regex_prefix, '\\4', db_entry)
      if (grepl('^V', db_date) | db_date == '') db_date = '0000-00-00'

      final_db_entries[[length(final_db_entries)+1]] = data.table(
         Primer1=fr_primers[1],
         Primer2=fr_primers[2],
         Primer1_sequence_5to3=f_seq,
         Primer2_sequence_3to5=r_seq,
         Hypervariable_region=hypervar_region,
         DB=sub('\\.nhr$', '', nhr_file),
         table=r_file,
         date=db_date,
         hr=hypervar_region, # This is used to iterate over so we update the other col
         hr_date=paste0(hypervar_region, '_', db_date)
      )

   }


   # Each hypervariable region should have a unique name, and we use that name
   # to call functions such as blast()
   # Automatically omit date from the latest versions ()
   final_db_table.dt = data.table::rbindlist(final_db_entries)
   print(head(final_db_table.dt))

   if (verbose) cat('.')
   final_db_table.dt = final_db_table.dt[order(Hypervariable_region, -date)]
   final_db_table.dt[
      , Hypervariable_region2 := na.omit(c(Hypervariable_region[1], hr_date[2:.N])), by=hr]
   final_db_table.dt[
      , Hypervariable_region := Hypervariable_region2
   ]
   final_db_table.dt[
      , Hypervariable_region2 := NULL
   ]
   final_db_table.dt[, c('hr', 'hr_date') := NULL]
   if (verbose) cat('OK.\n')

   # Overwrite old table
   # save old copy
   if (verbose) cat('* saving output file... ')
   old_table_filename = system.file('extdata', 'pcr_primers_table.txt', package=pname_l)
   if (backup_old) {
      new_table_filename = file.path(
         dirname(old_table_filename),
         paste0(
            'pcr_primers_table_', Sys.Date(), '_',
            gsub(':', '-', format(Sys.time(), '%X'), fixed=T),
            '_backup.txt'
         )
      )
      file.copy(old_table_filename, new_table_filename)
   }

   data.table::fwrite(final_db_table.dt, old_table_filename)
   rexmap_setoption('blast_dbs', final_db_table.dt)
   if (verbose) cat('OK.\n')
   # print(final_db_table.dt)

}




#' RExMap database from RExMapDB
#'
#' Update/replace database with RExMapDB latest version. This deletes the
#' current database files from find.package('rexmap')/inst/database, then
#' downloads the files from RExMapDB Github repository.
#'
#' @param delete_old Delete old files during an update.
#' @param verbose Display progress messeges.
#' @param use_api Use GitHub API. Warning: has 60 requests/hr hard limit.
#'
#' @export

download_database = function (delete_old=TRUE, verbose=TRUE, use_api=FALSE) {

   if (verbose) cat('RExMap database update\n')
   if (verbose) cat('  Check internet connection... ')
   if (!internet_connection()) stop('\nError: Internet connection unavailable.')
   if (verbose) cat('OK.\n')
   # local_db_path = file.path(find.package(pname_l), 'database')
   local_db_path = system.file('database', package='rexmap')
   if (verbose) cat('Local DB path:', local_db_path, '\n')
   curlpath = curl_path()
   if (verbose) cat('Curl:', curlpath, '\n')
   if (verbose) cat('Retrieving current database files:\n')
   local_db_files = dir(local_db_path, 'V[0-9].*\\.(?:txt|nin|nhr|nsq)',
                        full.names=TRUE)
   if (length(local_db_files) == 0) {
      if (verbose) cat('  No files found.\n')
      old_files.dt = data.table(filename=character())
   } else {
      old_files.dt = as.data.table(file.info(local_db_files, extra_cols=F))
      old_files.dt[, filename := basename(local_db_files)]
      # for (f in local_db_files) {
      #    file.remove(f)
      #    if (verbose) cat('  L deleted:', f, '\n')
      # }
   }
   if (verbose) cat('OK.\n')


   if (verbose) cat('Downloading database files:\n')

   if (use_api) {
      # Run curl | retrieve rexmapdb list of databases
      json_out = paste(system2(curlpath, c(
         '-L', 'https://api.github.com/repos/taolonglab/rexmapdb/contents/database/'
      ), stdout=T, stderr=F), collapse='')

      # Convert JSON output to data table
      db.dt = tryCatch(
         data.table::as.data.table(as.list(jsonlite::fromJSON(json_out)))[, c(1,4,11)],
         error = function (e) NULL)
      if (is.null(db.dt)) {
         if (verbose) cat('ERROR retrieving JSON from Github:', json_out, '\n')
         return()
      }
      db.dt = db.dt[grepl('V[0-9][-]?(V[0-9])?', name)]
      for (i in 1:nrow(db.dt)) {
         f = db.dt[i, 3]
         if (verbose) cat('   ', db.dt[i, name], '\n')
         # download.file(f, file.path(himap_database_path, basename(f)), extra=c(
         system2(curlpath, c(
            '-s',
            '-H', '"Accept: application/vnd.github.v3.raw"',
            '-L', f,
            '-o', file.path(local_db_path, db.dt[i, name])
         ), stdout=F, stderr=F)
      }
   } else {
      # Fall back to HTML extraction if API doesnt work
      db_list_out = system2(
         curlpath, c(
            '-s',
            'https://github.com/taolonglab/rexmapdb/tree/master/database'
         ),
         stdout=TRUE, stderr=FALSE
      )
      db_files_out = grep(
         'href="/taolonglab/rexmapdb/blob/master/database/V[0-9][^\\.]+\\.nsq"',
         db_list_out, value=T)
      db_file_prefix = 'https://github.com/taolonglab/rexmapdb/raw/master/database/'
      db_files = sub('.*href="/taolonglab/rexmapdb/blob/master/database/(V[0-9][^\\.]+?\\.nsq)".*',
                     '\\1', db_files_out)
      db_files = c(
         db_files,
         sub('\\.nsq$', '\\.nin', db_files),
         sub('\\.nsq$', '\\.nhr', db_files),
         sub('\\.nsq$', '\\.txt', db_files)
      )
      db_files_url = paste(db_file_prefix, db_files, sep='')
      for (i in 1:length(db_files)) {
         f = db_files[i]
         if (verbose) cat('    ', f, '\n')
         system2(curlpath, c(
            '-s',
            '-L', db_files_url[i],
            '-o ', file.path(local_db_path, f)
         ), stdout=F, stderr=F)
      }
      db.dt = data.table(name=basename(db_files))
   }
   db.dt[, Hypervariable_region := sub(
      '^(V[0-9][-]?[V]?[0-9]?)_.*$', '\\1', name)]
   db.dt[, Primer_pair := sub(
      '^V[0-9][-]?[V]?[0-9]?_([0-9]{1,4}F-[0-9]{1,4}[m]?R)_.*$', '\\1', name)]
   db.dt[, Database_name := sub('\\.[^\\.]+$', '', name)]

   if (verbose) cat('* OK.\n')

   # Delete old files that have not been updated
   if (verbose) cat('Deleting old files...')
   new_local_db_files = dir(local_db_path, 'V[0-9].*\\.(?:txt|nin|nhr|nsq)',
                            full.names=TRUE)
   new_files.dt = as.data.table(file.info(new_local_db_files, extra_cols=F))
   new_files.dt[, filename := basename(new_local_db_files)]
   files.dt = merge(old_files.dt, new_files.dt, by='filename', all.x=T)[
      !is.na(mtime.y) & !is.na(ctime.y) & !is.na(atime.y) & !is.na(size.y)
   ][mtime.x != mtime.y | size.x != size.y]
   if (nrow(files.dt) > 0) {
      if (delete_old) {
         fr = file.remove(files.dt[, filename])
      }
   }
   if (verbose) cat(' OK.\n')



   # Generate a reference table
   if (verbose) cat('  Downloading reference table...')
   if (use_api) {
      pcr_out = paste(trimws(system2(curlpath, c(
         '-s',
         '-H', '"Accept: application/vnd.github.v3.raw"',
         '-L', 'https://api.github.com/repos/taolonglab/rexmapdb/contents/primers/pcr_primers_table.txt'
         # '-o', system.file('extdata', 'pcr_primers_table.txt', package=pname_l)
      ), stdout=TRUE, stderr=FALSE)), collapse='\n')
   } else {
      pcr_out = system2(
         curlpath, c(
            '-s',
            '-L', 'https://raw.githubusercontent.com/taolonglab/rexmapdb/master/primers/pcr_primers_table.txt'
         ), stdout=T, stderr=F
      )
   }
   pcr.dt = tryCatch(
      data.table::fread(text=pcr_out),
      error = function (e) NULL)
   if (is.null(pcr.dt)) {
      if (verbose) cat('ERROR retrieving PCR primers table:', pcr_out, '\n')
      return()
   }
   pcr.dt[, Primer_pair := paste(Primer1, Primer2, sep='-')]
   pcr_db.dt = merge(
      pcr.dt,
      unique(
         db.dt[
            , .(Hypervariable_region, Primer_pair,
                DB=Database_name, table=paste(Database_name, '.txt', sep=''))
         ]
      ),
      by=c('Hypervariable_region', 'Primer_pair')
   )
   pcr_db.dt[, date := sub('^.*([0-9]{4}-[0-9]{2}-[0-9]{2}).*$', '\\1', DB)]
   if (verbose) cat('OK.\n')

   # Write pcr primers table
   if (verbose) cat('  Writing PCR primers table...')
   pcr_out_file = system.file('extdata', 'pcr_primers_table.txt',
                              package=pname_l)
   fwrite(pcr_db.dt, pcr_out_file, sep='\t')
   if (verbose) cat('OK.\n')
   # Reload rexmap options
   reload_blast_dbs()
   if (verbose) cat('Done.\n')

}




#' Update the RExMap database
#'
#' Updates the database from a RExMap GitHub repository with files pre-generated by
#' the RExMap authors. For manual database generation see:
#' https://www.github.com/taolonglab/himapdb .
#'
#' @param verbose Verbose output, print status message during execution.
#' @param source_repo Download RExMap DB from himap Github repository
#' (default: 'rexmap') or RExMap db repository ('rexmapdb').
#' @param update_ref_table Download and update reference table from the RExMap
#' repository. Typically this table can be generated automatically on update
#' when this is set to FALSE (default).
#'
#' @export
update_database = function (verbose=T, source_repo='rexmap', update_ref_table=FALSE) {
   if (verbose) cat('RExMap database update\n')
   if (verbose) cat('* check internet connection... ')
   if (!internet_connection()) stop('\nError: Internet connection unavailable.')
   if (verbose) cat('OK.\n')
   if (verbose) cat('* downloading database files:\n')
   himap_database_path = file.path(find.package(pname_l), 'inst', 'database')
   # Implements this method:
   # https://medium.com/@caludio/how-to-download-large-files-from-github-4863a2dbba3b
   curlpath = curl_path()
   if (source_repo == 'rexmap') {
      json_out = paste(system2(curlpath, c(
         '-L', 'https://api.github.com/repos/taolonglab/rexmap/contents/inst/database/'
      ), stdout=T, stderr=F), collapse='')
   } else if (source_repo == 'rexmapdb') {
      json_out = paste(system2(curlpath, c(
         '-L', 'https://api.github.com/repos/taolonglab/rexmapdb/contents/database/'
      ), stdout=T, stderr=F), collapse='')
   }

   db.dt = data.table::as.data.table(as.list(jsonlite::fromJSON(json_out)))[, c(1,4,11)]
   # Select only database files for hypervariable regions
   db.dt = db.dt[grepl('V[0-9][-]?(V[0-9])?', name)]
   for (i in 1:nrow(db.dt)) {
      f = db.dt[i, 3]
      if (verbose) cat('* -', db.dt[i, name], '\n')
      # download.file(f, file.path(himap_database_path, basename(f)), extra=c(
      system2(curlpath, c(
         '-s',
         '-H', '"Accept: application/vnd.github.v3.raw"',
         '-L', f,
         '-o', file.path(himap_database_path, db.dt[i, name])
      ), stdout=F, stderr=F)
   }
   if (verbose) cat('* OK.\n')
   if (update_ref_table) {
      if (verbose) cat('* downloading reference table...')
      system2(curlpath, c(
         '-s',
         '-H', '"Accept: application/vnd.github.v3.raw"',
         '-L', 'https://api.github.com/repos/taolonglab/rexmap/contents/inst/extdata/pcr_primers_table.txt',
         '-o', system.file('extdata', 'pcr_primers_table.txt', package=pname_l)
      ), stdout=F, stderr=F)
      if (verbose) cat('OK.\n')
   } else {
      update_pcr_primers_table(verbose=verbose)
   }

   if (verbose) cat('Done.\n')

}

