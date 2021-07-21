
#' Shortcut IF/ELSE function
#'
#' if \code{test} evaluates to TRUE, return \code{yes}. Otherwise return \code{no}.
#'
#' @export
#'
#' @examples
#' ie(1==2, '1 equals 2', '1 does not equal 2')
#' # Output: '1 does not equal 2'
ie = function(test, yes, no) {
  if (test) yes
  else no
}

#' Message reporting function
m = function (..., fill=TRUE, time_stamp=TRUE, verbose=TRUE) {
   if (verbose) {
      if (!time_stamp) {
         out_str = paste(..., collapse=' ')
      } else {
         out_str = paste(as.character(as.POSIXlt(Sys.time())), ' | ', ..., collapse=' ')
      }
      if (fill) {
         out_str = paste0(out_str, '\n')
      }
      cat(out_str)
   }
}

lu = function (x) length(unique(x))

#' Extract sample identifiers from FASTQ filenames
#'
#' @param filenames A character vector of filenames.
#' @param separator A single character used delimiting sample id from the
#' rest of FASTQ filename.
#'
#' @export
sampleids_from_filenames = function (filenames, separator='_') {
  # Get sample ids from FASTQ filenames
  sapply(strsplit(basename(filenames), separator, fixed=T), `[`, 1)
}

#' Read files from a folder satisfying a pattern.
#'
#' Shortcut for \code{sort(dir(path, pattern, full.names=T))}
#'
#' @param path Full path to the folder to list files in.
#' @param pattern Pattern for pattern matching. If '' then
#' list all files.
#'
#' @examples
#' read_files('~/data/diabimmune/fastq_tutorial', 'R1') # Get forward reads
#' read_files('~/data/diabimmune/fastq_tutorial', 'R2') # Get reverse reads
#'
#'
#' @export
read_files = function (path, pattern='') {
  sort(dir(path.expand(path), pattern, full.names=T))
}

#' Plot-friendly OSU labels
#'
#' Condense a character vector of strain names to a more readable shorter
#' string, consisting of species names with optional number of strains matched
#' for each species.
#'
#' This function expects input string (or a vector of strings) that consists of strain
#' names. By default, the strain names are separated by comma, and an underscore is used
#' as a whitespace symbol (these can be changed: see parameters `xsep`, `xws`).
#'
#'
#'
#' @param x (Required) A character vector of strain names.
#' @param xsep A character used for separating strain names in `x`.
#' Default: ','
#' @param xws A character used as a whitespace in `x`. Default: '_'
#' @param ws A character used as a whitespace in output (return value) of this
#' function. Default: ' ' (space)
#' @param sep A character string used as a species separator in output of this
#' function. Default: ', ' (comma followed by a space)
#' @param include_sp One of {'single', 'always', 'never'}. Controls display of
#' strains with unassigned species names (having 'sp.' after genus name). Single
#' keeps these names when this is the single strain in its genus. 'Always' keeps
#' these names always. 'Never' never keeps them in the output.
#' For example is the strain assignments are:
#'    'Blautia_sp._AF19-13LB,Blautia_sp._AF14-40,Fusicatenibacter_saccharivorans,
#'    Ruminococcus_sp._1,Ruminococcus_gnavus_2'
#' setting it to `single` (default) will return:
#'    'Blautia sp., Fusicatenibacter saccharivorans, Ruminococcus gnavus'
#' since these unnamed Blautia strains are the only strains in Blautia genus and named
#' strains are prioritized over unnamed strains for genera . If
#' this is set to never, then the 'sp.' strains are always ommited and the output
#' will show the best NAMED species only. Setting to always will return:
#'    'Blautia sp., Fusicatenibacter saccharivorans, Ruminococcus gnavus/sp.'
#' to indicate that there are more hits in Ruminococcus genus than the only named
#' one (gnavus).
#' Default: 'single'
#' @param decreasing Sort order of the output species. Default: TRUE sorts the
#' output starting with the species that had most strain hits.
#' @param show_count Display count of strains next to each species name in the
#' output. Default: TRUE
#' @param show_single If show_count == TRUE, this controls whether to display
#' number 1 for species that match only 1 strain. Default: FALSE (i.e. numbers
#' are only shown for species mapping to > 1 strain).
#' @param count_wrap A character vector of length 2 controling the string
#' format of the count output, Wraps the count in the two strings, i.e. 3 will
#' be shown as [3] by default.
#' @param group_genera A boolean controling the grouping of species belonging
#' to the same Genus. Default: TRUE.
#' @param group_genera_sep A character symbol used for separating different
#' species in the same Genus group.
#' @param trim_genera_len An integer controling the trimming length for long
#' Genus names. Any names longer than `trim_genera_len` will be trimmed to
#' `trim_general_len` followed by a `trim_symbol`.
#' Default: NULL (do not trim labels)
#' @param trim_species_len Same as `trim_genera_len` but for species names
#' instead of Genus names. Default: NULL (do not trim labels)
#' @param trim_symbol A character symbol used after trimmed names. Default: '.'
#' (A natural choice).
#' @param replace_bacterium A boolean controling the handling of strains with
#' unassigned species names that are not `sp.` but are `bacterium` or `cf.`. If
#' this is set to TRUE any names ending with `bacterium` and `cf.` are substituted
#' with `sp.` for consistency. Default: TRUE
#' @param multi_species_sub If there are too many species under a single
#' genus (greater than max_species), then use this string as replacement. First element
#' is the replacement if the hit is one species, second if there are more hits.
#' (default: `c(sp., spp.)``)
#' @param keep_family_strains Keep strains that do not even have genus assigned, such as
#' Ruminococcaceae sp. if there are other, strains with named genera.
#' @param keep_sp_strains Show full strain names for species names such as
#' Somegenus_sp._AB01. In this case full strain name is added to the species name
#' @param return_single_entry If the input string consists of a single strain name
#' should we just return that one and ignore other filters? Default: TRUE
#'
#'
#' @export
print_species = function (
   x, xsep=',', xws='_', ws=' ', sep=', ', include_sp='single', decreasing=TRUE,
   show_count=TRUE,
   show_single=FALSE,
   count_wrap=c('[', ']'),
   group_genera=TRUE,
   group_genera_sep='/',
   trim_genera_len=NULL,
   trim_species_len=NULL, trim_symbol='.',
   replace_bacterium=TRUE,
   max_species=NULL,
   multi_species_sub=c('sp.', 'spp.'),
   keep_family_strains=FALSE,
   keep_order_strains=FALSE,
   keep_class_strains=FALSE,
   keep_phylum_strains=FALSE,
   keep_sp_strains=FALSE,
   keep_sp_strains_symbol='_',
   return_single_entry=TRUE,
   debug=FALSE,
   show_progress=TRUE) {

   # x = input string of strain names
   #
   # * Options:
   # xsep = strain separator symbol in input (default: ,)
   # xws = whitespace symbol used in input (default: _)
   # ws = whitespace symbol used in output (default: [space])
   # include_sp = {'never', 'always', 'single'}
   #    if the list of species contains "sp." when to output that?
   #    if set to 'never', it will never show any strains with "sp.", i.e. no
   #      species name.
   #    if set to 'always', always keep "sp." strains in the output
   #    if set to 'single', the "sp." strains will be output only if that strain
   #      is the only strain in that genus.

   # trim_ = trim Genus and/or species names to specified lengths. Be careful:
   #    if this is too short, it can collapse multiple labels into one.

   # Load taxonomy and extract unique phylum, class, order, family names
   tax.dt = load_taxonomy()
   unique_phyla = tax.dt[, unique(phylum)]
   unique_classes = tax.dt[, unique(class)]
   unique_orders = tax.dt[, unique(order)]
   unique_families = tax.dt[, unique(family)]

   # Unnamed species labels; first entry will be the default output label
   # Second and third entries labels renamed to the first label (3 is hard-coded)
   # so keep that in mind when adding new unnamed species labels.
   sp_labels = c('sp\\.', 'bacterium', 'cf\\.')

   # replace_bacterium = replaces 'bacterium' species annotation with 'sp.'
   # x_out = c()

   x_out = sapply(1:length(x), function (i) {

      x_i = x[i]
      # Split string into individual strains
      x_strains = strsplit(x_i, xsep, fixed=T)[[1]]

      # Remove unnamed assignments like 'bacterium_LF3' is there are other
      # better assignments.
      unnamed_filter = !(x_strains %like% paste0('^bacterium', xws))
      if (length(x_strains[unnamed_filter]) > 0) {
         x_strains = x_strains[unnamed_filter]
      }

      # If we only have 1 assignment, return that one. The only change is replacing
      # the input whitespace symbol with the output whitespace symbol.
      if (length(x_strains) == 1 & return_single_entry) {
         if (!show_count) {
            count_regex = sub(']', '\\]',
                              sub('[', '\\[',
                                  paste0(count_wrap[1], '[0-9]+', count_wrap[2]),
                                  fixed=T
                              ), fixed=T)
            return(
               sub(paste0(ws, '+$'), '', sub(count_regex, '', gsub(xws, ws, x_strains[1], fixed=T)))
            )
         } else {
            return(gsub(xws, ws, x_strains[1], fixed=T))
         }
      }


      # If we filter out strains without strain name and species name, i.e. hits to
      # names ending with sp., cf., or bacterium, and we have other assignments left then
      # just use properly named ones.
      # sp_strain_filter = apply(
      #    sapply(
      #       sp_labels,
      #       function (slab) grepl(paste0(slab, '$'), x_strains)
      #    ),
      #    1,
      #    any
      # )

      sp_strain_filter = sapply(x_strains, function (xs_i) {
         # Check each strains vs each element of sp_labels.
         xs_i_unnamed = any(sapply(sp_labels,
                                   function (slab) grepl(paste0(slab, '$'), xs_i)))
         return(xs_i_unnamed)
      })

      if (length(x_strains[!sp_strain_filter]) > 0) {
         x_strains = x_strains[!sp_strain_filter]
      }

      # This changes all unnamed species names to 'sp.' for consistency.
      if (replace_bacterium) {
         for (j in 2:length(sp_labels)) {
            x_strains = sub(paste0(xws, sp_labels[j], '$'), paste0(xws, sp_labels[1]), x_strains)
            x_strains = sub(paste0(xws, sp_labels[j], xws), paste0(xws, sp_labels[1], xws), x_strains)
         }
      }

      # Find all unique genera and group by these
      #  vector with all genera for specific input
      genera = sub(paste0('^([^', xws, ']+)', xws, '.*'), '\\1', x_strains)
      #  species vector contains both actual genera and species

      # If keep_sp_strains is TRUE, then we generate "species" name for each "sp." strain
      # by adding the strain designation to sp as a species name
      if (!keep_sp_strains) {
         species = sub(paste0('^([^', xws, ']+)', xws, '([^', xws, ']+).*'),
                       paste0('\\1', ws, '\\2'), x_strains)
      } else {
         species = sub(paste0('^([^', xws, ']+)', xws, '([^', xws, ']+).*'),
                       paste0('\\1', ws, '\\2'), x_strains)
         species_sp = sub(paste0('^([^', xws, ']+)', xws, '(sp\\..*)'),
                          paste0('\\1', ws, '\\2'), x_strains)
         species = mapply(function (s1, s2) if (nchar(s2 > s1) & s2 %like% 'sp\\.') s2 else s1,
                          species, gsub(xws, ws, species_sp), USE.NAMES=F)
      }

      if (!keep_phylum_strains) {
         # This filter has value TRUE for strains that have an actual Genus name
         genera_as_phylum_filter = !(genera %in% unique_phyla)
         if (any(genera_as_phylum_filter)) {
            # Do we have any genera that is not an order level genus name?
            genera = genera[genera_as_phylum_filter]
            species = species[genera_as_phylum_filter]
         }
      }

      if (!keep_class_strains) {
         # This filter has value TRUE for strains that have an actual Genus name
         genera_as_class_filter = !(genera %in% unique_classes)
         if (any(genera_as_class_filter)) {
            # Do we have any genera that is not an order level genus name?
            genera = genera[genera_as_class_filter]
            species = species[genera_as_class_filter]
         }
      }

      if (!keep_order_strains) {
         # This filter has value TRUE for strains that have an actual Genus name
         genera_as_order_filter = !(genera %in% unique_orders)
         if (any(genera_as_order_filter)) {
            # Do we have any genera that is not an order level genus name?
            genera = genera[genera_as_order_filter]
            species = species[genera_as_order_filter]
         }
      }

      # Filter out family level genus names, if we have anything left. This will prioritize
      # strains named to proper Genus instead of strain where "Genus" is actually Family
      if (!keep_family_strains) {
         # This filter has value TRUE for strains that have an actual Genus name
         genera_as_family_filter = !(genera %in% unique_families)
         if (any(genera_as_family_filter)) {
            # Do we have any genera that is not a family level genus name?
            genera = genera[genera_as_family_filter]
            species = species[genera_as_family_filter]
         }
      }



      if (debug) {
         cat('Genera:\n')
         cat(genera)
         cat('\n')
         cat('Species:\n')
         cat(species)
         cat('\n')

      }
      # # This changes all unnamed species names to 'sp.' for consistency.
      # if (replace_bacterium) {
      #    species = sub(paste0(ws, 'bacterium'), paste0(ws, 'sp\\.'), species)
      #    species = sub(paste0(ws, 'cf\\.'), paste0(ws, 'sp\\.'), species)
      # }

      # if (keep_sp_strains) {
      #    sp_strain_filter = species %like% 'sp\\.$'
      #    # species[sp_strain_filter] = gsub(xws, ws, x_strains[sp_strain_filter], fixed=T)
      #    sp_strain_names = gsub(ws, keep_sp_strains_symbol, gsub(xws, keep_sp_strains_symbol, sub(
      #       '^[^_]+_sp\\.[_]?(.*)$',
      #       '\\1',
      #       x_strains[sp_strain_filter]
      #    ), fixed=T), fixed=T)
      #    species[sp_strain_filter] = paste(species[sp_strain_filter], sp_strain_names,
      #                                      sep=keep_sp_strains_symbol)
      # }


      # Pre-process genus/species names if they need trimming
      if (!is.null(trim_genera_len)) {
         if (debug) {
            cat('Trimming genera...\n')
         }
         if (trim_genera_len > 1) {
            # Trim using regex, cleaner this way and "obvious" when doing the species
            # (Genus species combos)

            # First pull up a filter of which genera actually need to be trimmed
            # and leave the rest as it is.
            #
            # The trimming works by taking trim_genera_len first characters of genus
            # name, then adding a trim_symbol (default .) afterwards. SO the total
            # genus name will have trim_genera_len+1 characters.
            trim_genera_filter = nchar(genera) > trim_genera_len+1
            genera[trim_genera_filter] = sub(
               paste0('^([^', ws, ']{', trim_genera_len, '}).*$'),
               '\\1',
               genera[trim_genera_filter]
            )
            species[trim_genera_filter] = sub(
               paste0('^([^', ws, ']{', trim_genera_len, '})[^', ws, ']+',
                      '(', ws, '[^', ws, ']+.*)$'),
               paste0('\\1', trim_symbol, '\\2'),
               species[trim_genera_filter]
            )

            # Very rarely genus name will end with a dot (saw it few times), in which case
            # appending a dot as a trim_symbol will make it two dots so cleanup those here.
            genera = gsub('\\.\\.', '\\.', genera)
            species = gsub('\\.\\.', '\\.', species)

         }
      }

      if (debug) {
         cat('Genera after trimming:\n')
         cat(genera)
         cat('\n')
         cat('Species after trimming:\n')
         cat(species)
         cat('\n')
      }

      if (!is.null(trim_species_len)) {
         if (trim_species_len > 1) {
            species = sub(
               paste0('^([^', ws, ']+', ws, ')([^', ws, ']{', trim_species_len, '})[^', ws, ']+$'),
               paste0('\\1\\2', trim_symbol),
               species
            )
            species = gsub('\\.\\.', '\\.', species)
         }
      }

      # Now process include_sp, by generating a T/F filter and applying it
      # If we NEVER include 'sp.' strains. This only makes sense if we have strain names
      # that are not "sp.". If all strains are sp. strains, then do nothing.
      if (include_sp == 'never') {
         sp_filter = !(species %like% 'sp\\.')  # TRUE for named species
         # if at least one element of sp_filter is TRUE, then we can take out sp. strains
         # if (all(!sp_filter)) { # Do we have ANY named species left?
         if (debug) {
            cat('sp_filter:\n')
            cat(sp_filter, '\n')
         }
         if (any(sp_filter)) { # Do we have ANY named species left?
            # sp_filter.list = list()
            # for (g in unique(genera)) {
            #    g_filter = (genera == g)
            #    if (length(unique(species[g_filter])) > 1) {
            #       sp_filter.list[[length(sp_filter.list)+1]] = !(
            #          (species %like% 'sp\\.') & (species %like% paste0('^', g))
            #       )
            #    } else {
            #       sp_filter.list[[length(sp_filter.list)+1]] = rep(TRUE, length(species))
            #    }
            # }
            # sp_filter = sapply(transpose(sp_filter.list), all)
         } else {
            sp_filter = rep(TRUE, length(species))
         }
      }

      # Are we filtering out all species? In this case fall back to single include_sp
      if (include_sp == 'single') {
         # For each genus, count number of species assignments. If that number is 1
         # and it's value is 'sp.' then it's a single species, so include it.
         # sp_filter = c()
         sp_filter = !(species %like% 'sp\\.')  # TRUE for named species
         if (any(sp_filter)) {
            sp_filter.list = list()
            for (g in unique(genera)) {
               g_filter = (genera == g)
               if (length(unique(species[g_filter])) > 1) {
                  sp_filter.list[[length(sp_filter.list)+1]] = !(
                     (species %like% 'sp\\.') & (species %like% paste0('^', g))
                  )
               } else {
                  sp_filter.list[[length(sp_filter.list)+1]] = rep(TRUE, length(species))
               }
            }
            sp_filter = sapply(transpose(sp_filter.list), all)
         } else {
            # All strains are sp. strains, so just return them all
            sp_filter = rep(TRUE, length(species))
         }
      }

      if (include_sp == 'always') {
         sp_filter = rep(TRUE, length(species))
      }


      # cat('sp_filter: ', paste(as.character(sp_filter), collapse=' '), '\n')
      # genera = genera[sp_filter]
      species_f = species[sp_filter]
      genera_f = genera[sp_filter]

      species.c = unique_sorted(species_f, decreasing=decreasing,
                                show_count=show_count, count_wrap=count_wrap,
                                show_single=show_single)

      # If group_genera == TRUE, group display by genus
      if (group_genera) {
         genera.c = sub(paste0('^([^', ws, ']+)', ws, '.*'), '\\1', species.c)
         genera_uniq.c = unique(genera.c)
         # We need to do this loop again after possible sp/genus omits from
         # unique_sorted()
         out.c = c()
         for (g in genera_uniq.c) {
            species_only_list = sub(
               paste0('^', g, ws), '', species.c[species.c %like% paste0('^', g, ws)]
            )
            if (!is.null(max_species)) {
               if (length(species_only_list) > max_species) {
                  if (length(species_only_list) == 1) {
                     species_only_list = multi_species_sub[1]
                  } else {
                     species_only_list = multi_species_sub[2]
                  }

               }
            }
            species_only = paste(
               species_only_list,
               collapse=group_genera_sep
            )
            out.c[length(out.c)+1] = paste(g, species_only, collapse=ws)
         }
         return(paste(out.c, collapse=sep))
      } else {
         return(paste(species.c, collapse=sep))
      }
   }) # end for x_i in x loop
   return(x_out)

}# , vectorize.args=c('x'), USE.NAMES=FALSE)


unique_sorted = function (v, decreasing=TRUE, show_count=TRUE,
                          count_wrap=c('[', ']'), show_single=FALSE) {
   # Sort vector of values by most common ones (unless decreasing==TRUE, in which
   # case just invert the output), optionally showing the count within count_wrap
   # If show_single==FALSE and show_count==TRUE, don't show [1] counts
   v_tab = sort(table(v), decreasing=decreasing)
   if (!show_count) {
      return(names(v_tab))
   } else {
      out = paste(
         names(v_tab), count_wrap[1], as.integer(v_tab), count_wrap[2], sep=''
      )
      if (!show_single) {
         out = gsub(paste0(count_wrap[1], 1, count_wrap[2]), '', out, fixed=T)
      }
      return(out)
   }
}

#' Strain vector to shorter strain string
#'
#' Condense a character vector of strain names to a more managable
#' shorter string.
#'
#' @param strains (Required) A character vector of strain names, with underscore _ instead of space
#' to separate genus, species and strain designations.
#' @param raw FALSE/TRUE If TRUE, then the output is just a list of concatenated strain
#' names (can be very large!). Default: FALSE.
#' @param deduplicate FALSE/TRUE If strain name occurs multiple times, count only unique names.
#' Default: TRUE.
#' @param nmax Integer. If the number of strains is below \code{nmax}, just list all the strain names.
#'
#' @export
print_strains = function (strains, raw=F, deduplicate=T,
                          nmax=rexmap_option('print_strains_nmax')) {
  # Input: vector of strains
  # Prints a single string with reduced list of strains
  # such that any non-_bacterium or _sp. strain is shown
  # on a species level if raw=FALSE
  # This is a much simpler function than print_species() and is used
  # for very basic OSU label output.

  if (deduplicate) {
     uniq_strains = unique(strains)
  } else {
     uniq_strains = strains
  }

  if (raw | length(uniq_strains) <= nmax) {
     return(paste(uniq_strains, collapse=','))
  }

  # Also add strains for species that occur only once!!
  if (length(strains)==1) {
     return(strains)
  } else {
    genuses = gsub('^[ ]*([^_]+)_.*', '\\1', uniq_strains)
    species = gsub('^[^_]+_([^_]+)[_]?.*', '\\1', uniq_strains)
    na_filter = grepl('^(sp\\.|bacterium)', species)
    # sp1 is the list of species that does not end with sp. or bacterium
    sp1 = paste(genuses[!na_filter], species[!na_filter], sep='_')
    # sp2 is the list of strains which do not have a species name (sp. or bacterium)
    sp2 = uniq_strains[na_filter]
    ft = as.table(sort(table(c(sp1, sp2)), decreasing=T))
    ft2 = ft[ft>1]
    # Get full strain names for species that occur only once
    sp_once = names(ft[ft==1])
    st_once = c()
    for (sp in sp_once) {
      st_once = c(st_once, grep(sp, uniq_strains, fixed=T, value=T))
    }
    out = c()
    if (length(ft2) > 0) {
      out = paste(paste0(paste(names(ft2), ft2, sep='_['), ']'), collapse=',')
    }
    if (length(out) > 0 & length(st_once) > 0) {
      out = paste(out, paste(st_once, collapse=','), sep=',')
    } else if (length(out) == 0 & length(st_once) > 0) {
      out = paste(st_once, collapse=',')
    }
    return(out)
    # return(gsub('_[1]', '', out, fixed=T))
  }
}

#' Summarize OSU labels
#'
#' Shorten sometimes long OSU labels which can list many strains or species
#' into a more managable string for printing.
#'
#'
#' @export
summarize_osu_label = Vectorize(function (x, show.sp=FALSE, sep=',') {
   x_split = strsplit(x, sep, fixed=T)[[1]]
   if (length(x_split) == 1) return(x)
   x_split = x_split[!grepl('_sp.', x_split, fixed=T)]
   if (length(x_split) == 0) {
      # There are no species without sp. annotation
      # Return genuses
      genuses = gsub('^[ ]*([^_]+)_.*', '\\1', strsplit(x, sep, fixed=T)[[1]])
      return()
   } else {
      x_split = print_strains(x_split, raw=F, nmax=1)
      return(x_split)
   }
}, 'x', USE.NAMES=F)


print_fixed_length_string = function (x, len=30) {
   if (nchar(x) >= len+3) {
      ndots=3
      left = floor((len-ndots)/2)
      right = ceiling((len-ndots)/2)
      return(paste0(substr(x, 1, left), '...', substr(x, nchar(x)-right+1, nchar(x))))
   } else if (nchar(x) == len+2) {
      ndots=2
      left = floor((len-ndots)/2)
      right = ceiling((len-ndots)/2)
      return(paste0(substr(x, 1, left), '..', substr(x, nchar(x)-right+1, nchar(x))))
   } else if (nchar(x) == len+1) {
      ndots=1
      left = floor((len-ndots)/2)
      right = ceiling((len-ndots)/2)
      return(paste0(substr(x, 1, left), '.', substr(x, nchar(x)-right+1, nchar(x))))
   } else {
      return(x)
   }
}

string_fixed_len = function (x, len=30) {
  if (is.na(x)) return(NA)
  if (nchar(x) > len) return(paste0(substr(x, 1, len-3), '...'))
  else return(x)
}

print.data.table2 = function (dt, width=rexmap_option('string_maxwidth'),
                             topn=rexmap_option('maxrows')) {
   dt2 = head(dt, topn)
   extra_rows = nrow(dt) - topn
   for (j in names(dt2)) {
      if (dt2[, class(get(j))] == 'character') {
         values = dt2[[j]]
         set(dt2, j = j, value = sapply(values, string_fixed_len, len=width))
      }
   }
   print.data.frame(dt2)
   if (extra_rows > 0) cat('...', paste0('(', nrow(dt), ' rows)'), fill=T)
}

#' Reverse complement sequence string
#'
#' Take a character string \code{seq_string} and reverse complement it using
#' Biostrings package function \code{reverseComplement}. Supports extended
#' nucleotide code, so good for reverse complementing PCR primer sequences.
#'
#' Shortcut for the monstrostity of:
#' \code{as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq_string)))}
#'
#' @export
reverse_complement = function (seq_string) {
   as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq_string)))
}

all_exist = function (files) all(file.exists(files))

#' Generate random sequence of given length
#'
#' @export
random_sequences = function(len=50, n=10, dictionary=c('A', 'C', 'G', 'T')) {
  do.call(paste0, replicate(len, sample(dictionary, n, TRUE), FALSE))
}


#' Write data table to file
#'
#' Like base function \code{\link{write.table}} but with normal defaults.
#'
#' @export
write_table = function (table, output, sep='\t', verbose=F, ...) {
  if (verbose) cat('* writing ', output, '...')
  write.table(table, output, sep=sep, quote=F, row.names=F, ...)
  if (verbose) cat('OK.', fill=T)
}

#' Filter and trim reads. Wrapper for dada2::filterAndTrim
#'
#' @export
filter_and_trim = function (fq_in, fq_out, maxEE=2,
                            multithread=rexmap_option('ncpu'),
                            verbose=rexmap_option('verbose'),
                            compress=FALSE,
                            ...) {
  dada2::filterAndTrim(fq_in, fq_out, maxEE=maxEE, multithread=multithread,
                       verbose=verbose, compress=compress, ...)
}


format_string = Vectorize(function (x, max_chars=50) {
   x_len = nchar(x)
   if (x_len == max_chars) {
      return(x)
   } else if (x_len > max_chars) {
      return(paste0(substring(x, 1, max_chars-3), '...'))
   } else {
      # x_len < max_width
      return(sprintf(paste0('% ', max_chars, 's'), x))
   }
}, vectorize.args=c('x'), SIMPLIFY=T, USE.NAMES=F)


datatable_to_string = function (dt, max_chars=50) {
   # Print data table to string; useful for logging
   columns = names(dt)
   column_classes = sapply(dt, class)
   # Find the maximum width of each column
   column_widths = sapply(columns, function (col) {
      min_width = min(
         c(max_chars,
           max(dt[, nchar(get(col))], nchar(col), na.rm=T)
         ), na.rm=T)
      return(min_width)
   })

   # Format individual columns, then transpose list so that the values are
   # listed by rows
   columns_values = transpose(lapply(
      columns,
      function (col) {
         column = c(col, dt[, as.character(get(col))])
         column_width = min(c(
            max_chars,
            max(nchar(column), na.rm=T)
         ), na.rm=T)
         format_string(column, max_chars=column_width)
      }
   ))

   # Collapse rows with space, columns with '\n'
   out = paste(unlist(lapply(columns_values, paste, collapse=' ')), collapse='\n ')
   return(paste0(' ', out))

}
