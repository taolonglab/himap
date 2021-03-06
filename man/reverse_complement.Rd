% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{reverse_complement}
\alias{reverse_complement}
\title{Reverse complement sequence string}
\usage{
reverse_complement(seq_string)
}
\description{
Take a character string \code{seq_string} and reverse complement it using
Biostrings package function \code{reverseComplement}. Supports extended
nucleotide code, so good for reverse complementing PCR primer sequences.
}
\details{
Shortcut for the monstrostity of:
\code{as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq_string)))}
}
