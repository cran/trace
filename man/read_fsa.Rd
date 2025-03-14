% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_fsa.R
\name{read_fsa}
\alias{read_fsa}
\title{Read fsa file}
\usage{
read_fsa(files)
}
\arguments{
\item{files}{a chr vector of fsa file names. For example, return all the fsa files in a directory with 'list.files("example_directory/", full.names = TRUE, pattern = ".fsa")'.}
}
\value{
A list of fragments_trace objects
}
\description{
Read fsa file into memory and create fragments_trace object
}
\details{
read_fsa is just a wrapper around \code{\link[seqinr:read.abif]{seqinr::read.abif()}} that reads the fsa file into memory and stores it inside a fragments_trace object. That enables you to use the next function \code{\link[=find_ladders]{find_ladders()}}.
}
\examples{

fsa_file <- read_fsa(system.file("abif/2_FAC321_0000205983_B02_004.fsa", package = "seqinr"))
plot_data_channels(fsa_file)

}
\seealso{
\code{\link[=find_ladders]{find_ladders()}}, \code{\link[=plot_data_channels]{plot_data_channels()}}
}
