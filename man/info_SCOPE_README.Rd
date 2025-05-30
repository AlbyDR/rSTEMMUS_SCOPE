% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/info_.R
\name{info_SCOPE_README}
\alias{info_SCOPE_README}
\title{information about the SCOPE inputs}
\usage{
info_SCOPE_README(
  patch = "D:/model/rSTEMMUS_SCOPE/",
  site_name = NA,
  run_name = NA
)
}
\arguments{
\item{patch}{the patch to the STEMMUS_SCOPE model directory, default \code{"D:/model/rSTEMMUS_SCOPE/"},}

\item{site_name, run_name}{the name of the location and the name of run. The last can be used to name runs with different model parameters or settings,}
}
\value{
a table with information about SCOPE parameters
}
\description{
\code{info_SCOPE_README} present the information about the SCOPE parameters from the readME spreadsheet of the input_data.xlsx file
}
\details{
This is a function show examples of typical input values ranges, parameters units, sensitivity of some
parameters per IBGP landcover class, among other information.
}
\examples{

\dontrun{
info_SCOPE_README(site_name = "DE-Hai", run_name = "ECdata_36")

formattable::formattable(info_SCOPE_README(site_name = "DE-Hai",
                         run_name = "ECdata_36")[[1]], align = c("l","l"))

formattable::formattable(info_SCOPE_README(site_name = "DE-Hai", run_name = "ECdata_36")[[2]],
                         align = c("l","c","c"), caption = "Table1. Example of LIDF values")

formattable::formattable(info_SCOPE_README(site_name = "DE-Hai", run_name = "ECdata_36")[[3]],
align = c("l","c","r","r","r","r","r","c","r","r"),
caption = "Table2. Examples of biochemical parameters for the model biochemical.m
Calibrated temperature sensitivity parameters for Vcmax and Resp")

formattable::formattable(info_SCOPE_README(site_name = "DE-Hai", run_name = "ECdata_36")[[4]],
align = c("l","c","c","r","r","r"), caption = "Table3. Typical ranges of input values")
}

}
\seealso{
Other information about parameters and settings: 
\code{\link{info_IGBP_classes}()},
\code{\link{info_SCOPE_options}()},
\code{\link{info_SCOPE_tsInputs}()},
\code{\link{info_STEMMUS_DefineConstants}()},
\code{\link{info_STEMMUS_ModelSettings}()},
\code{\link{info_STEMMUS_SoilConstants}()}
}
\concept{information about parameters and settings}
