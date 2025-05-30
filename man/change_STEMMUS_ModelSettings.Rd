% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_.R
\name{change_STEMMUS_ModelSettings}
\alias{change_STEMMUS_ModelSettings}
\title{to change the current STEMMUS Model settings}
\usage{
change_STEMMUS_ModelSettings(
  patch = "D:/model/rSTEMMUS_SCOPE/",
  site_name = NA,
  run_name = NA,
  R_depth = 350,
  J = 1,
  SWCC = 1,
  Hystrs = 0,
  Thmrlefc = 1,
  Soilairefc = 0,
  hThmrl = 1,
  W_Chg = 1,
  ThmrlCondCap = 1,
  ThermCond = 1,
  SSUR = "10^2",
  fc = 0.02,
  Tr = 20,
  T0 = 273.15,
  rwuef = 1,
  rroot = "1.5 * 1e-3",
  SFCC = 1,
  Tot_Depth = 500,
  Eqlspace = 0,
  NS = 1,
  NIT = 30,
  KT = 0,
  NL = 100
)
}
\arguments{
\item{patch}{the patch to the STEMMUS_SCOPE model directory, default "D:/model/rSTEMMUS_SCOPE/",}

\item{site_name, run_name}{the name of the location and the name of run will be changed in the runs folder. The last can be used to name runs with different model parameters or settings,}

\item{R_depth, J, SWCC, Hystrs, Thmrlefc, Soilairefc, hThmrl, W_Chg, ThmrlCondCap, ThermCond, SSUR, fc, Tr, T0, rwuef, rroot, SFCC, Tot_Depth, Eqlspace, NS, NIT, KT, NL}{defaut is c(350,1,1,0,1,0,1,1,1,1,"10^2",0.02,20,273.15,1,"1.5*1e-3",1,500,0,1,30,0,100)}
}
\value{
write the new STEMMUS settings
}
\description{
\code{change_STEMMUS_ModelSettings} change current model settings for STEMMUS
}
\details{
This function change current STEMMUS model settings @seealso \code{\link[=info_STEMMUS_ModelSettings]{info_STEMMUS_ModelSettings()}}.
\strong{important} differently from SCOPE options, any change in the settings will impact all
next runs and not a specific site_name/run_name run.
}
\examples{

\dontrun{
change_STEMMUS_ModelSettings(J = 1)

formattable::formattable(check_STEMMUS_ModelSettings(), align = c("l","c"))
}

}
\seealso{
Other check the current parameters and settings: 
\code{\link{change_SCOPE_constants}()},
\code{\link{change_SCOPE_options}()},
\code{\link{change_STEMMUS_DefineConstants}()},
\code{\link{change_STEMMUS_SoilConstants}()},
\code{\link{check_SCOPE_constants}()},
\code{\link{check_SCOPE_options}()},
\code{\link{check_SCOPE_tsInputs}()},
\code{\link{check_STEMMUS_DefineConstants}()},
\code{\link{check_STEMMUS_ForcingGlobal}()},
\code{\link{check_STEMMUS_ModelSettings}()},
\code{\link{check_STEMMUS_SoilConstants}()},
\code{\link{check_STEMMUS_SoilInitials}()},
\code{\link{check_STEMMUS_SoilProperties}()}
}
\concept{check the current parameters and settings}
