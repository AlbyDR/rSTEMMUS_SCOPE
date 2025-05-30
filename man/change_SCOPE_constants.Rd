% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_.R
\name{change_SCOPE_constants}
\alias{change_SCOPE_constants}
\title{to change the current SCOPE model default constants}
\usage{
change_SCOPE_constants(
  patch = "D:/model/rSTEMMUS_SCOPE/",
  site_name = NA,
  run_name = NA,
  LAT = NA,
  LON = NA,
  IGBP_veg_long = NA,
  hc = NA,
  veg_type_C = 3,
  Cab = 120,
  Cca = 30,
  Cdm = 0.015,
  Cw = 0.009,
  Cs = 0,
  Cant = 0,
  N = 1.4,
  rho_thermal = 0.01,
  tau_thermal = 0.01,
  Vcmo = 80,
  m = 9,
  BallBerry0 = 0.015,
  Type = 0,
  kV = 0.6396,
  Rdparam = 0.015,
  Tparam = c(0.2, 0.3, 288, 313, 328),
  Tyear = 9.74,
  beta = 0.507,
  kNPQs = 0,
  qLs = 1,
  stressfactor = 1,
  fqe = 0.01,
  spectrum = 1,
  rss = 500,
  rs_thermal = 0.06,
  cs = 1180,
  rhos = 1800,
  lambdas = 1.55,
  SMC = 25,
  BSMBrightness = 0.5,
  BSMlat = LAT,
  BSMlon = LON,
  LAI = 3,
  LIDFa = -0.35,
  LIDFb = -0.15,
  leafwidth = 0.1,
  z = 40,
  Rin = 600,
  Ta = 10,
  Rli = 300,
  p = 970,
  ea = 15,
  u = 2,
  Ca = 410,
  Oa = 209,
  zo = 2.5,
  d = 5.34,
  Cd = 0.3,
  rb = 10,
  CR = 0.35,
  CD1 = 20.6,
  Psicor = 0.2,
  CSSOIL = 0.01,
  rbs = 10,
  rwc = 0,
  startDOY = 0,
  endDOY = 366,
  timezn = 1,
  tts = 30,
  tto = 0,
  psi = 90,
  delHaV = 65330,
  delSV = 485,
  delHdV = 149250,
  delHaJ = 43540,
  delSJ = 495,
  delHdJ = 152040,
  delHaP = 53100,
  delSP = 490,
  delHdP = 150650,
  delHaR = 46390,
  delSR = 490,
  delHdR = 150650,
  delHaKc = 79430,
  delHaKo = 36380,
  delHaT = 37830,
  Q10 = 2,
  s1 = 0.3,
  s2 = 313.15,
  s3 = 0.2,
  s4 = 288.15,
  s5 = 1.3,
  s6 = 328.15
)
}
\arguments{
\item{patch}{the patch to the STEMMUS_SCOPE model directory, default "D:/model/rSTEMMUS_SCOPE/",}

\item{site_name, run_name}{the name of the location and the name of run. The last can be used to name runs with different model parameters or settings,}

\item{LAT, LON, startDOY, endDOY, timezn}{the latitude, longitude, timezone and initial and final ts doy of the year ,}

\item{IGBP_veg_long}{the main IGBP landcover class in the point of interest,}

\item{veg_type_C}{default is 3, alternately can be use 4}

\item{Cab, Cca, Cdm, Cw, Cs, Cant, N, rho_thermal, tau_thermal}{PROSPECT submodel parameters for the default values @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}

\item{Vcmo, m, BallBerry0, Type, kV, Rdparam, Tparam, Tyear, beta, kNPQs, qLs, stressfactor}{Leaf_Biochemical parameters for the default values @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}

\item{fqe}{Fluorescence parameters @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}

\item{spectrum, rss, rs_thermal, cs, rhos, lambdas, SMC, BSMBrightness, BSMlat, BSMlon}{Soil parameters @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}

\item{LAI, hc, LIDFa, LIDFb, leafwidth}{Canopy parameters @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}

\item{z, Rin, Rli, Ta, p, ea, u, Ca, Oa}{Meteo parameters @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}

\item{zo, d, Cd, rb, CR, CD1, Psicor, CSSOIL, rbs, rwc}{Aerodynamic parameters @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}

\item{tts, tto, psi}{Angles parameters @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}

\item{delHaV, delSV, delHdV, delHaJ, delSJ, delHdJ, delHaP, delSP, delHdP, delHaR, delSR, delHdR, delHaKc, delHaKo, delHaT, Q10, s1, s2, s3, s4, s5, s6}{Photosynthetic Temperature Dependence Functional parameters @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}}
}
\value{
write the new SCOPE model default constants in the input_data.xlsx file
}
\description{
\code{change_SCOPE_constants} change current SCOPE model default constants in the input_data.xlsx file
}
\details{
This function change current SCOPE model default constants selected to run, @seealso \code{\link[=check_SCOPE_constants]{check_SCOPE_constants()}}
}
\examples{

\dontrun{
change_SCOPE_constants(site_name = "DE-Hai",
                       run_name = "ECdata_36",
                       LAT = st_coordinates(EC_loc)[3,2],
                       LON = st_coordinates(EC_loc)[3,1],
                       IGBP_veg_long = IGBP_Modis_EC$class[3],
                       hc = hc_EC_points$hc[3],
                       veg_type_C = 3,
                       Cab = 120)

formattable::formattable(check_SCOPE_constants(site_name = "DE-Hai", run_name = "ECdata_36")[[1]],
                         align = c("l","c","c","l"))

}

}
\seealso{
Other check the current parameters and settings: 
\code{\link{change_SCOPE_options}()},
\code{\link{change_STEMMUS_DefineConstants}()},
\code{\link{change_STEMMUS_ModelSettings}()},
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
