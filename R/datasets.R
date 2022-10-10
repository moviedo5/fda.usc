#' aemet data
#' 
#' @description Series of daily summaries of 73 spanish weather stations selected for the
#' period 1980-2009.  The dataset contains geographic information of each
#' station and the average for the period 1980-2009 of daily temperature, daily
#' precipitation and daily wind speed.
#' @details  Meteorological State Agency of Spain (AEMET), \url{https://www.aemet.es/es/portada}. Government of Spain.\cr
#' It marks 36 UTF-8 string of names of stations and 3 UTF-8 string names of provinces through the function \code{\link{iconv}}.\cr 
#' 
#' @name aemet
#' @docType data
#' @format Elements of aemet:\cr \code{..$df:} Data frame with information of
#' each wheather station:
#' \itemize{
#' \item \code{ind:} Indicated weather station.
#' \item  \code{name:} Station Name.  36 marked UTF-8 strings.
#' \item  \code{province:}Province (region) of Spain. 36 marked UTF-8 strings
#' \item \code{altitude:} Altitude of the station (in meters).
#' \item   \code{year.ini:} Start year.
#' \item  \code{year.end:} End year.
#' \item  \code{longitude:} x geographic coordinate of the station (in decimal degrees).
#' \item  \code{latitude:} y geographic coordinate of the station (in decimal degrees).
#' }
#' The functional variables: 
#' \itemize{
#' \item  \code{...$temp}: mean curve of the average daily temperature
#' for the period 1980-2009 (in degrees Celsius, marked with UTF-8 string).
#' In leap years temperatures for February 28 and 29 were averaged.
#' \item \code{...$wind.speed}: mean curve of the average daily wind speed for the
#' period 1980-2009 (in m/s).
#' \item  \code{...$logprec}: mean curve of the log precipitation for
#' the period 1980-2009 (in log mm). Negligible precipitation 
#' (less than 1 tenth of mm) is replaced by \code{0.05} and no precipitation
#' (0.0 mm) is replaced by \code{0.01}.  Then the logarithm is applied.
#' }
#' @author Manuel Febrero Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @source The data were obtained from the FTP of AEMET in 2009.
#' @keywords datasets
#' @examples
#' \dontrun{
#' data(aemet)
#' names(aemet)
#' names(aemet$df)
#' class(aemet)<-c("ldata","list") # ldata object
#' lat <- ifelse(aemet$df$latitude<31,"red","blue")
#' plot(aemet,col=lat)
#' }
NULL

#' tecator data
#' 
#' Water, Fat and Protein content of meat samples
#' 
#' \code{absorp.fdata} absorbance data for 215 samples. The first 129 were
#' originally used as a training set endpoints the percentages of Fat, Water
#' and Protein.\cr for more details see tecator package
#' 
#' @name tecator
#' @docType data
#' @format The format is: \cr \code{..$absorp.fdata}: absorbance data.
#' \code{fdata} class object with: \cr \itemize{ \item \code{"data"}: Matrix of
#' class \code{fdata} with 215 curves (rows) discretized in 100 points or
#' argvals (columns).\cr \item \code{"argvals"}: 100 discretization points from
#' 850 to 1050mm \cr \item \code{"rangeval"}=(850,1050):
#' range(\code{"argvals"}) \item \code{"names"} list with: \code{main} an
#' overall title "Tecator data set", \code{xlab} title for \code{x} axis
#' "Wavelength (mm)" and \code{ylab} title for \code{y} axis "Absorbances". }
#' \code{..$y}: the percentages of Fat, Water and Protein.  The three contents
#' are determined by analytic chemistry.\cr
#' @author Manuel Febrero-Bande and Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' @keywords datasets
#' @examples
#' data(tecator)
#' names(tecator)
#' names(tecator$absorp.fdata)
#' names(tecator$y)
#' names(tecator$y)
#' class(tecator$absorp.fdata)
#' class(tecator$y)
#' dim(tecator$absorp.fdata)
#' dim(tecator$y)
#' 
NULL

#' poblenou data
#' 
#' NOx levels measured every hour by a control station in Poblenou in Barcelona
#' (Spain).
#' 
#' The dataset starts on 23 February and ends on 26 June, in 2005. We split the
#' whole sample of hourly measures in a dataset of functional trajectories of
#' 24 h observations (each curve represents the evolution of the levels in 1
#' day).\cr Twelve curves that contained missing data were eliminated.
#' 
#' @name poblenou
#' @docType data
#' @format The format is:\cr \code{..$nox}: \code{fdata} class object with: \cr
#' i.- \code{"data"}: Matrix with 115 curves (rows) discretized in 24 points or
#' argvals (columns).\cr ii.- \code{"argvals": 0:23}\cr iii.-
#' \code{"rangeval"=(0,23)}: range(\code{"argvals"}), \cr iv.- \code{"names"}
#' list with: \code{main} an overall title "NOx data set", \code{xlab} title
#' for \code{x} axis "Hours" and \code{ylab} title for \code{y} axis "NOx
#' (mglm^3)".\cr \cr \code{..$df}: Data Frame with (115x3) dimension.  \cr
#' "date" in the first column.\cr Second column ("day.week").  Factor levels:
#' "Monday" 1, "Tuesday" 2, "Wednesday" 3, "Thursday" 4, "Friday" 5, "Saturday"
#' 6 and "Sunday" 7.\cr Third column "day.festive".  Factor levels: "non
#' festive day" 0 and "festive day" 1.\cr
#' @author Febrero-Bande, M and Oviedo de la Fuente, Manuel
#' @references Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W.
#' (2008).  \emph{Outlier detection in functional data by depth measures with
#' application to identify abnormal NOx levels}. Environmetrics 19, 4, 331-345.
#' @source \url{https://mediambient.gencat.cat/ca/inici}
#' @keywords datasets
#' @examples
#' 
#' data(poblenou)
#' names(poblenou)
#' names(poblenou$nox) 
#' nox<-poblenou$nox
#' class(nox)
#' ind.weekend<-as.integer(poblenou$df[,"day.week"])>5
#' plot(nox,col=ind.weekend+1)
#' 
NULL



#' phoneme data
#' 
#' Phoneme curves
#' 
#' The following instructions have been used file: \cr
#' \url{https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/npfda-phondiscRS.txt}\cr
#' of \code{Phoneme dataset} file.
#' 
#' @name phoneme
#' @docType data
#' @format Elements of phoneme:\cr \code{..$learn}: learning sample of curves.
#' \code{fdata} class object with: i.- \code{"data"}: Matrix of class
#' \code{fdata} with 250 curves (rows) discretized in 150 points or argvals
#' (columns).\cr, ii.- \code{"argvals"}, iii.- \code{"rangeval"}:
#' range(\code{"argvals"}), iv.- \code{"names"} list with: \code{main} an
#' overall title "Phoneme learn", \code{xlab} title for \code{x} axis
#' "frequencies" and \code{ylab} title for \code{y} axis "log-periodograms".\cr
#' \cr \code{..$test}: testing sample of curves. \code{fdata} class object
#' with: i.- \code{"data"}: Matrix of class \code{fdata} with 250 curves (rows)
#' discretized in 150 points or argvals (columns).\cr, ii.- \code{"argvals"},
#' iii.- \code{"rangeval"}: range(\code{"argvals"}), iv.- \code{"names"} list
#' with: \code{main} an overall title "Phoneme learn", \code{xlab} title for
#' \code{x} axis "frequencies" and \code{ylab} title for \code{y} axis
#' "log-periodograms".\cr \cr \code{..$classlearn}:learning class numbers (as
#' factor). Factor levels: "sh" 1, "iy" 2, "dcl" 3, "aa" 4 and "ao" 5.\cr \cr
#' \code{..$classtest}: testing class numbers (as factor). Factor levels: "sh"
#' 1, "iy" 2, "dcl" 3, "aa" 4 and "ao" 5.\cr
#' @author Manuel Febrero-Bande and Manuel Oviedo de la Fuente
#' <manuel.oviedo@@udc.es>
#' @references Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}. Free
#' access on line at \url{https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/}
#' @source
#' \url{https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/npfda-datasets.html}
#' @keywords datasets
#' @examples
#' 
#' data(phoneme)
#' names(phoneme)
#' names(phoneme$learn)
#' class(phoneme$learn)
#' dim(phoneme$learn)
#' table(phoneme$classlearn)
#' 
NULL

#' Mithochondiral calcium overload (MCO) data set
#' 
#' The mithochondiral calcium overload (MCO) was measured in two groups
#' (control and treatment) every 10 seconds during an hour in isolated mouse
#' cardiac cells. In fact, due to technical reasons, the original experiment
#' [see Ruiz-Meana et al. (2000)] was performed twice, using both the "intact",
#' original cells and "permeabilized" cells (a condition related to the
#' mitochondrial membrane).
#' 
#' 
#' @name MCO
#' @docType data
#' @format Elements of MCO:\cr \code{..$intact}: \code{fdata} class object with
#' ``intact cells''curves,\cr \itemize{ \item \code{"data"}: Matrix of class
#' \code{fdata} with 89 intact cells curves (rows) measured every 10 seconds
#' during an hour in isolated mouse cardiac cell. \item \code{"argvals"}, 360
#' discretization points from seond 0 to 3590. \item \code{"rangeval"}:
#' range(\code{"argvals"}). \item \code{"names"} list with: \code{main} an
#' overall title "Control Intact Treatment", \code{xlab} title for \code{x}
#' axis "seconds" and \code{ylab} title for \code{y} axis "Ca". }
#' \code{..$classintact}: Factor levels of ``intact cells'' curves: "1" control
#' group and "2" treatment group.\cr
#' 
#' \code{..$permea}: \code{fdata} class object with ``permeabilized cells''
#' curves (whose membrane has been removed), \itemize{ \item \code{"data"}:
#' Matrix of class \code{fdata} with 90 permeabilizzed cells curves (rows)
#' measured every 10 seconds during an hour in isolated mouse cardiac cell.
#' \item \code{"argvals"}, 360 discretization points from seond 0 to 3590.
#' \item \code{"rangeval"}: range(\code{"argvals"}). \item \code{"names"} list
#' with: \code{main} an overall title "Control Intact Treatment", \code{xlab}
#' title for \code{x} axis "seconds" and \code{ylab} title for \code{y} axis
#' "Ca". } \code{..$classpermea}: Factor levels of ``permeabilized cells''
#' curves: "1" control group and "2" treatment group.\cr
#' @note The structure of the curves during the initial period (first 180
#' seconds) of the experiment shows a erratic behavior (not very relevant in
#' the experiment context) during this period.
#' @references
#' 
#' Ruiz--Meana M, Garcia-Dorado D, Pina P, Inserte J, Agullo L, Soler--Soler J.
#' Cariporide preserves mitochondrial proton gradient and delays ATP depletion
#' in cardiomyocytes during ischemic conditions. \emph{American Journal
#' Physiology Heart Circulatori Physiology}. 2003 Sep;285(3):H999--1006.
#' @keywords datasets
#' @examples
#' 
#' data(MCO)
#' names(MCO)
#' par(mfrow=c(1,2))
#' plot(MCO$intact,col=MCO$classintact)
#' plot(MCO$permea,col=MCO$classpermea)
#' 
NULL

