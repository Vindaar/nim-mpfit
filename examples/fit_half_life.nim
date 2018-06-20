import strutils
import seqmath
import sequtils
import strformat
import plotly
import ../src/mpfit_nim
import zero_functional
import chroma

const
  n = 1000
  filename = "data/half_life_muon.txt"

proc expH(p: seq[float], x: float): float =
  ## the function we'd like to fit. Any user defined function needs to
  ## be of the signature
  ## proc[T](p: seq[T], x: T): T
  ## i.e. conform to the FuncProto[T] type defined in mpfit_nim
  result = p[0] * exp(-p[1] * x)
  
proc parseHalfLifeData(filename: string): (seq[float], seq[float]) =
  ## Parse the input file. First create seq of tuple of floats
  ## then convert that to tuple of seq[float]
  let s = readFile(filename).splitLines --> filter('#' notin it and it.len > 0).
                                             map(it.splitWhitespace).
                                             map((it[0].parseFloat,
                                                  it[1].parseFloat))
  result[0] = s --> map(it[0])
  result[1] = s --> map(it[1])  
                                            
proc fitHalfLife(bins, counts, countsErr: seq[float]): (seq[float], mp_result) =
  ## the actual code which performs the fitting. Call the `fit` proc
  ## with the user defined function to be fitted as the first argument,
  ## the initial parameter guess as the second and finally x, y and y_err
  # start parameters
  let p = [1400.0, 1.0]
  # now just call fit
  let (pRes, res) = fit(expH, p, bins, counts, countsErr)
  echoResult(pRes, res = res)
  result = (pRes, res)
  echo &"The lifetime of the muon is ~ {1.0 / pRes[1]:.2f} µs"
                                            
proc plot(bins, counts, countsErr, xFit, yFit: seq[float]) =
  ## plot the data using plotly. The first plot is an ErrorBar plot
  ## the second a lines plot
  
  let d = Trace[float](mode: PlotMode.Markers, `type`: PlotType.ScatterGL,
                       xs: bins, ys: counts, name: "Muon lifetime measurement")
  d.marker = Marker[float](size: @[3.0])
  # add error bars for counts
  d.ys_err = newErrorBar(countsErr, color = Color(r: 0.5, g: 0.5, b: 0.5, a: 1.0))

  let dFit = Trace[float](mode: PlotMode.Lines, `type`: PlotType.ScatterGL,
                         xs: xFit, ys: yFit, name: "Fit to lifetime data")
  
  let layout = Layout(title: "Muon half life measurement", width: 1200, height: 800,
                      xaxis: Axis(title: "time / µs"),
                      yaxis: Axis(title: "# counts"), autosize: false)
  Plot[float](layout: layout, traces: @[d, dFit]).show()

when isMainModule:
  # first parse the data from the file
  let (bins, counts) = parseHalfLifeData(filename)
  # calculates errors: poisson errors on the counts
  let countsErr = counts.mapIt(sqrt(it))
  # perform the fit and echo results
  let (pRes, res) = fitHalfLife(bins, counts, countsErr)

  # calculate some datapoints for the plot of the fit
  let
    xFit = linspace(0, 10, 1000)
    # call fitted function for each value with final parameters
    yFit = xFit.mapIt(expH(pRes, it))

  # plot the data and the fit
  plot(bins, counts, countsErr, xFit, yFit)
  


