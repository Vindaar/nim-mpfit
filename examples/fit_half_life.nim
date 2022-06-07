import std / [strutils, sequtils, strformat]
import pkg / [zero_functional, seqmath]
import mpfit

const
  filename = "data/half_life_muon.txt"

func expH(p: seq[float], x: float): float =
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

when isMainModule:
  # first parse the data from the file
  let (bins, counts) = parseHalfLifeData(filename)
  # calculates errors: poisson errors on the counts
  let countsErr = counts.mapIt(sqrt(it))
  # perform the fit and echo results
  let (pRes, res) = fitHalfLife(bins, counts, countsErr)

  # plot the data and the fit
  import mpfit / plotting # import plotting convenience function
  plot(
    expH, pRes, # the function we fit and resulting fit params
    bins, counts, countsErr, # the input data & errors
    res, # the `mp_result` returned from the `fit` call
    xMin = 0.0, xMax = 10.0, # customize range of plot
    xlabel = "time / μs", ylabel = "# counts", title = "Muon half life measurement", # and labels
    outfile = "../media/muon_lifetime_measurement.png", # save as png
    verbose = false) # we set `verbose` to false, as we already `echoResult` manually
