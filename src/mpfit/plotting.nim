import std / [strutils, sequtils, math]
import ../mpfit
import ggplotnim, seqmath
export ggplotnim # so that user does not have to import it


# to work around generics bug...
let fnMin = f{`y` - `ey`}
let fnMax = f{`y` + `ey`}

proc plotImpl*[T](fn: FuncProto[T],
                  pRes: seq[float], xs, ys, ey: seq[float],
                  res: mp_result,
                  funcName: string,
                  outfile = "/tmp/plot_fit.pdf",
                  unicode = true,
                  xlabel = "", ylabel = "", title = "",
                  xMin = NaN, xMax = NaN,
                  xMargin = 0.05, nPoints = 500,
                  verbose = true,
                  precision = -1,
             ) =
  let resTxt = pretty(pRes, res, unicode = unicode, precision = precision, prefix = "")
  if verbose:
    echo "Fit results:"
    echo resTxt

  let (xLow, xHigh) = block:
    var xL = if classify(xMin) == fcNan: xs.min
             else: xMin
    var xH = if classify(xMax) == fcNan: xs.max
             else: xMax
    let xRange = abs(xH - xL)
    if classify(xMin) == fcNan and classify(xMax) == fcNaN:
      (xL - xRange * xMargin, xH + xRange * xMargin)
    else:
      (xL, xH)
  let xFit = linspace(xLow, xHigh, nPoints)
  let yFit = xFit.mapIt(fn(pRes, it))

  let title = if title.len > 0: title
              else: "Fit of user defined function: " & funcName
  let xlabel = if xlabel.len > 0: xlabel
               else: "x"
  let ylabel = if ylabel.len > 0: ylabel
               else: "y"
  ggplot(toDf({"x" : xs, "y" : ys, "ey" : ey}), aes("x", "y")) +
    geom_point() +
    geom_errorbar(aes = aes(yMin = fnMin, yMax = fnMax)) +
    geom_line(data = toDf(xFit, yFit), aes = aes("xFit", "yFit"), color = "orange") +
    margin(right = 8) +
    annotate(text = resTxt, left = 1.02, bottom = 0.5, font = font(family = "monospace")) +
    xlab(xlabel) + ylab(ylabel) +
    ggtitle(title) +
    ggsave(outfile, width = 800, height = 480)

template plot*(fn: untyped,
               pRes: seq[float], xs, ys, ey: seq[float],
               res: mp_result,
               outfile = "/tmp/plot_fit.pdf",
               unicode = true,
               xlabel = "", ylabel = "", title = "",
               xMin = NaN, xMax = NaN,
               xMargin = 0.05, nPoints = 500,
               verbose = true,
               precision = -1,
              ): untyped =
  plotImpl(fn, pRes, xs, ys, ey, res, astToStr(fn),
           outfile, unicode,
           xlabel, ylabel, title,
           xMin, xMax, xMargin, nPoints,
           verbose,
           precision)
