import std / [strutils, sequtils]
import ../mpfit
import ggplotnim


# to work around generics bug...
let fnMin = f{`y` - `ey`}
let fnMax = f{`y` + `ey`}

proc plotImpl*[T](fn: FuncProto[T],
                  pRes: seq[float], xs, ys, ey: seq[float],
                  res: mp_result,
                  funcName: string,
                  outfile = "/tmp/plot_fit.pdf",
                  unicode = true,
                  format = ffScientific,
                  precision = -1,
             ) =
  let resTxt = pretty(pRes, res, unicode = unicode, format = format, precision = precision, prefix = "")
  echo "Fit results:"
  echo resTxt

  let xFit = linspace(xs.min, xs.max, 500)
  let yFit = xFit.mapIt(fn(pRes, it))

  ggplot(toDf({"x" : xs, "y" : ys, "ey" : es}), aes("x", "y")) +
    geom_point() +
    geom_errorbar(aes = aes(yMin = fnMin, yMax = fnMax)) +
    geom_line(data = toDf(xFit, yFit), aes = aes("xFit", "yFit"), color = "orange") +
    margin(right = 8) +
    annotate(text = resTxt, left = 1.02, bottom = 0.5, font = font(family = "monospace")) +
    ggtitle("Fit of user defined function: " & funcName) +
    ggsave(outfile, width = 800, height = 480)

template plot*(fn: untyped,
               pRes: seq[float], xs, ys, ey: seq[float],
               res: mp_result,
               outfile = "/tmp/plot_fit.pdf",
               unicode = true,
               format = ffScientific,
               precision = 2
              ): untyped =
  plotImpl(fn, pRes, xs, ys, ey, res, astToStr(fn),
           outfile, unicode, format, precision)
