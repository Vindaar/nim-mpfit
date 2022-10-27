# Package

version       = "0.2.0"
author        = "Sebastian Schmidt"
description   = "Wrapper for the cMPFIT non-linear least squares fitting library (Levenberg-Marquardt)"
license       = "MIT"
srcDir        = "src"

# Dependencies

requires "nim >= 0.18.1"

import ospaths, strutils, strformat
const
  pkgName = "mpfit"
  orgFile = "docs" / (pkgName & ".org")
  rstFile = "docs" / (pkgName & ".rst")
  rstFileAuto = "docs" / (pkgName & "_autogen.rst")

proc basename(f: string): string =
  let (dir, name, ext) = f.splitFile
  result = name

proc removePrefix(f, prefix: string): string =
  result = f
  result.removePrefix(prefix)

# doc generation inspired by `strfmt`
task docs, "Generate HTML docs using the Org file":
  # https://github.com/jgm/pandoc/issues/4749
  exec "pandoc " & orgFile & " -o " & rstFile
  var files: seq[string]
  template walk(path: string, outf: untyped): untyped {.dirty.} =
    for filePath in listFiles(path):
      if filePath.endsWith(".nim"):
        let outfile = outf
        exec &"nim doc {outfile} {filePath}"
        files.add outfile.removePrefix("-o:")
  walk("src", "-o:index.html")
  walk("src" / pkgName, &"-o:{filePath.basename}.html")
  mvFile rstFile, rstFileAuto
  for f in files:
    let fname = f.basename & ".html"
    mvFile fname, "docs/" & $fname
