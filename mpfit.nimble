# Package

version       = "0.1.0"
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

# doc generation inspired by `strfmt`
task docs, "Generate HTML docs using the Org file":
  # https://github.com/jgm/pandoc/issues/4749
  exec "pandoc " & orgFile & " -o " & rstFile
  for filePath in listFiles("src"):
    if filePath.endsWith(".nim"):
      let outfile = "-o:docs/index.html"
      echo outfile
      echo filepath
      exec &"nim doc {outfile} {filePath}"
  for filePath in listFiles("src" / pkgName):
    if filePath.endsWith(".nim"):
      let outfile = &"-o:docs/{filePath.basename}.html"
      echo outfile
      echo filepath
      exec &"nim doc {outfile} {filePath}"

  mvFile rstFile, rstFileAuto
  #mvFile htmlFileNimDoc, htmlFileIndex
