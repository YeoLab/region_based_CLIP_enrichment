#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [combine_ReadsByLoc_files.pl]

inputs:

  readsByLocFiles:
    type: File[]
    inputBinding:
      position: 1
    label: "tabbed files containing number of reads per gene by location"
    doc: "tabbed files containing number of reads per gene by location"

  output:
    type: string
    label: "output file"
    doc: "output file"

outputs:
  outputFile:
    type: File
    outputBinding:
      glob: $(inputs.output)

stdout: $(inputs.output)
