#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [convert_ReadsByLoc_combined_significancecalls.pl]

inputs:

  combinedReadsByLocFile:
    type: File
    inputBinding:
      position: 1
    label: "tabbed files containing number of reads per gene by location"
    doc: "tabbed files containing number of reads per gene by location"

  clipMappedReadNumFile:
    type: File
    inputBinding:
      position: 2

  inputMappedReadNumFile:
    type: File
    inputBinding:
      position: 3

  l2fcWithPvalEnr:
    type: string
    inputBinding:
      position: 4

  l2fc:
    type: string
    inputBinding:
      position: 5

outputs:
  outputFile:
    type: File
    outputBinding:
      glob: $(inputs.l2fcWithPvalEnr)

  l2fcOutputFile:
    type: File
    outputBinding:
      glob: $(inputs.l2fc)