#!/usr/bin/env cwltool

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [count_reads_broadfeatures_frombamfi_PEmap.pl]

inputs:

  clipBamFile:
    type: File
    inputBinding:
      position: 1
    label: "IDR File"
    doc: "IDR File"
  gencodeGTFFile:
    type: File
    inputBinding:
      position: 2
    label: "gencode GTF file"
    doc: "gencode GTF file"
  gencodeTableBrowserFile:
    type: File
    inputBinding:
      position: 3
    label: "gencode parsed ucsc tableformat file"
    doc: "gencode parsed ucsc tableformat file"
  output:
    type: string
    inputBinding:
      position: 4
    label: "output file"
    doc: "output file"

outputs:
  outputFile:
    type: File
    outputBinding:
      glob: $(inputs.output)
