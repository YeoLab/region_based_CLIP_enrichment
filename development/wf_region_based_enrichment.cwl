#!/usr/bin/env cwltool

cwlVersion: v1.0
class: Workflow

requirements:
  - class: MultipleInputFeatureRequirement

inputs:

  clipBamFile:
    type: File
  inputBamFile:
    type: File

  gencodeGTFFile:
    type: File
  gencodeTableBrowserFile:
    type: File
  clipBroadFeatures:
    type: string
  inputBroadFeatures:
    type: string
  clipMappedReadNum:
    type: string
  inputMappedReadNum:
    type: string
  combinedOutput:
    type: string
  l2fcWithPvalEnr:
    type: string
  l2fc:
    type: string

outputs:
  clipBroadFeatureCountsFile:
    type: File
    outputSource: clipCountReadsBroadFeatures/outputFile

  inputBroadFeatureCountsFile:
    type: File
    outputSource: inputCountReadsBroadFeatures/outputFile

  combinedOutputFile:
    type: File
    outputSource: combineReadsByLocFiles/outputFile

  l2fcWithPvalEnrFile:
    type: File
    outputSource: significanceTest/outputFile
  l2fcFile:
    type: File
    outputSource: significanceTest/l2fcOutputFile

steps:
  clipCalculateMappedReadNum:
    run: calculate_readnum.cwl
    in:
      bamFile: clipBamFile
      output: clipMappedReadNum
    out:
      - outputFile

  inputCalculateMappedReadNum:
    run: calculate_readnum.cwl
    in:
      bamFile: inputBamFile
      output: inputMappedReadNum
    out:
      - outputFile

  clipCountReadsBroadFeatures:
    run: count_reads_broadfeatures_frombamfi_PEmap.cwl
    in:
      clipBamFile: clipBamFile
      gencodeGTFFile: gencodeGTFFile
      gencodeTableBrowserFile: gencodeTableBrowserFile
      output: clipBroadFeatures
    out:
      - outputFile

  inputCountReadsBroadFeatures:
    run: count_reads_broadfeatures_frombamfi_PEmap.cwl
    in:
      clipBamFile: inputBamFile
      gencodeGTFFile: gencodeGTFFile
      gencodeTableBrowserFile: gencodeTableBrowserFile
      output: inputBroadFeatures
    out:
      - outputFile

  combineReadsByLocFiles:
    run: combine_ReadsByLoc_files.cwl
    in:
      readsByLocFiles: [clipCountReadsBroadFeatures/outputFile, inputCountReadsBroadFeatures/outputFile]
      output: combinedOutput
    out:
      - outputFile

  significanceTest:
    run: convert_ReadsByLoc_combined_significancecalls.cwl
    in:
      combinedReadsByLocFile: combineReadsByLocFiles/outputFile
      clipMappedReadNumFile: clipCalculateMappedReadNum/outputFile
      inputMappedReadNumFile: inputCalculateMappedReadNum/outputFile
      l2fcWithPvalEnr: l2fcWithPvalEnr
      l2fc: l2fc
    out:
      - outputFile
      - l2fcOutputFile