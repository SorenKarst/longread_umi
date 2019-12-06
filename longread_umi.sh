#!/bin/bash
# DESCRIPTION
#    longread_umi: pipelines and tools for longread UMI processing.
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryan Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
#


### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [ name ...]

-- longread_umi: pipelines and tools for longread UMI processing.

where:
    -h   Show this help text.
    name Name of tool.
    ...  Commands for tool.

Pipelines:
* nanopore_pipeline:           Generate UMI consensus sequences from raw nanopore data with
                               terminal UMIs.
* qc_pipeline:                 Compare UMI consensus sequences with refence sequences.

Tools:
* check_primer:                Check primer positions in raw longread data. Important
                               input for UMI binning script.
* demultiplex_nb:              Post UMI consensus demultiplexing of using Nanopore 
                               native barcodes

For help with a specific tool type: longread_UMI <name> -h
"

### Terminal Arguments ---------------------------------------------------------


# Import user arguments
while getopts ':hz' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done
TOOL=$1
TOOL_ARG="${*:2}"

# Check missing arguments
if [ -z "${TOOL}" ]; then printf "\n No tool selected!\n\n"; echo "$USAGE"; exit 1; fi; 
if [ -z "${TOOL_ARG}" ]; then
  printf "\n No arguments provided to $TOOL!\n\n"
  TOOL_ARG="-h"
fi

# Paths
export LONGREAD_UMI_PATH="$(dirname "$(readlink -f "$0")")"


### Call tool or command ---------------------------------------------------------------

$LONGREAD_UMI_PATH/scripts/${TOOL}.sh $TOOL_ARG
