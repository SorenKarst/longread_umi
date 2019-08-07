#!/bin/bash
# DESCRIPTION
#    longread-UMI: pipelines and tools for longread UMI processing.
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    license	GNU General Public License
#
# To-do:


### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [ name ...]

-- longread_UMI: pipelines and tools for longread UMI processing.

where:
    -h   Show this help text.
    name Name of tool.
    ...  Commands for tool.

Pipelines:
* longread_UMI_pipeline:       Generate UMI consensus sequences from raw longread data with
                               terminal UMIs.
* longread_UMI_mockanalysis:   Compare UMI consensus sequences generated from the
                               ZymoBIOMICS MIcrobial community standard (D6305).

Commands:
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
PIPELINE_PATH=`readlink ${BASH_SOURCE[0]}`
export PIPELINE_PATH=${SCRIPT_PATH%/*}

### Call tool ------------------------------------------------------------------

$PIPELINE_PATH/scripts/$TOOL $TOOL_ARG
