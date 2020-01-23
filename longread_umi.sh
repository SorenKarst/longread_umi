#!/bin/bash
# DESCRIPTION
#    longread_umi: pipelines and tools for longread UMI processing.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#

### List tools ----------------------------------------------------------------
export LONGREAD_UMI_PATH="$(dirname "$(readlink -f "$0")")"

SCRIPT_LIST=$(
  find \
    $LONGREAD_UMI_PATH/scripts/ \
	-name '*.sh' \
	-print |\
  sed \
    -e '/install/d' \
	-e '/dependen/d' \
	-e 's|^.*/||' \
	-e 's/.sh$//' |\
  sort
  )

PIPES=$(
  echo "$SCRIPT_LIST" |\
  gawk '
    $0 ~ /pipeline/{print "  " $0}
  '
)

TOOLS=$(
  echo "$SCRIPT_LIST" |\
  gawk '
    $0 !~ /pipeline/{ print "  " $0}
  '
)

# Generate docs

if [ "$1" == "compile_docs" ]; then
  for SCRIPT in $SCRIPT_LIST; do
    DOC_OUT="$LONGREAD_UMI_PATH/docs/USAGE_${SCRIPT}.md"
    echo '```' > $DOC_OUT
    $LONGREAD_UMI_PATH/scripts/${SCRIPT}.sh -h >> $DOC_OUT
    echo '```' >> $DOC_OUT
  done
  exit 1
fi

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [ name ...]

-- longread_umi: pipelines and tools for longread UMI processing.

where:
    -h   Show this help text.
    name Name of tool or pipeline.
    ...  Commands for tool or pipeline.

Pipelines:

$(echo "$PIPES")

Tools:

$(echo "$TOOLS")

For help with a specific tool or pipeline:
longread_umi <name> -h
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

### Call tool or pipeline ---------------------------------------------------------------

$LONGREAD_UMI_PATH/scripts/${TOOL}.sh $TOOL_ARG