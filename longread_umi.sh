#!/bin/bash
# DESCRIPTION
#    longread_umi: pipelines and tools for longread UMI processing.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#

### Error handling ------------------------------------------------------------
#set -euo pipefail -o noclobber

### List tools ----------------------------------------------------------------
export LONGREAD_UMI_PATH="$(dirname "$(readlink -f "$0")")"

# Format list of scripts
SCRIPT_LIST=$(
  find \
    $LONGREAD_UMI_PATH/scripts/ \
    -name '*.sh' \
    -print |\
  sed \
    -e '/install_conda.sh/d' \
    -e '/install_dependencies.sh/d' \
    -e '/install_conda.sh/d' \
    -e '/dependencies.sh/d' \
    -e 's|^.*/||' \
    -e 's/.sh$//' |\
  sort
  )

# Call -h for all scripts
DOCS=""
for SCRIPT in $SCRIPT_LIST; do
DOCS="${DOCS}
\`\`\`
$($LONGREAD_UMI_PATH/scripts/${SCRIPT}.sh -h | sed -e "s|$LONGREAD_UMI_PATH|longread_umi|")
\`\`\`
  
  
"
done

# List piplines
PIPES=$(
  echo "$DOCS" |\
  gawk '
    $0 ~ /--.*pipeline/{
	    gsub("-- longread_umi", "  ", $0)
	    print $0
	}' |\
  column -t -s:
)

# List tools
TOOLS=$(
  echo "$DOCS" |\
  gawk '
    $0 ~ /--/ && $0 !~ /pipeline/{
      gsub("-- longread_umi", "  ", $0)
      print $0
    }
  ' |\
  column -t -s:
)

# Generate docs
if [ "$1" == "compile_docs" ]; then
  # Add docs for initiation script
  DOCS="
\`\`\`
$($LONGREAD_UMI_PATH/longread_umi.sh -h)
\`\`\`
  
  
${DOCS}
"
  
  # Define location in README.md
  LEAD='^## Usage$'
  TAIL='^## License$'
  
  echo "$DOCS" > $LONGREAD_UMI_PATH/docs/tmp.md
  sed -i \
  -e "/$LEAD/,/$TAIL/{ /$LEAD/{p; r $LONGREAD_UMI_PATH/docs/tmp.md
        }; /$TAIL/p; d }"  $LONGREAD_UMI_PATH/README.md
  rm $LONGREAD_UMI_PATH/docs/tmp.md
  exit 1
fi

### Description ----------------------------------------------------------------

USAGE="
-- longread_umi: pipelines and tools for longread UMI processing.

usage: $(basename "$0" .sh) [-h -v] ( name ...)

where:
    -h   Show this help text.
    -v   Show git version.
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
while getopts ':hzv' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    v) echo "$(git --git-dir ${LONGREAD_UMI_PATH}/.git describe --tag)"; exit 1;;
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