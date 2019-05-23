#!/bin/bash
# DESCRIPTION
#    Run longread-UMI-pipeline on test data. 
#    
# IMPLEMENTATION
#    author	Søren Karst (sorenkarst@gmail.com)
#               Ryans Ziels (ziels@mail.ubc.ca)
#    version	0.1.0
#    license	GNU General Public License

../longread_UMI_pipeline.sh -d test_reads.fq -s 10 -c 30 -t 1
