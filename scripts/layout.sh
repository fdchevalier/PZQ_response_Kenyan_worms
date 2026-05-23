#!/bin/bash
# Title:
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-03-20
# Modified in:
# Licence : GPL v3



#======#
# Aims #
#======#

aim="Generate a table listing plate position and well labels."



#==========#
# Versions #
#==========#

# v0.0 - 2021-03-20: creation

version=$(grep -i -m 1 "version" "$0" | cut -d ":" -f 2 | sed "s/^ *//g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -g|--grp name -w|--wells ranges -o|--out path -h|--help

Aim: $aim

Version: $version

Options:
    -g, --grp       group name corresponding to each well range (e.g.: Treated Control)
                        note: values must be space separated
    -w, --wells     well ranges (e.g.: 14-23 26-35)
                        note: ranges must be space separated
                              wells outside ranges are labelled \"blank\"
    -o, --out       path of the output file
    -h, --help      this message
    "
}


# Info message
function info {
    if [[ -t 1 ]]
    then
        echo -e "\e[32mInfo:\e[00m $1"
    else
        echo -e "Info: $1"
    fi
}


# Warning message
function warning {
    if [[ -t 1 ]]
    then
        echo -e "\e[33mWarning:\e[00m $1"
    else
        echo -e "Warning: $1"
    fi
}


# Error message
## usage: error "message" exit_code
## exit code optional (no exit allowing downstream steps)
function error {
    if [[ -t 1 ]]
    then
        echo -e "\e[31mError:\e[00m $1"
    else
        echo -e "Error: $1"
    fi

    if [[ -n $2 ]]
    then
        exit $2
    fi
}


# Clean up function for trap command
## Usage: clean_up file1 file2 ...
function clean_up {
    rm -rf $@
    exit 1
}



#===========#
# Variables #
#===========#

# Usage on empty command
[[ $# == 0 ]] && usage && exit 0

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -g|--grp     ) groups=("$2") ; shift 2
                        while [[ ! -z "$1" && $(echo "$1"\ | grep -qv "^-" ; echo $?) == 0 ]]
                        do
                            groups+=("$1")
                            shift
                        done ;;
        -w|--wells   ) wells=("$2") ; shift 2
                        while [[ ! -z "$1" && $(echo "$1"\ | grep -qv "^-" ; echo $?) == 0 ]]
                        do
                            wells+=("$1")
                            shift
                        done ;;
        -o|--out    ) output="$2" ; shift 2 ;;
        -h|--help   ) usage ; exit 0 ;;
        *           ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check the existence of mandatory options
[[ -z "$output" ]] && error "Output file is mandatory. Exiting..." 1
[[ -s "$output" ]] && error "Output file already exist. Exiting..." 1

# Sanity checks
[[ ${#groups[@]} != ${#wells[@]} ]] && error "Groups and well ranges should be of the same length. Exiting..." 1



#============#
# Processing #
#============#

set -euo pipefail

# Well numbering
mywells=$(eval echo {A..H}{01..12} | tr " " "\n")

mygroups=$(printf 'blank %.0s' {1..96} | tr " " "\n")

for i in ${!wells[@]}
do
    group=${groups[$i]}
    well=$(sed "s/-/,/" <<< "${wells[$i]}")

    mygroups=$(sed "${well}s/blank/$group/" <<< "$mygroups")

done

# Final table
echo -e "Well\tType" > "$output"
paste <(echo "$mywells") <(echo "$mygroups") >> "$output"

exit 0
