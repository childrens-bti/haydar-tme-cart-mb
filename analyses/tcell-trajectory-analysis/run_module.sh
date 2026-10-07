#!/bin/bash

set -e
set -o pipefail

module_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$module_dir"

# The refined T-cell object is produced by downstream-analyses script 06.
# Downloading the v5 data release is an alternative when that output is absent.

default_scripts=(
  "01-tcell-trajectory.Rmd"
  "02-tcell-cd4-cd8-annotation.Rmd"
  "03-tcell-cd4-cd8-slingshot.Rmd"
)

compartment="both"
scripts=()

while [ "$#" -gt 0 ]; do
  case "$1" in
    --compartment)
      if [ "$#" -lt 2 ]; then
        echo "--compartment requires CD4, CD8, or both" >&2
        exit 1
      fi
      compartment="$2"
      shift 2
      ;;
    --compartment=*)
      compartment="${1#*=}"
      shift
      ;;
    *)
      scripts+=("$1")
      shift
      ;;
  esac
done

case "$compartment" in
  CD4|cd4|CD4-like|cd4-like)
    compartment="CD4-like"
    ;;
  CD8|cd8|CD8-like|cd8-like)
    compartment="CD8-like"
    ;;
  both|Both|BOTH)
    compartment="both"
    ;;
  *)
    echo "Invalid compartment: $compartment. Use CD4, CD8, or both." >&2
    exit 1
    ;;
esac

if [ "${#scripts[@]}" -eq 0 ]; then
  scripts=("${default_scripts[@]}")
fi

for script in "${scripts[@]}"; do
  if [ ! -f "$script" ]; then
    echo "Trajectory script not found: $script" >&2
    exit 1
  fi

  echo "Rendering $script"
  script_name="$(basename "$script")"
  if [ "$script_name" = "04-tcell-cd4-cd8-tradeseq.Rmd" ]; then
    if [ "$compartment" = "CD4-like" ]; then
      report_file="${script_name%%-*}-tcell-cd4-tradeseq.html"
    elif [ "$compartment" = "CD8-like" ]; then
      report_file="${script_name%%-*}-tcell-cd8-tradeseq.html"
    else
      report_file="${script_name%.Rmd}.html"
    fi

    Rscript -e "rmarkdown::render('$script', params = list(compartment = '$compartment'), output_file = '$report_file')"
  else
    if [ "$compartment" != "both" ]; then
      echo "--compartment is supported only for analysis 04" >&2
      exit 1
    fi
    Rscript -e "rmarkdown::render('$script')"
  fi
done
