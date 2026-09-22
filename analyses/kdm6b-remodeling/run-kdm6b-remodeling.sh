#!/bin/bash
set -e
set -o pipefail

# The KDM6B workflow owns repository-relative data, figure, and result paths
# in the haydar-ad-hoc source repository, so render it from its submodule root.
module_dir="external/haydar-ad-hoc/analyses/kdm6b-remodeling"

if [ ! -d "$module_dir" ]; then
  echo "Missing KDM6B source module at $module_dir." >&2
  echo "Run: git submodule update --init --recursive" >&2
  exit 1
fi

cd "$module_dir"
bash run_module.sh
