#!/usr/bin/env bash
# Checks every documentation word budget against the table in AGENTS.md, which
# is the one place the budgets are listed (#287). Run from the repository root:
#   bash tools/check-word-budgets.sh
# A table this script cannot read fails rather than passes: a check that
# silently matches nothing would enforce nothing.
set -euo pipefail
export LC_ALL=C.UTF-8

fail=0
error() { echo "::error::$*"; fail=1; }

header='| file | its job | budget |'
grep -qxF "$header" AGENTS.md || {
  echo "::error::AGENTS.md has no line '$header'. If the table was reformatted, update tools/check-word-budgets.sh to read it."
  exit 1
}

# The table: the header, its separator, then every |-line up to the first that is not one.
rows=$(awk -v h="$header" '$0 == h {t = 1; getline; next} t && /^\|/ {print; next} t {exit}' AGENTS.md)

declare -A budget
while IFS= read -r row; do
  file=$(sed -nE 's/^\| `([^`]+)` \|.*$/\1/p' <<<"$row")
  last=$(sed -E 's/^.*\|[[:space:]]*([^|]*[^|[:space:]])[[:space:]]*\|[[:space:]]*$/\1/' <<<"$row")
  if [[ -z $file ]]; then
    error "AGENTS.md budget row has no backticked file name in its first cell: $row"
  elif [[ $last =~ ^[0-9]+$ ]]; then
    budget[$file]=$last
  elif [[ $last != none* ]]; then
    error "AGENTS.md budget row for $file ends in '$last', which is neither a number nor 'none': $row"
  fi
done <<<"$rows"

if (( ${#budget[@]} == 0 )); then
  error "AGENTS.md's budget table has no budgeted row."
fi

# A doc that states its own budget must be in the table, or deleting its row
# would stop it being checked.
for doc in *.md; do
  if grep -qE 'Keep this file under [0-9]+ words' "$doc" && [[ -z ${budget[$doc]:-} ]]; then
    error "$doc says 'Keep this file under N words' but has no budget row in AGENTS.md."
  fi
done

printf '%-18s %6s %6s %6s\n' file words budget spare
for file in $(printf '%s\n' "${!budget[@]}" | sort); do
  if [[ ! -f $file ]]; then
    error "AGENTS.md budgets $file, which does not exist."
    continue
  fi
  limit=${budget[$file]}
  words=$(wc -w < "$file" | tr -d ' ')
  printf '%-18s %6s %6s %6s\n' "$file" "$words" "$limit" "$((limit - words))"
  stated=$(grep -oE 'Keep this file under [0-9]+ words' "$file" | grep -oE '[0-9]+' | head -n1 || true)
  if [[ -z $stated ]]; then
    error "$file has no 'Keep this file under N words' line; add one saying $limit."
  elif [[ $stated != "$limit" ]]; then
    error "$file says 'Keep this file under $stated words' but AGENTS.md budgets it at $limit."
  fi
  if (( words > limit )); then
    error "$file has $words words, $((words - limit)) over its budget of $limit. Something comes out before something goes in."
  fi
done

exit $fail
