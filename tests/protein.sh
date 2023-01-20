#!/bin/bash

# Strict mode, exit on error
set -euo pipefail
IFS=$'\n\t'

tmpfile="$(mktemp /tmp/evcouplings.XXXXXX)"
tmpfile2="$(mktemp /tmp/evcouplings.XXXXXX)"
echo "$tmpfile2"
../bin/plmc -c "$tmpfile" --save-weights "$tmpfile2" -le 16.0 -lh 0.01 -m 100 -g -f DYR_ECOLI ../example/protein/DHFR.a2m

# Compare couplings lines except for last decimal digit
comp_value=0
diff --brief -L "$tmpfile" -L DHFR.couplings \
    <(while IFS= read -r line; do echo "${line%?}"; done < "$tmpfile") \
    <(while IFS= read -r line; do echo "${line%?}"; done < DHFR.couplings) ||
    comp_value=$?

weights_comp=0
diff --brief "$tmpfile2" DHFR_weights_raw.txt || weights_comp=$?

if [[ weights_comp -eq 1 ]]
then
    echo "Warning: Computed sequence weights differ from expected values"
    echo "tmpfile: $tmpfile2"
else
    echo "."
    rm "$tmpfile2"
fi

if [[ $comp_value -eq 1 ]]
then
    echo "Error: Couplings differ (using computed sequence weights)"
    echo "tmpfile: $tmpfile"
    exit 1
else
    echo "."
    rm "$tmpfile"
fi


# Test with weights (if true above)
echo "Testing weight loading"
tmpfile="$(mktemp /tmp/evcouplings.XXXXXX)"
tmpfile2="$(mktemp /tmp/evcouplings.XXXXXX)"
../bin/plmc -c "$tmpfile" -le 16.0 -lh 0.01 -m 100 -g -f DYR_ECOLI -w DHFR_weights_raw.txt --save-weights "$tmpfile2" ../example/protein/DHFR.a2m

# Compare couplings lines except for last decimal digit
comp_value=0
diff --brief -L "$tmpfile" -L DHFR.couplings \
    <(while IFS= read -r line; do echo "${line%?}"; done < "$tmpfile") \
    <(while IFS= read -r line; do echo "${line%?}"; done < DHFR.couplings) ||
    comp_value=$?

weights_comp=0
diff --brief "$tmpfile2" DHFR_weights_raw.txt || weights_comp=$?

if [[ weights_comp -eq 1 ]]
then
    echo "Warning: Loaded sequence weights differ from expected values"
    echo "tmpfile: $tmpfile2"
else
    echo "."
    rm "$tmpfile2"
fi

if [[ $comp_value -eq 1 ]]
then
    echo "Error: Couplings differ (using loaded sequence weights)"
    exit 1
else
    echo "."
fi
rm "$tmpfile"
