#!/bin/sh
# Checks the repository before a commit or release. Exits with an error if:
#   1. any file is larger than 50 MB (GitHub warns at 50 MB and refuses 100 MB)
#   2. any input or derived data file would be committed
#   3. any file contains a machine-specific path or personal identifier
#   4. any CSV header contains record-level PCCF+ or address identifiers
# Run from the repository root: sh tools/check_release.sh
cd "$(dirname "$0")/.." || exit 1
fail=0
# files that would be committed (all files when the folder is not yet a git repository)
if [ -d .git ]; then FILES=$(git ls-files --cached --others --exclude-standard); else FILES=$(find . -type f ! -path './.git/*' | sed 's|^\./||'); fi

big=$(echo "$FILES" | tr '\n' '\0' | xargs -0 -I{} find {} -size +50M 2>/dev/null)
[ -n "$big" ] && { echo "FAIL files over 50 MB:"; echo "$big"; fail=1; }

tracked=$(echo "$FILES" | grep '^data/' | grep -v -e '/README.md$' -e '^data/raw/IOR-2021/' -e '^data/modis_request.json$')
[ -n "$tracked" ] && { echo "FAIL data files would be committed:"; echo "$tracked"; fail=1; }

paths=$(echo "$FILES" | grep -v '^tools/check_release.sh$' | tr '\n' '\0' | xargs -0 grep -Il -e '/Users/' -e 'Dropbox' -e 'C:\\\\Users' 2>/dev/null)
[ -n "$paths" ] && { echo "FAIL machine-specific paths in:"; echo "$paths"; fail=1; }

# record identifiers or coordinates in any published CSV; aggregate tables are short
ID='(^|,)"?(ADDR_GUID|LOC_GUID|ref_DAUID|pccf_DAUID|pccf_PCODE|PCODE|BG_X|BG_Y|BF_REPPOINT_X|BF_REPPOINT_Y|pccf_LAT|pccf_LONG|pccf_DMT|pccf_SACtype)"?(,|$)'
for f in $(find outputs correction -name '*.csv' 2>/dev/null); do
  head -1 "$f" | grep -q -i -E "$ID" && { echo "FAIL record-level columns in $f"; fail=1; }
  [ "$(wc -l < "$f")" -gt 50000 ] && { echo "FAIL more than 50,000 rows (not an aggregate table?) in $f"; fail=1; }
done
# postal-code tables (NAR only): must not carry any PCCF+ field
for f in $(find outputs correction -name '*.csv.gz' 2>/dev/null); do
  gzip -dc "$f" | head -1 | grep -q -i -E 'pccf|DAUID|DMT|SAC|ADDR_GUID' && { echo "FAIL PCCF+ or record fields in $f"; fail=1; }
done

[ $fail -eq 0 ] && echo "check_release: OK" || echo "check_release: FAILED"
exit $fail
