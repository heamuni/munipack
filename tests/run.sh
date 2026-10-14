#!/bin/bash
# Runs the tests tests/t_*.sh (or the ones given as arguments), each in its
# own empty directory, against the `munipack` on PATH. Exit status 0 when
# all pass. In the Docker build: docker build --target test .
#
#   tests/run.sh                 all tests
#   tests/run.sh t_cli t_phcal   selected tests
#   KEEP=1 tests/run.sh ...      keep the working directories

here=$(cd "$(dirname "$0")" && pwd)
cd "$here" || exit 1
if [ $# -gt 0 ]; then
    tests=()
    for t in "$@"; do tests+=("${t%.sh}.sh"); done
else
    tests=(t_*.sh)
fi

munipack --version | head -1
passed=0; failed=()
for t in "${tests[@]}"; do
    work=$(mktemp -d "${TMPDIR:-/tmp}/munipack-${t%.sh}.XXXX")
    start=$SECONDS; keep=$KEEP
    echo "=== ${t%.sh}"
    if (cd "$work" && bash -c ". '$here/lib.sh'; . '$here/$t'"); then
        echo "PASS ${t%.sh} ($((SECONDS - start)) s)"
        passed=$((passed + 1))
    else
        echo "FAIL ${t%.sh} ($((SECONDS - start)) s), files in $work"
        failed+=("${t%.sh}")
        keep=1
    fi
    [ -n "$keep" ] || rm -rf "$work"
done

echo "=== $passed passed, ${#failed[@]} failed ${failed[*]}"
[ ${#failed[@]} -eq 0 ]
