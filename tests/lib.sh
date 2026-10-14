# Helpers for the tests in tests/t_*.sh (sourced by run.sh).
#
# Every test runs in its own empty directory with `munipack` on PATH and
# stops at the first failed check. Checks print one line each, so the log
# of a failed run shows what passed before it.

set -e -o pipefail

fail() { echo "FAIL: $*" >&2; exit 1; }
ok() { echo "  ok: $*"; }

# check NAME VALUE LO HI: LO <= VALUE <= HI (numbers, any format awk reads)
check() {
    awk -v v="$2" -v lo="$3" -v hi="$4" \
        'BEGIN { exit !(v == v + 0 && lo <= v + 0 && v + 0 <= hi) }' \
        || fail "$1 = $2, expected $3 .. $4"
    ok "$1 = $2 ($3 .. $4)"
}

# expect_fail NAME CMD...: the command exits with status 1 (an error, not a crash)
expect_fail() {
    local name=$1; shift
    local st=0
    "$@" > /dev/null 2>&1 || st=$?
    [ "$st" -eq 1 ] || fail "$name: exit status $st, expected 1"
    ok "$name: exit status 1"
}

# key FILE[EXT] KEY: value of a header keyword ('' when missing)
key() { munipack fits -K "$2" --value "$1" 2> /dev/null | awk "NR == 1"; }

# nrows FILE[EXT]: rows of a table
nrows() { key "$1" NAXIS2; }

# colnum FILE[EXT] NAME: 1-based index of a table column
colnum() {
    munipack fits -lh "$1" | awk -v n="$2" -F "[=/']" \
        '/^TTYPE/ { gsub(/ /, "", $3); if ($3 == n) { sub("TTYPE", "", $1); print $1 + 0; exit } }'
}

# column FILE[EXT] NAME...: the columns as whitespace-separated text
column() {
    local f=$1; shift
    local idx=()
    for c in "$@"; do
        local i; i=$(colnum "$f" "$c")
        [ -n "$i" ] || fail "no column $c in $f"
        idx+=("$i")
    done
    munipack fits -lt "$f" | awk -v cols="${idx[*]}" \
        'BEGIN { n = split(cols, c, " ") }
         { s = ""; for (k = 1; k <= n; k++) s = s (k > 1 ? " " : "") $(c[k]); print s }'
}

# imstats FILE [SECTION]: "mean std n" of image pixels (SECTION as x1:x2,y1:y2)
imstats() {
    local f=$1
    [ -n "$2" ] && f="$1[$2]"
    munipack fits -li "$f" | awk '{ s += $3; ss += $3 * $3; n++ }
        END { m = s / n; v = (ss - n * m * m) / (n - 1)
              printf "%.10g %.10g %d\n", m, sqrt(v > 0 ? v : 0), n }'
}

# imsum FILE: checksum of the pixel values (bit-identical images)
imsum() { munipack fits -li "$1" | md5sum | cut -d" " -f1; }

# mkcat FILE.fits N SEED [RA DEC RADIUS MAGMIN MAGMAX MINSEP]
# A synthetic catalogue: N stars at random positions within RADIUS [deg] of
# RA, DEC, uniform in V from MAGMIN to MAGMAX, B - V uniform in 0 .. 1.5,
# no two stars closer than MINSEP [deg]. Columns RAJ2000, DEJ2000 [deg],
# pmRA, pmDE [mas/yr, zero], Vmag, e_Vmag, Bmag, e_Bmag (errors 0.001 mag:
# the frames of artificial render the magnitudes exactly).
mkcat() {
    local out=$1 n=$2 seed=$3 ra=${4:-180} dec=${5:-0} rad=${6:-0.08}
    local m1=${7:-10} m2=${8:-15} sep=${9:-0.004}
    local lst=${out%.fits}.lst
    awk -v n="$n" -v seed="$seed" -v ra0="$ra" -v dec0="$dec" -v r="$rad" \
        -v m1="$m1" -v m2="$m2" -v sep="$sep" '
    # Park-Miller generator: reproducible (mawk ignores srand(seed))
    function rnd() { s = (s * 16807) % 2147483647; return s / 2147483647 }
    BEGIN {
        s = seed % 2147483646 + 1; pi = atan2(0, -1); cd = cos(dec0 * pi / 180)
        k = 0; tries = 0
        while (k < n && tries < 100 * n) {
            tries++
            dx = (2 * rnd() - 1) * r; dy = (2 * rnd() - 1) * r
            if (dx * dx + dy * dy > r * r) continue
            near = 0
            for (j = 1; j <= k; j++)
                if ((dx - x[j])^2 + (dy - y[j])^2 < sep * sep) { near = 1; break }
            if (near) continue
            k++; x[k] = dx; y[k] = dy
        }
        print "# BEGIN HDU 0"
        print "SIMPLE  =                    T"
        print "BITPIX  =                    8"
        print "NAXIS   =                    0"
        print "EXTEND  =                    T"
        print "END"
        print "# END HDU 0"
        print "# BEGIN HDU 1"
        print "XTENSION= '\''BINTABLE'\''"
        print "BITPIX  =                    8"
        print "NAXIS   =                    2"
        print "NAXIS1  =                   48"
        printf "NAXIS2  = %20d\n", k
        print "PCOUNT  =                    0"
        print "GCOUNT  =                    1"
        print "TFIELDS =                    8"
        print "TTYPE1  = '\''RAJ2000'\''"
        print "TFORM1  = '\''D'\''"
        print "TUNIT1  = '\''deg'\''"
        print "TTYPE2  = '\''DEJ2000'\''"
        print "TFORM2  = '\''D'\''"
        print "TUNIT2  = '\''deg'\''"
        print "TTYPE3  = '\''pmRA'\''"
        print "TFORM3  = '\''D'\''"
        print "TUNIT3  = '\''mas/yr'\''"
        print "TTYPE4  = '\''pmDE'\''"
        print "TFORM4  = '\''D'\''"
        print "TUNIT4  = '\''mas/yr'\''"
        print "TTYPE5  = '\''Vmag'\''"
        print "TFORM5  = '\''E'\''"
        print "TTYPE6  = '\''e_Vmag'\''"
        print "TFORM6  = '\''E'\''"
        print "TTYPE7  = '\''Bmag'\''"
        print "TFORM7  = '\''E'\''"
        print "TTYPE8  = '\''e_Bmag'\''"
        print "TFORM8  = '\''E'\''"
        print "EXTNAME = '\''SYNTHETIC'\''"
        print "END"
        for (j = 1; j <= k; j++) {
            v = m1 + (m2 - m1) * rnd(); bv = 1.5 * rnd()
            printf "%.10f %.10f 0 0 %.4f 0.001 %.4f 0.001\n", ra0 + x[j] / cd, dec0 + y[j], v, v + bv
        }
        print "# END HDU 1"
    }' > "$lst"
    munipack fits --restore "$lst" > /dev/null
    [ "$(nrows "$out[1]")" = "$n" ] || fail "mkcat: $(nrows "$out[1]") of $n stars"
}

# The columns of mkcat for artificial and astrometry.
CATCOLS="--col-ra RAJ2000 --col-dec DEJ2000 --col-pm-ra pmRA --col-pm-dec pmDE --col-mag Vmag"

# match CAT.fits[EXT] CALCOLS FILE[EXT] COLS TOL: join two tables by RA, Dec
# (the first two columns of each) within TOL [deg]; prints the remaining
# columns of the catalogue then of the table for every matched pair.
match() {
    local tol=$5
    { column "$1" $2 | sed 's/^/C /'; column "$3" $4 | sed 's/^/T /'; } | awk -v tol="$tol" '
        $1 == "C" { nc++; cra[nc] = $2; cde[nc] = $3; $1 = $2 = $3 = ""; crest[nc] = $0; next }
        { pi = atan2(0, -1); best = 0; bd = tol * tol
          for (i = 1; i <= nc; i++) {
              dx = ($2 - cra[i]) * cos($3 * pi / 180); dy = $3 - cde[i]
              d = dx * dx + dy * dy
              if (d < bd) { bd = d; best = i }
          }
          if (best) { $1 = $2 = $3 = ""; print crest[best], $0 } }' \
        | awk '{ $1 = $1; print }'
}

# stats < numbers: "median robust-sd n" (robust sd = 1.4826 MAD)
stats() {
    sort -g | awk '{ v[++n] = $1 }
        END { if (n == 0) { print "nan nan 0"; exit }
              m = (n % 2) ? v[(n + 1) / 2] : (v[n / 2] + v[n / 2 + 1]) / 2
              for (i = 1; i <= n; i++) a[i] = (v[i] > m) ? v[i] - m : m - v[i]
              k = asort_(a, n)
              printf "%.6g %.6g %d\n", m, 1.4826 * k, n }
        function asort_(a, n,   i, j, t) {
            for (i = 2; i <= n; i++) { t = a[i]; for (j = i - 1; j > 0 && a[j] > t; j--) a[j + 1] = a[j]; a[j + 1] = t }
            return (n % 2) ? a[(n + 1) / 2] : (a[n / 2] + a[n / 2 + 1]) / 2 }'
}
