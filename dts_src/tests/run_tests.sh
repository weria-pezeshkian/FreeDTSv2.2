#!/bin/bash
# Usage: run_tests.sh [test_name...]   (default: every tests/test_*.cpp)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/.."
BUILD="$HERE/build"
OMP="${OMP:-/opt/homebrew/opt/libomp}"
CXX="${CXX:-clang++}"
FLAGS="-O2 -std=c++17 -Xpreprocessor -fopenmp -I$OMP/include"
mkdir -p "$BUILD"

newest_header=$(ls -t "$SRC"/*.h "$SRC"/Registry/*.h 2>/dev/null | head -1)
objs=()
pids=()
for cpp in "$SRC"/*.cpp; do
    name=$(basename "$cpp" .cpp)
    [ "$name" = "FreeDTS" ] && continue
    obj="$BUILD/$name.o"
    objs+=("$obj")
    if [ ! -f "$obj" ] || [ "$cpp" -nt "$obj" ] || [ "$newest_header" -nt "$obj" ]; then
        echo "compiling $name"
        $CXX -c $FLAGS "$cpp" -o "$obj" &
        pids+=($!)
        while [ "$(jobs -r | wc -l)" -ge 8 ]; do sleep 0.2; done
    fi
done
for pid in "${pids[@]}"; do wait "$pid" || exit 1; done

tests=("$@")
if [ ${#tests[@]} -eq 0 ]; then
    for t in "$HERE"/test_*.cpp; do tests+=("$(basename "$t" .cpp)"); done
fi

status=0
for t in "${tests[@]}"; do
    $CXX $FLAGS "$HERE/$t.cpp" "${objs[@]}" -L"$OMP/lib" -lomp -o "$BUILD/$t"
    tmp=$(mktemp -d)
    cp "$HERE"/data/* "$tmp"/
    if (cd "$tmp" && "$BUILD/$t"); then echo "PASS $t"; else echo "FAIL $t"; status=1; fi
    rm -rf "$tmp"
done
exit $status
