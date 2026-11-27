#!/usr/bin/env bash
set -eu
root=$(cd "$(dirname "$0")" && pwd)
if [[ $# == 0 || $1 == --help ]]; then
    echo "Usage: $0 CASE.cfg [key=value ...]"
    exit $(( $# == 0 ? 2 : 0 ))
fi
run=$(python3 "$root/scripts/prepare_case.py" "$@")
printf 'Run directory: %s\n' "$run"
solver=$(sed -n 's/^solver=//p' "$run/effective.cfg")
finish() {
    code=$?
    printf '%s\n' "$code" > "$run/exit-code"
    if [[ $code != 0 ]]; then tail -n 20 "$run/run.log" >&2; fi
}
trap finish EXIT
{
    printf 'launch '
    LC_ALL=C printf '%q ' "$root/run_case.sh" "$@"
    printf '\nbuild '
    LC_ALL=C printf '%q ' "$root/build.sh" "$solver" "$run/bin"
    LC_ALL=C printf '\nexecute cd %q && %q\n' "$run/work" "$run/bin/$solver"
    LC_ALL=C printf 'FC=%q FFLAGS=%q LDFLAGS=%q SDKROOT=%q\n' "${FC:-gfortran}" "${FFLAGS:--O2 -ffree-line-length-none}" "${LDFLAGS:-}" "${SDKROOT:-}"
    "${FC:-gfortran}" --version | head -n 1
    if git -C "$root" rev-parse HEAD >/dev/null 2>&1; then
        git -C "$root" rev-parse HEAD
        git -C "$root" status --short --untracked-files=normal
    fi
} > "$run/command.log"
"$root/build.sh" "$solver" "$run/bin" > "$run/run.log" 2>&1
(cd "$run/work" && "$run/bin/$solver") >> "$run/run.log" 2>&1
