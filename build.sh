#!/usr/bin/env bash
set -e
root=$(cd "$(dirname "$0")" && pwd)

targets="poisson_p1 poisson_p2 poisson_p3 advec_p1 advec_p2 advec_p3 convecdiff_supg laplace_p1 laplace_p2 laplace_test mesh_p2 mesh_p3 fields_advection fields_diffusion grid_generator"
sources() {
    local degree
    local elements=("$root/src/triangle_quadrature.f90" "$root/src/triangle_elements.f90")
    case "$1" in
        poisson_p[123])
            files=("${elements[@]}" "$root/src/poisson_solver.f90" "$root/poisson_isopara/$1/pgm/poisson.f90") ;;
        advec_p[123])
            degree=${1#advec_}
            files=("${elements[@]}" "$root/src/supg_numerics.f90" "$root/src/supg_io.f90" "$root/src/supg_advection.f90"
                   "$root/adv_supg/$1/pgm/advection_$degree.f90") ;;
        convecdiff_supg)
            files=("${elements[@]}" "$root/src/supg_numerics.f90" "$root/src/supg_io.f90" "$root/convecdiff_supg/pgm/addif.f90") ;;
        laplace_p1) files=("$root/laplace/pgm/laplace_boundary.f90" "$root/laplace/pgm/potential_FEM_new.f90") ;;
        laplace_p2) files=("$root/laplace/pgm/laplace_boundary.f90" "$root/laplace/pgm/Potential_FEM_quadratic.f90") ;;
        laplace_test) files=("$root/laplace/pgm/laplace_boundary.f90" "$root/laplace/pgm/potential_test.f90") ;;
        mesh_p[23])
            degree=${1#mesh_}
            files=("$root/ｍeshconv_p1_to_$degree/mesh_p1_to_$degree.f90") ;;
        test_triangle) files=("${elements[@]}" "$root/adv_supg/tests/test_triangle.f90") ;;
        fields_advection) files=("$root/adv_supg/advec_p1/mkdat/makdata.f90") ;;
        fields_diffusion) files=("$root/convecdiff_supg/mkdat/makdata.f90") ;;
        grid_generator) files=("$root/tools/gr90.f90") ;;
        *) echo "Unknown target: $1. Available: $targets" >&2; exit 2 ;;
    esac
}

if [[ ${1:-} == --list ]]; then printf '%s\n' $targets; exit 0; fi
if (( $# > 2 )); then echo "Usage: $0 [target|all] [build-directory]" >&2; exit 2; fi
target=${1:-all}
dest=${2:-"$root/build"}
mkdir -p "$dest"
dest=$(cd "$dest" && pwd)
temp=$(mktemp -d "$dest/.compile.XXXXXXXX")
trap 'rm -rf "$temp"' EXIT
cd "$temp"
fc=${FC:-gfortran}
read -r -a flags <<< "${FFLAGS:--O2 -ffree-line-length-none}"
read -r -a ldflags <<< "${LDFLAGS:-}"
sdk=()
if [[ $(uname -s) == Darwin ]]; then
    sdkroot=${SDKROOT:-/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk}
    if [[ ! -d "$sdkroot" ]]; then sdkroot=$(xcrun --sdk macosx --show-sdk-path); fi
    sdk=(-isysroot "$sdkroot")
fi
if [[ $target == all ]]; then read -r -a selected <<< "$targets"; else selected=("$target"); fi
for name in "${selected[@]}"; do
    sources "$name"
    mkdir "$name"
    compile=("$fc" "${flags[@]}" "${sdk[@]}" "${files[@]}" "${ldflags[@]}" -o solver)
    printf '%q ' "${compile[@]}"
    printf '\n'
    (cd "$name" && "${compile[@]}")
    mv "$name/solver" "$dest/$name"
    printf 'Built %s\n' "$dest/$name"
done
