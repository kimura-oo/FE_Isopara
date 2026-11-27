"""Exercise polynomial quadrature, prescribed nonzero boundaries, and zero initial mass."""

import math
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]


def test_transport(binary, degree, value, work, boundary_only=False):
    # Supplied meshes have area four and both boundary and interior degrees of freedom.
    divisions = {1: 20, 2: 10, 3: 5}[degree]
    mesh = ROOT / f"adv_supg/advec_p{degree}/data/Corn/mes_{divisions}/cube_{divisions}.mes"
    nnode = int(mesh.read_text().splitlines()[0].split()[0])
    ids = [int(row.split()[1]) for row in mesh.with_suffix(".bc").read_text().splitlines()[1:] if row.strip()]
    case = work / f"p{degree}-{value}-{boundary_only}"
    case.mkdir()
    (case / "settings").write_text("-1 4\n1e-3\n1e-14 200\n")
    (case / "boundary").write_text(str(len(ids)) + "\n" +
                                  "".join(f"{i} {node} {value}\n" for i, node in enumerate(ids, 1)))
    initial = 0 if boundary_only else value
    (case / "initial").write_text("".join(f"{node} {initial}\n" for node in range(1, nnode+1)))
    (case / "file.txt").write_text("\n".join((
        "settings", str(mesh), str(mesh.with_suffix(".uv")), "boundary", "initial",
        "field", "step_", "profile", "mass")) + "\n")
    run = subprocess.run([binary, "--config", case / "file.txt", "--steps", "2"], cwd=work,
                         text=True, capture_output=True)
    if run.returncode:
        raise AssertionError(run.stdout + run.stderr)
    initial_mass = float(next(s.split()[-1] for s in run.stdout.splitlines() if "Initial mass integral" in s))
    if not boundary_only:
        assert math.isclose(initial_mass, 4*value, abs_tol=1e-11)
    history = (case / "mass").read_text().splitlines()
    if value == 0:
        assert history.pop(0).startswith("#")
    assert len(history) == 2
    for step, row in enumerate(history, 1):
        step_id, diagnostic = row.split()
        assert int(step_id) == step
        assert math.isfinite(float(diagnostic))
        if not boundary_only:
            assert math.isclose(float(diagnostic), 100 if value else 0, abs_tol=1e-11)
    rows = (case / "field").read_text().splitlines()
    assert len(rows) == 3*(nnode+1)
    fixed = set(ids)
    free_values = []
    for step in range(3):
        start = step*(nnode+1)
        assert int(rows[start].split()[0]) == step
        free = []
        for row in rows[start+1:start+nnode+1]:
            node, field = row.split()
            field = float(field)
            assert math.isfinite(field)
            if not boundary_only or int(node) in fixed:
                assert math.isclose(field, value, abs_tol=1e-11)
            else:
                free.append(field)
        free_values.append(free)
    if boundary_only:
        assert free_values[0] != free_values[1] != free_values[2]


def test_time_refinement(binary, work):
    # Fixed mesh isolates the time integrator when successive numerical solutions are compared.
    case = work / "refinement"
    case.mkdir()
    mesh = ROOT / "adv_supg/advec_p1/data/Corn/mes_20/cube_20.mes"
    nnode = int(mesh.read_text().splitlines()[0].split()[0])
    (case / "initial").write_text("".join(f"{i} 0\n" for i in range(1, nnode+1)))
    (case / "velocity").write_text("".join(f"{i} 0.5 0\n" for i in range(1, nnode+1)))
    (case / "file.txt").write_text("\n".join((
        "settings", str(mesh), "velocity", str(mesh.with_suffix(".bc")), "initial", "field", "step_")) + "\n")
    fields = []
    for steps in (2, 4, 8):
        (case / "settings").write_text(f"-1 {steps}\n{0.2/steps}\n1e-22 1000\n0.02\n")
        result = subprocess.run([binary], cwd=case, capture_output=True, text=True)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)
        rows = (case / "field").read_text().splitlines()[-nnode:]
        fields.append([float(row.split()[1]) for row in rows])
    changes = [math.sqrt(sum((a-b)**2 for a,b in zip(coarse, fine))/nnode)
               for coarse, fine in zip(fields, fields[1:])]
    assert 3.3 < changes[0]/changes[1] < 4.8, changes


def main():
    binaries = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        triangle_test = binaries / "test_triangle"
        if not triangle_test.is_file():
            subprocess.run([ROOT / "build.sh", "test_triangle", work], check=True)
            triangle_test = work / "test_triangle"
        subprocess.run([triangle_test], cwd=work, check=True)
        for degree in (1, 2, 3):
            for value in (0, 2):
                test_transport(binaries / f"advec_p{degree}", degree, value, work)
            test_transport(binaries / f"advec_p{degree}", degree, 2, work, boundary_only=True)
        test_time_refinement(binaries / "convecdiff_supg", work)
    print("Triangle rules, SUPG transport/boundary cases and time refinement passed.")


if __name__ == "__main__":
    main()
