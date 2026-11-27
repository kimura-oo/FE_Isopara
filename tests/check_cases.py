#!/usr/bin/env python3
"""Exercise the public FEM entry points and the numerical regression suites."""
import math
import os
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def check_converter(order, binary, tmp):
    source = ROOT / "poisson_isopara/poisson_p1/dat/gene/mesh/cube_test_2.mes"
    target = tmp / f"elevated_p{order}.mes"
    subprocess.run([binary, source, target], check=True)
    lines = target.read_text().splitlines()
    nnode, nelm = map(int, lines[0].split())
    assert (nnode, nelm) == ((25 if order == 2 else 49), 8)
    cells = [list(map(int, row.split()[1:])) for row in lines[1+nnode:]]
    shared = {}
    for cell in cells:
        for side in range(3):
            a, b = cell[side], cell[(side+1) % 3]
            start = 3+(order-1)*side
            ids = cell[start:start+order-1]
            if a > b:
                ids = list(reversed(ids))
            key = tuple(sorted((a,b)))
            assert shared.setdefault(key, ids) == ids
    assert len(shared) == 16


with tempfile.TemporaryDirectory(prefix="fem-cases-") as temp:
    tmp = Path(temp)
    binaries = tmp / "bin"
    binaries.mkdir()
    for config in sorted((ROOT / "cases").glob("*.cfg")):
        name = config.stem
        log = subprocess.check_output([str(ROOT / "run_case.sh"), str(config), "output_root=" + str(tmp)],
                                      cwd=tmp, text=True)
        run = Path(log.strip().split("Run directory: ")[-1])
        assert (run / "exit-code").read_text().strip() == "0"
        assert (run / "effective.cfg").stat().st_size and (run / "command.log").stat().st_size
        os.symlink(run / "bin" / name, binaries / name)
        if name.startswith("poisson") or name.startswith("laplace"):
            text = (run / "res/FEM.res").read_text()
            assert text.strip() and "nan" not in text.lower() and "inf" not in text.lower()
        if name.startswith("laplace"):
            mesh = (run / "inputs/mesh.mes").read_text().splitlines()
            count = int(mesh[0].split()[0])
            xs = [float(row.split()[1]) for row in mesh[1:count+1]]
            phi = [float(row.split()[1].replace("D","E")) for row in text.splitlines()[:count]]
            if name == "laplace_p1":
                expected = [10-x/2 for x in xs]
            elif name == "laplace_p2":
                expected = xs
            else:
                expected = [2*x-x*x/2 for x in xs]
            assert max(abs(u-v) for u,v in zip(phi,expected)) < 1e-6, name
        if name.startswith("advec") or name == "convecdiff_supg":
            vtk = sorted((run / "res").glob("*.vtk"))
            assert len(vtk) == 3  # initial state, then both completed steps
            assert "nan" not in vtk[-1].read_text().lower()
        if name == "poisson_p2":
            lines = (run / "res/out.vtk").read_text().splitlines()
            start = next(i for i, line in enumerate(lines) if line.startswith("POINTS ")) + 1
            points = [list(map(float, s.split())) for s in lines[start:start+25]]
            start = lines.index("LOOKUP_TABLE default") + 1
            values = [float(s) for s in lines[start:start+25]]
            error = max(abs(v - x*x - y*y) for v, (x,y,_) in zip(values, points))
            assert math.isfinite(error) and error < 1e-8, error
            print(f"P2 quadratic maximum nodal error: {error:.3g}")
        print("Passed", name, flush=True)
    subprocess.run(["python3", ROOT / "poisson_isopara/tests/check_poisson.py",
                    *(binaries / f"poisson_p{i}" for i in (1,2,3))], check=True)
    subprocess.run(["python3", ROOT / "adv_supg/tests/test_supg.py", binaries], check=True)
    for degree in (2,3):
        name = f"mesh_p{degree}"
        subprocess.run([ROOT / "build.sh", name, binaries], check=True)
        check_converter(degree, binaries / name, tmp)
    print("Passed P2 and P3 mesh connectivity")
