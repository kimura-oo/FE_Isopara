"""Check all three Poisson executables against polynomial solutions on two triangles."""
from pathlib import Path
import subprocess
import sys
import tempfile


def make_mesh(degree, directory):
    # Use integer lattice coordinates to identify a shared node exactly.
    nodes = {}
    cells = []
    for corners in (((0, 0), (degree, 0), (degree, degree)),
                    ((0, 0), (degree, degree), (0, degree))):
        local = list(corners)
        for start, end in zip(corners, corners[1:] + corners[:1]):
            for k in range(1, degree):
                local.append(tuple(((degree-k)*a+k*b)//degree for a, b in zip(start, end)))
        if degree == 3:
            local.append(tuple(sum(axis)//3 for axis in zip(*corners)))
        cells.append([nodes.setdefault(key, len(nodes)+1) for key in local])
    rows = [f'{len(nodes)} {len(cells)}']
    # Reverse the file records while preserving IDs to test the mesh reader.
    rows.extend(f'{idx} {x/degree:.16g} {y/degree:.16g}'
                for (x, y), idx in reversed(list(nodes.items())))
    rows.extend(f'{idx} ' + ' '.join(map(str, cell)) for idx, cell in enumerate(cells, 1))
    (directory/'mesh').write_text('\n'.join(rows)+'\n')
    return nodes


def write_boundary(directory, fixed, flux_edges):
    rows = [str(len(fixed))]
    rows.extend(f'{k} {node} {value:.16g}' for k, (node, value) in enumerate(fixed, 1))
    if flux_edges:
        rows.append(f'{len(flux_edges)} NEUMANN')
        for k, (edge, flux) in enumerate(flux_edges, 1):
            rows.append(f'{k} ' + ' '.join(map(str, edge)) + f' {flux}')
    (directory/'boundary').write_text('\n'.join(rows)+'\n')


def solution(directory, binary, args, expected):
    subprocess.run([str(binary), *args], cwd=directory, check=True)
    lines = (directory/'solution.vtk').read_text().splitlines()
    start = lines.index('LOOKUP_TABLE default') + 1
    values = [float(line) for line in lines[start:start+len(expected)]]
    error = max(abs(values[idx-1]-value) for idx, value in expected.items())
    if not error < 1e-10:
        raise AssertionError(f'{binary.name}: max nodal error {error:g}')


def check(degree, binary, directory):
    (directory/'file_FEM.dat').write_text('mesh\nboundary\nsolution.res\nsolution.vtk\nx.res\ny.res\n')
    nodes = make_mesh(degree, directory)
    (directory/'poisson.nml').write_text("&poisson source='sine' /\n")
    # u=2+x. Prescribe the left boundary and integrate outward flux on
    # the right with its local orientation reversed. The two other sides have zero flux.
    fixed = [(idx, 2) for (x, _), idx in nodes.items() if x == 0]
    right = [nodes[degree, y] for y in [degree, 0, *range(degree-1, 0, -1)]]
    write_boundary(directory, fixed, [(right, 1)])
    solution(directory, binary, ['--source', 'zero'], {idx: 2+x/degree for (x, _), idx in nodes.items()})

    # P2 and P3 reproduce u=x^2+y^2 with -Delta u=-4 on affine triangles.
    if degree >= 2:
        expected = {idx: (x*x+y*y)/degree**2 for (x, y), idx in nodes.items()}
        fixed = [(idx, expected[idx]) for (x, y), idx in nodes.items()
                 if x in (0, degree) or y in (0, degree)]
        write_boundary(directory, fixed, [])
        (directory/'poisson.nml').write_text("&poisson source='quadratic' /\n")
        solution(directory, binary, [], expected)
    # The diagonal shared by both triangles is not an exterior Neumann edge.
    diagonal = [nodes[k, k] for k in [0, degree, *range(1,degree)]]
    write_boundary(directory, fixed, [(diagonal, 1)])
    rejected = subprocess.run([str(binary)], cwd=directory, capture_output=True, text=True)
    if rejected.returncode == 0 or 'exterior edge' not in rejected.stderr:
        raise AssertionError(f'{binary.name}: invalid internal boundary accepted or unexpected error: {rejected.stderr}')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        raise SystemExit('usage: check_poisson.py P1_BINARY P2_BINARY P3_BINARY')
    for degree, name in enumerate(sys.argv[1:], 1):
        with tempfile.TemporaryDirectory() as tmp:
            check(degree, Path(name).resolve(), Path(tmp))
    print('P1/P2/P3 Poisson patch and boundary tests passed')
