"""Stage the original mesh, boundary and field formats for one FEM case."""
from case_config import prepare


def stage(case):
    solver = case.value("solver")
    mesh = case.input_file("mesh", "mesh.mes")
    boundary = case.input_file("boundary", "boundary.bc")
    (case.run / "res").mkdir()
    if solver in ("poisson_p1", "poisson_p2", "poisson_p3"):
        source = case.value("source", "quadratic" if solver == "poisson_p2" else "sine")
        if source not in ("zero", "quadratic", "sine"):
            raise ValueError("source must be zero, quadratic or sine")
        (case.work / "poisson.nml").write_text(f"&poisson source='{source}' /\n")
        paths = [mesh, boundary, "../res/FEM.res", "../res/out.vtk",
                 "../res/outfile_x05.res", "../res/outfile_y05.res"]
        name = "file_FEM.dat"
    elif solver in ("laplace_p1", "laplace_p2", "laplace_test"):
        paths = [mesh, boundary, "../res/FEM.res", "../res/out2.res", "../res/integration.res"]
        name = "file_FEM.dat"
    elif solver in ("advec_p1", "advec_p2", "advec_p3", "convecdiff_supg"):
        velocity = case.input_file("velocity", "velocity.uv")
        initial = case.input_file("initial", "initial.ini")
        steps = case.number("steps", int)
        dt = case.number("dt", float)
        tolerance = case.number("tolerance", float, "1e-12")
        iterations = case.number("max_iterations", int, 1000)
        if steps < 1 or dt <= 0 or tolerance <= 0 or iterations < 1:
            raise ValueError("steps, dt, tolerance and max_iterations must be positive")
        records = [f"-1 {steps}", str(dt), f"{tolerance} {iterations}"]
        paths = ["inp.dat", mesh, velocity, boundary, initial, "../res/field.res", "../res/advection_"]
        if solver == "convecdiff_supg":
            alpha = case.number("diffusivity", float)
            if alpha < 0:
                raise ValueError("diffusivity must be nonnegative")
            records.append(str(alpha))
        else:
            sample_step = case.number("sample_step", int, 1000 if solver == "advec_p1" else 50)
            if sample_step < 0:
                raise ValueError("sample_step must be nonnegative")
            records.append(f"sample_step={sample_step}")
            paths += ["../res/profile.res", "../res/conservation.res"]
        (case.work / "inp.dat").write_text("\n".join(records) + "\n")
        name = "file.txt"
    else:
        raise ValueError(f"unknown solver: {solver}")
    (case.work / name).write_text("\n".join(paths) + "\n")


if __name__ == "__main__":
    prepare(stage)
