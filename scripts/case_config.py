"""Read one case and preserve the resolved configuration beside its outputs."""
import math
import shutil
import sys
import tempfile
from pathlib import Path


def assignments(lines):
    result = {}
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, sep, value = line.partition("=")
        key, value = key.strip(), value.strip()
        if not sep or not key.isidentifier() or not value or "\n" in value or "\r" in value:
            raise ValueError(f"expected key=value, got {line!r}")
        if key in result:
            raise ValueError(f"duplicate setting: {key}")
        result[key] = value
    return result


class Case:
    def __init__(self, filename, overrides):
        self.base = filename.resolve().parent
        self.options = assignments(filename.read_text(encoding="utf-8").splitlines())
        self.options.update(assignments(overrides))
        self.effective = {}
        repo = Path(__file__).resolve().parents[1]
        output = (self.base / self.options.pop("output_root", str(repo / "runs"))).resolve()
        output.mkdir(parents=True, exist_ok=True)
        self.effective["output_root"] = str(output)
        self.run = Path(tempfile.mkdtemp(prefix=filename.stem + "-", dir=output))
        self.work = self.run / "work"
        self.work.mkdir()
        (self.run / "inputs").mkdir()

    def value(self, key, default=None):
        value = self.options.pop(key, default)
        if value is None:
            raise ValueError(f"missing setting: {key}")
        self.effective[key] = str(value)
        return value

    def number(self, key, convert, default=None):
        value = str(self.value(key, default))
        try:
            number = convert(value.lower().replace("d", "e"))
        except ValueError:
            raise ValueError(f"invalid {key}: {value}") from None
        if not math.isfinite(number):
            raise ValueError(f"{key} must be finite")
        return number

    def input_file(self, key, name):
        source = self.base / self.value(key)
        target = Path("inputs") / name
        shutil.copyfile(source, self.run / target)
        self.effective[key] = str(target)
        return "../" + str(target)


def prepare(stage):
    case = None
    try:
        case = Case(Path(sys.argv[1]), sys.argv[2:])
        stage(case)
        if case.options:
            raise ValueError("unknown setting(s): " + ", ".join(sorted(case.options)))
        text = "".join(f"{k}={v}\n" for k, v in case.effective.items())
        (case.run / "effective.cfg").write_text(text, encoding="utf-8")
        print(case.run)
    except (OSError, ValueError) as exc:
        message = f"Case setup failed: {exc}"
        if case:
            (case.run / "run.log").write_text(message + "\n")
            (case.run / "exit-code").write_text("2\n")
            print(f"Run directory: {case.run}", file=sys.stderr)
        print(message, file=sys.stderr)
        sys.exit(2)
