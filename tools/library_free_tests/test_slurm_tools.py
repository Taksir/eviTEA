"""The helper of tools/slurm/record_pass.sbatch (call_args.py).   Run: python -B tools/library_free_tests/test_slurm_tools.py

S1  the recorded `evitea call` command comes back with --te-assets, --output and --work-dir replaced (both `--opt value` and `--opt=value`),
    a missing one added, everything else in its place.
S2  the command line of a command that is not `call` is refused.
S3  the script writes NUL-separated words, reading a manifest in the layout `evitea call` writes (architecture_id, schema_version,
    values.command_line); a top-level command_line is also read; a manifest with none is refused.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools" / "slurm"))
import call_args as C  # noqa: E402

FAIL, NCHK = [], [0]


def check(name, ok, detail=""):
    NCHK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


new = {"--te-assets": "/new/te", "--output": "/new/out", "--work-dir": "/new/work"}
old = ["evitea", "call", "--msbwt-index", "/m/idx", "--sample-id", "S", "--te-assets", "/old/te", "--output=/old/out", "--threads", "4",
       "--reference-assets", "/ref"]
got = C.call_args(old, new)
check("S1 replaced in place, missing one added", got == ["--msbwt-index", "/m/idx", "--sample-id", "S", "--te-assets", "/new/te", "--output", "/new/out",
                                                        "--threads", "4", "--reference-assets", "/ref", "--work-dir", "/new/work"], got)
try:
    C.call_args(["evitea", "library-free", "--output", "x"], new)
    check("S2 a command that is not call is refused", False)
except SystemExit:
    check("S2 a command that is not call is refused", True)
d = Path(tempfile.mkdtemp(prefix="callargs_"))
SCRIPT = str(HERE.parents[1] / "tools" / "slurm" / "call_args.py")


def run(document, name):
    (d / name).write_text(json.dumps(document))
    return subprocess.run([sys.executable, "-B", SCRIPT, str(d / name), "/new/te", "/new/out", "/new/work"], capture_output=True)


real = {"architecture_id": "evitea.run-manifest", "schema_version": "evitea.run-manifest.v1", "values": {"evitea_version": "x", "command_line": old}}
out = run(real, "real.json")
check("S3 NUL-separated words from a manifest in the written layout", out.returncode == 0 and out.stdout.endswith(b"\0")
      and out.stdout.decode().split("\0")[:-1] == got, out)
out = run({"command_line": old}, "flat.json")
check("S3 a top-level command_line is read too", out.returncode == 0 and out.stdout.decode().split("\0")[:-1] == got, out)
out = run({"architecture_id": "a", "schema_version": "b", "values": {}}, "none.json")
check("S3 a manifest with no command_line is refused", out.returncode != 0 and b"no command_line" in out.stderr, out)
print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
sys.exit(1 if FAIL else 0)
