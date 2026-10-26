"""Print the `evitea call` arguments recorded in run_manifest.json, with --te-assets, --output and --work-dir
replaced, as NUL-separated words (for `mapfile -d ''`).   usage: call_args.py <run_manifest.json> <te-assets> <output> <work-dir>"""
import json
import sys

SET = ("--te-assets", "--output", "--work-dir")


def recorded_command_line(document):
    """Return the command line of a run_manifest.json (under "values" in a real manifest: architecture_id, schema_version, values)."""
    for holder in (document.get("values"), document):
        if isinstance(holder, dict) and isinstance(holder.get("command_line"), list):
            return holder["command_line"]
    raise SystemExit("run_manifest.json holds no command_line")


def call_args(command_line, values):
    words = list(command_line)
    if "call" not in words:
        raise SystemExit("no `call` in the recorded command line")
    words = words[words.index("call") + 1:]
    out, seen, i = [], set(), 0
    while i < len(words):
        w = words[i]
        name = w.split("=", 1)[0]
        if name in values:
            seen.add(name)
            out += [name, values[name]]
            i += 1 if "=" in w else 2
        else:
            out.append(w)
            i += 1
    for name in SET:
        if name not in seen:
            out += [name, values[name]]
    return out


if __name__ == "__main__":
    manifest, te, output, work = sys.argv[1:5]
    with open(manifest, encoding="utf-8") as handle:
        command_line = recorded_command_line(json.load(handle))
    sys.stdout.write("".join(w + "\0" for w in call_args(command_line, dict(zip(SET, (te, output, work))))))
