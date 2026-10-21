#!/usr/bin/env python3
"""Documentation coverage and error check for TROLL (run by `make docs-check`).

Reads Doxygen's XML output and warning log. Prints the share of documented
entities by kind and by file. Exits with status 1 if Doxygen reported a
documentation *error* (malformed command, comment attached to nothing,
unknown parameter name...), or if coverage falls below --min.
Missing documentation (undocumented entities, missing @param) is reported
but does not fail the check.
"""
import argparse, collections, glob, os, re, sys
import xml.etree.ElementTree as ET

# Doxygen warnings about *missing* documentation: reported, not errors
MISSING = re.compile(r"is not documented|are not \(all\) documented|are not documented|return type of member .* is not documented")
# Doxygen 1.9.x misparses the (valid) nested Markdown lists of NEWS.md; GitHub, pandoc and
# Doxygen >= 1.14 render them correctly. This parser message is ignored for NEWS.md only.
KNOWN_PARSER_ISSUES = re.compile(r"NEWS\.md:\d+: warning: Invalid list item found")


def text(e):
    return "".join(e.itertext()).strip() if e is not None else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("docs_dir", help="Doxygen output directory (contains xml/ and doxygen-warnings.txt)")
    ap.add_argument("--min", type=float, default=0.0, help="minimum overall coverage in percent (default 0: report only)")
    a = ap.parse_args()

    kinds = collections.defaultdict(lambda: [0, 0])
    files = collections.defaultdict(lambda: [0, 0])
    params = [0, 0]  # functions with parameters, of which all documented with @param
    seen = set()
    for f in glob.glob(os.path.join(a.docs_dir, "xml", "*.xml")):
        if f.endswith("index.xml"):
            continue
        for m in ET.parse(f).getroot().iter("memberdef"):
            kind = m.get("kind")
            if kind not in ("function", "variable", "define", "typedef", "enum") or m.get("id") in seen:
                continue
            seen.add(m.get("id"))
            loc = m.find("location")
            path = loc.get("file") if loc is not None else "?"
            documented = bool(text(m.find("briefdescription")) or text(m.find("detaileddescription"))
                              or text(m.find("inbodydescription")))
            kinds[kind][0] += 1
            kinds[kind][1] += documented
            files[path][0] += 1
            files[path][1] += documented
            if kind == "function":
                # declarations in headers often leave parameters unnamed: use the definition name too
                plist = [p for p in m.findall("param") if text(p.find("declname")) or text(p.find("defname"))]
                if plist:
                    params[0] += 1
                    dd = m.find("detaileddescription")
                    if dd is not None and len(list(dd.iter("parameteritem"))) >= len(plist):
                        params[1] += 1
    total = sum(v[0] for v in kinds.values())
    documented = sum(v[1] for v in kinds.values())
    if total == 0:
        print("No Doxygen XML found in", a.docs_dir)
        return 1
    print("Documentation coverage (entities with a description):")
    for k, (t, d) in sorted(kinds.items(), key=lambda x: -x[1][0]):
        print(f"  {k:9s} {d:4d} / {t:<4d} {100 * d / t:5.1f}%")
    print(f"  {'all':9s} {documented:4d} / {total:<4d} {100 * documented / total:5.1f}%")
    print(f"  functions with parameters documented with @param: {params[1]} / {params[0]}")
    print("Lowest coverage by file:")
    for path, (t, d) in sorted(files.items(), key=lambda x: x[1][1] / x[1][0])[:8]:
        print(f"  {path:45s} {d:4d} / {t:<4d} {100 * d / t:5.1f}%")

    errors = []
    log = os.path.join(a.docs_dir, "doxygen-warnings.txt")
    if os.path.exists(log):
        for line in open(log, errors="replace"):
            if "warning:" in line and not MISSING.search(line) and not KNOWN_PARSER_ISSUES.search(line):
                errors.append(line.rstrip())
    print(f"Doxygen documentation errors: {len(errors)}")
    for e in errors:
        print("  " + e)

    status = 0
    if errors:
        status = 1
    if 100 * documented / total < a.min:
        print(f"Coverage {100 * documented / total:.1f}% is below the minimum {a.min}%")
        status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
