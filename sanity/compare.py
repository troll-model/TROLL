#!/usr/bin/env python3
"""Compare two directories of TROLL outputs.

Default mode is bit-exact (every line identical), which is the right check for
pure refactoring. With --rtol/--atol, numeric fields are compared with a
tolerance, which is the right check after deliberate numerical changes
(compiler upgrade, float -> double, reordered sums...).

Lines matching any --ignore regex are dropped before comparison
(by default: timing, and the random-seed and version lines, which older references lack).

Exit status: 0 if everything matches, 1 otherwise.
"""
import argparse
import math
import os
import re
import sys

DEFAULT_IGNORE = [r"Average computation time", r"^Random seed :", r"^TROLL version :"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reference", help="directory with reference outputs")
    p.add_argument("candidate", help="directory with outputs to check")
    p.add_argument("--rtol", type=float, default=0.0, help="relative tolerance for numeric fields (default 0 = exact)")
    p.add_argument("--atol", type=float, default=0.0, help="absolute tolerance for numeric fields (default 0)")
    p.add_argument("--ignore", action="append", default=None,
                   help="regex of lines to ignore (repeatable; replaces the defaults)")
    p.add_argument("--max-report", type=int, default=5, help="differing lines to print per file")
    p.add_argument("--include", default=None,
                   help="regex: compare only files whose name matches (e.g. stand-level aggregates)")
    return p.parse_args()


def load(path, ignore):
    with open(path, "r", errors="replace") as fh:
        return [line.rstrip("\n") for line in fh if not any(r.search(line) for r in ignore)]


def to_float(token):
    try:
        return float(token)
    except ValueError:
        return None


def fields_close(a, b, rtol, atol):
    """Return (equal, relative_difference) for two tokens."""
    if a == b:
        return True, 0.0
    fa, fb = to_float(a), to_float(b)
    if fa is None or fb is None:
        return False, math.inf
    if math.isnan(fa) and math.isnan(fb):
        return True, 0.0
    diff = abs(fa - fb)
    scale = max(abs(fa), abs(fb))
    rel = diff / scale if scale > 0 else 0.0
    return diff <= atol + rtol * scale, rel


def compare_file(ref_path, cand_path, rtol, atol, ignore, max_report):
    ref, cand = load(ref_path, ignore), load(cand_path, ignore)
    problems = []
    max_rel = 0.0
    if len(ref) != len(cand):
        problems.append(f"line count differs: reference {len(ref)}, candidate {len(cand)}")
    exact = rtol == 0.0 and atol == 0.0
    for n, (lr, lc) in enumerate(zip(ref, cand), start=1):
        if lr == lc:
            continue
        if exact:
            ok = False
        else:
            tr, tc = lr.split(), lc.split()
            ok = len(tr) == len(tc)
            if ok:
                for a, b in zip(tr, tc):
                    equal, rel = fields_close(a, b, rtol, atol)
                    max_rel = max(max_rel, rel)
                    ok = ok and equal
        if not ok:
            problems.append(f"line {n}:\n      ref: {lr[:200]}\n      new: {lc[:200]}")
            if len(problems) >= max_report:
                problems.append("...")
                break
    return problems, max_rel


def main():
    args = parse_args()
    ignore = [re.compile(r) for r in (args.ignore if args.ignore is not None else DEFAULT_IGNORE)]
    ref_files = set(os.listdir(args.reference))
    cand_files = set(os.listdir(args.candidate))
    if args.include:
        keep = re.compile(args.include)
        ref_files = {f for f in ref_files if keep.search(f)}
        cand_files = {f for f in cand_files if keep.search(f)}
        if not ref_files:
            print(f"No reference file matches --include {args.include!r}")
            return 1
    failed = 0
    for name in sorted(ref_files - cand_files):
        print(f"{name} [MISSING in candidate]")
        failed += 1
    for name in sorted(cand_files - ref_files):
        print(f"{name} [MISSING in reference]")
        failed += 1
    common = sorted(ref_files & cand_files)
    for name in common:
        problems, max_rel = compare_file(os.path.join(args.reference, name), os.path.join(args.candidate, name),
                                         args.rtol, args.atol, ignore, args.max_report)
        note = f" (max rel diff {max_rel:.3g})" if max_rel > 0 and max_rel != math.inf else ""
        if problems:
            failed += 1
            print(f"{name} [DIFF]{note}")
            for pr in problems:
                print(f"    {pr}")
        else:
            print(f"{name} [OK]{note}")
    print("-----------------------------------")
    mode = "exact" if args.rtol == 0 and args.atol == 0 else f"rtol={args.rtol:g}, atol={args.atol:g}"
    if failed:
        print(f"{failed} of {len(common)} files differ ({mode}).")
        return 1
    print(f"All {len(common)} files match ({mode}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
