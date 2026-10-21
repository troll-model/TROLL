#!/usr/bin/env python3
"""Prepare the pkgdown sources of the TROLL website from the repository documents.

Called by site/build.sh. Generates (never edit these by hand; they are overwritten):
  site/index.md                 from README.md
  site/NEWS.md                  from NEWS.md, with one heading per version (pkgdown format)
  site/vignettes/<page>.Rmd     from docs/*.md, CONTRIBUTING.md, sanity/README.md,
                                TERMS_OF_REFERENCE.md
Links between these documents become links between site pages; links to any
other file of the repository point to that file on GitHub.
"""
import os
import re
import sys

REPO_URL = "https://github.com/troll-model/TROLL/blob/dev/"

# repository document -> page of the site, relative to the site root
PAGES = {
    "README.md": "index.html",
    "NEWS.md": "news/index.html",
    "docs/inputs.md": "articles/inputs.html",
    "docs/outputs.md": "articles/outputs.html",
    "docs/equations.md": "articles/equations.html",
    "docs/open-questions.md": "articles/open-questions.html",
    "CONTRIBUTING.md": "articles/contributing.html",
    "sanity/README.md": "articles/testing.html",
    "TERMS_OF_REFERENCE.md": "articles/terms-of-reference.html",
}
ARTICLES = {  # source -> (vignette name, title)
    "docs/inputs.md": ("inputs", "Input files"),
    "docs/outputs.md": ("outputs", "Output files"),
    "docs/equations.md": ("equations", "Equations and their implementation"),
    "docs/open-questions.md": ("open-questions", "Open questions for the authors"),
    "CONTRIBUTING.md": ("contributing", "Contributing to TROLL"),
    "sanity/README.md": ("testing", "Regression tests"),
    "TERMS_OF_REFERENCE.md": ("terms-of-reference", "Terms of Reference"),
}


def rewrite_links(text, source, depth):
    """Rewrite relative Markdown links of document `source`; depth = folder depth of the page."""
    up = "../" * depth
    base = os.path.dirname(source)

    def repl(m):
        label, target = m.group(1), m.group(2)
        if re.match(r"^(https?:|mailto:|#)", target):
            return m.group(0)
        path, _, anchor = target.partition("#")
        resolved = os.path.normpath(os.path.join(base, path)).replace("\\", "/")
        if resolved in PAGES:
            href = up + PAGES[resolved] + ("#" + anchor if anchor else "")
        else:
            href = REPO_URL + resolved + ("#" + anchor if anchor else "")
        return f"[{label}]({href})"

    return re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", repl, text)


def strip_title(text):
    """Remove the first level-1 heading (pkgdown shows the YAML title instead)."""
    return re.sub(r"\A\s*# [^\n]*\n", "", text, count=1)


def news_to_pkgdown(text):
    """'-   TROLL x.y.z: summary' bullets -> '# TROLL x.y.z' headings, items dedented."""
    out = []
    for line in text.split("\n"):
        m = re.match(r"^-\s+TROLL (\d+\.\d+\.\d+):\s*(.*)$", line)
        if m:
            out += ["", f"# TROLL {m.group(1)}", "", f"*{m.group(2).strip()}*", ""]
        elif re.match(r"^\s{2,}-\s+", line):
            out.append("- " + re.sub(r"^\s+-\s+", "", line))
        else:
            out.append(line)
    return "\n".join(out).strip() + "\n"


def main():
    root, site = sys.argv[1], sys.argv[2]
    read = lambda p: open(os.path.join(root, p), encoding="utf-8").read()
    os.makedirs(os.path.join(site, "vignettes"), exist_ok=True)

    with open(os.path.join(site, "index.md"), "w", encoding="utf-8") as f:
        f.write(rewrite_links(read("README.md"), "README.md", 0))
    with open(os.path.join(site, "NEWS.md"), "w", encoding="utf-8") as f:
        f.write(news_to_pkgdown(rewrite_links(read("NEWS.md"), "NEWS.md", 1)))
    for src, (name, title) in ARTICLES.items():
        body = strip_title(rewrite_links(read(src), src, 1))
        with open(os.path.join(site, "vignettes", name + ".Rmd"), "w", encoding="utf-8") as f:
            f.write(f'---\ntitle: "{title}"\n---\n\n<!-- Generated from {src} by site/prepare.py: edit the source instead. -->\n\n{body}')
    print(f"prepared index.md, NEWS.md and {len(ARTICLES)} articles")


if __name__ == "__main__":
    main()
