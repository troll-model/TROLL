#!/usr/bin/env python3
"""Turn Doxygen's XML output into API reference pages of the TROLL website.

Called by site/build.sh after Doxygen has run. Writes site/vignettes/api-*.Rmd (generated,
ignored by git): one page per main part of the code, with the same look, navigation and search
as the rest of the site. Every documented function, variable and macro appears on one page;
the script fails if any would be left out.

Usage: doxygen_api.py XML_DIR SITE_DIR
"""
import os
import re
import sys
import xml.etree.ElementTree as ET

REPO = "https://github.com/troll-model/TROLL/blob/dev/"

CONTEXT_PARTS = ["Context", "FileIO", "InputBuffers", "OutputConfig", "ModelOptions", "Grid", "TimeState",
                 "Climate", "LookupTables", "Soil", "Diagnostics", "PointCloud", "RNGState", "SpeciesState",
                 "Intraspecific", "CrownGeometry", "SimParams", "SimFields"]
SMALL_TYPES = ["leafFluxes", "RespGPPResult", "ParamSpec"]
CONSTANT_FILES = ["include/constants.hpp", "include/troll_defines.hpp", "include/kernel_precision.hpp"]
FILE_ORDER = ["src/troll.cpp", "src/simulation.cpp", "src/memory.cpp", "src/input.cpp", "src/output.cpp",
              "src/crown.cpp", "src/lookup_tables.cpp", "src/internal/lookup_tables_cache.cpp",
              "src/internal/param_registry.cpp", "include/quad.hpp", "include/troll.hpp",
              "include/internal/param_registry.hpp", "include/internal/lookup_tables_cache.hpp"]


# ----------------------------------------------------------------------------- text helpers
def inline(e):
    """Doxygen inline markup -> Markdown text (one line)."""
    if e is None:
        return ""
    out = [e.text or ""]
    for c in e:
        tag = c.tag
        if tag == "computeroutput":
            out.append("`" + "".join(c.itertext()).strip() + "`")
        elif tag in ("bold",):
            out.append("**" + inline(c).strip() + "**")
        elif tag in ("emphasis",):
            out.append("*" + inline(c).strip() + "*")
        elif tag == "ulink":
            out.append(f"[{inline(c).strip()}]({c.get('url')})")
        elif tag == "linebreak":
            out.append(" ")
        elif tag in ("parameterlist", "simplesect", "itemizedlist", "orderedlist", "programlisting"):
            pass  # handled by block()
        else:
            out.append(inline(c))
        out.append(c.tail or "")
    return re.sub(r"\s+", " ", "".join(out))


def block(desc):
    """brief/detailed description -> (markdown paragraphs, params, returns)."""
    paras, params, returns = [], [], ""
    if desc is None:
        return paras, params, returns
    for para in desc.iter("para"):
        text = inline(para).strip()
        if text:
            paras.append(text)
        for lst in para.findall("itemizedlist"):
            items = ["- " + inline(li.find("para")).strip() for li in lst.findall("listitem") if li.find("para") is not None]
            if items:
                paras.append("\n".join(items))
        for pl in para.findall("parameterlist"):
            for item in pl.findall("parameteritem"):
                names = [inline(n).strip() for n in item.iter("parametername")]
                d = " ".join(inline(p).strip() for p in item.find("parameterdescription").findall("para"))
                params.append((", ".join(names), d))
        for ss in para.findall("simplesect"):
            if ss.get("kind") == "return":
                returns = " ".join(inline(p).strip() for p in ss.findall("para"))
    # remove paragraphs that only repeat nested list text
    return [p for p in paras if p], params, returns


def cell(s):
    return (s or "").replace("|", "\\|").replace("\n", " ").strip()


def location(m):
    loc = m.find("location")
    if loc is None:
        return "", ""
    f = loc.get("bodyfile") or loc.get("file") or ""
    line = loc.get("bodystart") or loc.get("line") or ""
    f = re.sub(r"^.*?/(src|include|tests)/", r"\1/", f)
    return f, line


def source_link(m):
    f, line = location(m)
    if not f:
        return ""
    anchor = f"#L{line}" if line and line != "-1" else ""
    return f"[`{f}`{':' + line if anchor else ''}]({REPO}{f}{anchor})"


def slug(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


# ----------------------------------------------------------------------------- renderers
def render_function(m, level="###"):
    name = m.findtext("name")
    sig = f"{(m.findtext('type') or '').strip()} {m.findtext('qualifiedname') or name}{m.findtext('argsstring') or ''}".strip()
    paras, params, returns = block(m.find("briefdescription"))
    p2, params2, returns2 = block(m.find("detaileddescription"))
    paras += p2; params += params2; returns = returns or returns2
    out = [f"{level} {name}", "", "```cpp", re.sub(r"\s+", " ", sig), "```", ""]
    out += [p + "\n" for p in paras]
    if params:
        out += ["| Parameter | Description |", "| --- | --- |"]
        out += [f"| `{cell(n)}` | {cell(d)} |" for n, d in params]
        out.append("")
    if returns:
        out += [f"**Returns:** {returns}", ""]
    link = source_link(m)
    if link:
        out += [f"Source: {link}", ""]
    return out


def member_table(members, with_type=True):
    rows = ["| Member | Type | Description |" if with_type else "| Member | Description |",
            "| --- | --- | --- |" if with_type else "| --- | --- |"]
    for m in members:
        paras, _, _ = block(m.find("briefdescription"))
        p2, _, _ = block(m.find("detaileddescription"))
        desc = " ".join(paras + p2)
        name = m.findtext("name") + (m.findtext("argsstring") or "")
        if with_type:
            rows.append(f"| `{cell(name)}` | `{cell((m.findtext('type') or '').strip())}` | {cell(desc)} |")
        else:
            rows.append(f"| `{cell(name)}` | {cell(desc)} |")
    return rows


def compound_doc(c):
    paras, _, _ = block(c.find("briefdescription"))
    p2, _, _ = block(c.find("detaileddescription"))
    return paras + p2


def write(site, name, title, intro, body):
    path = os.path.join(site, "vignettes", f"{name}.Rmd")
    text = f'---\ntitle: "{title}"\n---\n\n<!-- Generated from the Doxygen XML by site/doxygen_api.py: document the C++ code instead. -->\n\n'
    text += intro.strip() + "\n\n" + "\n".join(body).rstrip() + "\n"
    open(path, "w", encoding="utf-8").write(text)


# ----------------------------------------------------------------------------- main
def main():
    xml_dir, site = sys.argv[1], sys.argv[2]
    os.makedirs(os.path.join(site, "vignettes"), exist_ok=True)
    idx = ET.parse(os.path.join(xml_dir, "index.xml")).getroot()
    compounds = {}
    for c in idx.findall("compound"):
        f = os.path.join(xml_dir, c.get("refid") + ".xml")
        if os.path.exists(f):
            for cd in ET.parse(f).getroot().iter("compounddef"):
                compounds[(cd.get("kind"), cd.findtext("compoundname"))] = cd

    covered = set()

    def members(cd, kinds=("function", "variable", "define", "enum", "typedef")):
        ms = []
        for sec in cd.findall("sectiondef"):
            for m in sec.findall("memberdef"):
                if m.get("kind") in kinds:
                    ms.append(m)
        return ms

    # -- classes Tree and Species
    for cls, page, title, intro in [
        ("Tree", "api-tree", "API: Tree",
         "A `Tree` occupies one site of the grid; an empty site is a `Tree` with age 0. Methods are listed first, then data members. Units are given in the descriptions."),
        ("Species", "api-species", "API: Species",
         "A `Species` holds the traits read from the species file and per-species output sums."),
    ]:
        cd = compounds[("class", cls)]
        ms = members(cd)
        funcs = [m for m in ms if m.get("kind") == "function"]
        vars_ = [m for m in ms if m.get("kind") == "variable"]
        body = [p + "\n" for p in compound_doc(cd)]
        body += ["## Methods", ""]
        for m in funcs:
            body += render_function(m)
        body += ["## Data members", ""] + member_table(vars_) + [""]
        covered.update(m.get("id") for m in ms)
        write(site, page, title, intro, body)

    # -- Context and its parts
    body = []
    for name in CONTEXT_PARTS:
        cd = compounds.get(("struct", name))
        if cd is None:
            continue
        ms = members(cd)
        body += [f"## {name}", ""] + [p + "\n" for p in compound_doc(cd)] + member_table(ms) + [""]
        covered.update(m.get("id") for m in ms)
    write(site, "api-context", "API: Context (model state)",
          "`Context` holds all the state of a simulation and is passed by reference to the model functions (it replaces the global variables of TROLL 4.0). Each section below is one of its parts.",
          body)

    # -- small types
    body = []
    for name in SMALL_TYPES:
        cd = compounds.get(("struct", name))
        if cd is None:
            continue
        ms = members(cd)
        funcs = [m for m in ms if m.get("kind") == "function"]
        others = [m for m in ms if m.get("kind") != "function"]
        body += [f"## {name}", ""] + [p + "\n" for p in compound_doc(cd)]
        if others:
            body += member_table(others) + [""]
        for m in funcs:
            body += render_function(m, "###")
        covered.update(m.get("id") for m in ms)
    write(site, "api-types", "API: small types", "Small structures used to pass results between functions.", body)

    # -- free functions, grouped by the file where they are defined
    files = {k[1]: cd for k, cd in compounds.items() if k[0] == "file"}
    by_file = {}
    seen = set()
    for fname, cd in files.items():
        for m in members(cd, ("function",)):
            if m.get("id") in seen:
                continue
            seen.add(m.get("id"))
            f, _ = location(m)
            by_file.setdefault(f, []).append(m)
    body = []
    order = [f for f in FILE_ORDER if f in by_file] + sorted(f for f in by_file if f not in FILE_ORDER)
    for f in order:
        ms = sorted(by_file[f], key=lambda m: int((location(m)[1] or "0") if (location(m)[1] or "0").lstrip("-").isdigit() else 0))
        body += [f"## {f}", ""]
        for m in ms:
            body += render_function(m)
        covered.update(m.get("id") for m in ms)
    write(site, "api-functions", "API: functions",
          "Free functions (not members of a class), grouped by the file that defines them, in order of appearance. Tree methods are on the [Tree page](api-tree.html).",
          body)

    # -- constants, options, macros, global variables
    body = []
    for fname, cd in files.items():
        ms = [m for m in members(cd, ("define", "variable", "enum", "typedef")) if m.get("id") not in covered]
        if not ms:
            continue
        path = location(ms[0])[0] or fname
        body += [f"## {path}", ""] + [p + "\n" for p in compound_doc(cd)] + member_table(ms, with_type=False) + [""]
        covered.update(m.get("id") for m in ms)
    write(site, "api-constants", "API: constants, options and globals",
          "Physical constants, compile-time model options (see the README for how to switch them), precision macros of the leaf kernel, and the few remaining global variables.",
          body)

    # -- overview
    pages = [("api-tree", "Tree", "the tree class: growth, physiology, mortality, outputs"),
             ("api-species", "Species", "species traits and per-species sums"),
             ("api-context", "Context (model state)", "all state of a simulation, in 17 parts"),
             ("api-functions", "Functions", "free functions by source file: simulation loop, inputs, outputs, crowns, lookup tables"),
             ("api-constants", "Constants, options and globals", "physical constants, compile-time options, macros"),
             ("api-types", "Small types", "`leafFluxes`, `RespGPPResult`, `ParamSpec`")]
    body = ["| Page | Contents |", "| --- | --- |"] + [f"| [{t}]({n}.html) | {d} |" for n, t, d in pages]
    body += ["", "The [full Doxygen documentation](../api/index.html) adds the source code browser, file lists and cross-references."]
    write(site, "api-index", "API reference",
          "Reference of the C++ code, generated from its Doxygen comments at each build of the website. Every function lists its parameters with their units; every data member has a description.",
          body)

    # -- completeness check: every member of the XML must be on a page
    all_ids = set()
    for k, cd in compounds.items():
        if k[0] in ("class", "struct", "file"):
            all_ids.update(m.get("id") for m in members(cd))
    missing = all_ids - covered - seen
    if missing:
        print(f"ERROR: {len(missing)} members not placed on any page, e.g. {sorted(missing)[:5]}")
        return 1
    print(f"API pages written: 7 pages, {len(covered | seen)} members")
    return 0


if __name__ == "__main__":
    sys.exit(main())
