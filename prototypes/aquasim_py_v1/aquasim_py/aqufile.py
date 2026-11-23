"""Reader for AQUASIM .aqu files (nested {..} brace records, SaveStart/SaveEnd in fileio.cpp).

parse(text) -> nested Python lists of strings, mirroring the C++ Load() call order.
Each component class (e.g. SEDCOMP::Load) then maps positional fields to attributes.
"""
from __future__ import annotations


def parse(text: str):
    text = text.replace("\r", "").replace("\n", "")  # AQUASIM wraps lines at fixed width
    stack, cur, buf, i = [], [], [], 0
    depth0 = None
    for ch in text:
        if ch == "{":
            stack.append(cur); cur = []; buf = []
        elif ch == "}":
            if buf:
                cur.append("".join(buf)); buf = []
            elif cur == [] :
                pass
            node = cur
            cur = stack.pop()
            # a brace holding only one string is a scalar field
            cur.append(node[0] if len(node) == 1 and isinstance(node[0], str) else node)
        else:
            buf.append(ch)
    return cur


def load(path: str):
    with open(path, "r", encoding="latin-1") as f:
        return parse(f.read())


def find_section(tree, name: str):
    """Return the record following a tag like 'COMPSYS', searching recursively."""
    if isinstance(tree, list):
        for k, item in enumerate(tree):
            if item == name and k + 1 < len(tree):
                return tree[k + 1]
            r = find_section(item, name)
            if r is not None:
                return r
    return None
