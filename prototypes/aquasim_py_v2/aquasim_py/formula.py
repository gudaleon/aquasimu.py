"""AQUASIM formula/variable evaluator (port of formvar.cpp, othervar.cpp, aqvar.cpp subset).

Supports CONSTVAR, FORMVAR (+-*/^, functions, comparisons, if-then-else-endif),
REALLISTVAR (LINEAR interpolation, constant extrapolation), VARLISTVAR (LINEAR),
PROGVAR (t, z, theta(EPSL), calcnum, zone), STATEVAR. Values may be numpy arrays
(one entry per grid point) so a whole column is evaluated in one pass.
"""
import re
import numpy as np

_FUNCS = {"exp": np.exp, "log": np.log, "ln": np.log, "log10": np.log10, "sqrt": np.sqrt,
          "sin": np.sin, "cos": np.cos, "tan": np.tan, "abs": np.abs, "min": np.minimum,
          "max": np.maximum, "_where": np.where, "np": np}


def translate(expr: str) -> str:
    s = " " + expr.strip() + " "
    s = s.replace("^", "**")
    s = re.sub(r"\band\b", "&", s); s = re.sub(r"\bor\b", "|", s); s = re.sub(r"\bnot\b", "~", s)
    s = s.replace("!=", "!=").replace("<>", "!=")
    pat = re.compile(r"\bif\b((?:(?!\bif\b).)*?)\bthen\b((?:(?!\bif\b).)*?)\belse\b((?:(?!\bif\b).)*?)\bendif\b", re.S)
    while True:
        s2 = pat.sub(lambda m: "_where((%s),(%s),(%s))" % (m.group(1), m.group(2), m.group(3)), s)
        if s2 == s:
            break
        s = s2
    return s.strip()


class VarSys:
    def __init__(self, records):
        """records: list of (type, fields) from the .aqu VARSYS section."""
        self.defs = {}
        for typ, f in records:
            self.defs[f[1]] = (typ, f)
        self.code = {}
        for name, (typ, f) in self.defs.items():
            if typ == "FORMVAR":
                self.code[name] = compile(translate(f[4]), name, "eval")
        self.prog = {}       # PROGVAR type -> name
        for name, (typ, f) in self.defs.items():
            if typ == "PROGVAR":
                self.prog[f[4]] = name
        self._exprcache = {}

    def compile_expr(self, expr):
        c = self._exprcache.get(expr)
        if c is None:
            c = self._exprcache[expr] = compile(translate(expr), expr, "eval")
        return c

    def env(self, state, prog):
        return _Env(self, state, prog)


class _Env(dict):
    def __init__(self, vs, state, prog):
        super().__init__()
        self.vs = vs
        self.update(state)
        for ptype, val in prog.items():
            if ptype in vs.prog:
                self[vs.prog[ptype]] = val

    def __missing__(self, name):
        d = self.vs.defs.get(name)
        if d is None:
            if name in _FUNCS:
                return _FUNCS[name]
            raise NameError(name)
        typ, f = d
        if typ == "CONSTVAR":
            v = float(f[4])
        elif typ == "FORMVAR":
            v = eval(self.vs.code[name], _FUNCS, self)
        elif typ == "REALLISTVAR":
            arg = self[f[4]]
            xy = np.array(f[14], float).reshape(-1, 2)
            v = np.interp(arg, xy[:, 0], xy[:, 1])
        elif typ == "VARLISTVAR":
            arg = self[f[4]]
            pairs = f[7]
            xs = np.array(pairs[0::2], float)
            vals = [self[n] for n in pairs[1::2]]
            arg_a = np.asarray(arg, float)
            if arg_a.ndim == 0 and float(arg) in xs:
                v = vals[list(xs).index(float(arg))]
            else:
                k = np.clip(np.searchsorted(xs, arg_a) - 1, 0, len(xs) - 2)
                w = (arg_a - xs[k]) / (xs[k + 1] - xs[k])
                w = np.clip(w, 0, 1)
                v = (1 - w) * np.choose(k, vals) + w * np.choose(k + 1, vals)
        elif typ == "STATEVAR":
            v = 0.0
        else:
            raise NameError("%s (%s)" % (name, typ))
        self[name] = v
        return v

    def ev(self, expr):
        return eval(self.vs.compile_expr(expr), _FUNCS, self)
