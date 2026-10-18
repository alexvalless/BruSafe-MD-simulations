"""Execute a .ipynb in-process, without a Jupyter kernel.

Why this exists: the environment these notebooks were developed in could not bind a
socket, so no kernel could start. On a normal machine `jupyter lab` is the better
tool and this harness is only useful for unattended or batch runs.

Cells execute in __main__ so that functions defined in a cell stay picklable for
PyMC's multiprocessing workers. Every figure — whether closed by plt.show() or left
open at the end of a cell — is written to figs/<tag>_cellNN_k.png. Per-cell timing
and output go to handoff/log_<tag>.txt, and a machine-readable summary with the
error list to handoff/run_<tag>.json.

    cd notebooks
    python ../tools/exec_nb.py BruSafe_01_M1_expression_IPTG.ipynb run01

Recommended environment for reproducible timing:

    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    PYTENSOR_FLAGS="base_compiledir=$PWD/.pytensor_cache" python tools/exec_nb.py ...
"""
import sys as _sys, os as _os, io as _io, ast as _ast, json as _json
import time as _time, traceback as _tb, contextlib as _ctx
import matplotlib as _mpl
_mpl.use("Agg")
import matplotlib.pyplot as _plt
import nbformat as _nbf


def _run(_src, _tag):
    _os.makedirs("figs", exist_ok=True); _os.makedirs("handoff", exist_ok=True)
    _nb = _nbf.read(_src, as_version=4)
    _cells = [(i, c.source) for i, c in enumerate(_nb.cells) if c.cell_type == "code"]
    _st = {"cell": -1, "k": 0}; _figs = []

    def _save_open():
        for num in _plt.get_fignums():
            fig = _plt.figure(num)
            if not fig.get_axes():
                _plt.close(fig); continue
            _st["k"] += 1
            fn = f"figs/{_tag}_cell{_st['cell']:02d}_{_st['k']}.png"
            fig.savefig(fn, dpi=150, bbox_inches="tight"); _figs.append(fn); _plt.close(fig)

    _plt.show = lambda *a, **k: _save_open()

    def _display(*objs):
        for o in objs:
            print(o.to_string() if hasattr(o, "to_string") else repr(o))

    _main = vars(_sys.modules["__main__"])
    _main["display"] = _display; _main["get_ipython"] = lambda: None

    def _patch_smc():
        """PyMC's SMC spawns workers; fork keeps the __main__ namespace they need."""
        pm = _main.get("pm")
        if pm is None or getattr(pm.sample_smc, "_fork_patched", False): return
        _orig = pm.sample_smc
        def _smc(*a, **k):
            k.setdefault("mp_ctx", "fork"); return _orig(*a, **k)
        _smc._fork_patched = True; pm.sample_smc = _smc

    _log, _status, _t_all = _io.StringIO(), [], _time.time()
    for _i, _code in _cells:
        _st["cell"], _st["k"] = _i, 0
        _t0, _buf, _err = _time.time(), _io.StringIO(), None
        try:
            _tree = _ast.parse(_code)
            _body, _tail = _tree.body[:-1], _tree.body[-1:] if _tree.body else []
            with _ctx.redirect_stdout(_buf), _ctx.redirect_stderr(_buf):
                if _body:
                    exec(compile(_ast.Module(body=_body, type_ignores=[]), f"<cell {_i}>", "exec"), _main)
                if _tail:
                    _node = _tail[0]
                    if isinstance(_node, _ast.Expr):   # echo the last expression, as IPython does
                        _val = eval(compile(_ast.Expression(_node.value), f"<cell {_i}>", "eval"), _main)
                        if _val is not None: _display(_val)
                    else:
                        exec(compile(_ast.Module(body=_tail, type_ignores=[]), f"<cell {_i}>", "exec"), _main)
                _save_open()
        except Exception as e:
            _err = f"{type(e).__name__}: {e}"; _buf.write("\n" + _tb.format_exc()[-2000:])
        _patch_smc()
        _dt = round(_time.time() - _t0, 1)
        _status.append({"cell": _i, "wall_s": _dt, "error": _err})
        _log.write(f"\n===== cell {_i}  ({_dt} s){'  ERROR: ' + _err if _err else ''} =====\n")
        _log.write(_buf.getvalue())
        print(f"cell {_i:>3}  {_dt:>8.1f}s  {'ERR ' + _err if _err else 'ok'}", flush=True)

    with open(f"handoff/log_{_tag}.txt", "w") as fh: fh.write(_log.getvalue())
    _s = {"notebook": _src, "tag": _tag, "wall_s": round(_time.time() - _t_all, 1),
          "n_cells": len(_status), "n_errors": sum(1 for s in _status if s["error"]),
          "errors": [s for s in _status if s["error"]], "n_figures": len(_figs), "figures": _figs}
    with open(f"handoff/run_{_tag}.json", "w") as fh: _json.dump(_s, fh, indent=1)
    print(_json.dumps({k: v for k, v in _s.items() if k != "figures"}, indent=1))


if __name__ == "__main__":
    if len(_sys.argv) != 3:
        raise SystemExit("usage: python exec_nb.py <notebook.ipynb> <tag>")
    _run(_sys.argv[1], _sys.argv[2])
