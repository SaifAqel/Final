import logging, functools, time
TRACE_LEVEL_NUM = 5
logging.addLevelName(TRACE_LEVEL_NUM, "TRACE")

class _TraceLogger(logging.Logger):
    def trace(self, msg, *a, **k):
        if self.isEnabledFor(TRACE_LEVEL_NUM):
            self._log(TRACE_LEVEL_NUM, msg, a, **k)
logging.setLoggerClass(_TraceLogger)

_FMT = "%(asctime)s | %(levelname)s | %(name)s | stage=%(stage)s step=%(step)s | %(message)s"
_DATE = "%Y-%m-%d %H:%M:%S"

class _Stage(logging.Filter):
    def filter(self, r):
        if not hasattr(r,"stage"): r.stage="-"
        if not hasattr(r,"step"): r.step="-"
        return True

def setup_logging(level: int | str = logging.INFO):
    if isinstance(level,str):
        lvl = logging.getLevelName(level.upper())
        level = lvl if isinstance(lvl,int) else logging.INFO
    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(level)
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter(_FMT,_DATE))
    h.addFilter(_Stage())
    root.addHandler(h)

def _fmt(v):
    try:    return f"{v:.6g~P}"  # pint Quantity pretty
    except: return repr(v)

def trace_calls(name: str|None=None, values: bool=False):
    def _wrap(fn):
        qual = name or f"{fn.__module__}.{fn.__qualname__}"
        log = logging.getLogger(qual)
        @functools.wraps(fn)
        def _inner(*a, **k):
            stage = k.get("stage","-")
            log.trace("enter", extra={"stage": stage, "step": fn.__name__})
            if values:
                arg_s = ", ".join([*map(_fmt, a),
                                   *[f"{kk}={_fmt(v)}" for kk,v in k.items()]])
                log.trace(f"args: {arg_s}", extra={"stage": stage, "step": fn.__name__})
            t0=time.perf_counter()
            try:
                out=fn(*a, **k)
                dt=(time.perf_counter()-t0)*1000
                if values:
                    log.trace(f"ret: {_fmt(out)}", extra={"stage": stage, "step": fn.__name__})
                log.trace(f"exit ok in {dt:.2f} ms", extra={"stage": stage, "step": fn.__name__})
                return out
            except Exception as e:
                log.exception(f"exit err: {e}", extra={"stage": stage, "step": fn.__name__})
                raise
        return _inner
    return _wrap

