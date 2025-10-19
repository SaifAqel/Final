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

def trace_calls(name: str|None=None):
    def _wrap(fn):
        qual = name or f"{fn.__module__}.{fn.__qualname__}"
        log = logging.getLogger(qual)
        @functools.wraps(fn)
        def _inner(*a, **k):
            log.trace("enter", extra={"stage": k.get("stage","-"), "step": fn.__name__})
            t0=time.perf_counter()
            try:
                out=fn(*a, **k)
                log.trace(f"exit ok in {(time.perf_counter()-t0)*1000:.2f} ms",
                          extra={"stage": k.get("stage","-"), "step": fn.__name__})
                return out
            except Exception as e:
                log.exception(f"exit err: {e}", extra={"stage": k.get("stage","-"), "step": fn.__name__})
                raise
        return _inner
    return _wrap
