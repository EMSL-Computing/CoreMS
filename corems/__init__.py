__author__ = "Yuri E. Corilo"
__version__ = "4.0.1"
import time
import os
import sys
import hashlib

# Package documentation for pdoc: README plus the full install how-to so the
# landing page uses the same styling as the rest of the API site.
_pkg_dir = os.path.dirname(__file__)
_repo_root = os.path.join(_pkg_dir, "..")
_fallback_doc = (
    "CoreMS: A comprehensive mass spectrometry framework for software "
    "development and data analysis of small molecules analysis."
)

_doc_parts = []
for _rel in (
    "README.md",
    os.path.join("docs", "user", "installation.md"),
):
    _path = os.path.join(_repo_root, _rel)
    if not os.path.exists(_path):
        continue
    try:
        with open(_path, "r", encoding="utf-8") as _fh:
            _text = _fh.read().strip()
        if _text:
            _doc_parts.append(_text)
    except Exception as e:
        print(f"Warning: Could not read {_rel} for package docs. Error: {e}")

__doc__ = "\n\n---\n\n".join(_doc_parts) if _doc_parts else _fallback_doc


def timeit(print_time=True):
    def decorator(method):
        def timed(*args, **kw):
            # Extract print_time from kwargs if provided
            local_print_time = kw.pop('print_time', print_time)
            ts = time.time()
            result = method(*args, **kw)
            te = time.time()
            if "log_time" in kw:
                name = kw.get("log_name", method.__name__.upper())
                kw["log_time"][name] = int((te - ts) * 1000)
            elif local_print_time:
                print("%r  %2.2f ms" % (method.__name__, (te - ts) * 1000))
            return result
        return timed
    return decorator


class SuppressPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def corems_md5(fname):
    bytes_io = fname.open("rb").read()

    md5_returned = hashlib.sha256(bytes_io).hexdigest()

    return "{}:{}".format("sha256", md5_returned)
