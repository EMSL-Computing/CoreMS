__author__ = "Yuri E. Corilo"
__version__ = "3.5.0"
import time
import os
import sys
import hashlib

# Get the path to the README file
readme_path = os.path.join(os.path.dirname(__file__), "..", "README.md")

# Read the contents of the README file if it exists
if os.path.exists(readme_path):
    try:
        with open(readme_path, "r", encoding="utf-8") as readme_file:
            __doc__ = readme_file.read()
    except Exception as e:
        __doc__ = "CoreMS: A comprehensive mass spectrometry framework for software development and data analysis of small molecules analysis."
        print(f"Warning: Could not read README.md file. Error: {e}")
else:
    __doc__ = "CoreMS: A comprehensive mass spectrometry framework for software development and data analysis of small molecules analysis."



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
