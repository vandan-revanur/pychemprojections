"""
pytest_runner.py – tiny shim that lets Bazel's py_test execute pytest.

Bazel invokes this file as __main__.  It forwards sys.argv (which Bazel
populates with the test file paths and any --args) straight to pytest and
exits with pytest's return code.
"""

import sys
import pytest

if __name__ == "__main__":
    sys.exit(pytest.main(sys.argv[1:]))

