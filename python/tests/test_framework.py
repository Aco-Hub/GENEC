"""Test framework helpers matching the Fortran test_framework.f90 interface."""


def assert_approx(actual, expected, tol, name=""):
    diff = abs(actual - expected)
    assert diff <= tol, f"FAIL: {name} got={actual:.12e} expected={expected:.12e} tol={tol:.2e} diff={diff:.2e}"


def assert_rel(actual, expected, reltol, name=""):
    denom = max(abs(expected), 1.0e-30)
    rel = abs(actual - expected) / denom
    assert rel <= reltol, f"FAIL: {name} got={actual:.12e} expected={expected:.12e} reltol={reltol:.2e} rel={rel:.2e}"


def assert_equal_int(actual, expected, name=""):
    assert actual == expected, f"FAIL: {name} got {actual}, expected {expected}"


def assert_true(condition, name=""):
    assert condition, f"FAIL: {name}"
