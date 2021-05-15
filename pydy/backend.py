try:
    from sympy.core import backend as sm
    USE_SYMENGINE = sm.USE_SYMENGINE
except:
    import sympy as sm      # noqa: F401
    USE_SYMENGINE = False
