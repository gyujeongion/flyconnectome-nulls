"""Build-time preprocessing: convert Unicode super/subscript runs to pandoc ^...^ / ~...~ so any text font can render them."""
import re, sys
SUP = dict(zip("⁰¹²³⁴⁵⁶⁷⁸⁹⁻ᵀ", "0123456789−T")); SUB = dict(zip("₀₁₂₃₄₅₆₇₈₉", "0123456789"))
t = open(sys.argv[1], encoding="utf-8").read()
t = re.sub("[" + "".join(SUP) + "]+", lambda m: "^" + "".join(SUP[c] for c in m.group(0)) + "^", t)
t = re.sub("[" + "".join(SUB) + "]+", lambda m: "~" + "".join(SUB[c] for c in m.group(0)) + "~", t)
open(sys.argv[2], "w", encoding="utf-8").write(t)
