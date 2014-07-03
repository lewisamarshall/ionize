"""Represent ions in aqueous solution."""

from ion import ion
from solution import solution
# from load_ion import load_ion

if __name__ == "__main__":
    a = ion('a', [-1], [5], [-10])
    print a
    s = solution([a], [0.05])
    pass
