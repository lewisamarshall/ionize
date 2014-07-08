"""Represent ions in aqueous solution."""

from Ion import Ion
from Solution import Solution
# from load_ion import load_ion

if __name__ == "__main__":
    a = Ion('A', [-1], [5], [-10])
    b = Ion('B', [1], [8], [10])
    print a
    # s = Solution([a, b], [0.05, 0.03])
    pass
