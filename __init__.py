"""Represent ions in aqueous solution."""

from Ion import Ion
from Solution import Solution
# from load_ion import load_ion

if __name__ == "__main__":
    a = Ion('A', [-1], [-2], [-10])
    b = Ion('B', [1], [7], [10])
    s = Solution([a, b], [0.03, 0.06])
    print s.I
    print s.pH
    print s.conductivity()
