from setuptools import setup, find_packages

# Read long description from readme.md.
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except:
    long_description = None

# Read version from package.
from ionize.__version__ import __version__

setup(name='ionize',
      version=__version__,
      author='Lewis A. Marshall',
      author_email='lewis.a.marshall@gmail.com',
      url="http://lewisamarshall.github.io/ionize/",
      classifiers=[
          "Programming Language :: Python",
          "Development Status :: 3 - Alpha",
          "Environment :: Console",
          "Intended Audience :: Science/Research",
          "Operating System :: OS Independent",
          "Topic :: Software Development :: Libraries :: Python Modules",
          "Topic :: Scientific/Engineering :: Chemistry",
          ],
      use_2to3=False,
      license='LICENSE',
      description='A package for calculating electrolyte properties.',
      long_description=long_description,
      packages=find_packages(),
      requires=['numpy', 'scipy', 'biopython', 'click'],
      package_data={'ionize': ['Database/ion_data.json']},
      entry_points={'console_scripts': ['ionize = ionize.__main__:cli']},
      test_suite="ionize.tests",
      )
