from setuptools import setup, find_packages

setup(name='ionize',
      version='0.2.0c',
      author='Lewis A. Marshall',
      author_email='lewis.a.marshall@gmail.com',
      url="http://lewisamarshall.github.io/ionize/",
      classifiers=[
          "Programming Language :: Python",
          "Development Status :: 3 - Alpha",
          "Environment :: Console",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
          "Operating System :: OS Independent",
          "Topic :: Software Development :: Libraries :: Python Modules",
          "Topic :: Scientific/Engineering :: Chemistry",
          ],
      use_2to3 = True,
      license='LICENSE',
      description='A package for calculating electrolyte properties.',
      long_description=open('README.txt').read(),
      packages=find_packages(),
      requires=['numpy', 'scipy'],
      package_data={'ionize': ['ions_db.json',
                               'nightingale_data/*.txt']
                    }
      )
