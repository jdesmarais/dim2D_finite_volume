from distutils.core import setup

with open('README.txt', 'r') as readme:
    README_TEXT = readme.read()

setup(name='fortranfile',
      version='0.2.1',
      author='Neil Martinsen-Burrell',
      author_email='neilmartinsenburrell@gmail.com',
      url='http://launchpad.net/fortranfile',
      description='Read and write unformatted fortran output files',
      long_description=README_TEXT,
      classifiers=['License :: OSI Approved :: MIT License',
          'Programming Language :: Python',
          'Programming Language :: Fortran',
          'Operating System :: OS Independent',
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering',
          ],
      py_modules=['fortranfile'],
      )
