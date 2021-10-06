# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#

from setuptools import setup, find_packages
import re
import sys

url = "https://github.com/tomchartier/SHERIFS"

README = """
Python Toolkit for the construction of fault based models
Copyright (C) 2017-2021 Thomas Chartier
"""


def get_version():
    version_re = r"^__version__\s+=\s+['\"]([^'\"]*)['\"]"
    version = None

    package_init = 'sherifs/__init__.py'
    for line in open(package_init, 'r'):
        version_match = re.search(version_re, line, re.M)
        if version_match:
            version = version_match.group(1)
            break
    else:
        sys.exit('__version__ variable not found in %s' % package_init)

    return version


version = get_version()

setup(
    name='sherifs',
    version=version,
    description=README,
    url=url,
    packages=find_packages(exclude=['tests', 'tests.*']),
    # Minimal requirements, for a complete list see requirements.txt
    install_requires=[
        'geojson',
    ],
    python_requires='>=3.8',
    author='Thomas Chartier',
    author_email='thomas.chartier@globalquakemodel.org',
    maintainer='Thomas Chartier',
    maintainer_email='thomas.chartier@globalquakemodel.org',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2.1',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ],
    namespace_packages=['sherifs'],
    keywords="seismic hazard",
    license="LGPL2.1",
    platforms=["any"],
    package_data={"sherifs": ["README.md", "LICENSE"]},
    include_package_data=True,
    zip_safe=False,
)
