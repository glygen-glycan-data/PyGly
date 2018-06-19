#!/usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

setup(name="PyGly",
      version="1.1.0",
      description="Python package for reading, writing, and manipulating glycans",
      packages=['pygly'],
      )
