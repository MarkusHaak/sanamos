#!/usr/bin/env python

# Big parts taken from:
# https://github.com/joerick/python-ctypes-package-sample
# https://github.com/himbeles/ctypes-example?tab=readme-ov-file

import os
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

class build_ext(_build_ext):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_ext_filename(self, ext_name):
        return ext_name + ".so"

class bdist_wheel_abi_none(_bdist_wheel):
    def finalize_options(self):
        _bdist_wheel.finalize_options(self)
        self.root_is_pure = False

    def get_tag(self):
        python, abi, plat = _bdist_wheel.get_tag(self)
        return "py3", "none", plat

setup(
    name="suffix_array",
    ext_modules=[
        Extension(
            name="suffix_array.sa",
            sources=["src/suffix_array/sa.c"],
            extra_compile_args=["--std=gnu11", "-O3", "-march=native", "-lrt", "-fopenmp", "-lm", "-funsigned-char"],
            extra_link_args=["-fPIC", "-shared", "-fopenmp"]
        ),
    ],
    cmdclass={"build_ext": build_ext, "bdist_wheel": bdist_wheel_abi_none},
)