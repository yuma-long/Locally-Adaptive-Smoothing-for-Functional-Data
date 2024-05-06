import sys, os
import setuptools
from setuptools import setup
from setuptools.command.build_ext import build_ext
from glob import glob
from pybind11.setup_helpers import Pybind11Extension

class get_pybind_include(object):
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11

        return pybind11.get_include(self.user)


include_dirs = [get_pybind_include(), get_pybind_include(user=True), "cpp/include"]

ext_modules = [
    Pybind11Extension(
        "tv_denoiser",
        sources=["cpp/pybind/flsa_pybind.cpp"] + glob("cpp/core/flsa.cpp"),
        include_dirs=include_dirs,
        language="c++",
    )
]


def is_options_supported(compiler, option):
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".cpp") as f:
        f.write("int main (int argc, char **argv) { return 0;}")
        try:
            compiler.compile([f.name], extra_postargs=[option])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def choose_compiler_version(compiler):
    options = ["-std=c++20", "-std=c++17", "-std=c++14", "-std=c++11"]

    for option in options:
        if is_options_supported(compiler, option):
            return option

    raise RuntimeError("Compiler not found: at least C++11 support is needed")


class BuildExt(build_ext):
    c_opts = {
        "msvc": ["/EHsc", "/O2"],
        "unix": ["-O3"],
    }
    l_opts = {
        "msvc": [],
        "unix": [],
    }

    if sys.platform == "darwin":
        darwin_opts = [
            "-mmacosx-version-min=10.8",
        ]
        c_opts["unix"] = darwin_opts + c_opts["unix"]
        l_opts["unix"] = darwin_opts + l_opts["unix"]

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == "unix":
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(choose_compiler_version(self.compiler))
            if is_options_supported(self.compiler, "-fvisibility=hidden"):
                opts.append("-fvisibility=hidden")
        elif ct == "msvc":
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        build_ext.build_extensions(self)


setup(
    name="tv_denoiser",
    version="1.0.0",
    packages=setuptools.find_packages(),
    author="Yuma Takeda",
    author_email="utklav1511@gmail.com",
    ext_modules=ext_modules,
    install_requires=["pybind11>=2.9"],
    cmdclass={"build_ext": BuildExt},
    zip_safe=False,
)
