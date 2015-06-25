try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

from Cython.Distutils import build_ext
import pkg_resources
import platform
import re
import subprocess


def is_clang(bin):
    proc = subprocess.Popen([bin, '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    output = '\n'.join([stdout, stderr])
    return not re.search(r'clang', output) is None

class my_build_ext( build_ext ):
    def build_extensions(self):
        binary = self.compiler.compiler[0]
        if is_clang(binary):
            for e in self.extensions:
                e.extra_compile_args.append('-stdlib=libc++')
                if platform.system() == 'Darwin':
                    e.extra_compile_args.append('-mmacosx-version-min=10.7')
        build_ext.build_extensions(self)


compile_args = ['-std=c++1y']

data_dir = pkg_resources.resource_filename("autowrap", "data_files")

ext = Extension("bpp",
                sources = ['bpp.pyx',
                           'src/Alignment.cpp',
                           'src/ModelFactory.cpp',
                           'src/SiteContainerBuilder.cpp'],
                language="c++",
                include_dirs = [data_dir],
                libraries=['bpp-core', 'bpp-seq', 'bpp-phyl'],
                extra_compile_args=compile_args,
               )

setup(cmdclass={'build_ext':my_build_ext},
      name="bpp",
      author='Kevin Gori',
      author_email='kgori@ebi.ac.uk',
      description='Pairwise distances by maximum likelihood',
      url='https://github.com/kgori/bpp.git',
      version="0.0.13",
      scripts=['bin/pairdist', 'bin/simulate'],
      ext_modules = [ext],
      install_requires=[
        'autowrap',
        'cython',
      ],
     )
