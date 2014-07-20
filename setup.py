from distutils.core import setup, Extension

from Cython.Distutils import build_ext
import pkg_resources

data_dir = pkg_resources.resource_filename("autowrap", "data_files")

ext = Extension("pairdists",
                sources = ['pairdists.pyx',
                           'src/Distance.cpp',
                           'src/ModelFactory.cpp',
                           'src/SiteContainerBuilder.cpp'],
                language="c++",
                include_dirs = [data_dir],
                libraries=['bpp-core', 'bpp-seq', 'bpp-phyl'],
                extra_compile_args=['-std=c++11'],
               )

setup(cmdclass={'build_ext':build_ext},
      name="pairdists",
      version="0.0.1",
      ext_modules = [ext]
     )