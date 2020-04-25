import setuptools
#from numpy.distutils.core import Extension
from numpy.distutils.core import setup
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

vnumber="0.2.4"

#libraries=[('toolkit', dict(sources=['toolkit.f90']))]
#extensions=[Extension(name='LovelacePM.toolkit', sources=['LovelacePM/toolkit.f90'], language='f90'), \
#    Extension(name='LovelacePM.fdyn', sources=['LovelacePM/fdyn.f90'], language='f90')]

setup(
    name="LovelacePM",
    version=vnumber,
    author="Pedro de Almeida Secchi",
    author_email="pedrosecchimail@gmail.com",
    description="Python based, open source vortex ring panel method code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pedrosecchi67/LovelacePM",
    packages=['LovelacePM'],
    package_data={'':['*.dat', 'VERSION_NOTES.txt']},
    #libraries=libraries,
    #ext_modules=extensions,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy', 'func-timeout', 'scipy', 'matplotlib', 'cloudpickle', 'pyqtgraph', 'PyOpenGL', 'LoveUpdate'],
    python_requires='>=3.6',
)

try:
    import LoveUpdate as lupd
    lupd.release_note_report(fname='VERSION_NOTES.txt', fdir=os.path.dirname(__file__), program='LovelacePM', version='v'+vnumber)
except:
    print('LoveUpdate not yet available for release notes')