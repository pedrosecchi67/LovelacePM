import setuptools
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


#libraries=[('toolkit', dict(sources=['toolkit.f90']))]
extensions=[Extension(name='LovelacePM.toolkit', sources=['LovelacePM/toolkit.f90'], language='f90'), \
    Extension(name='LovelacePM.fdyn', sources=['LovelacePM/fdyn.f90'], language='f90')]

setup(
    name="LovelacePM",
    version="0.1.5",
    author="Pedro de Almeida Secchi",
    author_email="pedrosecchimail@gmail.com",
    description="Python based, open source vortex ring panel method code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pedrosecchi67/LovelacePM",
    packages=['LovelacePM'],
    #libraries=libraries,
    ext_modules=extensions,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy', 'scipy', 'matplotlib', 'LoveUpdate'],
    python_requires='>=3.6',
)
