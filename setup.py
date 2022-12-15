from setuptools import setup, find_packages

setup(
    name='CVM_Toolkit',
    version='0.1.0',
    author='Sayan Samanta',
    author_email='sayan_samanta@brown.edu',
    packages=find_packages(),
    scripts=['scripts/sro_correction.py',],
    license='LICENSE.txt',
    description='A CVM Optimizer',
    long_description=open('README.md').read(),
    install_requires=[
        "ase >= 3.22.1",
        "numpy >= 1.23.5",
        "scipy >= 1.9.0",
        "sympy >= 1.11",
        "matplotlib >= 3.5.3",
        "pandas >= 1.4.3",
    ],
)
