import os
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='bio_assembly_refinement',
    version='0.3.0',
    description='Assembly refinement tools, mostly useful for (but not limited to) pacbio bacterial assembly',
    long_description=read('README.md'),
    packages = find_packages(),
    author='Nishadi De Silva',
    author_email='nds@sanger.ac.uk',
    url='https://github.com/nds/bio_assembly_refinement',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    install_requires=['nose >= 1.3', 'pyfastaq >= 3.5.0'],
    license='GPLv3',
)
