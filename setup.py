from setuptools import setup

setup(
    name='mixemt',
    version='0.1',
    description='',
    url='http://github.com/svohr/mixemt',
    author='Samuel H. Vohr',
    author_email='svohr@soe.ucsc.edu',
    license='MIT',
    packages=['mixemt'],
    package_data={'mixemt': ['ref/rCRS.mtDNA.fa', 'ref/RSRS.mtDNA.fa',
                             'phylotree/mtDNA_tree_Build_16.csv',
                             'phylotree/mtDNA_tree_Build_17.csv',
                             'ref/README.txt', 'phylotree/README.md']},
    scripts=['bin/mixemt'],
    install_requires=['numpy', 'scipy', 'pysam', 'biopython']
    )
