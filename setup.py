from setuptools import setup

setup(
    name='mixemt',
    version='0.1',
    description='A tool for the phylogenetic interpretation of '
                'human mtDNA sequence mixtures',
    url='http://github.com/svohr/mixemt',
    author='Samuel H. Vohr',
    author_email='svohr@soe.ucsc.edu',
    license='MIT',
    packages=['mixemt'],
    package_data={'mixemt': ['ref/rCRS.mtDNA.fa', 'ref/RSRS.mtDNA.fa',
                             'ref/rCRS.mtDNA.fa.fai', 'ref/RSRS.mtDNA.fa.fai',
                             'phylotree/mtDNA_tree_Build_16.csv',
                             'phylotree/mtDNA_tree_Build_17.csv',
                             'ref/README.txt', 'phylotree/README.md']},
    scripts=['bin/mixemt'],
    install_requires=['numpy', 'scipy', 'pysam', 'biopython', 'numba', 'ttd']
    )
