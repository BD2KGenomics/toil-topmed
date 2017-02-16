from version import version, required_versions
from setuptools import find_packages, setup


kwargs = dict(
    name='toil-topmed',
    version=version,
    description="UC Santa Cruz Computational Genomics Lab's Toil-based TOPMed pipeline",
    author='UCSC Computational Genomics Lab',
    author_email='cgl-toil@googlegroups.com',
    url="https://github.com/BD2KGenomics/toil-topmed",
    install_requires=[x + y for x, y in required_versions.iteritems()],
    tests_require=['pytest==2.8.3'],
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={
        'console_scripts': ['toil-topmed = toil_topmed.topmed_cgl_pipeline:main']})


setup(**kwargs)
