from setuptools import setup, find_packages
import versioneer

setup(
    name='ngs-destruct',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Destruct is a tool for joint prediction of rearrangement breakpoints from single or multiple tumour samples',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    url='http://bitbucket.org/dranew/destruct',
    download_url='https://bitbucket.org/dranew/destruct/get/v{}.tar.gz'.format(versioneer.get_version()),
    keywords=['scientific', 'sequence analysis', 'cancer'],
    classifiers=[],
    entry_points={'console_scripts': ['destruct = destruct.ui.main:main']},
    package_data={'destruct': ['data/*_chr_map.tsv']},
    setup_requires=[
        'numpy',
        'setuptools',
    ],
    install_requires=[
        'numpy',
        'pandas',
        'pypeliner',
        'pygenes',
        'pyyaml',
        'matplotlib',
        'seaborn',
        'pysam',
        'blossomv',
    ],
    dependency_links=[
        'https://github.com/amcpherson/blossomv/archive/v2.04_r4.tar.gz#egg=blossomv',
    ],
)

