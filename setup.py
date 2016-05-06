from setuptools import setup, find_packages
import versioneer

setup(
    name='destruct',
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
    scripts=[
        'destruct/scripts/destruct_run.py',
        'destruct/scripts/destruct_create_ref_data.py',
    ],
    package_data={'destruct': ['data/*_chr_map.tsv']},
)

