from setuptools import setup, find_packages

setup(
    name='destruct',
    version='0.3.1',
    packages=find_packages(),
    scripts=[
        'destruct/run_destruct.py',
        'destruct/create_ref_data.py',
    ],
    package_data={'destruct': ['data/*_chr_map.tsv']},
)
