"""The setup module for SequenceBouncer"""

# Always prefer setuptools over distutils
from setuptools import setup, find_namespace_packages
import pathlib

# Get the long description from the README file
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')


setup(
    name='sequence_bouncer',
    version='1.23.1',
    description='A setuptools fork for SequenceBouncer',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Cory Dunn',
    author_email='cory.dunn@helsinki.fi',
    maintainer='Patmanidis Stefanos',
    maintainer_email='stefanpatman91@gmail.com',
    package_dir={'': 'src'},
    packages=find_namespace_packages(
        include=('itaxotools*',),
        where='src',
    ),
    python_requires='>=3.9, <3.10',
    install_requires=[
        'biopython==1.78',
        'matplotlib==3.4.2',
        'numpy==1.20.2',
        'pandas==1.2.4',
        ],
    extras_require={
        'dev': ['pytest>=6.2.5'],
    },
    entry_points={
        'console_scripts': [
            'SequenceBouncer = itaxotools.sequence_bouncer.__main__:main',
        ]
    },
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',  # noqa
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
    include_package_data=True,
)
