from setuptools import setup, find_packages

setup(
    name='shiba2sashimi',
    version='v0.1.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
    	"numpy>=1.18.0",
		"matplotlib>=3.1.0",
		"scipy>=1.4.0",
		"pysam>=0.15.0"
    ],
    entry_points={
        'console_scripts': [
            'shiba2sashimi=shiba2sashimi.main:main',
        ],
    },
    author='Naoto Kubota',
    author_email='the.owner.of.the.dream123@gmail.com',
    description='A utility for creating sashimi plot from Shiba output',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/NaotoKubota/shiba2sashimi',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)