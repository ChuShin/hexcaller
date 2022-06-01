from setuptools import setup, find_packages

setup(
    name='hexcaller',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'deeptools',
        'Click'
    ],
    entry_points={
        'console_scripts': [
            'hexcaller = hexcaller.src.hexcaller:call_cnv',
        ],
    },
)
