from setuptools import setup

setup(
    name='paleopy',
    version='1.0.0',
    description='Calculating Signal and Background Expectations for Paleo Detectors.',
    author='Thomas Edwards and Bradley Kavanagh',
    author_email='tedwards2412@gmail.nl',
    url = 'https://github.com/tedwards2412/paleo_detectors',
    packages=['paleopy'],
    package_data={'paleopy': [] },
    long_description=open('README.md').read(),
)