from setuptools import setup

setup(
    name='PYLIQ',
    version='1.0.0',
    description='Free soil liquefaction analysis software from SPT data',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    include_package_data=True,
    author='Jorge Ordóñez',
    author_email='jorgeordonezr@gmail.com',
    url='https://github.com/jorge243/PYLIQ',
    packages=['PYLIQ'],
    license=['MIT'],
)