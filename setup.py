"""
dropletbuilder
Droplet on graphene builder
"""
from setuptools import setup

short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:]),


setup(
    # Self-descriptive entries which should always be present
    name='dropletbuilder',
    author='Felix Tiet',
    author_email='ftiet21@gmail.com',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    packages=['dropletbuilder', "dropletbuilder.tests"],

    # Optional include package data to ship with your package
    # Comment out this line to prevent the files from being packaged with your software
    # Extend/modify the list to include/exclude other items as need be
    package_data={'dropletbuilder': ["data/*.dat"]
                  },

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # author_email='me@place.org',      # Author email
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

    # entry_points={
    #     'mbuild.plugins':[
    #         "GrapheneDroplet = dropletbuilder.dropletbuilder:GrapheneDroplet",
    #     ],
    # }

)
