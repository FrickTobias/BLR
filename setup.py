import sys
from setuptools import setup, find_packages

with open("README.md") as f:
    long_description = f.read()

setup(
    name="blr",
    author="Tobias Frick",
    url="https://github.com/TobiasFrick/BLR/",
    description="Barcoded long reads pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="",
    python_requires=">=3.6",
    package_dir={"": "src"},
    packages=find_packages("src"),
    entry_points={"console_scripts": ["blr = blr.__main__:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
#       "License :: ...",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ]
)
