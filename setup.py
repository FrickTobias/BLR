from setuptools import setup, find_namespace_packages

with open("README.md") as f:
    long_description = f.read()

setup(
    name="blr",
    author="Tobias Frick",
    url="https://github.com/TobiasFrick/BLR/",
    description="Barcoded long reads pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.6",
    install_requires=[
        "pysam",
        "dnaio",
        "tqdm",
        "snakemake",
        "importlib_resources; python_version<'3.7'",
    ],
    extras_require={
        "dev": ["flake8"],
    },
    package_dir={"": "src"},
    packages=find_namespace_packages("src"),
    package_data={"blr": ["Snakefile", "rules/*.smk", "config.schema.yaml", "blr.yaml"]},
    entry_points={"console_scripts": ["blr = blr.__main__:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ]
)
