"""Set up the pdb2cif package."""
import setuptools


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name="go2pdb",
    version="1.0.0",
    description=(
        "This code identifies Protein Data Bank (PDB) entries with specific "
        "Gene Ontology (GO) annotations."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires=">=3.8",
    author="Nathan Baker",
    author_email="nathanandrewbaker@gmail.com",
    packages=setuptools.find_packages(),
    install_requires=[
        "requests",
        "biopython",
        "pandas",
        "openpyxl",
        "numpy",
        "mmcif_pdbx",
    ],
    tests_require=["pytest"],
    entry_points={"console_scripts": ["go2pdb=go2pdb.__main__:main"]},
    keywords="science chemistry biophysics biochemistry",
    url="https://github.com/Electrostatics/go2pdb",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
