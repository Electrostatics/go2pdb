"""Set up the pdb2cif package."""
import setuptools


setuptools.setup(
    name="go2pdb",
    version="0.0.1",
    description=(
        "This code downloads and analyzes Protein Data Bank entries with specific GO annotations."
    ),
    python_requires=">=3.8",
    author="Nathan Baker",
    author_email="nathan.baker@pnnl.gov",
    packages=setuptools.find_packages(),
    # package_data={"": ["./config.ini"]},
    install_requires=["requests", "biopython", "pandas", "openpyxl", "numpy", "mmcif_pdbx"],
    tests_require=["pytest"],
    # entry_points={"console_scripts": ["mvalue=osmolytes.main:console"]},
    keywords="science chemistry biophysics biochemistry",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Common Public License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
