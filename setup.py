from setuptools import setup, find_packages
from os import path

root_dir = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(root_dir, "README.org"), encoding="utf-8") as f:
    long_description = f.read()

with open(path.join(root_dir, "VERSION")) as version_file:
    version = version_file.read().strip()

setup(
    name="gemma",
    version=version,
    description="Genome-wide efficient 'exact' mixed-model analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/genetics-statistics/gemmalib",
    author="Gemma authors",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    keywords="gemma gwas genomics variants lmm linear-models",
    packages=find_packages(),
    python_requires=">=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, !=3.5.*, <4",
    install_requires=["numpy",],
    extras_require={"dev": ["flake8", "black", "pytest"]},
    include_package_data=True,
    entry_points={"console_scripts": ["gemma = gemma2.__main__:main"]},
    project_urls={  # Optional
        "Bug Reports": "https://github.com/genetics-statistics/gemmalib/issues",
        "Source": "https://github.com/genetics-statistics/gemmalib",
    },
)
