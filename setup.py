import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bayesian_mirt",
    version="0.0.0a0",
    author="jiguang_li",
    author_email="jiguang@chicagobooth.edu",
    description="Exploring Bayesian MIRT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/smulla/Bayesian-MIRT",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: Other/Proprietary License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3"
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.7",
    install_requires=[
        # PyMC3 dropped support for python3.6, but we can make it work with version-specific dependencies
        "numpy",
        "pandas",
    ]
)
