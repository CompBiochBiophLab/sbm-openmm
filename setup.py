import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sbmOpenMM", 
    version="0.0.1",
    author="Martin Floor, Kengjie Li",
    author_email="martinfloor@gmail.com",
    description="An OpenMM package for Structure Based Model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CompBiochBiophLab/sbm-openmm",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
    'openmm>=7.0',
    'numpy>=1.15',
    ],
)
