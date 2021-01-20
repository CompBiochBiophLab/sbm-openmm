import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sbmOpenMM", 
    version="1.0.0",
    author="Martin Floor, Kengjie Li",
    author_email="martinfloor@gmail.com",
    description="An OpenMM package for simulating protein Structure Based Models.",
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
    'numpy>=1.15',
    ],
)
