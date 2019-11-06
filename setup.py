import setuptools

setuptools.setup(
    name="sbmOpenMM", 
    version="0.0.1",
    author="Martin Floor, Kengjie Li",
    author_email="martinfloor@gmail.com",
    description="An OpenMM package for Structure Based Model",
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
