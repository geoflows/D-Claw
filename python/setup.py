#! /usr/bin/env python


import setuptools

setuptools.setup(
    name="dclaw",
    version=4.0,
    author="",
    author_email="",
    description="",
    long_description="",
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=2.7",
    install_requires=open("requirements.txt", "r").read().splitlines(),
)
