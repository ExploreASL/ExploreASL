from setuptools import setup, find_packages

CLASSIFIERS = [
    "License :: OSI Approved :: MIT License",
    "Intended Audience :: Developers",
    "Intended Audience :: Healthcare Industry",
    "Intended Audience :: Science/Research",
    "Development Status :: Unstable",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Medical Science Apps."]

with open("README.md") as reader:
    long_description = reader.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='ExploreASL_GUI',
    version='0.0.1',
    python_requires=">=3.7",
    packages=find_packages(),
    package_dir={"": "ExploreASL_GUI"},
    url='https://github.com/MauricePasternak/ExploreASL_GUI',
    license='GPLv3',
    author='Maurice Pasternak',
    author_email='maurice.pasternak@utoronto.ca',
    description='A graphical interface with the ExploreASL program',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=required,
    package_data={
        "JSON_LOGIC": ["ExecutorTranslators.json", "GraphingParameters.json"],
        "media": ["*.png", "*.gif", "*.svg", "*.ico"],
        "External": ["DCM2NIIX/*.exe"]
    },
    include_package_data=True,
    classifiers=CLASSIFIERS
)
