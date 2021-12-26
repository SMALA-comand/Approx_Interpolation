from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="approxinterpol",
    packages=find_packages(),
    version="0.0.3",
    description='Ready-made Approximation and interpolation',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Mark Kozlov",
    author_email="mark.k.2012@yandex.ru",
    url="https://github.com/SMALA-comand/Approx_Interpolation",
    install_requires=['scipy', 'numpy', 'sympy'],
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)
