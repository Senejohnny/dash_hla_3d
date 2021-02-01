
from setuptools import setup, find_packages
setup(
    name = 'hla3d',
    version="0.2.0",
    author="Danial Senejohnny",
    author_email="d.senejohnny@gmail.com",
    description="HLA Epitope 3D visualisation",
    packages = find_packages(),
    install_requires=[
        'pandas',
    ],
)
