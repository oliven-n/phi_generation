from setuptools import setup, find_packages

setup(
    name='your_package_name',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'langchain',
        'ipykernel',
        # Other dependencies...
    ],
    # ... other metadata
)