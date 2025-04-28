from setuptools import setup, find_packages

setup(
    name='phi_gen',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'langchain-community',  # Replace langchain with langchain-community
        'langchain-openai',     # New package
        'openai',               # New package
        'chromadb',             # New package
        'python-dotenv',        # New package
        'pytest',               # New package
        'pandas',               # New package
        'ipykernel',            # Already included
        'tabulate',
        'anndata'             # New package
    ],
    # ... other metadata
)