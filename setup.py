from setuptools import setup, find_packages

packages = find_packages()

setup(
    name='oligo_nn',
    version='0.1',
    packages=packages,
    install_requires=['numpy', 'matplotlib', 'pandas', 'torch', 'torchvision', 'torchsummary', 'biopython',
                      'import-ipynb', 'scikit-learn']
)
