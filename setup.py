import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gisansexplorer", # Replace with your own username
    version="0.0.20",
    author="Juan M. Carmona Loaiza",
    author_email="j.carmona.loaiza@fz-juelich.de",
    description="Nicos data reduction and visualisation tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/juanmcloaiza/gisansexplorer",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    entry_points={
    'console_scripts':
        ['gisansexplorer=gisansexplorer:entry_point']
                },
    python_requires='>=3.6',
    data_files=[('resources', ['gisansexplorer/resources/Icon.png', 'gisansexplorer/resources/show_test_map.npy'])]#,
)
