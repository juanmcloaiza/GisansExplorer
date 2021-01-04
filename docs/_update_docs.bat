REM for more info: https://sphinx-rtd-tutorial.readthedocs.io/en/latest/build-the-docs.html
./make.bat clean
sphinx-apidoc -o ./source ../gisansexplorer
./make.bat html
