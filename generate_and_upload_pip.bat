
pip uninstall gisansexplorer
rm -r build
rm -r dist
rm -r .\gisansexplorer.egg-info

python setup.py sdist bdist_wheel
python -m twine upload dist/*

pip install gisansexplorer