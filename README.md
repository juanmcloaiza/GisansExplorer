# GisansExplorer

Simple **GUI** that interactively helps in the reduction of [NICOS](https://nicos-controls.org/) [GISANS](http://www.gisaxs.de/theory.html) data.

![Screenshot](./Screenshot.png)

# Install:
```
 $ git clone https://github.com/juanmcloaiza/GisansExplorer.git
 $ cd NeutronDataViewer
 ```
 
 [**Use Python 3.7**](https://realpython.com/installing-python/) from here on:
 
 Linux / Mac:
 ```
 $ python -m venv ./GisansExplorerEnv
 $ source ./GisansExplorerEnv/bin/activate
```

Windows:
```
 $ python -m venv .\GisansExplorerEnv
 $ python -m venv --system-site-packages GisansExplorerEnv
 $ python -m venv --upgrade GisansExplorerEnv
 $ .\GisansExplorerEnv\Scripts\activate
```

# Run:
```
 $ python GisansExplorer.py
```

**N.B.** Should requirements not be satisfied, run the following in the terminal within the GisansExplorerEnv [virtual environment](https://docs.python.org/3/tutorial/venv.html) created above:
```
 (GisansExplorerEnv) $ pip install --upgrade pip
 (GisansExplorerEnv) $ pip install -r requirements.txt
 ```
