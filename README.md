lenrmc
======

A small set of scripts for simulating the transmission of photons of various energies
through different media and their rate of detection by different types of detector.

Right now the functionality is basic and the results no doubt inaccurate or incorrect.

# Setup

To get set up, run the following commands:

```
% git clone https://github.com/emwalker/lenrmc.git
% cd lenrmc
% virtualenv .python
New python executable in .python/bin/python
Installing setuptools, pip...done.
% .python/bin/pip install -r requirements.txt
% .python/bin/nosetests -w test
..............
----------------------------------------------------------------------
Ran 14 tests in 0.015s

OK
```

# Running

Use the following command to execute the script:

```
% .python/bin/python scripts/photons.py
```
