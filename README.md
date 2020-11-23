reactions
======

Explore possible nuclear reactions resulting from the collision of two parent nuclides:

```
%  python scripts/calc.py "H+Li" --simple
d + 6Li → 2·4He + 22373 keV
d + 6Li → ɣ + 8Be + 22281 keV
t + 6Li → ɣ + 9Be + 17688 keV
p + 7Li → 2·4He + 17346 keV
p + 7Li → ɣ + 8Be + 17254 keV
t + 7Li → ɣ + 10Be + 17249 keV
d + 7Li → ɣ + 9Be + 16694 keV
t + 6Li → n + 2·4He + 16116 keV
t + 6Li → n + 8Be + 16024 keV
t + 6Li → 4He + 5He + 15381 keV
d + 7Li → n + 2·4He + 15122 keV
d + 7Li → n + 8Be + 15030 keV
d + 7Li → 4He + 5He + 14387 keV
t + 7Li → n + 9Be + 10437 keV
t + 7Li → 4He + 6He + 9840 keV
t + 7Li → 2·n + 8Be + 8773 keV
t + 7Li → n + 4He + 5He + 8130 keV
t + 7Li → 2·5He + 7395 keV
p + 6Li → ɣ + 7Be + 5607 keV
d + 6Li → p + 7Li + 5027 keV
p + 6Li → 3He + 4He + 4020 keV
d + 6Li → n + 7Be + 3382 keV
d + 6Li → p + t + 4He + 2559 keV
d + 6Li → n + 3He + 4He + 1795 keV
d + 6Li → 3He + 5He + 1060 keV
t + 6Li → d + 7Li + 994 keV
t + 6Li → p + 8Li + 802 keV
d + 6Li → t + 5Li + 593 keV
```


### Requirements

To get this working, you'll need:
* Python 3
* [Virtualenv wrapper](https://virtualenvwrapper.readthedocs.org/en/latest/)

### Setup

Run the following commands to get set up:

```
% git clone https://github.com/emwalker/reactions.git
% cd reactions
% source /usr/local/bin/virtualenvwrapper.sh
% mkvirtualenv reactions --python /usr/local/bin/python3
% pip install -r requirements.txt
% nosetests
...............................................................
----------------------------------------------------------------------
Ran 63 tests in 0.211s

OK
```

### Running

Activate the virtual environment to start:
```
% workon reactions
% python scripts/calc.py "H+Li"
...
```

Use `deactivate` to leave virtual environment:
```
% deactivate
```
