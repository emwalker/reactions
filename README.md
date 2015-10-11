lenrmc
======

This project is intended to hold a set of tools for those interested in exploring different
reaction pathways.  Currently `scripts/reactions.py` is the only tool available.

# Requirements

The following components are needed to get this up and running:
* Python 3
* [Virtualenv wrapper](https://virtualenvwrapper.readthedocs.org/en/latest/)

# Setup

To get set up, run the following commands:

```
% git clone https://github.com/emwalker/lenrmc.git
% cd lenrmc
% mkvirtualenv reactions --python /usr/local/bin/python3
% pip install -r requirements.txt
% nosetests
...............................................................
----------------------------------------------------------------------
Ran 63 tests in 0.211s

OK
$
```

# Running

Use the following command to execute the script:

```
% # switch to the virtual environment used for this project
% workon reactions
% python scripts/reactions.py "H+Li"
d + 6Li → p + 7Li + 5027 keV                            in nature, n-transfer         ✗ 6Li [L15],   ✗ 7Li [L15]
d + 6Li → p + t + 4He + 2559 keV                        n-transfer, t, α, →β-         ✗ 6Li [L15]
d + 6Li → t + 5Li + 593 keV                             n-transfer, t, →p, →β-        ✗ 6Li [L15]
d + 6Li → 2·4He + 22373 keV                             in nature, α                  ✗ 6Li [L15]
p + 7Li → 2·4He + 17346 keV                             in nature, α
d + 7Li → 4He + 5He + 14387 keV                         α, →n
p + 6Li → 3He + 4He + 4020 keV                          in nature, α                  ✗ 6Li [L15]
d + 6Li → 3He + 5He + 1060 keV                          →n                            ✗ 6Li [L15]
d + 6Li → ɣ + 8Be + 22281 keV                           ɣ, →α                         ✗ 6Li [L15]
p + 7Li → ɣ + 8Be + 17254 keV                           ɣ, →α
d + 7Li → ɣ + 9Be + 16694 keV                           in nature, ɣ
p + 6Li → ɣ + 7Be + 5607 keV                            ɣ, →ε                         ✗ 6Li [L15]
d + 7Li → n + 2·4He + 15122 keV                         n, α, →β-
d + 7Li → n + 8Be + 15030 keV                           n, →α, →β-
d + 6Li → n + 7Be + 3382 keV                            n, →β-, →ε                    ✗ 6Li [L15]
d + 6Li → n + 3He + 4He + 1795 keV                      n, α, →β-                     ✗ 6Li [L15]

[L15] 2015 Lugano E-Cat test by Levi et al.
%
% python scripts/reactions.py "H+Li" --ascii --simple

d + 6Li => p + 7Li + 5027 keV
d + 6Li => p + t + 4He + 2559 keV
d + 6Li => t + 5Li + 593 keV
d + 6Li => 2*4He + 22373 keV
d + 6Li => gamma + 8Be + 22281 keV
p + 7Li => 2*4He + 17346 keV
p + 7Li => gamma + 8Be + 17254 keV
d + 7Li => gamma + 9Be + 16694 keV
d + 7Li => 4He + 5He + 14387 keV
p + 6Li => gamma + 7Be + 5607 keV
p + 6Li => 3He + 4He + 4020 keV
d + 6Li => 3He + 5He + 1060 keV
d + 7Li => n + 2*4He + 15122 keV
d + 7Li => n + 8Be + 15030 keV
d + 6Li => n + 7Be + 3382 keV
d + 6Li => n + 3He + 4He + 1795 keV
% # Leave the virtual environment used for this project
% deactivate
```
