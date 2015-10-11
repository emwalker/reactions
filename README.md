lenrmc
======

Explore different reaction pathways by specifying the parent isotopes or elements:

```
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
%
```


# Requirements

The following components are required:
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
%
```

# Running

Before doing anything else, make sure you are in the virtual environment that is used
with this project:
```
% workon reactions
```

## Examples

In this example, extended information is provided about possible reactions between stable
isotopes of lithium and stable isotopes of hydrogen:

```
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
```

To print out information for several kinds of reaction involving different combinations
of parent isotopes, use quotation marks and separate the parents with commas:
```
% python scripts/reactions.py "Li+Ni, H+Ni" --spins
6Li + 61Ni → p + 4He + 62Ni + 6898 keV                  in nature, n-transfer, α  1+, 3/2-             0+, 0+, 1/2+             ✗ 6Li [L15],  ✗ 61Ni [S05],  ✓ 62Ni [L15]
d + 58Ni → p + 59Ni + 6775 keV                          n-transfer, →β+           0+, 1+               1/2+, 3/2-
d + 60Ni → p + 61Ni + 5596 keV                          in nature, n-transfer     0+, 1+               1/2+, 3/2-              ✗ 60Ni [S05],  ✗ 61Ni [L15],  ✓ 61Ni [S05]
6Li + 58Ni → p + 4He + 59Ni + 5301 keV                  n-transfer, α, →β+        0+, 1+               0+, 1/2+, 3/2-           ✗ 6Li [L15]
6Li + 61Ni → 5Li + 62Ni + 4931 keV                      n-transfer, →p            1+, 3/2-             0+, 3/2-                 ✗ 6Li [L15],  ✗ 61Ni [S05],  ✓ 62Ni [L15]
d + 62Ni → p + 63Ni + 4613 keV                          n-transfer, →β-           0+, 1+               1/2+, 1/2-              ✗ 62Ni [L15]
6Li + 60Ni → p + 4He + 61Ni + 4122 keV                  in nature, n-transfer, α  0+, 1+               0+, 1/2+, 3/2-           ✗ 6Li [L15],  ✗ 60Ni [S05],  ✗ 61Ni [L15],  ✓ 61Ni [S05]
d + 64Ni → p + 65Ni + 3873 keV                          n-transfer, →β-           0+, 1+               1/2+, 5/2-
7Li + 61Ni → 6Li + 62Ni + 3345 keV                      in nature, n-transfer     3/2-, 3/2-           0+, 1+                  ✗ 61Ni [S05],   ✓ 6Li [L15],  ✓ 62Ni [L15]
6Li + 58Ni → 5Li + 59Ni + 3335 keV                      n-transfer, →p, →β+       0+, 1+               3/2-, 3/2-               ✗ 6Li [L15]
6Li + 62Ni → p + 4He + 63Ni + 3139 keV                  n-transfer, α, →β-        0+, 1+               0+, 1/2+, 1/2-           ✗ 6Li [L15],  ✗ 62Ni [L15]
6Li + 64Ni → p + 4He + 65Ni + 2400 keV                  n-transfer, α, →β-        0+, 1+               0+, 1/2+, 5/2-           ✗ 6Li [L15]
6Li + 60Ni → 5Li + 61Ni + 2156 keV                      n-transfer, →p            0+, 1+               3/2-, 3/2-               ✗ 6Li [L15],  ✗ 60Ni [S05],  ✗ 61Ni [L15],  ✓ 61Ni [S05]
7Li + 61Ni → d + 4He + 62Ni + 1871 keV                  in nature, n-transfer, α  3/2-, 3/2-           0+, 0+, 1+              ✗ 61Ni [S05],  ✓ 62Ni [L15]
7Li + 58Ni → 6Li + 59Ni + 1748 keV                      n-transfer, →β+           0+, 3/2-             1+, 3/2-                 ✓ 6Li [L15]
6Li + 62Ni → 5Li + 63Ni + 1173 keV                      n-transfer, →p, →β-       0+, 1+               1/2-, 3/2-               ✗ 6Li [L15],  ✗ 62Ni [L15]
7Li + 60Ni → 6Li + 61Ni + 569 keV                       in nature, n-transfer     0+, 3/2-             1+, 3/2-                ✗ 60Ni [S05],   ✓ 6Li [L15],  ✗ 61Ni [L15],  ✓ 61Ni [S05]
6Li + 64Ni → 5Li + 65Ni + 434 keV                       n-transfer, →p, →β-       0+, 1+               3/2-, 5/2-               ✗ 6Li [L15]
6Li + 64Ni → n + p + 68Zn + 1635 keV                    n, →β-                    0+, 1+               0+, 1/2+, 1/2+           ✗ 6Li [L15]

[many reactions omitted ...]

d + 58Ni → n + 59Cu + 1194 keV                          n, →β+, →β-               0+, 1+               1/2+, 3/2-
6Li + 60Ni → n + 4He + 61Cu + 1102 keV                  n, α, →β+, →β-            0+, 1+               0+, 1/2+, 3/2-           ✗ 6Li [L15],  ✗ 60Ni [S05]
7Li + 64Ni → 2·n + 69Ga + 994 keV                       n, →β-                    0+, 3/2-             1/2+, 1/2+, 3/2-
6Li + 62Ni → n + p + 66Zn + 880 keV                     n, →β-                    0+, 1+               0+, 1/2+, 1/2+           ✗ 6Li [L15],  ✗ 62Ni [L15]
7Li + 62Ni → n + p + 67Zn + 681 keV                     n, →β-                    0+, 3/2-             1/2+, 1/2+, 5/2-        ✗ 62Ni [L15]
7Li + 61Ni → n + 12C + 55Mn + 326 keV                   n, →β-                    3/2-, 3/2-           0+, 1/2+, 5/2-          ✗ 61Ni [S05]
6Li + 60Ni → n + p + 64Zn + 258 keV                     n, →β-                    0+, 1+               0+, 1/2+, 1/2+           ✗ 6Li [L15],  ✗ 60Ni [S05]

[L15] 2015 Lugano E-Cat test by Levi et al.
[M96] 1996 Miley and Patterson Infinite Energy reprint
[O97] 1997 Ohmori Fusion Technology transmutation study
[S04] 2004 Savvatimova Ti isotope study
[S05] 2005 Savvatimova isotope study
[U15] 2015 Urutskoev paper in Peter Gluck's blog
[V03] 2003 Violante et al. Ni-hydride plasmon study
%
```

## Leaving the virtual environment
To leave the Python virtualenv environment you set up, run this command:
```
% deactivate
```
