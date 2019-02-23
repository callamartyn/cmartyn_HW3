## HW3 Skeleton

[![Build
Status](https://travis-ci.org/callamartyn/cmartyn_HW3)](https://travis-ci.org/callamartyn/cmartyn_HW3)

## structure

SW.py contains all of the relevant functions I used to solve Part 1 of the homework.
optimize.py contains my optimization function.


## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `hw3skeleton/__main__.py`) can be run as
follows

```
python -m hw2skeleton -P data
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.
