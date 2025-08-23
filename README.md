# pychemprojections

[![License](https://img.shields.io/badge/license-BSD-green)](LICENSE.txt)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)

pychemprojections is a python library for visualizing various 2D projections of molecules.

The following projections are implemented:
* Fischer Projection
* Newman Projection
* WedgeDash Projection

## Note

### Fischer Projection

The Fischer projection drawing only supports SMILES which contain one or more chiral carbon atom(s).

It does not support certain SMILES of compounds where cyclic rings are formed with any of chiral carbon atoms present in the main chain.
Therefore, if you would like to draw the Fischer projection of sugars, use the open chain form isomeric SMILES of the compound.

For example, to draw the Fischer projection of Glucose, use the open chain isomeric SMILES which is:

`C([C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O`

instead of the closed chain SMILES:

`C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O`

## Installing
Install the library using pip:

```
pip install git+https://github.com/vandan-revanur/pychemprojections
```

## Examples
Example notebooks for all the projections implemented are available in [examples](examples) directory.

## Contributing
Please follow the [contributing guide](CONTRIBUTING.md) to understand how to contribute to the repo.

## Authors
* **Vandan Revanur**

## License

This project is licensed under the BSD License - see the [LICENSE.txt](LICENSE.txt) file for details
