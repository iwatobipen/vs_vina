# vs_vina

## requirements

- [AutoDock Vina](https://anaconda.org/bioconda/autodock-vina)
- [ADFR](https://ccsb.scripps.edu/projects/docking/)
- [pymol](https://pymol.org/2/)
- [rdkit](http://rdkit.org/)
- [pdb-tools](http://www.bonvinlab.org/pdb-tools/)


## run selfdocling workflow

``` shell
# to run following command, all packages which are described in requirements should be installed
$ git clone gh repo clone iwatobipen/vs_vina
$ cd pdbtool
$ pdb_fetch 1iep
$ python auto_dock_wf.py 1iep.pdb --ligand STI
```


## To Do

- [x] Write example code for self docking
- [ ] Write example code for cross docking
- [ ] add option for vina setting
- [ ] analysis docking score