for f in *.pdbqt; do b=`basename $f`; echo Processing ligand $b; mkdir -p bs3-1; vina --config donf_bs3-1.txt --ligand $f --out bs3-1/$f --log bs3-1/$f.txt; done
