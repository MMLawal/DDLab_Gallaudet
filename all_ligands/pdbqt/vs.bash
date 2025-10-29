for f in *.pdbqt; do b=`basename $f`; echo Processing ligand $b; mkdir -p pbdplk1; vina --config conf_bs1.txt --ligand $f --out pbdplk1/$f --log pbdplk1/$f.txt; done
