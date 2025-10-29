#!/bin/bash
antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -j 5 -at gaff -c bcc -nc 0
parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod
tleap -s -f tleap.in 

Nvidia-smi pmon -c 