# Pymol script for representing the sectors of 1IGS.

delete all
load 1IGS.pdb, main
hide all
bg_color white
show cartoon, (chain A)

color white

sele sector_1, (resi 33,35,42,43,,,,,48,50,52,55,60,65,66,70,77,80,88,99,103,111,113,116,118,119,141,142,149,155,158,162,166,169,172,173,174,177,190,196,199,201,207,213,218,227,228,229,236,,246)& (chain A)
color red, sector_1



sele sector_2, (resi 33,45,,,,51,53,59,,73,,81,82,83,93,95,100,109,110,113,117,118,121,123,124,128,131,135,145,148,152,157,159,160,161,165,178,180,182,184,185,189,202,208,210,211,218,231,233,234,)& (chain A)
color blue, sector_2

