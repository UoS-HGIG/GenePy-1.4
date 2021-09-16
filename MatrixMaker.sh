#!/bin/bash
module load conda/py2-latest

cat $1/* > temp1
echo "CAT file created"
#python /dics/em3e17/software/generate_final_matrix.py temp1 temp2
python /mainfs/hgig/private/software/GENEPY_v1.3/utils/generate_final_matrix.py  temp1 $2

echo "MATRIX GENERATED"
#python /dics/em3e17/software/normalise_matrices.py temp2 temp3
#echo "MATRIX NORMALISED"
#python /dics/em3e17/software/gdi_scale_iri.py temp3 $2
#echo "MATRIX GDI SCALED"
rm temp1
#rm temp1 temp2 temp3

