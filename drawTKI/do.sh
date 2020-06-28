which mkexe.sh

mkdir -p output

mkexe.sh drawTKI -lstyle -I../style  -L../style || return 1

for kPiZero in 0 1
do
    ./drawTKI ${kPiZero} > seedrawTKI${kPiZero}.log
done


