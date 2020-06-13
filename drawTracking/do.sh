which compile_exe.sh

mkdir -p output

compile_exe.sh drawTracking -lstyle -I../style  -L../style || return 1

for kPiZero in 0 1
do 
    for kProton in 0 1
    do 
        ./drawTracking $kPiZero $kProton > seedrawTracking${kPiZero}${kProton}.log
    done
done

