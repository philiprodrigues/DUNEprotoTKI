which mkexe.sh

mkdir -p output

mkexe.sh drawTracking -lstyle -I../style  -L../style || return 1

for kPiZero in 0 1
do 
    for kTrackingProton in 0 1
    do 
        ./drawTracking $kPiZero $kTrackingProton > seedrawTracking${kPiZero}${kTrackingProton}.log
    done
done

