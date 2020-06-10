which compile_exe.sh

mkdir -p output

compile_exe.sh drawTracking -lstyle -I../style  -L../style || return 1

for ii in 0 1
do 
    for jj in 0 1
    do 
        ./drawTracking $ii $jj > seedrawTracking${ii}${jj}.log
    done
done

