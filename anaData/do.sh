which compile_exe.sh


#####

mkdir -p output

compile_exe.sh anaData -lstyle -I../style  -L../style || return 1

#exit

for ii in 0 1
do 
    for jj in 0 1
    do
        ./anaData $ii $jj >seeana${ii}${jj}.log 
    done
done
