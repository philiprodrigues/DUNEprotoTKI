code=$1

rm -f lib${code}.so  *.o *Dict.* 

incpath="-I$(root-config --incdir) -I$(pwd)"
echo $incpath

echo 
comm="g++ $incpath -fPIC -c ${code}.cxx -o ${code}.o"
echo $comm
eval $comm || exit 1

rootver=$(root-config --version | awk -F\. '{print $1}')
if [ ${rootver} == 5 ]
then
    rcopt="-c"
elif [ ${rootver} == 6 ]
then 
    rcopt=""
else
    echo unknown root version
    exit 1
fi

comm="rootcint -f ${code}Dict.cxx ${rcopt} $incpath  ${code}.h ${code}LinkDef.h"
echo $comm
eval $comm  || exit 1

comm="g++ $incpath  -fPIC -c ${code}Dict.cxx -o ${code}Dict.o"
echo $comm
eval $comm || exit 1


comm="g++ -shared -O3 -Wall -Werror $incpath $(root-config --libs) ${code}.o ${code}Dict.o -o lib${code}.so"
echo $comm
eval $comm
    
rm -f  ${code}.o ${code}Dict.h ${code}Dict.cxx ${code}Dict.o
#${code}Dict_rdict.pcm #this is needed!

echo
echo "=========================     ${code} done!     ========================="
echo

#ln -s $(pwd)/lib${code}.so ../lib/
