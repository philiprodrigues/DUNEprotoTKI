which compile_exe.sh

mkdir -p output

compile_exe.sh drawTKI -lstyle -I../style  -L../style || return 1

./drawTKI 0 > seedrawTKI0.log
./drawTKI 1 > seedrawTKI1.log
