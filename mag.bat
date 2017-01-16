
cd Magneto_c

gcc -m64 -I"D:/R_installs_and_docs/R-3.1.3/include" -O2 -Wall -std=gnu99 -mtune=core2 -c Rmagneto.c -o Rmagneto.o

gcc -m64 -shared -s -static-libgcc -o Rmagneto.dll Rmagneto.o -LD:/R_installs_and_docs/R-3.1.3/bin/x64 -lR