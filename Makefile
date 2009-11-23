CC = gcc -g -O3
all: fftw-staple

fftw-staple.o: fftw-staple.c rectangle.h
	$(CC) -c fftw-staple.c

fftw-staple: fftw-staple.o parse_geometry.o
	$(CC) -o fftw-staple   -L/usr/lib -lImlib2 -lfreetype -lz -L/usr/X11R6/lib -lX11 -lXext -ldl -lm -lfftw3 -lm fftw-staple.o parse_geometry.o
