all: fftw-staple
fftw-staple: fftw-staple.c
	gcc -g -o fftw-staple   -L/usr/lib -lImlib2 -lfreetype -lz -L/usr/X11R6/lib -lX11 -lXext -ldl -lm -lfftw3 -lm fftw-staple.c
