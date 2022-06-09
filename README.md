# StochasticProcessCouplingProject
Set of tools to analyze general coupling of stochastic processes 


These are some simple programs on the study of coupled random walks. All the figures in the article "On the synchronization of coupled random walks and applications to social systems" are reproducible here and many more features are available such as implementation in exotic geometries and movie making through a set of .png images. 

The contents are free and any modification is permitted. I just ask that your modifications to be available to the general public as well. 


This project is written in FORTRAN95 language and o run it you should have the appropriate compiler. I run it with GNU FORTRAN compiler: gfortran, gcc and g++. Download the main archive, uncomment in any text editor the routine you want to run (there are instructions and comments in the text) go to the directory where this archive is saved to compile it. 

I compile such files with:

gfortran -O2 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace Main.f95 -o MainExec 

But one can compile it with just 

gfortran Main.f95

Afterward, if you get no problems, type: 

./MainExec

If you run like me, or type:

./a.out 

if you run the second option


That's all! Have lots of fun!  

