# Crodinger
High order NLSE solvers in C

This is a port of my SchrodingerMAT library. Work in progress, more details to be added.

Important Note: I have had to rewrite virtually all my code using an extremely ugly approach so that I can easily run it in double/extended/quad precision with minimal modification. The code is, thus, very hard to read and debug. This is NOT my normal coding style and I do not encourage anyone to do that in general, I simply had no other option due to time constraints and I didn't want to maintain 3 different copies of each file.

As a summary, here's how it works.
double --> double precision.
long double --> extended precision (if supported on machine and by compiler, C standard does not guarantee it).
__float128 --> quad precision using GCC 4.6 or higher.

Essentially, I wrote a bunch of macros:
MYTYPE --> double/long double/__float128, the precision of the data
MYFFTW(func) --> This converts fftw functions into the version matching the selected precision
e.g. MYFFTW(_malloc) --> fftwl_malloc if long double is chosen.
D(num) --> this converts floating point literals to make sure they are correct.
e.g. D(2.213) = 2.213L is long double is chosen.

C doesn't support overloading so I didn't have much of an option here. I'm even considering a full rewrite in C++ because this is just ugly.
