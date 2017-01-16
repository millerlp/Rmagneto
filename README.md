Rmagneto
R-compatible port of Merlin Oz's C functions to implement magnetometer and accelerometer ellipsoid corrections

	```
	Rmagneto.c, a program to calculate bias corrections for
	magnetometer and accelerometer data using an algorithm derived from
	Li, Q. and J.G. Griffiths (2004). Least squares ellipsoid specific fitting. Geometric Modeling and Processing, 2004. Proceedings, IEEE
	
	as implemented by Merlin Oz 2013.
	https://sites.google.com/site/sailboatinstruments1/home
	Copyright (C) 2013 www.sailboatinstruments.blogspot.com
  
  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
  associated documentation files (the "Software"), to deal in the Software without restriction, including 
  without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
  copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the 
  following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
```

Modified for use with R 3.1.3 by Luke Miller July 2015

The files Rmagneto.c, and Rmagneto.h should be sufficient to compile a Windows dll file Rmagneto.dll for R which 
will implement a single function to be called from R, 'calibrate'. 

// For the calibrate function, from R we will expect to receive pointers to 3 input vectors
// X, Y, Z, which for now will be double precision floating point values.
// norm will be a floating point value for the local magnetic field total norm (or 
// acceleration norm of 1000 milli-g). 
// b will be an empty 1x3 vector to hold the offset corrections output data
// Xcorr will be an empty 1x3 vector to hold the coefficients for the X-axis 
// scaling + soft iron correction
// Ycorr will be an empty 1x3 vector to hold the coefficients for the Y-axis 
// scaling + soft iron correction
// Zcorr will be an empty 1x3 vector to hold the coefficients for the Z-axis 
// scaling + soft iron correction
