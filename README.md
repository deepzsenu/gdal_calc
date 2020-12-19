	#!/usr/bin/env python
	# -*- coding: utf-8 -*-
#******************************************************************************
#
#  Project:  GDAL
#### Purpose:  Command line raster calculator with numpy syntax
#### Author:   Chris Yesson, chris.yesson@ioz.ac.uk
	#
	#******************************************************************************
- 0	### Copyright (c) 2010, Chris Yesson <chris.yesson@ioz.ac.uk>
- 1	### Copyright (c) 2010-2011, Even Rouault <even dot rouault at mines-paris dot org>
- 2	####  Copyright (c) 2016, Piers Titus van der Torren <pierstitus@gmail.com>
- 3	#
- 4	#####  Permission 
  is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

	##  The above copyright notice and this permission notice shall be included
	###  in all copies or substantial portions of the Software.
-  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
	  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
	 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
	 DEALINGS IN THE SOFTWARE.
   # ******************************************************************************
	
################################################################
# Command line raster calculator with numpy syntax. Use any basic arithmetic supported by numpy arrays such as +-*\ along with logical operators such as >.  Note that all files must have the same dimensions, but no projection checking is performed.  Use gdal_calc.py --help for list of options.

# example 1 - add two files together
# gdal_calc.py -A input1.tif -B input2.tif --outfile=result.tif --calc="A+B"

# example 2 - average of two layers
# gdal_calc.py -A input.tif -B input2.tif --outfile=result.tif --calc="(A+B)/2"
# example 3 - set values of zero and below to null
# gdal_calc.py -A input.tif --outfile=result.tif --calc="A*(A>0)" --NoDataValue=0
################################################################
