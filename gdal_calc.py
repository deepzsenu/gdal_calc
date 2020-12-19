	################################################################
from optparse import OptionParser, Values
import os
import os.path
import sys
import numpy

from osgeo import gdal
from osgeo import gdalnumeric

# create alphabetic list for storing input layers
AlphaList=["A","B","C","D","E","F","G","H","I","J","K","L","M",
           "N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
	
# set up some default nodatavalues for each datatype
DefaultNDVLookup={'Byte':255, 'UInt16':65535, 'Int16':-32767, 'UInt32':4294967293, 'Int32':-2147483647, 'Float32':3.402823466E+38, 'Float64':1.7976931348623158E+308}
	
def DoesDriverHandleExtension(drv, ext):
    exts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
    return exts is not None and exts.lower().find(ext.lower()) >= 0

def GetExtension(filename):
    ext = os.path.splitext(filename)[1]
    if ext.startswith('.'):
        ext = ext[1:]
    return ext
	
def GetOutputDriversFor(filename):
    drv_list = []
    ext = GetExtension(filename)
    for i in range(gdal.GetDriverCount()):
        drv = gdal.GetDriver(i)
        if (drv.GetMetadataItem(gdal.DCAP_CREATE) is not None or \
            drv.GetMetadataItem(gdal.DCAP_CREATECOPY) is not None) and \
           drv.GetMetadataItem(gdal.DCAP_RASTER) is not None:
            if len(ext) > 0 and DoesDriverHandleExtension(drv, ext):
                drv_list.append( drv.ShortName )
            else:
                prefix = drv.GetMetadataItem(gdal.DMD_CONNECTION_PREFIX)
                if prefix is not None and filename.lower().startswith(prefix.lower()):
                   drv_list.append( drv.ShortName )

    # GMT is registered before netCDF for opening reasons, but we want
    # netCDF to be used by default for output.
    if ext.lower() == 'nc' and len(drv_list) == 0 and \
       drv_list[0].upper() == 'GMT' and drv_list[1].upper() == 'NETCDF':
           drv_list = [ 'NETCDF', 'GMT' ]

    return drv_list
    
    
def GetOutputDriverFor(filename):
    drv_list = GetOutputDriversFor(filename)
    if len(drv_list) == 0:
      ext = GetExtension(filename)
	      if len(ext) == 0:
	            return 'GTiff'
	       else:
            raise Exception("Cannot guess driver for %s" % filename)
    elif len(drv_list) > 1:
        print("Several drivers matching %s extension. Using %s" % (ext, drv_list[0]))
    return drv_list[0]

################################################################
def doit(opts, args):

    if opts.debug:
        print("gdal_calc.py starting calculation %s" %(opts.calc))

    # set up global namespace for eval with all functions of gdalnumeric
    global_namespace = dict([(key, getattr(gdalnumeric, key))
        for key in dir(gdalnumeric) if not key.startswith('__')])

    if not opts.calc:
        raise Exception("No calculation provided.")
    elif not opts.outF:
        raise Exception("No output file provided.")

    if opts.format is None:
        opts.format = GetOutputDriverFor(opts.outF)

    ################################################################
    # fetch details of input layers
    ################################################################

    # set up some lists to store data for each band
    myFiles=[]
    myBands=[]
    myAlphaList=[]
    myDataType=[]
    myDataTypeNum=[]
    myNDV=[]
    DimensionsCheck=None
1
    # loop through input files - checking dimensions
141	    for myI, myF in opts.input_files.items():
142	        if not myI.endswith("_band"):
143	            # check if we have asked for a specific band...
144	            if "%s_band" % myI in opts.input_files:
145	                myBand = opts.input_files["%s_band" % myI]
146	            else:
147	                myBand = 1
148	
149	            myFile = gdal.Open(myF, gdal.GA_ReadOnly)
150	            if not myFile:
151	                raise IOError("No such file or directory: '%s'" % myF)
152	
153	            myFiles.append(myFile)
154	            myBands.append(myBand)
155	            myAlphaList.append(myI)
156	            myDataType.append(gdal.GetDataTypeName(myFile.GetRasterBand(myBand).DataType))
157	            myDataTypeNum.append(myFile.GetRasterBand(myBand).DataType)
158	            myNDV.append(myFile.GetRasterBand(myBand).GetNoDataValue())
159	            # check that the dimensions of each layer are the same
160	            if DimensionsCheck:
161	                if DimensionsCheck != [myFile.RasterXSize, myFile.RasterYSize]:
162	                    raise Exception("Error! Dimensions of file %s (%i, %i) are different from other files (%i, %i).  Cannot proceed" % \
163	                            (myF, myFile.RasterXSize, myFile.RasterYSize, DimensionsCheck[0], DimensionsCheck[1]))
164	            else:
165	                DimensionsCheck = [myFile.RasterXSize, myFile.RasterYSize]
166	
167	            if opts.debug:
168	                print("file %s: %s, dimensions: %s, %s, type: %s" %(myI,myF,DimensionsCheck[0],DimensionsCheck[1],myDataType[-1]))
169	
170	    # process allBands option
171	    allBandsIndex=None
172	    allBandsCount=1
173	    if opts.allBands:
174	        try:
175	            allBandsIndex=myAlphaList.index(opts.allBands)
176	        except ValueError:
177	            raise Exception("Error! allBands option was given but Band %s not found.  Cannot proceed" % (opts.allBands))
178	        allBandsCount=myFiles[allBandsIndex].RasterCount
179	        if allBandsCount <= 1:
180	            allBandsIndex=None
181	
182	    ################################################################
183	    # set up output file
184	    ################################################################
185	
186	    # open output file exists
187	    if os.path.isfile(opts.outF) and not opts.overwrite:
188	        if allBandsIndex is not None:
189	            raise Exception("Error! allBands option was given but Output file exists, must use --overwrite option!")
190	        if opts.debug:
191	            print("Output file %s exists - filling in results into file" %(opts.outF))
192	        myOut=gdal.Open(opts.outF, gdal.GA_Update)
193	        if [myOut.RasterXSize,myOut.RasterYSize] != DimensionsCheck:
194	            raise Exception("Error! Output exists, but is the wrong size.  Use the --overwrite option to automatically overwrite the existing file")
195	        myOutB=myOut.GetRasterBand(1)
196	        myOutNDV=myOutB.GetNoDataValue()
197	        myOutType=gdal.GetDataTypeName(myOutB.DataType)
198	
199	    else:
200	        # remove existing file and regenerate
201	        if os.path.isfile(opts.outF):
202	            os.remove(opts.outF)
203	        # create a new file
204	        if opts.debug:
205	            print("Generating output file %s" %(opts.outF))
206	
207	        # find data type to use
208	        if not opts.type:
209	            # use the largest type of the input files
210	            myOutType=gdal.GetDataTypeName(max(myDataTypeNum))
211	        else:
212	            myOutType=opts.type
213	
214	        # create file
215	        myOutDrv = gdal.GetDriverByName(opts.format)
216	        myOut = myOutDrv.Create(
217	            opts.outF, DimensionsCheck[0], DimensionsCheck[1], allBandsCount,
218	            gdal.GetDataTypeByName(myOutType), opts.creation_options)
219	
220	        # set output geo info based on first input layer
221	        myOut.SetGeoTransform(myFiles[0].GetGeoTransform())
222	        myOut.SetProjection(myFiles[0].GetProjection())
223	
224	        if opts.NoDataValue!=None:
225	            myOutNDV=opts.NoDataValue
226	        else:
227	            myOutNDV=DefaultNDVLookup[myOutType]
228	
229	        for i in range(1,allBandsCount+1):
230	            myOutB=myOut.GetRasterBand(i)
231	            myOutB.SetNoDataValue(myOutNDV)
232	            # write to band
233	            myOutB=None
234	
235	    if opts.debug:
236	        print("output file: %s, dimensions: %s, %s, type: %s" %(opts.outF,myOut.RasterXSize,myOut.RasterYSize,myOutType))
237	
238	    ################################################################
239	    # find block size to chop grids into bite-sized chunks
240	    ################################################################
241	
242	    # use the block size of the first layer to read efficiently
243	    myBlockSize=myFiles[0].GetRasterBand(myBands[0]).GetBlockSize();
244	    # store these numbers in variables that may change later
245	    nXValid = myBlockSize[0]
246	    nYValid = myBlockSize[1]
247	    # find total x and y blocks to be read
248	    nXBlocks = (int)((DimensionsCheck[0] + myBlockSize[0] - 1) / myBlockSize[0]);
249	    nYBlocks = (int)((DimensionsCheck[1] + myBlockSize[1] - 1) / myBlockSize[1]);
250	    myBufSize = myBlockSize[0]*myBlockSize[1]
251	
252	    if opts.debug:
253	        print("using blocksize %s x %s" %(myBlockSize[0], myBlockSize[1]))
254	
255	    # variables for displaying progress
256	    ProgressCt=-1
257	    ProgressMk=-1
258	    ProgressEnd=nXBlocks*nYBlocks*allBandsCount
259	
260	    ################################################################
261	    # start looping through each band in allBandsCount
262	    ################################################################
263	
264	    for bandNo in range(1,allBandsCount+1):
265	
266	        ################################################################
267	        # start looping through blocks of data
268	        ################################################################
269	
270	        # loop through X-lines
271	        for X in range(0,nXBlocks):
272	
273	            # in the rare (impossible?) case that the blocks don't fit perfectly
274	            # change the block size of the final piece
275	            if X==nXBlocks-1:
276	                nXValid = DimensionsCheck[0] - X * myBlockSize[0]
277	                myBufSize = nXValid*nYValid
278	
279	            # find X offset
280	            myX=X*myBlockSize[0]
281	
282	            # reset buffer size for start of Y loop
283	            nYValid = myBlockSize[1]
284	            myBufSize = nXValid*nYValid
285	
286	            # loop through Y lines
287	            for Y in range(0,nYBlocks):
288	                ProgressCt+=1
289	                if 10*ProgressCt/ProgressEnd%10!=ProgressMk and not opts.quiet:
290	                    ProgressMk=10*ProgressCt/ProgressEnd%10
291	                    from sys import version_info
292	                    if version_info >= (3,0,0):
293	                        exec('print("%d.." % (10*ProgressMk), end=" ")')
294	                    else:
295	                        exec('print 10*ProgressMk, "..",')
296	
297	                # change the block size of the final piece
298	                if Y==nYBlocks-1:
299	                    nYValid = DimensionsCheck[1] - Y * myBlockSize[1]
300	                    myBufSize = nXValid*nYValid
301	
302	                # find Y offset
303	                myY=Y*myBlockSize[1]
304	
305	                # create empty buffer to mark where nodata occurs
306	                myNDVs = None
307	
308	                # make local namespace for calculation
309	                local_namespace = {}
310	
311	                # fetch data for each input layer
312	                for i,Alpha in enumerate(myAlphaList):
313	
314	                    # populate lettered arrays with values
315	                    if allBandsIndex is not None and allBandsIndex==i:
316	                        myBandNo=bandNo
317	                    else:
318	                        myBandNo=myBands[i]
319	                    myval=gdalnumeric.BandReadAsArray(myFiles[i].GetRasterBand(myBandNo),
320	                                          xoff=myX, yoff=myY,
321	                                          win_xsize=nXValid, win_ysize=nYValid)
322	
323	                    # fill in nodata values
324	                    if myNDV[i] is not None:
325	                        if myNDVs is None:
326	                            myNDVs = numpy.zeros(myBufSize)
327	                            myNDVs.shape=(nYValid,nXValid)
328	                        myNDVs=1*numpy.logical_or(myNDVs==1, myval==myNDV[i])
329	
330	                    # add an array of values for this block to the eval namespace
331	                    local_namespace[Alpha] = myval
332	                    myval=None
333	
334	
335	                # try the calculation on the array blocks
336	                try:
337	                    myResult = eval(opts.calc, global_namespace, local_namespace)
338	                except:
339	                    print("evaluation of calculation %s failed" %(opts.calc))
340	                    raise
341	
342	                # Propagate nodata values (set nodata cells to zero
343	                # then add nodata value to these cells).
344	                if myNDVs is not None:
345	                    myResult = ((1*(myNDVs==0))*myResult) + (myOutNDV*myNDVs)
346	                elif not isinstance(myResult, numpy.ndarray):
347	                    myResult = numpy.ones( (nYValid,nXValid) ) * myResult
348	
349	                # write data block to the output file
350	                myOutB=myOut.GetRasterBand(bandNo)
351	                gdalnumeric.BandWriteArray(myOutB, myResult, xoff=myX, yoff=myY)
352	
353	    if not opts.quiet:
354	        print("100 - Done")
355	    #print("Finished - Results written to %s" %opts.outF)
356	
357	    return
358	
359	################################################################
360	def Calc(calc, outfile, NoDataValue=None, type=None, format=None, creation_options=[], allBands='', overwrite=False, debug=False, quiet=False, **input_files):
361	    """ Perform raster calculations with numpy syntax.
362	    Use any basic arithmetic supported by numpy arrays such as +-*\ along with logical
363	    operators such as >. Note that all files must have the same dimensions, but no projection checking is performed.
364	
365	    Keyword arguments:
366	        [A-Z]: input files
367	        [A_band - Z_band]: band to use for respective input file
368	
369	    Examples:
370	    add two files together:
371	        Calc("A+B", A="input1.tif", B="input2.tif", outfile="result.tif")
372	
373	    average of two layers:
374	        Calc(calc="(A+B)/2", A="input1.tif", B="input2.tif", outfile="result.tif")
375	
376	    set values of zero and below to null:
377	        Calc(calc="A*(A>0)", A="input.tif", A_Band=2, outfile="result.tif", NoDataValue=0)
378	    """
379	    opts = Values()
380	    opts.input_files = input_files
381	    opts.calc = calc
382	    opts.outF = outfile
383	    opts.NoDataValue = NoDataValue
384	    opts.type = type
385	    opts.format = format
386	    opts.creation_options = creation_options
387	    opts.allBands = allBands
388	    opts.overwrite = overwrite
389	    opts.debug = debug
390	    opts.quiet = quiet
391	
392	    doit(opts, None)
393	
394	def store_input_file(option, opt_str, value, parser):
395	    if not hasattr(parser.values, 'input_files'):
396	        parser.values.input_files = {}
397	    parser.values.input_files[opt_str.lstrip('-')] = value
398	
399	def main():
400	    usage = """usage: %prog --calc=expression --outfile=out_filename [-A filename]
401	                    [--A_band=n] [-B...-Z filename] [other_options]"""
402	    parser = OptionParser(usage)
403	
404	    # define options
405	    parser.add_option("--calc", dest="calc", help="calculation in gdalnumeric syntax using +-/* or any numpy array functions (i.e. log10())", metavar="expression")
406	    # limit the input file options to the ones in the argument list
407	    given_args = set([a[1] for a in sys.argv if a[1:2] in AlphaList] + ['A'])
408	    for myAlpha in given_args:
409	        parser.add_option("-%s" % myAlpha, action="callback", callback=store_input_file, type=str, help="input gdal raster file, you can use any letter (A-Z)", metavar='filename')
410	        parser.add_option("--%s_band" % myAlpha, action="callback", callback=store_input_file, type=int, help="number of raster band for file %s (default 1)" % myAlpha, metavar='n')
411	
412	    parser.add_option("--outfile", dest="outF", help="output file to generate or fill", metavar="filename")
413	    parser.add_option("--NoDataValue", dest="NoDataValue", type=float, help="output nodata value (default datatype specific value)", metavar="value")
414	    parser.add_option("--type", dest="type", help="output datatype, must be one of %s" % list(DefaultNDVLookup.keys()), metavar="datatype")
415	    parser.add_option("--format", dest="format", help="GDAL format for output file", metavar="gdal_format")
416	    parser.add_option(
417	        "--creation-option", "--co", dest="creation_options", default=[], action="append",
418	        help="Passes a creation option to the output format driver. Multiple "
419	        "options may be listed. See format specific documentation for legal "
420	        "creation options for each format.", metavar="option")
421	    parser.add_option("--allBands", dest="allBands", default="", help="process all bands of given raster (A-Z)", metavar="[A-Z]")
422	    parser.add_option("--overwrite", dest="overwrite", action="store_true", help="overwrite output file if it already exists")
423	    parser.add_option("--debug", dest="debug", action="store_true", help="print debugging information")
424	    parser.add_option("--quiet", dest="quiet", action="store_true", help="suppress progress messages")
425	
426	    (opts, args) = parser.parse_args()
427	    if not hasattr(opts, "input_files"):
428	        opts.input_files = {}
429	
430	    if len(sys.argv) == 1:
431	        parser.print_help()
432	        sys.exit(1)
433	    elif not opts.calc:
434	        print("No calculation provided. Nothing to do!")
435	        parser.print_help()
436	        sys.exit(1)
437	    elif not opts.outF:
438	        print("No output file provided. Cannot proceed.")
439	        parser.print_help()
440	        sys.exit(1)
441	    else:
442	        try:
443	            doit(opts, args)
444	        except IOError as e:
445	            print(e)
446	            sys.exit(1)
447	
448	
449	if __name__ == "__main__":
450	    main()
