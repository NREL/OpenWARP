import matlab.engine
import sys, getopt
#import numpy as np

## GDFMesh(fname,tX,CG,nfobj)
# fname = filename of GDF file in the format 'file.GDF' 
# tX =  translations, tx(1) 
# CG =  position of gravity center, CG(1,3)
# nfobj = target number of panels for AquaPlus Mesh. nfobj(1)

# inputFile
fname=''
tx=''
CG=''
nfobj=''

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"f:t:g:p")
 
###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    	if o == '-f':
        	fname=a
    	elif o == '-t':
        	tX=a
    	elif o == '-g':
        	CG= [0,0,0]
		CG=a
		#CG = np.array([0., 0., 0.])
   	elif o == '-p':
		nfobj=a
        else:
        	print("Usage: %s -i input -o output" % sys.argv[0])
 
# Display input and output file name passed as the args
print ("Input file : %s" % (fname) )
print  fname
eng = matlab.engine.start_matlab()
print "start_matlab"
#eng.GDFMesh(flap-meshed-quad.gdf,1,(1,3),1)
print "GDFMesh0"
fun= eng.GDFMesh(fname,tX,CG,nfobj)
print fun
print "GDFMesh1"

def GDFMesh(fname,tX,CG,nfobj):
	print  fname
	eng = matlab.engine.start_matlab()
	#eng.GDFMesh(flap-meshed-quad.gdf,1,(1,3),1)
	eng.GDFMesh(fname,tX,CG,nfobj)
	return

