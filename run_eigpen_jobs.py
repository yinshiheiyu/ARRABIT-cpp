#!/usr/bin/python
# create job files
import glob
import os 
import math 
import datetime
import shutil
import sys, getopt

def main(argv):
   suffix = '' 
   nev = 0 
   tol = "1e-2"
   ChebyD = 0
   deflate = 0
   post = 1

   synatx = 'run_eigpen_jobs.py -s <suffix> -n <nev> -t <tol> -d <deflate> -c <ChebyD>'
  
   try:
       opts, args = getopt.getopt(argv,"hs:n:t:d:c:p:",["suffix=","nev=",\
                    "tol=","deflate","ChebyD","post"])
   except getopt.GetoptError:
      print synatx
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print synatx 
         sys.exit()
      elif opt in ("-s", "--suffix"):
         suffix = '_'+arg
      elif opt in ("-n", "--nev"):
         nev = arg
      elif opt in ("-t", "--tol"):
         tol = arg
      elif opt in ("-d", "--deflate"):
         deflate = arg
      elif opt in ("-c", "--ChebyD"):
         ChebyD = arg
      elif opt in ("-p", "--pos"):
         post = arg
      else:
         assert False, "unhandled option"

   # suffix is appended to the filenames
   print synatx
   print 'suffix is: ', suffix, ', nev: ', nev, ', tol: ', tol, ', ChebyD: ', ChebyD 
   
   matpath = '/scratch/scratchdirs/cyang/matbin/'
   namelist = ['Andrews', 'C60', 'c_65', 'cfd1', \
   	    'finance', 'Ga10As10H30', 'Ga3As3H12', 'OPF3754',\
               'shallow_water1', 'Si10H16', 'Si5H12', 'SiO', 'wathen100'] 
   
   nevlist = [600, 200, 500, 700, 700, 1000, 600, 200, 800, 200, 200, 400, 300]


   #stol = "%3.0e" %tol
   suffix = suffix + 'k' + str(nev) + 't' + tol

   if deflate > 0:
     suffix = suffix + '_d' + str(deflate)

   suffix = suffix + '_p' + str(post)

   if ChebyD > 0:
     suffix = suffix + '_Che' + str(ChebyD)

   now = datetime.datetime.now()
   dname = str(now)[:10]+suffix
   #dname = str(now)[:10]+'-'+str(now.hour)+':'+str(now.minute)+suffix  
   src = os.getcwd()
   out = 'OUTPUT'
   det1 = os.path.join(src,out)
   det2 = os.path.join(det1,dname)
   dets = os.path.join(det2,'script')
   
   # create output directory
   if not os.path.isdir(det1):
       os.mkdir(det1)
   
   if not os.path.isdir(det2):
       os.mkdir(det2)
   
   if not os.path.isdir(dets):
       os.mkdir(dets)
   
   print "running jobs..."
   
   # change work directory 
   os.chdir(dets) 
   for idx in range(len(namelist)):
   #for idx in range(2):
       name = namelist[idx]
       names = namelist[idx]+suffix
   
       # first create the parameter file for eigpen
       inpfile = name+'.input'
       fo = open(inpfile, 'w')
       fo.write('name     ' + '\'' + matpath+name+'.fhb\'\n')

       # the number of eigenvalues
       if nev == 0: 
           fo.write('nev      ' + str(nevlist[idx]) + '\n')
           fo.write('kp       ' + str(int(0.1*nevlist[idx])) + '\n')
       else:
           fo.write('nev      ' + str(nev) + '\n')
           fo.write('kp       ' + str(int(0.1*float(nev))) + '\n')
 
       fo.write('tol      ' + str(tol)+'\n')
       fo.write('maxit    ' + '10000\n')
       fo.write('shift    ' + '10\n')
       fo.write('print    ' + '1\n')
       fo.write('restart  ' + '2\n')
       fo.write('post   ' + str(post)+'\n')
       fo.write('refine   ' + '0\n')
       fo.write('deflate  ' + str(deflate)+'\n')
       fo.write('ChebyD   ' + str(ChebyD)+'\n')
   
       fo.close()
   
       # then create batch script file
       pbsfile = 'run_'+name+'.pbs'
       
       fo = open(pbsfile, 'w')
       fo.write('# may change debug to regular if runtime exceeds 30 min\n')
       #fo.write('#PBS -q debug\n')
       fo.write('#PBS -q regular\n')
       fo.write('#PBS -N '+name+'\n')
       fo.write('#PBS -o '+name+'.log\n')
       fo.write('#PBS -e '+name+'.err\n')
       fo.write('#PBS -V\n')
       fo.write('#PBS -l mppwidth=24\n')
       fo.write('#PBS -l walltime=6:00:00\n\n\n')
       fo.write('cd $PBS_O_WORKDIR\n')
       fo.write('setenv OMP_NUM_THREADS 24\n')
       fo.write('aprun -n 1 -d 24 '+os.path.join(src,'driveEigpen ')\
          +inpfile + ' > ' + os.path.join(det2,name) + '.out\n' )
       fo.close()

       command = 'qsub '+pbsfile
       print command 
       os.system(command)

   # change work directory 
   os.chdir(src) 

if __name__ == "__main__":
   main(sys.argv[1:])
