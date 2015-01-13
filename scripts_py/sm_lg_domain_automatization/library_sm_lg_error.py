#!/usr/bin/python

'''
@description
useful functions to generate the error files that are computed
as the relative difference b/w the small and large domain
simulations
'''


import os          #for os based functions (cp,mv,mkdir...)
import sys         #for interactign with the system (printing error msg)
import subprocess  #for running processes in shell
import shutil      #for additional os based functions (rm -rf)
import fnmatch

# compile the executable able to create the error files
def generate_error_exe():
    '''
    @description
    compile the executable able to create the error files
    '''

    
    #1) verify that the fortran file to be compiled exists
    exePath = os.getenv('augeanstables')
    exePath+= '/scripts_fortran/error_computation'
    exePath+= '/compute_error_file'
    if(not os.path.isfile(exePath+'.f')):
        print 'library_sm_lg_error'
        print 'generate_error_exe'
        sys.exit('the fortran file '+exePath+' does not exist')


    #2) modify the system environment variable to
    #   compile the executable with optimizations
    os.putenv('AUGEANSTABLES_PROFILE','false')
    
    
    #3) compile the fortran file
    print 'generating the executable to compute the error: ...'
    cmd = 'cd '+os.path.dirname(exePath)+' && '
    cmd+= 'make cleanall > /dev/null && '
    cmd+= 'make '+os.path.basename(exePath)+' 1> /dev/null 2>&1 && '
    cmd+= 'make clean > /dev/null'
    subprocess.call(cmd, shell=True)

    #4) check whether the executable has been generated
    if(not os.path.isfile(exePath)):
        print 'library_sm_lg_error'
        print 'generate_error_exe'
        sys.exit('***error when generating the error executable***')

    print 'generating the executable to compute the error: done'
    return exePath


# compare two data files and generate an error file
def compare_data(exePath,
                 dataPath_sm_domain,
                 dataPath_lg_domain,
                 errorPath):
    '''
    @description
    compare two data files and generate an error file
    '''
    
    
    #1) generate the error file using the executable
    cmd = exePath
    cmd+=' -s '+dataPath_sm_domain
    cmd+=' -l '+dataPath_lg_domain
    cmd+=' -o '+errorPath

    subprocess.call(cmd, shell=True)


    #2) verify that the output error file has been generated
    if(not os.path.isfile(errorPath)):
        print 'library_sm_lg_error'
        print 'compare_data'
        sys.exit('error when generating the error file '+errorPath)

    return


# get the total number of data files in folder
def get_nb_data_files(cdir):

    nb_data_files = 0

    for file in os.listdir(cdir):
        if fnmatch.fnmatch(file, 'data[0-9]*.nc'):
            
            step = int(file.split('data')[1].split('.nc')[0])
            if(step>nb_data_files):
                nb_data_files = step

    return nb_data_files

    
# generate error files by comparing the files saved in the
# folders storing the small domain and the large domain
# simulation results
def compare_folders(exePath,
                    dataDir_sm_domain,
                    dataDir_lg_domain):
    '''
    @description
    generate error files by comparing the files saved in the
    folders storing the small domain and the large domain
    simulation results
    '''


    #1) verify that the folders exist for the small and large
    #   domain simulations
    if(not os.path.isdir(dataDir_sm_domain)):
        print 'library_sm_lg_error'
        print 'compare_folders'
        sys.exit('folder '+dataDir_sm_domain+'does not exist')

    if(not os.path.isdir(dataDir_lg_domain)):
        print 'library_sm_lg_error'
        print 'compare_folders'
        sys.exit('folder '+dataDir_lg_domain+'does not exist')


    #2) create the directory for the error files
    errorDir = os.path.abspath(os.path.join(dataDir_sm_domain,'..','error'))
    
    if(os.path.isdir(errorDir)):
        shutil.rmtree(errorDir)

    os.mkdir(errorDir)


    #3) get the number of data files in small and large domain
    nbFiles_sm_domain = get_nb_data_files(dataDir_sm_domain)
    nbFiles_lg_domain = get_nb_data_files(dataDir_lg_domain)

    if(nbFiles_sm_domain != nbFiles_lg_domain):
        print 'library_sm_lg_error'
        print 'compare_folders'
        sys.exit('different number of files in small and large folders')
        
    print 'generating error files: ...'
    print str(nbFiles_sm_domain)+' file to be processed'

    #4) generate the error files for matching files
    for i in range(0,min(nbFiles_sm_domain,nbFiles_lg_domain)+1):
        
        dataPath_sm_domain = dataDir_sm_domain+'/data'+str(i)+'.nc'
        dataPath_lg_domain = dataDir_lg_domain+'/data'+str(i)+'.nc'
        errorPath          = errorDir+'/error'+str(i)+'.nc'

        compare_data(exePath,
                     dataPath_sm_domain,
                     dataPath_lg_domain,
                     errorPath)

    print 'generating error files: done'


if __name__=='__main__':

    #test: generate_error_exe
    exePath = generate_error_exe()
    print 'errorPath: ', exePath


    #test: compare_data
    mainDir = '/home/jdesmarais/projects'
    mainDir+= '/dim2d_0.995_0.1'

    dataPath_sm_domain = mainDir+'/sm_domain/data0.nc'
    dataPath_lg_domain = mainDir+'/lg_domain/data0.nc'
    errorPath= mainDir+'/error0.nc'

    compare_data(exePath,
                 dataPath_sm_domain,
                 dataPath_lg_domain,
                 errorPath)


    #test: compare_folders
    dataDir_sm_domain = mainDir+'/sm_domain'
    dataDir_lg_domain = mainDir+'/lg_domain'
    
    compare_folders(exePath,
                    dataDir_sm_domain,
                    dataDir_lg_domain)
