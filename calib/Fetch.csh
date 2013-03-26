#!/bin/tcsh
#=============================================================================
#+
# NAME:
#   Fetch.csh
#
# PURPOSE:
#   Download all the calibration data you need to run the Pangloss 
#   example.
#
# COMMENTS:
#
# USAGE:
#   Fetch.csh
#
# INPUTS:
#
# OPTIONAL INPUTS:
#   -h --help
#   -v --verbose
#   -x --clobber                      Overwrite data already owned.
#
# OUTPUTS:
#   Millennium/kappa_example.fits     Ray-traced convergence map
#   Millennium/catalog_example.txt    Catalog of galaxy and halo properties
#
#   SHMR/HaloMassRedshift.catalog     Halo mass,z catalog to allow
#                                      empirical halo mass function to 
#                                      be constructed. This is needed by
#                                      the M*-Mh relation code.
#
# DEPENDENCIES:
#
#   wget
#
# BUGS:
#  
# REVISION HISTORY:
#   2013-03-21  started Marshall & Collett (Oxford)
#-
#=======================================================================

unset noclobber

# Set defaults:

set help = 0
set vb = 0
set klobber = 0
set urls = ()

# Parse command line:

while ( $#argv > 0 )
   switch ($argv[1])
   case -h:           #  print help
      set help = 1
      shift argv
      breaksw
   case --{help}:  
      set help = 1
      shift argv
      breaksw
   case -x:           #  clobber
      set klobber = 1
      shift argv
      breaksw
   case --{clobber}:  
      set klobber = 1
      shift argv
      breaksw
   case *:
      shift argv
      breaksw
   endsw
end

#-----------------------------------------------------------------------

if ($help) then
  more $0
  goto FINISH
endif

set BACK = `echo $cwd`

set website = "http://www.ast.cam.ac.uk/~tcollett/Pangloss/calib"

  echo "${0:t}: Downloading data files from Tom's website:"
  echo "${0:t}:   $website"
  echo "${0:t}: into current working directory:"
  echo "${0:t}:   $BACK"
  if ($klobber) echo "${0:t}: Clobbering contents of existing directories!"

#-----------------------------------------------------------------------

set targets = (\
    Millennium/kappa_example.fits \
    Millennium/catalog_example.txt \
    SHMR/HaloMassRedshiftCatalog.pickle \
    ../example/example0_catalog_lightcone.txt\
    ../example/example1_catalog_lightcone.txt\
    ../example/example2_catalog_lightcone.txt\
    ../example/example3_catalog_lightcone.txt\
    ../example/example4_catalog_lightcone.txt\
    ../example/example5_catalog_lightcone.txt\
    ../example/example6_catalog_lightcone.txt\
    ../example/example7_catalog_lightcone.txt\
    ../example/example8_catalog_lightcone.txt\
    ../example/example9_catalog_lightcone.txt\
    )

foreach file ( $targets )

  set dir = $file:h
  mkdir -p $dir
  chdir $dir

  if ($klobber) then

    \rm -rf $file:t
    set done = 0

  else if (-e $file:t) then
  
    echo "${0:t}: $file exists, skipping..."
    set done = 0
    
  else
  
    set now = `date '+%Y%m%d-%H%M%S'`
    set logfile = ".wget.$file:t.$now.log"

    set url = "$website/$file"

    wget "$url" \
      -O $file:t \
      -e robots=off \
      >& $logfile
      
    set done = 1

  endif
    
  chdir $BACK
  
  if ($done) then
    echo "${0:t}: log stored in $dir/$logfile"
    echo "${0:t}: result:"
    du -sh $file
  endif

end

FINISH:

#=======================================================================
