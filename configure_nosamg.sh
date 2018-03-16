export mkGetfemInc=/home/daniele/software/getfem/getfem/include/
export mkGetfemLib=/home/daniele/software/getfem/getfem/lib/
#export mkBoostInc=
export mkQhullLib=/home/daniele/software/getfem/qhull/lib/
export PATH=${mkGetfemInc}/../bin/:$PATH
export LD_LIBRARY_PATH=${mkQhullLib}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${mkGetfemLib}:${LD_LIBRARY_PATH}

export WITH_SAMG=0
