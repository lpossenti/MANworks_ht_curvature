export mkGetfemInc=/home/luca/Documenti/getfem5.2/include/
export mkGetfemLib=/home/luca/Documenti/getfem5.2/lib/
#export mkBoostInc=
export mkQhullLib=/home/luca/Documenti/getfem5.2/lib/
export PATH=${mkGetfemInc}/../bin/:$PATH
export LD_LIBRARY_PATH=${mkQhullLib}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${mkGetfemLib}:${LD_LIBRARY_PATH}

export WITH_SAMG=0
