source /u/sw/etc/profile
module load gcc-glibc/5
module load getfem
module load qhull
module load boost 
export mkVtkInc=/opt/VTK/include
export mkVtkLib=/opt/VTK/lib
export mkVtkHome=/opt/VTK
export SAMG=/opt/lib/samg
#
export LD_LIBRARY_PATH=$SAMG:$LD_LIBRARY_PATH
export SVD_LICENSE_FILE=@nisserver.mate.polimi.iexport SVD_LICENSE_FILE=@nisserver.mate.polimi.it
# maximum number of threads
export OMP_NUM_THREADS=17
