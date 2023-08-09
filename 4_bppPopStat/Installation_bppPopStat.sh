################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
###########################################################################
#   Steps for installation of BppPopStat on Linux-HPC and macOS           #
#                                                                         #
# Files needed:                                                           #
#             * Git installed                                             #
#                                                                         #
###########################################################################
#
#
#
######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

###############
#  Linux-HPC  #
###############
# BppPopStat installation
cd /home/pereira/software/BPP_danilo

module load cmake/3.10.2 
module load gcc/6.2.1
module load perl/5.30.1

git clone https://github.com/BioPP/bpp-core.git
cd bpp-core
cmake -DCMAKE_INSTALL_PREFIX=/home/pereira/software/BPP_danilo/
make install

git clone https://github.com/BioPP/bpp-seq.git
cd bpp-seq
cmake -DCMAKE_INSTALL_PREFIX=/home/pereira/software/BPP_danilo/
make install

git clone https://github.com/BioPP/bpp-phyl.git
cd bpp-phyl
cmake -DCMAKE_INSTALL_PREFIX=/home/pereira/software/BPP_danilo/
make install

git clone https://github.com/BioPP/bpp-popgen.git
cd bpp-popgen
cmake -DCMAKE_INSTALL_PREFIX=/home/pereira/software/BPP_danilo/
make install

git clone https://github.com/BioPP/bppsuite.git
cd bppsuite
cmake -DCMAKE_INSTALL_PREFIX=/home/pereira/software/BPP_danilo/ ./
make # compile
make install # move files to the installation directory (this will create a $bpp_dir/bin/ directory)

# since the libraries were installed in another place then the default (here it was /home/pereira/software/BPP_danilo), the path to each library needs to be added to the .bash_profile permanently, or exported before usage. My .bash_profile is checked every time I open wallace, and I added into it the following:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pereira/software/BPP_danilo:/home/pereira/software/BPP_danilo/bpp-phyl/src:/home/pereira/software/BPP_danilo/bpp-popgen/src:/home/pereira/software/BPP_danilo/bpp-seq/src:/home/pereira/software/BPP_danilo/bpp-core/src


###########
#  macOS  #
###########

# local mac
git clone https://github.com/BioPP/bpp-core.git
cd bpp-core
cmake -DCMAKE_INSTALL_PREFIX=/Users/danilo/Documents/Git_Repository/BPP_danilo/
make install

git clone https://github.com/BioPP/bpp-seq.git
cd bpp-seq
cmake -DCMAKE_INSTALL_PREFIX=/Users/danilo/Documents/Git_Repository/BPP_danilo/
make install

git clone https://github.com/BioPP/bpp-phyl.git
cd bpp-phyl
cmake -DCMAKE_INSTALL_PREFIX=/Users/danilo/Documents/Git_Repository/BPP_danilo/
make install

git clone https://github.com/BioPP/bpp-popgen.git
cd bpp-popgen
cmake -DCMAKE_INSTALL_PREFIX=/Users/danilo/Documents/Git_Repository/BPP_danilo/
make install

git clone https://github.com/BioPP/bppsuite.git
cd bppsuite
cmake -DCMAKE_INSTALL_PREFIX=/Users/danilo/Documents/Git_Repository/BPP_danilo/ ./
make # compile
make install # move files to the installation directory (this will create a $bpp_dir/bin/ directory)

# since the libraries were installed in another place then the default (here it was /Users/danilo/Documents/Git_Repository/BPP_danilo/), the path to each library needs to be added to the .bash_profile permanently, or exported before usage. My .bash_profile is checked every time I open wallace, and I added into it the following:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/danilo/Documents/Git_Repository/BPP_danilo/:/Users/danilo/Documents/Git_Repository/BPP_danilo/bpp-phyl/src://Users/danilo/Documents/Git_Repository/BPP_danilo/bpp-popgen/src:/Users/danilo/Documents/Git_Repository/BPP_danilo/bpp-seq/src:/Users/danilo/Documents/Git_Repository/BPP_danilo/bpp-core/src

git clone https://github.com/BioPP/grapes
cd grapes
cmake -DCMAKE_INSTALL_PREFIX=/Users/danilo/Documents/Git_Repository/BPP_danilo/
make install

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################
