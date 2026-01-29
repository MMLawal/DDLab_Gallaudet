#start with plumed, if you have previously installed Amber, uninstall it by going to 'build directory' and excuting 'make clean && make distclean' then rm -rf amberfolder

cd $HOME
wget https://github.com/plumed/plumed2/releases/download/v2.10.0/plumed-2.10.0.tgz
tar -xvzf plumed-2.10.0.tgz
cd plumed-2.10.0

./configure --prefix=$HOME/plumed --enable-modules=all --enable-mpi
make -j4
make install

#add to .bashrc
export PATH=$HOME/plumed/bin:$PATH
export LD_LIBRARY_PATH=$HOME/plumed/lib:$LD_LIBRARY_PATH

#Amber24 installation
Note: https://ambermd.org/GetAmber.php named it as pmemd24, you may want to rename (amber24) for convenience
put thetar.bz2 in you home directory

cd $HOME
tar -xvjf Amber24.tar.bz2
cd amber24_src/build
edit run_cmake MPI=TRUE, DPLUMED for plumed

  cmake $AMBER_PREFIX/amber24_src \
    -DCMAKE_INSTALL_PREFIX=$AMBER_PREFIX/amber24 \
    -DCOMPILER=GNU  \
    -DMPI=TRUE -DCUDA=TRUE -DINSTALL_TESTS=TRUE \
    -DDOWNLOAD_MINICONDA=FALSE -DBUILD_PYTHON=FALSE \
    -DBUILD_PERL=FALSE -DBUILD_GUI=FALSE \
    -DPMEMD_ONLY=TRUE -DCHECK_UPDATES=FALSE \
    -DPLUMED_LIBRARY=/u/username/plumed/lib/libplumed.so \
    -DPLUMED_KERNEL_LIBRARY=/u/username/plumed/lib/libplumedKernel.so \
    -DPLUMED_INCLUDES=/u/username/plumed/include \
    -DPLUMED_WORKS=/u/username/plumed \
    2>&1 | tee  cmake.log

execute ./run_cmake, then
make install

put 'source /u/username/amber24/amber.sh' in ~/.bashrc and execute source ~/.bashrc
