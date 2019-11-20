FROM ubuntu

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install gcc g++ gfortran wget cpio zlib1g-dev && \
  cd /tmp && \
  wget -q http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/15275/l_mkl_2019.3.199.tgz && \
  tar -xzf l_mkl_2019.3.199.tgz && \
  cd l_mkl_2019.3.199 && \
  sed -i 's/ACCEPT_EULA=decline/ACCEPT_EULA=accept/g' silent.cfg && \
  sed -i 's/ARCH_SELECTED=ALL/ARCH_SELECTED=INTEL64/g' silent.cfg && \
  sed -i 's/COMPONENTS=DEFAULTS/COMPONENTS=;intel-comp-l-all-vars__noarch;intel-comp-nomcu-vars__noarch;intel-openmp__x86_64;intel-tbb-libs__x86_64;intel-mkl-common__noarch;intel-mkl-installer-license__noarch;intel-mkl-core__x86_64;intel-mkl-core-rt__x86_64;intel-mkl-doc__noarch;intel-mkl-doc-ps__noarch;intel-mkl-gnu__x86_64;intel-mkl-gnu-rt__x86_64;intel-mkl-common-ps__noarch;intel-mkl-core-ps__x86_64;intel-mkl-common-c__noarch;intel-mkl-core-c__x86_64;intel-mkl-common-c-ps__noarch;intel-mkl-tbb__x86_64;intel-mkl-tbb-rt__x86_64;intel-mkl-gnu-c__x86_64;intel-mkl-common-f__noarch;intel-mkl-core-f__x86_64;intel-mkl-gnu-f-rt__x86_64;intel-mkl-gnu-f__x86_64;intel-mkl-f95-common__noarch;intel-mkl-f__x86_64;intel-mkl-psxe__noarch;intel-psxe-common__noarch;intel-psxe-common-doc__noarch;intel-compxe-pset/g' silent.cfg && \
  ./install.sh -s silent.cfg && \
  cd .. && rm -rf * && \
  rm -rf /opt/intel/.*.log /opt/intel/compilers_and_libraries_2019.3.199/licensing && \
  echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/intel.conf && \
  ldconfig && \
  echo "source /opt/intel/mkl/bin/mklvars.sh intel64" >> /etc/bash.bashrc

ENV LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2019.3.199/linux/tbb/lib/intel64_lin/gcc4.7:/opt/intel/compilers_and_libraries_2019.3.199/linux/compiler/lib/intel64_lin:/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin
ENV CPATH=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/include
ENV NLSPATH=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin/locale/%l_%t/%N
ENV LIBRARY_PATH=/opt/intel/compilers_and_libraries_2019.3.199/linux/tbb/lib/intel64_lin/gcc4.7:/opt/intel/compilers_and_libraries_2019.3.199/linux/compiler/lib/intel64_lin:/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin
ENV MKLROOT=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV PKG_CONFIG_PATH=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/bin/pkgconfig

RUN wget -q https://dl.bintray.com/boostorg/release/1.71.0/source/boost_1_71_0.tar.gz && \
  tar -xzf boost_1_71_0.tar.gz && \
  cd boost_1_71_0 && \
  ./bootstrap.sh && \
  ./b2 install
  
COPY src /GEM/src/

RUN apt-get -y install make && \
  cd /GEM/src/ && \
  env && \
  pwd && \
  ls -l && \
  make && \
  mv /GEM/src/GEM /GEM/GEM

