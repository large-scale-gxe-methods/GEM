FROM robbyjo/ubuntu-mkl:18.04-2019.1

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install make zlib1g-dev

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

RUN cd /GEM/src/ && \
  env && \
  pwd && \
  ls -l && \
  make && \
  mv /GEM/src/GEM /GEM/GEM

