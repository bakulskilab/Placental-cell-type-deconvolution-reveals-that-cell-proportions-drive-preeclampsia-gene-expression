# https://github.com/statgen/popscle/issues/21 Had to manually edit a file to avoid make command failing at 6%
# http://www.htslib.org/download/ for htslib installation instruction, --prefix=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies/ used to install locally
# https://github.com/samtools/htslib/releases/ Had to install version 1.10.2 because >1.11 breaks popscle compilation

cmake -DHTS_INCLUDE_DIRS=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies/include/ \
  -DHTS_LIBRARIES=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies/lib/libhts.a ..