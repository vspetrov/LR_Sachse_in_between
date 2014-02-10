import os
env = Environment(ENV = os.environ)
env['CC']='icl'

env.Append(CPPPATH = ['/workspace/opencv/build/include'])
env.Append(LIBPATH = ['/workspace/opencv/build/x64/vc11/lib'])
env.Append(LIBS = ['opencv_core248','opencv_highgui248','opencv_imgproc248'])
sources = [
'LR_cell.cpp',
'LR_lattice.cpp',
'Sachse_fibrolast.cpp',
'main.cpp'
]

headers = [
'LR_cell.h',
'LR_lattice.h',
'Sachse_fibroblast.h'
]

env.Append(CCFLAGS="-O3")
env.Program(target="lr", source=sources)