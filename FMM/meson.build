project('FMM', 'c')

m_dep = meson.get_compiler('c').find_library('m', required : false)

fmm_lib = static_library(
  'fmm',
  ['test_artUpdate.c',
   'SoSFunction.c', 'eik_grid.c', 'files_methods.c',
   'linAlg.c', 'opti_method.c', 'priority_queue.c', 'coord.c',
   'faces.c', 'facets.c', 'neighbors.c', 'triMesh_2D.c', 'fmm_2d.c', 'path.c'])

executable('test_artUpdate', fmm_lib, dependencies : m_dep)

executable('test_artUpdateOriginal', fmm_lib, dependencies : m_dep)
