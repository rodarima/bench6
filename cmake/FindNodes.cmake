include(GNUInstallDirs)

find_library(NODES_LIBRARY NAMES nodes)
find_path(NODES_INCLUDE_DIR nodes.h)
find_file(NODES_WRAPPER NAMES nodes-main-wrapper.o PATH_SUFFIXES "lib")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Nodes DEFAULT_MSG
  NODES_LIBRARY NODES_INCLUDE_DIR NODES_WRAPPER)

if(NOT NODES_FOUND)
  return()
endif()

if(NOT TARGET Nodes::nodes)
  add_library(Nodes::nodes SHARED IMPORTED)
  set_target_properties(Nodes::nodes PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${NODES_INCLUDE_DIR}"
    IMPORTED_LOCATION ${NODES_LIBRARY})
endif()

if(NOT TARGET Nodes::wrapper)
  add_library(Nodes::wrapper STATIC IMPORTED)
  set_target_properties(Nodes::wrapper PROPERTIES
    IMPORTED_LOCATION ${NODES_WRAPPER})
  target_compile_options(Nodes::wrapper INTERFACE "-fompss-2")
  target_link_libraries(Nodes::wrapper INTERFACE Nodes::nodes)
endif()
