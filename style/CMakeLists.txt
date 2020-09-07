root_generate_dictionary(styleDict style.h LINKDEF styleLinkDef.h)

add_library(style SHARED styleDict.cxx style.cxx)
target_include_directories(style PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(style PUBLIC ROOT::Hist ROOT::Tree ROOT::Physics ROOT::Graf ROOT::Gui)

install(TARGETS style DESTINATION "${CMAKE_INSTALL_LIBDIR}")
