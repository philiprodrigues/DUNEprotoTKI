add_executable(anaData anaData.cxx)
target_link_libraries(anaData ROOT::Core style)

install(TARGETS anaData DESTINATION "${CMAKE_INSTALL_BINDIR}")
