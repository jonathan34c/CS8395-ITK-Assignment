cmake_minimum_required(VERSION 3.16)
project (RegistrationITK)

find_package (ITK REQUIRED)
include (${ITK_USE_FILE})

add_executable (RegistrationITK RegistrationITK.cxx)

target_link_libraries (RegistrationITK ${ITK_LIBRARIES})
