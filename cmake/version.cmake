# Simple target for printing Project version
add_custom_target(print_version
                  COMMAND ${CMAKE_COMMAND} -E echo "PROJECT_NAME    = ${PROJECT_NAME}"
                  COMMAND ${CMAKE_COMMAND} -E echo "PROJECT_DATE    = ${${PROJECT_NAME}_DATE}"
                  COMMAND ${CMAKE_COMMAND} -E echo "PROJECT_VERSION = ${PROJECT_VERSION}")