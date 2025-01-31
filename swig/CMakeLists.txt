if(NOT SWIG)
  return()
endif()

find_package(Python3 REQUIRED)  # Trouver Python3

set(PYTHON_SCRIPT  ${PROJECT_SOURCE_DIR}/swig/classanalysis.py)
set(PYTHON_SCRIPT2 ${PROJECT_SOURCE_DIR}/swig/toString.py)

set(FILE_TO_READ ${PROJECT_SOURCE_DIR}/swig/swig_inc.i)

set(SWIG_PACKAGE_ROOT_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
set(SWIG_PACKAGE_DESTINATION_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>)
set(GENERATED_FILE1 ${CMAKE_CURRENT_BINARY_DIR}/newobject.i)
set(GENERATED_FILE2 ${CMAKE_CURRENT_BINARY_DIR}/toString.i)



# Définition d'une cible qui exécute le script Python une seule fois
add_custom_command(
    OUTPUT ${GENERATED_FILE1}
    COMMAND ${Python3_EXECUTABLE} ${PYTHON_SCRIPT} ${FILE_TO_READ}
            ${GENERATED_FILE1}
    DEPENDS ${SOURCE_DIR} ${PYTHON_SCRIPT}
    COMMENT "Generating ${GENERATED_FILE1} from ${SOURCE_FILE}"
    VERBATIM
)

add_custom_command(
    OUTPUT  ${GENERATED_FILE2}
    COMMAND ${Python3_EXECUTABLE} ${PYTHON_SCRIPT2} 
            ${GENERATED_FILE2}
    DEPENDS ${SOURCE_DIR} ${PYTHON_SCRIPT2}
    COMMENT "Generating ${GENERATED_FILE2} from ${SOURCE_FILE}"
    VERBATIM
)

# Création d'une cible dédiée pour la génération du fichier
add_custom_target(generate_output DEPENDS ${GENERATED_FILE1} ${GENERATED_FILE2})

set(CMAKE_SWIG_FLAGS "-I${CMAKE_BINARY_DIR}")