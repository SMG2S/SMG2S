find_package(Doxygen REQUIRED)
find_package(Python REQUIRED)

message(${Python_EXECUTABLE})

execute_process(COMMAND ${Python_EXECUTABLE} -c "import sphinx;import breathe"  RESULT_VARIABLE sphinx_ret)

if (sphinx_ret EQUAL "0")
  set(Sphinx_FOUND 1)
else()
  set(Sphinx_FOUND 0)
endif()

if(Sphinx_FOUND)
  message("Sphinx_FOUND")
endif()

set(PYTHON_DEP_STRING "Please install Python package breathe (at least version 4.18.0) and sphinx_rtd_theme.")

execute_process(COMMAND ${Python_EXECUTABLE} -c "from packaging import version; import sys; import breathe; sys.exit(version.parse(\"4.18.0\") > version.parse(breathe.__version__))" RESULT_VARIABLE ret)
if (ret)
  message(FATAL_ERROR "Missing Documentation dependency breathe. ${PYTHON_DEP_STRING}")
endif()
execute_process(COMMAND ${Python_EXECUTABLE} -c "import sphinx_rtd_theme" RESULT_VARIABLE ret)
if (ret)
  message(FATAL_ERROR "Missing Documentation dependency sphinx_rtd_theme. ${PYTHON_DEP_STRING}")
endif()


set(DOXYGEN_INPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
message("DOXYGEN_INPUT_DIRECTORY = ${DOXYGEN_INPUT_DIRECTORY}")
set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
message("DOXYGEN_OUTPUT_DIRECTORY = ${DOXYGEN_OUTPUT_DIRECTORY}")
set(DOXYGEN_HTML_OUTPUT ${CMAKE_BINARY_DIR}/sphinx/html/doxygen)
message("DOXYGEN_HTML_OUTPUT = ${DOXYGEN_HTML_OUTPUT}")

set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIRECTORY}/xml/index.xml)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_XML YES)
set(DOXYGEN_ENABLE_PREPROCESSING YES)
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
set(DOXYGEN_INTERNAL_DOCS NO)
set(DOXYGEN_EXCLUDE_PATTERNS ${PROJECT_SOURCE_DIR}/*utils*)
set(DOXYGEN_EXTRACT_PRIVATE YES)

make_directory(${DOXYGEN_HTML_OUTPUT})

doxygen_add_docs(Doxygen
    ${DOXYGEN_INPUT_DIRECTORY}
    DOXYGEN_PROJECT_NAME SMG2S
    DOXYGEN_PROJECT_NUMBER 1.2.0
    COMMENT "Generating docs")

set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR})
set(SPHINX_BUILD ${CMAKE_BINARY_DIR}/sphinx/html)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

if (UNIX)
  set(separator ":")
else ()
  set(separator ";")
endif()

set(RST_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.rst
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py)

add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different ${RST_FILES} ${CMAKE_CURRENT_BINARY_DIR}    
                   COMMAND  ${Python_EXECUTABLE} -m sphinx -b html
                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS
                   ${RST_FILES}
                   Doxygen
                   MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
                   COMMENT "Generating HTML documentation with Sphinx into {CMAKE_BINARY_DIR}/sphinx/html/")

add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})
