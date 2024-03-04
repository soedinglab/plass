include_directories(commons)
add_subdirectory(assembler)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)
add_subdirectory(util)

add_executable(plass
        ${commons_source_files}
        ${plass_assembler_source_files}
        ${plass_workflow_source_files}
        ${util_source_files}
        LocalCommandDeclarations.h
        plass.cpp
        )
mmseqs_setup_derived_target(plass)
target_link_libraries(plass version kerasify flash)

add_dependencies(plass local-generated)

install(TARGETS plass DESTINATION bin)


add_executable(penguin
        ${commons_source_files}
        ${penguin_assembler_source_files}
        ${penguin_workflow_source_files}
        ${util_source_files}
        LocalCommandDeclarations.h
        penguin.cpp
        )
mmseqs_setup_derived_target(penguin)
target_link_libraries(penguin version kerasify flash)

add_dependencies(penguin local-generated)

install(TARGETS penguin DESTINATION bin)