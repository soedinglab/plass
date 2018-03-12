include(MMseqsResourceCompiler)

set(COMPILED_RESOURCES
        assembler.sh
        hybridassembler.sh
        nuclassembler.sh
        predict_coding_acc9260_56x96.model
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(local-generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})
