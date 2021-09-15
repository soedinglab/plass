include(MMseqsResourceCompiler)

set(COMPILED_RESOURCES
        assemble.sh
        nuclassemble.sh
        guidedNuclAssemble.sh
        predict_coding_acc9260_56x96.model
        predict_coding_acc9623_57x32x64.model
        predict_coding_acc9540_57x32x64.model
        predict_coding_acc9642_57x32x64.model
        predict_coding_acc9598_57x32x64.model
        predict_coding_acc9743_57x32x64.model
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(local-generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})
