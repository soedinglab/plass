add_library(kerasify OBJECT
        keras_model.cpp
        keras_model.h
        )
set_target_properties(kerasify PROPERTIES COMPILE_FLAGS ${MMSEQS_CXX_FLAGS} LINK_FLAGS ${MMSEQS_CXX_FLAGS})
