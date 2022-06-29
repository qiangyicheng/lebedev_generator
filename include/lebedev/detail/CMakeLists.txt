if(ENABLE_WOLFRAM_SCRIPT)
    ##rule_point_expand.h
    set(CURRENT_TARGET_NAME rule_point_expand)
    set(CURRENT_GENERATOR_DIR ${PROJECT_SOURCE_DIR}/generator/${CURRENT_TARGET_NAME})
    add_custom_command(
        OUTPUT 
            ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TARGET_NAME}.h
        COMMAND 
            wolframscript -file 
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.wls 
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.h.template
            ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TARGET_NAME}.h
        DEPENDS 
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.wls 
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.h.template
        VERBATIM
    )
    add_custom_target(${CURRENT_TARGET_NAME}
        DEPENDS 
            ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TARGET_NAME}.h
            ${WOLFRAM_SCRIPT_TARGET_LIST}
    )
    append_wolfram_script_target_list(${CURRENT_TARGET_NAME})
    
    
    ##rule_point_data.h
    set(CURRENT_TARGET_NAME rule_point_data)
    set(CURRENT_GENERATOR_DIR ${PROJECT_SOURCE_DIR}/generator/${CURRENT_TARGET_NAME})
    set(CURRENT_SUPPORT_DIR ${PROJECT_SOURCE_DIR}/support/${CURRENT_TARGET_NAME})
    add_custom_command(
        OUTPUT 
            ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TARGET_NAME}.h
        COMMAND 
            wolframscript -file 
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.wls
            ${CURRENT_SUPPORT_DIR}/sphere_lebedev_rule.cpp
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.h.template
            ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TARGET_NAME}.h
        DEPENDS
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.wls
            ${CURRENT_SUPPORT_DIR}/sphere_lebedev_rule.cpp
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.h.template
        VERBATIM
    )
    add_custom_target(${CURRENT_TARGET_NAME}
        DEPENDS 
            ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_TARGET_NAME}.h
            ${WOLFRAM_SCRIPT_TARGET_LIST}
    )
    append_wolfram_script_target_list(${CURRENT_TARGET_NAME})

    ##
    add_dependencies(${PROJECT_NAME} rule_point_expand rule_point_data)
endif()

target_sources(${PROJECT_NAME}
    INTERFACE
    point_type.h
    rule_info.h
    #${CMAKE_CURRENT_SOURCE_DIR}/rule_point_data.h ## not included, since the dependency is set by add_dependencies before (which has a private scope)
    #${CMAKE_CURRENT_SOURCE_DIR}/rule_point_expand.h ## not included, since the dependency is set by add_dependencies before (which has a private scope)
    rule_point_summary.h
)