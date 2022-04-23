if(ENABLE_WOLFRAM_SCRIPT)
    ##rule_point_expand.h
    add_custom_command(
        OUTPUT 
            rule_point_expand.h
        COMMAND 
            wolframscript -file 
            ${PROJECT_SOURCE_DIR}/generator/rule_point_expand/rule_point_expand.wls 
            ${PROJECT_SOURCE_DIR}/generator/rule_point_expand/rule_point_expand.h.template
            ${CMAKE_CURRENT_SOURCE_DIR}/rule_point_expand.h
        DEPENDS 
            ${PROJECT_SOURCE_DIR}/generator/rule_point_expand/rule_point_expand.wls 
            ${PROJECT_SOURCE_DIR}/generator/rule_point_expand/rule_point_expand.h.template
        VERBATIM
    )
    add_custom_target(rule_point_expand
        DEPENDS 
            rule_point_expand.h
            ${WOLFRAM_SCRIPT_TARGET_LIST}
    )
    append_wolfram_script_target_list(rule_point_expand)
    ##rule_point_data.h
    add_custom_command(
        OUTPUT 
            rule_point_data.h
        COMMAND 
            wolframscript -file 
            ${PROJECT_SOURCE_DIR}/generator/rule_point_data/rule_point_data.wls
            ${PROJECT_SOURCE_DIR}/support/rule_point_data/sphere_lebedev_rule.cpp
            ${PROJECT_SOURCE_DIR}/generator/rule_point_data/rule_point_data.h.template
            ${CMAKE_CURRENT_SOURCE_DIR}/rule_point_data.h
        DEPENDS
            rule_point_expand
            ${PROJECT_SOURCE_DIR}/generator/rule_point_data/rule_point_data.wls
            ${PROJECT_SOURCE_DIR}/support/rule_point_data/sphere_lebedev_rule.cpp
            ${PROJECT_SOURCE_DIR}/generator/rule_point_data/rule_point_data.h.template
        VERBATIM
    )
    add_custom_target(rule_point_data
        DEPENDS 
            rule_point_data.h
            ${WOLFRAM_SCRIPT_TARGET_LIST}
    )
    append_wolfram_script_target_list(rule_point_data)

    ##
    add_dependencies(${PROJECT_NAME} rule_point_expand rule_point_data)
endif()

target_sources(${PROJECT_NAME}
    INTERFACE
    point_type.h
    rule_info.h
    rule_point_data.h
    rule_point_expand.h
    rule_point_summary.h
)