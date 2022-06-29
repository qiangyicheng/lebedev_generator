set(GROUP_LIST 
    2D009
    3D225 3D229 3D230
)

string(REPLACE ";" " " GROUP_LIST_SPACE_SEPRATED_STRING "${GROUP_LIST}")

macro(add_prefix_surfix prefix surfix list_name)
   set(${list_name}_TMP)
   foreach(l IN LISTS ${list_name})
      list(APPEND ${list_name}_TMP ${prefix}${l}${surfix} )
   endforeach()

   set(${list_name} ${${list_name}_TMP})
   unset(${list_name}_TMP)
endmacro(add_prefix_surfix)

set(GROUP_INFO_FILENAMES_H ${GROUP_LIST})
set(GROUP_INFO_FILENAMES_CUH ${GROUP_LIST})
set(GROUP_INFO_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/group)
set(GROUP_INFO_SURFIX .h)
add_prefix_surfix(${GROUP_INFO_PREFIX} ${GROUP_INFO_SURFIX} GROUP_INFO_FILENAMES_H)
set(GROUP_INFO_SURFIX .cuh)
add_prefix_surfix(${GROUP_INFO_PREFIX} ${GROUP_INFO_SURFIX} GROUP_INFO_FILENAMES_CUH)

if(ENABLE_WOLFRAM_SCRIPT)
    ##group*.cuh group*.h
    set(CURRENT_TARGET_NAME group_info)
    set(CURRENT_GENERATOR_DIR ${PROJECT_SOURCE_DIR}/generator/${CURRENT_TARGET_NAME})
    add_custom_command(
        OUTPUT 
            ${GROUP_INFO_FILENAMES_CUH}
            ${GROUP_INFO_FILENAMES_H}
        COMMAND 
            wolframscript -file 
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.wls 
            ${CURRENT_GENERATOR_DIR}/groupXXX${GROUP_INFO_SURFIX}.template
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${GROUP_LIST_SPACE_SEPRATED_STRING}
        DEPENDS 
            ${CURRENT_GENERATOR_DIR}/GroupInfoGenerator.wl
            ${CURRENT_GENERATOR_DIR}/${CURRENT_TARGET_NAME}.wls 
            ${CURRENT_GENERATOR_DIR}/groupXXX${GROUP_INFO_SURFIX}.template
    VERBATIM
    )
    add_custom_target(${CURRENT_TARGET_NAME}
        DEPENDS 
            ${GROUP_INFO_FILENAMES_CUH}
            ${GROUP_INFO_FILENAMES_H}
            ${WOLFRAM_SCRIPT_TARGET_LIST}
    )
    append_wolfram_script_target_list(${CURRENT_TARGET_NAME})

    ##
    add_dependencies(${PROJECT_NAME} group_info)
else()
    ##when wolframscropt engine is not found, headers must be allready generated
    set(CURRENT_TARGET_NAME group_info)
    add_library(${CURRENT_TARGET_NAME} INTERFACE
        ${GROUP_INFO_FILENAMES_CUH}
        ${GROUP_INFO_FILENAMES_H}
    )
    add_dependencies(${PROJECT_NAME} group_info)
endif()

target_sources(${PROJECT_NAME}
    INTERFACE
    device_helper.h
    #${GROUP_INFO_FILENAMES_CUH} ## not included, since the dependency is set by add_dependencies before (which has a private scope)
    #${GROUP_INFO_FILENAMES_H} ## not included, since the dependency is set by add_dependencies before (which has a private scope)
)