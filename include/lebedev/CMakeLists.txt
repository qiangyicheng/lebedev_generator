add_subdirectory(detail)
add_subdirectory(group_info)

target_sources(${PROJECT_NAME}
    INTERFACE
    lebedev_info.h
)