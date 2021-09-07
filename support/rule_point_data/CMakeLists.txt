target_sources(support 
    PUBLIC 
    sphere_lebedev_rule.hpp 
    PRIVATE
    sphere_lebedev_rule.cpp
)

target_include_directories(support 
    PUBLIC .
)