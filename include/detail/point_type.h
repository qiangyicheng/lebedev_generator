//This file defines the types of different Lebedev points.
//The number OPTRN denotes the number of unique equivalent points under the octahedra group 

#pragma once

namespace lebedev{
    namespace detail{
        enum class LEBEDEV_POINT_TYPE{
            OPTRN6_001 = 0,
            OPTRN12_0AA = 1,
            OPTRN8_AAA = 2,
            OPTRN24_AAB = 3,
            OPTRN24_AB0 = 4,
            OPTRN48_ABC = 5,
            OPTRN0_EMPTY = 6
        };
    }
}