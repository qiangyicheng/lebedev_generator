// This file defines the types of different Lebedev points.
// The number OPTRN denotes the number of unique equivalent points under the octahedra group

#pragma once

namespace lebedev
{
    namespace detail
    {
        enum class LEBEDEV_POINT_TYPE
        {
            OPTRN6_00C = 0,
            OPTRN12_0BB = 1,
            OPTRN8_AAA = 2,
            OPTRN24_AAC = 3,
            OPTRN24_AB0 = 4,
            OPTRN48_ABC = 5,
            OPTRN0_EMPTY = 6
        };
    }
}