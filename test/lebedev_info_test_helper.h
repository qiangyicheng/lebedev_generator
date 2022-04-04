#pragma once

#include <type_traits>

#include "lebedev_info.h"

template <lebedev::detail::LEBEDEV_POINT_TYPE Type>
auto get_reference_genoh_results(const lebedev::detail::c_array<double, 3> p)
{
    using namespace lebedev::detail;
    constexpr size_t size = point_type_multiplicity(Type);
    constexpr size_t code = static_cast<size_t>(Type) + 1;

    c_array<c_array<double, 3>, size> ans;

    double a, b;
    if constexpr (Type == LEBEDEV_POINT_TYPE::OPTRN6_00C)
    {
        a = p[2];
        b = 0;
    }
    else if constexpr (Type == LEBEDEV_POINT_TYPE::OPTRN12_0BB)
    {
        a = p[1];
        b = 0;
    }
    else if constexpr (Type == LEBEDEV_POINT_TYPE::OPTRN8_AAA)
    {
        a = p[0];
        b = 0;
    }
    else if constexpr (Type == LEBEDEV_POINT_TYPE::OPTRN24_AAC)
    {
        a = p[0];
        b = p[2];
    }
    else if constexpr (Type == LEBEDEV_POINT_TYPE::OPTRN24_AB0)
    {
        a = p[0];
        b = p[1];
    }
    else if constexpr (Type == LEBEDEV_POINT_TYPE::OPTRN48_ABC)
    {
        a = p[0];
        b = p[1];
    }

    double x[size];
    double y[size];
    double z[size];
    double w[size];

    gen_oh(code, a, b, 0, x, y, z, w);

    for (int itr = 0; itr < size; ++itr)
    {
        ans[itr][0] = x[itr];
        ans[itr][1] = y[itr];
        ans[itr][2] = z[itr];
    }
    return ans;
}
