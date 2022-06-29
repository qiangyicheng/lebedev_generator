#pragma once

#include <type_traits>

#include "qutility/c_array.h"

namespace lebedev
{
    namespace detail
    {
        using qutility::c_array::c_array;
        using qutility::c_array::inner_product;

        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto oau_size = inner_product(LebedevInfo::point_type_count_, GroupInfo::fields_required_table_);

        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto get_oau_positions()
        {
            c_array<c_array<double, 3>, oau_size<LebedevInfo, GroupInfo>> ans{};

            size_t count = 0;
            for (size_t itr_l = 0; itr_l < LebedevInfo::point_type_total_; ++itr_l)
            {
                auto point_type = LebedevInfo::point_type_list_[itr_l];
                auto point_shift = LebedevInfo::point_shift_list_[itr_l];
                auto point_type_index = static_cast<size_t>(point_type);
                for (size_t itr_g = 0; itr_g < GroupInfo::fields_required_table_[point_type_index]; ++itr_g)
                {
                    auto field_index = GroupInfo::field_indexes_table_[point_type_index][itr_g];
                    ans[count][0] = LebedevInfo::PointDataFull::pos_list_[point_shift + field_index][0];
                    ans[count][1] = LebedevInfo::PointDataFull::pos_list_[point_shift + field_index][1];
                    ans[count][2] = LebedevInfo::PointDataFull::pos_list_[point_shift + field_index][2];
                    ++count;
                }
            }

            return ans;
        }
        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto oau_positions = get_oau_positions<LebedevInfo, GroupInfo>();

        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto get_oau_types()
        {
            using EleT = std::remove_cvref_t<decltype(LebedevInfo::point_type_list_[0])>;
            c_array<EleT, oau_size<LebedevInfo, GroupInfo>> ans{};

            size_t count = 0;
            for (size_t itr_l = 0; itr_l < LebedevInfo::point_type_total_; ++itr_l)
            {
                auto point_type = LebedevInfo::point_type_list_[itr_l];
                auto point_type_index = static_cast<size_t>(point_type);
                for (size_t itr_g = 0; itr_g < GroupInfo::fields_required_table_[point_type_index]; ++itr_g)
                {
                    ans[count] = point_type;
                    ++count;
                }
            }

            return ans;
        }
        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto oau_types = get_oau_types<LebedevInfo, GroupInfo>();

        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto get_oau_weights()
        {
            c_array<double, oau_size<LebedevInfo, GroupInfo>> ans{};

            size_t count = 0;
            for (size_t itr_l = 0; itr_l < LebedevInfo::point_type_total_; ++itr_l)
            {
                auto point_type = LebedevInfo::point_type_list_[itr_l];
                auto point_type_index = static_cast<size_t>(point_type);
                for (size_t itr_g = 0; itr_g < GroupInfo::fields_required_table_[point_type_index]; ++itr_g)
                {
                    ans[count] = LebedevInfo::PointData::weight_list_[itr_l];
                    ++count;
                }
            }

            return ans;
        }
        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto oau_weights = get_oau_weights<LebedevInfo, GroupInfo>();

        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto get_oau_submultiplicities()
        {
            using EleT = std::remove_cvref_t<decltype(GroupInfo::field_submultiplicities_table_[0][0])>;
            c_array<EleT, oau_size<LebedevInfo, GroupInfo>> ans{};

            size_t count = 0;
            for (size_t itr_l = 0; itr_l < LebedevInfo::point_type_total_; ++itr_l)
            {
                auto point_type = LebedevInfo::point_type_list_[itr_l];
                auto point_type_index = static_cast<size_t>(point_type);
                for (size_t itr_g = 0; itr_g < GroupInfo::fields_required_table_[point_type_index]; ++itr_g)
                {
                    ans[count] = GroupInfo::field_submultiplicities_table_[point_type_index][itr_g];
                    ++count;
                }
            }

            return ans;
        }
        template <typename LebedevInfo, typename GroupInfo>
        constexpr auto oau_submultiplicities = get_oau_submultiplicities<LebedevInfo, GroupInfo>();
    }

    using detail::oau_positions;
    using detail::oau_types;
    using detail::oau_weights;
    using detail::oau_submultiplicities;

}