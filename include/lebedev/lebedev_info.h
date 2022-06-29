#pragma once

#include "detail/point_type.h"
#include "detail/rule_info.h"
#include "detail/rule_point_data.h"
#include "detail/rule_point_expand.h"
#include "detail/rule_point_summary.h"

namespace lebedev
{
    using detail::c_array;
    using detail::LEBEDEV_POINT_TYPE;

    //Note that Rule starts from 1
    template <size_t Rule>
    struct LebedevInfo
    {
    public:
        static constexpr bool if_available_ = detail::if_available(Rule);
        static constexpr size_t precision_ = detail::precision(Rule);
        static constexpr size_t available_index_ = detail::available_index(Rule);
        static constexpr size_t total_num_of_points_ = detail::point_count(Rule);
        static constexpr c_array<size_t, 6> point_type_count_ = detail::point_type_count(Rule);
        static constexpr c_array<size_t, 6> point_type_accum_count_ = detail::point_type_accum_count(Rule);
        static constexpr c_array<size_t, 6> point_type_multiplicity_ = detail::point_type_multiplicity();

        static constexpr size_t point_type_total_ = point_type_accum_count_[5];

        using PointData = detail::RulePointData<Rule>;
        using PointDataFull = detail::RulePointDataFull<Rule>;

        template <LEBEDEV_POINT_TYPE Type, size_t index,
                  typename = typename std::enable_if_t<std::greater<size_t>()(point_type_count_[static_cast<size_t>(Type)], index)>>
        struct UniquePoint
        {
            static constexpr size_t type_id_ = static_cast<size_t>(Type);
            static constexpr size_t point_id_ = index + (type_id_ > 0 ? point_type_accum_count_[type_id_ - 1] : 0);
            static constexpr c_array<double, 3> pos_ = PointData::pos_list_[point_id_];
            static constexpr double weight_ = PointData::weight_list_[point_id_];
        };

    private:
        // obtain the type of corresponding unique Lebdev point, namely equivalent points under octahedra group is considered as one
        // index should be less than point_type_total_
        static constexpr LEBEDEV_POINT_TYPE point_type(size_t index)
        {
            if (index < point_type_accum_count_[0])
                return LEBEDEV_POINT_TYPE::OPTRN6_00C;
            if (index < point_type_accum_count_[1])
                return LEBEDEV_POINT_TYPE::OPTRN12_0BB;
            if (index < point_type_accum_count_[2])
                return LEBEDEV_POINT_TYPE::OPTRN8_AAA;
            if (index < point_type_accum_count_[3])
                return LEBEDEV_POINT_TYPE::OPTRN24_AAC;
            if (index < point_type_accum_count_[4])
                return LEBEDEV_POINT_TYPE::OPTRN24_AB0;
            if (index < point_type_accum_count_[5])
                return LEBEDEV_POINT_TYPE::OPTRN48_ABC;
            return LEBEDEV_POINT_TYPE::OPTRN0_EMPTY;
        }

        struct EmptyPointInfo
        {
            static constexpr LEBEDEV_POINT_TYPE point_type_ = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY;
            static constexpr size_t point_multiplicity_ = 0;
            static constexpr size_t point_shift_ = 0;
        };

        // obtain the infomation of corresponding unique Lebdev point, namely equivalent points under octahedra group is considered as one
        template <size_t Index>
        struct PointInfo
        {
            static_assert(Index < point_type_total_, "Index must be smaller than point_type_total_");
            using LastPointInfo = typename qutility::traits::static_if<Index == 0, EmptyPointInfo, PointInfo<Index - 1>>::type;
            static constexpr LEBEDEV_POINT_TYPE point_type_ = point_type(Index);
            static constexpr size_t point_multiplicity_ = point_type_multiplicity_[static_cast<size_t>(point_type_)];
            static constexpr size_t point_shift_ = LastPointInfo::point_shift_ + LastPointInfo::point_multiplicity_;
        };

        // Here the dual {} is not used since this will cause error when N==0
        template <size_t N, size_t... Is>
        static constexpr c_array<LEBEDEV_POINT_TYPE, N> point_type_list_impl(qutility::c_array::seq<Is...>)
        {
            return {PointInfo<Is>::point_type_...};
        };

        template <size_t N, size_t... Is>
        static constexpr c_array<size_t, N> point_multiplicity_list_impl(qutility::c_array::seq<Is...>)
        {
            return {PointInfo<Is>::point_multiplicity_...};
        };

        template <size_t N, size_t... Is>
        static constexpr c_array<size_t, N> point_shift_list_impl(qutility::c_array::seq<Is...>)
        {
            return {PointInfo<Is>::point_shift_...};
        };

    public:
        static constexpr auto point_type_list_ = point_type_list_impl<point_type_total_>(qutility::c_array::gen_seq<point_type_total_>{});
        static constexpr auto point_multiplicity_list_ = point_multiplicity_list_impl<point_type_total_>(qutility::c_array::gen_seq<point_type_total_>{});
        static constexpr auto point_shift_list_ = point_shift_list_impl<point_type_total_>(qutility::c_array::gen_seq<point_type_total_>{});
    };
}