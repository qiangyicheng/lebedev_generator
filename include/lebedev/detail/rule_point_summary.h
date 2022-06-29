#pragma once

#include "rule_point_expand.h"

namespace lebedev
{
    namespace detail
    {
        using qutility::c_array::c_array;

        template <size_t Rule>
        struct RulePointDataFull
        {
            using RulePointDataT = RulePointData<Rule>;

        private:
            using SeqT = qutility::c_array::gen_seq<RulePointDataT::point_type_total_>;
            template <size_t... Is>
            constexpr static auto type_list_impl(qutility::c_array::seq<Is...>)
            {
                return qutility::c_array::join((type_expand<RulePointDataT::type_list_[Is]>())...);
            }
            template <size_t... Is>
            constexpr static auto weight_list_impl(qutility::c_array::seq<Is...>)
            {
                return qutility::c_array::join((weight_expand<RulePointDataT::type_list_[Is]>(RulePointDataT::weight_list_[Is]))...);
            }
            template <size_t... Is>
            constexpr static auto pos_list_impl(qutility::c_array::seq<Is...>)
            {
                return qutility::c_array::join((point_expand<RulePointDataT::type_list_[Is]>(RulePointDataT::pos_list_[Is]))...);
            }

        public:
            constexpr static auto type_list_ = type_list_impl(SeqT());
            constexpr static auto weight_list_ = weight_list_impl(SeqT());
            constexpr static auto pos_list_ = pos_list_impl(SeqT());
        };

    }
}