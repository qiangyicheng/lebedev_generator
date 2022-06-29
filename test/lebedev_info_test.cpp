// c++ headers
#include <type_traits>
#include <cmath>
#include <iomanip>

// gtest headers
#include <gtest/gtest.h>

// current project headers
#include "lebedev/lebedev_info.h"
#include "lebedev/group_info/group2D009.h"
#include "lebedev/group_info/group3D230.h"
#include "lebedev/orientational_asymmetric_unit.h"

// support headers
#include "sphere_lebedev_rule.hpp"

// helpers
#include "lebedev_info_test_helper.h"

TEST(CArrayTest, BasicCArrayCreation)
{
    using qutility::c_array::c_array;
    constexpr c_array<int, 5> test_arr{1, 2, 3, 4, 5};
    EXPECT_EQ(test_arr, test_arr);

    constexpr c_array<int, 5> test_arr_2{2, 3, 3, 3, 3};
    EXPECT_NE(test_arr, test_arr_2);

    constexpr c_array<int, 7> test_arr_3{1, 2, 3, 4, 5};
    EXPECT_NE(test_arr, test_arr_3);
}

TEST(CArrayTest, ArrayJoining)
{
    using qutility::c_array::c_array;
    constexpr c_array<int, 2> test_arr_0{1, 2};
    constexpr c_array<int, 3> test_arr_1{3, 4, 5};
    constexpr c_array<int, 3> test_arr_2{2, 3, 4};
    constexpr c_array<int, 8> test_arr_3{1, 2, 3, 4, 5, 2, 3, 4};

    constexpr auto ans = qutility::c_array::join(test_arr_0, test_arr_1, test_arr_2);

    EXPECT_EQ(test_arr_3, ans);
}

TEST(LebedevInfoTest, ValidRule15)
{
    using lebedev::LebedevInfo;
    using qutility::c_array::c_array;
    using LebedevRule15 = LebedevInfo<15>;

    EXPECT_EQ(LebedevRule15::precision_, 31) << "Incorrect Precision";
    EXPECT_EQ(LebedevRule15::if_available_, true) << "This rule should be valid";

    EXPECT_EQ(LebedevRule15::available_index_, 15);

    LebedevRule15::point_multiplicity_list_;
    LebedevRule15::point_shift_list_;

    c_array<size_t, 6> point_type_accum_count_{{1, 1, 2, 8, 10, 13}};
    EXPECT_EQ(LebedevRule15::point_type_accum_count_, point_type_accum_count_);

    c_array<size_t, 6> point_type_count{{1, 0, 1, 6, 2, 3}};
    EXPECT_EQ(LebedevRule15::point_type_count_, point_type_count);

    EXPECT_EQ(LebedevRule15::point_type_list_, LebedevRule15::PointData::type_list_);

    LebedevRule15::point_type_multiplicity_;
    LebedevRule15::point_type_total_;
    EXPECT_EQ(LebedevRule15::total_num_of_points_, 350);
}

TEST(LebedevInfoTest, GenOH00C)
{
    using namespace lebedev::detail;

    constexpr LEBEDEV_POINT_TYPE Type = LEBEDEV_POINT_TYPE::OPTRN6_00C;
    constexpr size_t pick_rule = 10;
    constexpr size_t pick_index = 0;

    constexpr c_array<double, 3> p = lebedev::LebedevInfo<pick_rule>::UniquePoint<Type, pick_index>::pos_;

    constexpr auto points = point_expand<Type>(p);
    auto ref_points = get_reference_genoh_results<Type>(p);

    for (int i = 0; i < point_type_multiplicity(Type); ++i)
    {
        EXPECT_DOUBLE_EQ(ref_points[i][0], points[i][0]);
        EXPECT_DOUBLE_EQ(ref_points[i][1], points[i][1]);
        EXPECT_DOUBLE_EQ(ref_points[i][2], points[i][2]);
    }
    // EXPECT_EQ(ref_points, points);
}

TEST(LebedevInfoTest, GenOH0BB)
{
    using namespace lebedev::detail;

    constexpr LEBEDEV_POINT_TYPE Type = LEBEDEV_POINT_TYPE::OPTRN12_0BB;
    constexpr size_t pick_rule = 65;
    constexpr size_t pick_index = 0;

    constexpr c_array<double, 3> p = lebedev::LebedevInfo<pick_rule>::UniquePoint<Type, pick_index>::pos_;

    constexpr auto points = point_expand<Type>(p);
    auto ref_points = get_reference_genoh_results<Type>(p);

    for (int i = 0; i < point_type_multiplicity(Type); ++i)
    {
        EXPECT_DOUBLE_EQ(ref_points[i][0], points[i][0]);
        EXPECT_DOUBLE_EQ(ref_points[i][1], points[i][1]);
        EXPECT_DOUBLE_EQ(ref_points[i][2], points[i][2]);
    }
    // EXPECT_EQ(ref_points, points);
}

TEST(LebedevInfoTest, GenOHAAA)
{
    using namespace lebedev::detail;

    constexpr LEBEDEV_POINT_TYPE Type = LEBEDEV_POINT_TYPE::OPTRN8_AAA;
    constexpr size_t pick_rule = 65;
    constexpr size_t pick_index = 0;

    constexpr c_array<double, 3> p = lebedev::LebedevInfo<pick_rule>::UniquePoint<Type, pick_index>::pos_;

    constexpr auto points = point_expand<Type>(p);
    auto ref_points = get_reference_genoh_results<Type>(p);

    for (int i = 0; i < point_type_multiplicity(Type); ++i)
    {
        EXPECT_DOUBLE_EQ(ref_points[i][0], points[i][0]);
        EXPECT_DOUBLE_EQ(ref_points[i][1], points[i][1]);
        EXPECT_DOUBLE_EQ(ref_points[i][2], points[i][2]);
    }
    // EXPECT_EQ(ref_points, points);
}

TEST(LebedevInfoTest, GenOHAAC)
{
    using namespace lebedev::detail;

    constexpr LEBEDEV_POINT_TYPE Type = LEBEDEV_POINT_TYPE::OPTRN24_AAC;
    constexpr size_t pick_rule = 65;
    constexpr size_t pick_index = 0;

    constexpr c_array<double, 3> p = lebedev::LebedevInfo<pick_rule>::UniquePoint<Type, pick_index>::pos_;

    constexpr auto points = point_expand<Type>(p);
    auto ref_points = get_reference_genoh_results<Type>(p);

    for (int i = 0; i < point_type_multiplicity(Type); ++i)
    {
        EXPECT_DOUBLE_EQ(ref_points[i][0], points[i][0]);
        EXPECT_DOUBLE_EQ(ref_points[i][1], points[i][1]);
        EXPECT_DOUBLE_EQ(ref_points[i][2], points[i][2]);
    }
    // EXPECT_EQ(ref_points, points);
}

TEST(LebedevInfoTest, GenOHAB0)
{
    using namespace lebedev::detail;

    constexpr LEBEDEV_POINT_TYPE Type = LEBEDEV_POINT_TYPE::OPTRN24_AB0;
    constexpr size_t pick_rule = 65;
    constexpr size_t pick_index = 0;

    constexpr c_array<double, 3> p = lebedev::LebedevInfo<pick_rule>::UniquePoint<Type, pick_index>::pos_;

    constexpr auto points = point_expand<Type>(p);
    auto ref_points = get_reference_genoh_results<Type>(p);

    for (int i = 0; i < point_type_multiplicity(Type); ++i)
    {
        EXPECT_DOUBLE_EQ(ref_points[i][0], points[i][0]);
        EXPECT_DOUBLE_EQ(ref_points[i][1], points[i][1]);
        EXPECT_DOUBLE_EQ(ref_points[i][2], points[i][2]);
    }
    // EXPECT_EQ(ref_points, points);
}

TEST(LebedevInfoTest, GenOHABC)
{
    using namespace lebedev::detail;

    constexpr LEBEDEV_POINT_TYPE Type = LEBEDEV_POINT_TYPE::OPTRN48_ABC;
    constexpr size_t pick_rule = 65;
    constexpr size_t pick_index = 0;

    constexpr c_array<double, 3> p = lebedev::LebedevInfo<pick_rule>::UniquePoint<Type, pick_index>::pos_;

    constexpr auto points = point_expand<Type>(p);
    auto ref_points = get_reference_genoh_results<Type>(p);

    for (int i = 0; i < point_type_multiplicity(Type); ++i)
    {
        EXPECT_DOUBLE_EQ(ref_points[i][0], points[i][0]);
        EXPECT_DOUBLE_EQ(ref_points[i][1], points[i][1]);
        EXPECT_DOUBLE_EQ(ref_points[i][2], points[i][2]);
    }
    // EXPECT_EQ(ref_points, points);
}

TEST(LebedevInfoTest, PointDataFull)
{
    using namespace lebedev::detail;
    constexpr size_t rule = 23;
    using RulePointDataFullT = RulePointDataFull<rule>;

    EXPECT_EQ(RulePointDataFullT::pos_list_.size(), point_count(rule));
}

TEST(GroupInfoTest, GroupInfo2D009)
{
    using GroupInfoT = lebedev::group2D009::host::GroupInfo<1, 32, 32, int>;

    constexpr auto field_indexes_table = GroupInfoT::field_indexes_table_;
    constexpr auto field_submultiplicities_table = GroupInfoT::field_submultiplicities_table_;
}

TEST(GroupInfoTest, GroupInfo3D230)
{
    using GroupInfoT = lebedev::group3D230::host::GroupInfo<32, int>;

    constexpr auto field_indexes_table = GroupInfoT::field_indexes_table_;
    constexpr auto field_submultiplicities_table = GroupInfoT::field_submultiplicities_table_;
}

TEST(OrientationalAsymmetricUnitTest, LebedevInfo15GroupInfo3D230)
{
    using LebedevInfoT = lebedev::LebedevInfo<15>;
    using GroupInfoT = lebedev::group3D230::host::GroupInfo<32, int>;

    LebedevInfoT::point_type_list_;
    LebedevInfoT::point_shift_list_;
    LebedevInfoT::PointDataFull::pos_list_;
    GroupInfoT::fields_required_table_;
    GroupInfoT::field_indexes_table_;
    constexpr auto oau_positions = lebedev::oau_positions<LebedevInfoT, GroupInfoT>;
    constexpr auto oau_types = lebedev::oau_types<LebedevInfoT, GroupInfoT>;
    constexpr auto oau_weights = lebedev::oau_weights<LebedevInfoT, GroupInfoT>;
    constexpr auto oau_submultiplicities = lebedev::oau_submultiplicities<LebedevInfoT, GroupInfoT>;
    constexpr auto total_weight = qutility::c_array::inner_product(oau_weights, oau_submultiplicities);
    EXPECT_DOUBLE_EQ(total_weight, 1.0);
}

TEST(OrientationalAsymmetricUnitTest, LebedevInfo15GroupInfo2D009)
{
    using LebedevInfoT = lebedev::LebedevInfo<15>;
    using GroupInfoT = lebedev::group2D009::host::GroupInfo<1, 32, 64, int>;

    LebedevInfoT::point_type_list_;
    LebedevInfoT::point_shift_list_;
    LebedevInfoT::PointDataFull::pos_list_;
    GroupInfoT::fields_required_table_;
    GroupInfoT::field_indexes_table_;
    constexpr auto oau_positions = lebedev::oau_positions<LebedevInfoT, GroupInfoT>;
    constexpr auto oau_types = lebedev::oau_types<LebedevInfoT, GroupInfoT>;
    constexpr auto oau_weights = lebedev::oau_weights<LebedevInfoT, GroupInfoT>;
    constexpr auto oau_submultiplicities = lebedev::oau_submultiplicities<LebedevInfoT, GroupInfoT>;
    constexpr auto total_weight = qutility::c_array::inner_product(oau_weights, oau_submultiplicities);
    EXPECT_DOUBLE_EQ(total_weight, 1.0);
}
