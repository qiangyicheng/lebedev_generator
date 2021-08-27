#include <gtest/gtest.h>

#include "lebedev_info.h"

TEST(CArrayTest, BasicCArrayCreation) {
  using qutility::c_array::c_array;
  constexpr c_array<int,5> test_arr{1,2,3,4,5};
  EXPECT_EQ(test_arr, test_arr);

  constexpr c_array<int,5> test_arr_2{2,3,3,3,3};
  EXPECT_NE(test_arr, test_arr_2);

  constexpr c_array<int,7> test_arr_3{1,2,3,4,5};
  EXPECT_NE(test_arr, test_arr_3);
}

TEST(LebedevInfoTest, ValidRule15) {
  using qutility::c_array::c_array;
  using lebedev::LebedevInfo;
  using LebedevRule15=LebedevInfo<15>;

  EXPECT_EQ(LebedevRule15::precision_, 31)<<"Incorrect Precision";
  EXPECT_EQ(LebedevRule15::if_available_, true)<<"This rule should be valid";

  EXPECT_EQ(LebedevRule15::available_index_, 15);

  LebedevRule15::point_multiplicity_list_;
  LebedevRule15::point_shift_list_;

  c_array<size_t, 6> point_type_accum_count_{{1,1,2,8,10,13}};
  EXPECT_EQ(LebedevRule15::point_type_accum_count_,point_type_accum_count_);

  c_array<size_t, 6> point_type_count{{1,0,1,6,2,3}};
  EXPECT_EQ(LebedevRule15::point_type_count_,point_type_count);
  
  EXPECT_EQ(LebedevRule15::point_type_list_, LebedevRule15::PointData::type_list);

  LebedevRule15::point_type_multiplicity_;
  LebedevRule15::point_type_total_;
  EXPECT_EQ(LebedevRule15::total_num_of_points_,350);

}