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

TEST(LebedevInfoTest, ValidRules) {
  using lebedev::lebedev::LebedevInfo;
  using LebedevRule15=LebedevInfo<15>;

  EXPECT_EQ(LebedevRule15::precision_, 31)<<"Incorrect Precision";
  EXPECT_EQ(LebedevRule15::if_available_, true);

}