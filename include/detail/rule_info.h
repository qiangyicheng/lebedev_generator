#pragma once

#include <type_traits>

#include "qutility/c_array.h"
#include "qutility/traits.h"

#include "point_type.h"

namespace lebedev {
	namespace detail {
		using qutility::c_array::c_array;

		struct LebedevRuleInfo {
			static constexpr size_t rule_max_ = 65;
			static constexpr size_t rule_avaliable_ = 32;
			static constexpr c_array<size_t, 6> point_type_multiplicity_ = { {
				6,   12,    8,   24,   24,   48
			} };
			static constexpr c_array<size_t, rule_max_> available_table_ = { {
				1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
				1,    1,    1,    1,    1,    0,    1,    0,    0,    1,
				0,    0,    1,    0,    0,    1,    0,    0,    1,    0,
				0,    1,    0,    0,    1,    0,    0,    1,    0,    0,
				1,    0,    0,    1,    0,    0,    1,    0,    0,    1,
				0,    0,    1,    0,    0,    1,    0,    0,    1,    0,
				0,    1,    0,    0,    1
			} };
			static constexpr c_array<size_t, rule_max_> available_index_table_ = { {
				1,    2,    3,    4,    5,    6,    7,    8,    9,   10,
			   11,   12,   13,   14,   15,    0,   16,    0,    0,   17,
				0,    0,   18,    0,    0,   19,    0,    0,   20,    0,
				0,   21,    0,    0,   22,    0,    0,   23,    0,    0,
			   24,    0,    0,   25,    0,    0,   26,    0,    0,   27,
				0,    0,   28,    0,    0,   29,    0,    0,   30,    0,
				0,   31,    0,    0,   32
			} };
			static constexpr c_array<size_t, rule_max_> point_count_table_ = { {
				6,   14,   26,   38,   50,   74,   86,  110,  146,  170,
			  194,  230,  266,  302,  350,  386,  434,  482,  530,  590,
			  650,  698,  770,  830,  890,  974, 1046, 1118, 1202, 1274,
			 1358, 1454, 1538, 1622, 1730, 1814, 1910, 2030, 2126, 2222,
			 2354, 2450, 2558, 2702, 2810, 2930, 3074, 3182, 3314, 3470,
			 3590, 3722, 3890, 4010, 4154, 4334, 4466, 4610, 4802, 4934,
			 5090, 5294, 5438, 5606, 5810
			} };
			static constexpr c_array<size_t, rule_max_> precision_table_ = { {
				3,    5,    7,    9,   11,   13,   15,   17,   19,   21,
			   23,   25,   27,   29,   31,   33,   35,   37,   39,   41,
			   43,   45,   47,   49,   51,   53,   55,   57,   59,   61,
			   63,   65,   67,   69,   71,   73,   75,   77,   79,   81,
			   83,   85,   87,   89,   91,   93,   95,   97,   99,  101,
			  103,  105,  107,  109,  111,  113,  115,  117,  119,  121,
			  123,  125,  127,  129,  131
			} };
			static constexpr c_array<c_array<size_t, 6>, rule_avaliable_> point_type_count_table_ = { {
			   {{1,    0,    0,    0,    0,    0}},
			   {{1,    0,    1,    0,    0,    0}},
			   {{1,    1,    1,    0,    0,    0}},
			   {{1,    0,    1,    0,    1,    0}},
			   {{1,    1,    1,    1,    0,    0}},
			   {{1,    1,    1,    1,    1,    0}},
			   {{1,    0,    1,    2,    1,    0}},
			   {{1,    0,    1,    3,    1,    0}},
			   {{1,    1,    1,    3,    0,    1}},
			   {{1,    1,    1,    3,    1,    1}},
			   {{1,    1,    1,    4,    1,    1}},
			   {{1,    0,    1,    5,    2,    1}},
			   {{1,    1,    1,    5,    1,    2}},
			   {{1,    0,    1,    6,    2,    2}},
			   {{1,    0,    1,    6,    2,    3}},
			   {{1,    1,    1,    7,    2,    4}},
			   {{1,    0,    1,    9,    3,    6}},
			   {{1,    1,    1,   10,    3,    9}},
			   {{1,    0,    1,   12,    4,   12}},
			   {{1,    1,    1,   13,    4,   16}},
			   {{1,    0,    1,   15,    5,   20}},
			   {{1,    1,    1,   16,    5,   25}},
			   {{1,    0,    1,   18,    6,   30}},
			   {{1,    1,    1,   19,    6,   36}},
			   {{1,    0,    1,   21,    7,   42}},
			   {{1,    1,    1,   22,    7,   49}},
			   {{1,    0,    1,   24,    8,   56}},
			   {{1,    1,    1,   25,    8,   64}},
			   {{1,    0,    1,   27,    9,   72}},
			   {{1,    1,    1,   28,    9,   81}},
			   {{1,    0,    1,   30,   10,   90}},
			   {{1,    1,    1,   31,   10,  100}}
			} };
			static constexpr c_array<c_array<size_t, 6>, rule_avaliable_> point_type_accum_count_table_ = { {
			   {{1,    1,    1,    1,    1,    1}},
			   {{1,    1,    2,    2,    2,    2}},
			   {{1,    2,    3,    3,    3,    3}},
			   {{1,    1,    2,    2,    3,    3}},
			   {{1,    2,    3,    4,    4,    4}},
			   {{1,    2,    3,    4,    5,    5}},
			   {{1,    1,    2,    4,    5,    5}},
			   {{1,    1,    2,    5,    6,    6}},
			   {{1,    2,    3,    6,    6,    7}},
			   {{1,    2,    3,    6,    7,    8}},
			   {{1,    2,    3,    7,    8,    9}},
			   {{1,    1,    2,    7,    9,   10}},
			   {{1,    2,    3,    8,    9,   11}},
			   {{1,    1,    2,    8,   10,   12}},
			   {{1,    1,    2,    8,   10,   13}},
			   {{1,    2,    3,   10,   12,   16}},
			   {{1,    1,    2,   11,   14,   20}},
			   {{1,    2,    3,   13,   16,   25}},
			   {{1,    1,    2,   14,   18,   30}},
			   {{1,    2,    3,   16,   20,   36}},
			   {{1,    1,    2,   17,   22,   42}},
			   {{1,    2,    3,   19,   24,   49}},
			   {{1,    1,    2,   20,   26,   56}},
			   {{1,    2,    3,   22,   28,   64}},
			   {{1,    1,    2,   23,   30,   72}},
			   {{1,    2,    3,   25,   32,   81}},
			   {{1,    1,    2,   26,   34,   90}},
			   {{1,    2,    3,   28,   36,  100}},
			   {{1,    1,    2,   29,   38,  110}},
			   {{1,    2,    3,   31,   40,  121}},
			   {{1,    1,    2,   32,   42,  132}},
			   {{1,    2,    3,   34,   44,  144}}
			} };
		};

		inline constexpr size_t point_type_multiplicity(LEBEDEV_POINT_TYPE PType) {
			auto index=static_cast<size_t>(PType);
			if ( index<6 ) return LebedevRuleInfo::point_type_multiplicity_[static_cast<size_t>(PType)];
			else return 0;
		}

        inline constexpr bool if_available(size_t rule) {
			if (rule < 1) return false;
			if (rule > LebedevRuleInfo::rule_max_) return false;
			return LebedevRuleInfo::available_table_[rule - 1];
		}

		inline constexpr size_t precision(size_t rule) {
			return LebedevRuleInfo::precision_table_[rule - 1];
		}

		inline constexpr size_t available_index(size_t rule) {
			return LebedevRuleInfo::available_index_table_[rule - 1];
		}

		inline constexpr size_t point_count(size_t rule) {
			return LebedevRuleInfo::point_count_table_[rule - 1];
		}

        inline constexpr c_array<size_t, 6> point_type_count(size_t rule) {
			return !if_available(rule) ? c_array<size_t, 6>{} : LebedevRuleInfo::point_type_count_table_[available_index(rule) - 1];
		}

        inline constexpr c_array<size_t, 6> point_type_accum_count(size_t rule) {
			return !if_available(rule) ? c_array<size_t, 6>{} : LebedevRuleInfo::point_type_accum_count_table_[available_index(rule) - 1];
		}

        inline constexpr c_array<size_t, 6> point_type_multiplicity() {
			return LebedevRuleInfo::point_type_multiplicity_;
		}
	}
}
