#pragma once



#include "../lebedev_info.h"

namespace lebedev
{
    namespace group3D230
    {
        namespace host
        {
            using lebedev::c_array;
            using lebedev::LEBEDEV_POINT_TYPE;

            namespace detail
            {
                using lebedev::detail::point_type_multiplicity;
                using qutility::c_array::maximum;
                using qutility::c_array::add_padding;

                template <typename SignedT = int>
                constexpr SignedT operation_number = 96;
                template <typename SignedT = int>
                constexpr c_array<SignedT, 6> fields_required_table = {{1, 1, 1, 1, 1, 1}};
                template <typename SignedT = int>
                constexpr c_array<SignedT, 6> extract_num_table = 
                {{
                    (SignedT)point_type_multiplicity(LEBEDEV_POINT_TYPE::OPTRN6_00C),
                    (SignedT)point_type_multiplicity(LEBEDEV_POINT_TYPE::OPTRN12_0BB),
                    (SignedT)point_type_multiplicity(LEBEDEV_POINT_TYPE::OPTRN8_AAA),
                    (SignedT)point_type_multiplicity(LEBEDEV_POINT_TYPE::OPTRN24_AAC),
                    (SignedT)point_type_multiplicity(LEBEDEV_POINT_TYPE::OPTRN24_AB0),
                    (SignedT)point_type_multiplicity(LEBEDEV_POINT_TYPE::OPTRN48_ABC)
                }};
                template <typename SignedT = int>
                constexpr c_array<SignedT, 6> reconstruct_num_table = {{96, 96, 96, 96, 96, 96}};
                template <size_t NX, typename SignedT = int>
                constexpr SignedT lengthx = NX/4;
                template <size_t NY, typename SignedT = int>
                constexpr SignedT lengthy = NY/4;
                template <size_t NZ, typename SignedT = int>
                constexpr SignedT lengthz = NZ/4;
                template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
                constexpr SignedT fields_required = fields_required_table<SignedT>[static_cast<int>(PType)];
                template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
                constexpr SignedT extract_num = extract_num_table<SignedT>[static_cast<int>(PType)];
                template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
                constexpr SignedT reconstruct_num = reconstruct_num_table<SignedT>[static_cast<int>(PType)];
                template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
                constexpr c_array<SignedT, fields_required<PType, SignedT>> field_indexes;
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN6_00C, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN6_00C, SignedT> = {{0}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN12_0BB, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN12_0BB, SignedT> = {{0}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN8_AAA, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN8_AAA, SignedT> = {{0}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN24_AAC, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN24_AAC, SignedT> = {{0}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN24_AB0, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN24_AB0, SignedT> = {{0}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN48_ABC, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN48_ABC, SignedT> = {{0}};
                template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
                constexpr c_array<SignedT, fields_required<PType, SignedT>> field_submultiplicities;
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN6_00C, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN6_00C, SignedT> = {{6}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN12_0BB, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN12_0BB, SignedT> = {{12}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN8_AAA, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN8_AAA, SignedT> = {{8}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN24_AAC, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN24_AAC, SignedT> = {{24}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN24_AB0, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN24_AB0, SignedT> = {{24}};
                template <typename SignedT>
                constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN48_ABC, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN48_ABC, SignedT> = {{48}};


                namespace kernel
                {
                    /// <summary>
                    /// This function determines whether an element is contained in the selected group3D230 asymmetric unit.
                    /// The asymmetric unit is selected as -1/8 <= y <= 1/8 && -1/8 <= z <= 1/8 && 0 <= x <= 1/4 && Max[-y, y, -z, z] <= x.
                    /// x falls in the range of [0, 1/4].
                    /// y falls in the range of [-1/8, 1/8].
                    /// z falls in the range of [-1/8, 1/8].
                    /// Note that these restriction is never checked.
                    /// </summary>
                    template <size_t NXYZ, typename SignedT = int>
                    inline constexpr bool element(SignedT x, SignedT y, SignedT z)
                    {
                        constexpr SignedT NXYZ8 = (SignedT)(NXYZ / 8);
                        return x >= NXYZ8 ? true :
                            (y >= -x - 1 && y <= x && z >= -x - 1 && z <= x);
                    }

                    /// <summary>
                    /// This function calculates the shift of a group3D230 asymmetric unit.
                    /// The asymmetric unit is selected as -1/8 <= y <= 1/8 && -1/8 <= z <= 1/8 && 0 <= x <= 1/4 && Max[-y, y, -z, z] <= x.
                    /// x falls in the range of [0, 1/4].
                    /// y falls in the range of [-1/8, 1/8].
                    /// z falls in the range of [-1/8, 1/8].
                    /// Note that these restriction is never checked.
                    /// </summary>
                    template <size_t NXYZ, typename SignedT = int>
                    inline constexpr SignedT index_xyz(SignedT x, SignedT y, SignedT z)
                    {
                        constexpr SignedT NXYZ8 = (SignedT)(NXYZ / 8);
                        constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
                        constexpr SignedT half_count = (2 * NXYZ8 * (1 + NXYZ8) * (2 * NXYZ8 + 1)) / 3;
                        return (x < NXYZ8) ?
                            (2 * x * (1 + x) * (2 * x + 1)) / 3 + 2 * (y + x + 1) * (x + 1) + (z + x + 1) :
                            half_count + (x - NXYZ8) * NXYZ4 * NXYZ4 + (y + NXYZ8) * NXYZ4 + (z + NXYZ8);
                    }

                    /// <summary>
                    /// This function calculates the number of the sample points in a group3D230 asymmetric unit.
                    /// The asymmetric unit is selected as -1/8 <= y <= 1/8 && -1/8 <= z <= 1/8 && 0 <= x <= 1/4 && Max[-y, y, -z, z] <= x.
                    /// </summary>
                    template <size_t NXYZ, typename SignedT = int>
                    inline constexpr SignedT total_elements()
                    {
                        constexpr SignedT NXYZ8 = (SignedT)(NXYZ / 8);
                        constexpr SignedT NXYZ4 = (SignedT)(NXYZ / 4);
                        return index_xyz<NXYZ, SignedT>(NXYZ4 - 1, NXYZ8 - 1, NXYZ8 - 1 + 1/*one-past-last*/);
                    }

                    /// <summary>
                    /// This function return the relative positions of the symmetry operations in lebedev grids
                    /// In detail, the data at lebedev point i in spatial-compressed format, 
                    /// is equivalent to the data at extract_field_pos<i,...> in orientational-compressed format
                    /// Ideally this should be consteval function which is not currently supported by CUDA 11.x
                    /// </summary>
                    template <int LebedevID, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
                    inline constexpr SignedT extract_field_pos()
                    {
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
                            constexpr SignedT list[6] =
                            {
                                  0,   0,   0,   0,   0,   0
                            };
                            return list[LebedevID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
                            constexpr SignedT list[12] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[LebedevID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
                            constexpr SignedT list[8] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[LebedevID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
                            constexpr SignedT list[24] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[LebedevID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
                            constexpr SignedT list[24] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[LebedevID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
                            constexpr SignedT list[48] = 
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[LebedevID];
                        }
                        return -1;
                    }

                    /// <summary>
                    /// This function return the relative positions of the symmetry operations in lebedev grids
                    /// In detail, the operation i, mapping data in spatial-compressed format at reconstruct_lebedev_point_pos<i,...>() 
                    /// to the data of asymunit No. i in orientational-compressed format in reconstruct_field_pos<i,...>()
                    /// Ideally this should be consteval function which is not currently supported by CUDA 11.x
                    /// </summary>
                    template <int OptrnID, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
                    inline constexpr SignedT reconstruct_lebedev_point_pos()
                    {
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
                            constexpr SignedT list[96] =
                            {
                                  0,   0,   1,   1,   2,   2,   3,   3,   4,   4,   5,   5, 
                                  1,   1,   0,   0,   4,   4,   5,   5,   2,   2,   3,   3, 
                                  1,   1,   0,   0,   3,   3,   2,   2,   5,   5,   4,   4, 
                                  0,   0,   1,   1,   5,   5,   4,   4,   3,   3,   2,   2, 
                                  0,   0,   1,   1,   2,   2,   3,   3,   4,   4,   5,   5, 
                                  1,   1,   0,   0,   4,   4,   5,   5,   2,   2,   3,   3, 
                                  1,   1,   0,   0,   3,   3,   2,   2,   5,   5,   4,   4, 
                                  0,   0,   1,   1,   5,   5,   4,   4,   3,   3,   2,   2
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
                            constexpr SignedT list[96] =
                            {
                                  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11, 
                                  5,   7,   4,   6,   2,   0,   3,   1,   8,  10,   9,  11, 
                                  3,   2,   1,   0,   7,   6,   5,   4,  11,  10,   9,   8, 
                                  6,   4,   7,   5,   1,   3,   0,   2,  11,   9,  10,   8, 
                                  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11, 
                                  5,   7,   4,   6,   2,   0,   3,   1,   8,  10,   9,  11, 
                                  3,   2,   1,   0,   7,   6,   5,   4,  11,  10,   9,   8, 
                                  6,   4,   7,   5,   1,   3,   0,   2,  11,   9,  10,   8
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
                            constexpr SignedT list[96] =
                            {
                                  0,   1,   2,   3,   0,   3,   1,   2,   0,   2,   3,   1, 
                                  4,   5,   6,   7,   4,   7,   5,   6,   4,   6,   7,   5, 
                                  5,   4,   6,   7,   5,   7,   4,   6,   5,   6,   7,   4, 
                                  1,   0,   2,   3,   1,   3,   0,   2,   1,   2,   3,   0, 
                                  0,   1,   2,   3,   0,   3,   1,   2,   0,   2,   3,   1, 
                                  4,   5,   6,   7,   4,   7,   5,   6,   4,   6,   7,   5, 
                                  5,   4,   6,   7,   5,   7,   4,   6,   5,   6,   7,   4, 
                                  1,   0,   2,   3,   1,   3,   0,   2,   1,   2,   3,   0
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
                            constexpr SignedT list[96] =
                            {
                                  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11, 
                                 12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23, 
                                 13,  12,  14,  15,  23,  22,  20,  21,  18,  19,  17,  16, 
                                  1,   0,   2,   3,  11,  10,   8,   9,   6,   7,   5,   4, 
                                  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11, 
                                 12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23, 
                                 13,  12,  14,  15,  23,  22,  20,  21,  18,  19,  17,  16, 
                                  1,   0,   2,   3,  11,  10,   8,   9,   6,   7,   5,   4
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
                            constexpr SignedT list[96] =
                            {
                                  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11, 
                                 12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23, 
                                  1,   0,   3,   2,   5,   4,   7,   6,   9,   8,  11,  10, 
                                 13,  12,  15,  14,  17,  16,  19,  18,  21,  20,  23,  22, 
                                  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11, 
                                 12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23, 
                                  1,   0,   3,   2,   5,   4,   7,   6,   9,   8,  11,  10, 
                                 13,  12,  15,  14,  17,  16,  19,  18,  21,  20,  23,  22
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
                            constexpr SignedT list[96] = 
                            {
                                  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11, 
                                 12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23, 
                                 24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35, 
                                 36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47, 
                                  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11, 
                                 12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23, 
                                 24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35, 
                                 36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47
                            };
                            return list[OptrnID];
                        }
                        return -1;
                    }

                    /// <summary>
                    /// This function return the relative positions of the symmetry operations in lebedev grids
                    /// In detail, the operation i, mapping data in spatial-compressed format at reconstruct_lebedev_point_pos<i,...>() 
                    /// to the data of asymunit No. i in orientational-compressed format in reconstruct_field_pos<i,...>()
                    /// Ideally this should be consteval function which is not currently supported by CUDA 11.x
                    /// </summary>
                    template <int OptrnID, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
                    inline constexpr SignedT reconstruct_field_pos()
                    {
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
                            constexpr SignedT list[96] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
                            constexpr SignedT list[96] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
                            constexpr SignedT list[96] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
                            constexpr SignedT list[96] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
                            constexpr SignedT list[96] =
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[OptrnID];
                        }
                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
                            constexpr SignedT list[96] = 
                            {
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
                            };
                            return list[OptrnID];
                        }
                        return -1;
                    }

                    /// <summary>
                    /// Note that this kernel requires that outer should be the number of lebedev points.
                    /// However, this is not checked in the kernel. Any block exceeds this constraint ***will return without doing anything***
                    /// </summary>
                    template<size_t NXYZ, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
                    void extract(double* dst_lebedev, const double* src_real) {
                        constexpr SignedT Nx = NXYZ;
                        constexpr SignedT Ny = NXYZ;
                        constexpr SignedT Nz = NXYZ;

                        constexpr SignedT NXYZ8 = NXYZ / 8;
                        static_assert(NXYZ == NXYZ8 * 8, "template parameter NXYZ must be divided by 8");

                        constexpr SignedT Sx = 0;
                        constexpr SignedT Sy = -NXYZ8;
                        constexpr SignedT Sz = -NXYZ8;

                        constexpr SignedT NAsymunit = total_elements<NXYZ, SignedT>();

                        for(SignedT outer = 0; outer < extract_num<PType, SignedT>; ++outer)
                        {
                            for(SignedT x = Sx + 0; x < Sx + lengthx<Nx, SignedT>; x += 1)
                            {
                                for(SignedT y = Sy + 0; y < Sy + lengthy<Ny, SignedT>; y += 1)
                                {
                                    for(SignedT z = Sz + 0; z < Sz + lengthz<Nz, SignedT>; z += 1)
                                    {                               
                                        SignedT Asymunit_index = index_xyz<NXYZ, SignedT>(x, y, z);
                                        SignedT real_index = 0;
                                        
    #define GROUP3D230_INDEX_CASES(LebedevID, x, y, z) \
                                        {\
                                            constexpr SignedT field_shift = extract_field_pos<LebedevID, PType, SignedT>() * Nx * Ny * Nz;\
                                            real_index = field_shift + (((x) % Nx) * Ny + ((y) % Ny)) * Nz + ((z) % Nz); break;\
                                        }

                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 3:GROUP3D230_INDEX_CASES(3, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 7:GROUP3D230_INDEX_CASES(7, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 24:GROUP3D230_INDEX_CASES(24, (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 25:GROUP3D230_INDEX_CASES(25, (4*NXYZ8 + x), (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 26:GROUP3D230_INDEX_CASES(26, (x), (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 27:GROUP3D230_INDEX_CASES(27, (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 28:GROUP3D230_INDEX_CASES(28, (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 29:GROUP3D230_INDEX_CASES(29, (4*NXYZ8 + y), (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 30:GROUP3D230_INDEX_CASES(30, (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 31:GROUP3D230_INDEX_CASES(31, (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z), (x));
                                            case 32:GROUP3D230_INDEX_CASES(32, (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 33:GROUP3D230_INDEX_CASES(33, (4*NXYZ8 + z), (x), (-1 + 4*NXYZ8 - y));
                                            case 34:GROUP3D230_INDEX_CASES(34, (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 35:GROUP3D230_INDEX_CASES(35, (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x), (8*NXYZ8 + y));
                                            case 36:GROUP3D230_INDEX_CASES(36, (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 37:GROUP3D230_INDEX_CASES(37, (6*NXYZ8 + y), (6*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 38:GROUP3D230_INDEX_CASES(38, (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 39:GROUP3D230_INDEX_CASES(39, (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 40:GROUP3D230_INDEX_CASES(40, (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 41:GROUP3D230_INDEX_CASES(41, (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 42:GROUP3D230_INDEX_CASES(42, (6*NXYZ8 + x), (6*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 43:GROUP3D230_INDEX_CASES(43, (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 44:GROUP3D230_INDEX_CASES(44, (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 45:GROUP3D230_INDEX_CASES(45, (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 46:GROUP3D230_INDEX_CASES(46, (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 47:GROUP3D230_INDEX_CASES(47, (6*NXYZ8 + z), (6*NXYZ8 + y), (6*NXYZ8 + x));
                                            default: return;
                                            }
                                        }
    #ifdef GROUP3D230_INDEX_CASES
    #undef GROUP3D230_INDEX_CASES
    #endif // GROUP3D230_INDEX_CASES
                                        if (element<NXYZ>(x, y, z)) {
                                            dst_lebedev[outer * NAsymunit + Asymunit_index] = src_real[real_index];
                                        }
                                    }
                                }
                            }
                        }

                    }
                    
                    /// <summary>
                    /// Note that this kernel requires that outer should be the number of operations multiples by number of leading lebedev points.
                    /// However, this is not checked in the kernel. Any block exceeds this constraint ***will return without doing anything***
                    /// </summary>
                    template<size_t NXYZ, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
                    void reconstruct(double* dst_real, const double* src_lebedev) {
                        constexpr SignedT Nx = NXYZ;
                        constexpr SignedT Ny = NXYZ;
                        constexpr SignedT Nz = NXYZ;

                        constexpr SignedT NXYZ8 = NXYZ / 8;
                        static_assert(NXYZ == NXYZ8 * 8, "template parameter NXYZ must be divided by 8");

                        constexpr SignedT Sx = 0;
                        constexpr SignedT Sy = -NXYZ8;
                        constexpr SignedT Sz = -NXYZ8;

                        constexpr SignedT NAsymunit = total_elements<NXYZ, SignedT>();

                        for(SignedT outer = 0; outer < reconstruct_num<PType, SignedT>; ++outer)
                        {
                            for(SignedT x = Sx + 0; x < Sx + lengthx<Nx, SignedT>; x += 1)
                            {
                                for(SignedT y = Sy + 0; y < Sy + lengthy<Ny, SignedT>; y += 1)
                                {
                                    for(SignedT z = Sz + 0; z < Sz + lengthz<Nz, SignedT>; z += 1)
                                    {                               
                                        SignedT Asymunit_index = index_xyz<NXYZ, SignedT>(x, y, z);
                                        SignedT lebdev_pos = 0;
                                        SignedT real_index = 0;

    #define GROUP3D230_INDEX_CASES(OptrnID, x, y, z) \
                                        {\
                                            constexpr SignedT temppos = reconstruct_lebedev_point_pos<OptrnID, PType, SignedT>();\
                                            lebdev_pos = temppos;\
                                            constexpr SignedT field_shift = reconstruct_field_pos<OptrnID, PType, SignedT>() * Nx * Ny * Nz;\
                                            real_index = field_shift + (((x) % Nx) * Ny + ((y) % Ny)) * Nz + ((z) % Nz); break;\
                                        }

                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 24:GROUP3D230_INDEX_CASES(24, (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 25:GROUP3D230_INDEX_CASES(25, (4*NXYZ8 + x), (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 26:GROUP3D230_INDEX_CASES(26, (x), (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 27:GROUP3D230_INDEX_CASES(27, (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 28:GROUP3D230_INDEX_CASES(28, (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 29:GROUP3D230_INDEX_CASES(29, (4*NXYZ8 + y), (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 30:GROUP3D230_INDEX_CASES(30, (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 31:GROUP3D230_INDEX_CASES(31, (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z), (x));
                                            case 32:GROUP3D230_INDEX_CASES(32, (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 33:GROUP3D230_INDEX_CASES(33, (4*NXYZ8 + z), (x), (-1 + 4*NXYZ8 - y));
                                            case 34:GROUP3D230_INDEX_CASES(34, (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 35:GROUP3D230_INDEX_CASES(35, (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x), (8*NXYZ8 + y));
                                            case 36:GROUP3D230_INDEX_CASES(36, (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 37:GROUP3D230_INDEX_CASES(37, (6*NXYZ8 + y), (6*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 38:GROUP3D230_INDEX_CASES(38, (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 39:GROUP3D230_INDEX_CASES(39, (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 40:GROUP3D230_INDEX_CASES(40, (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 41:GROUP3D230_INDEX_CASES(41, (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 42:GROUP3D230_INDEX_CASES(42, (6*NXYZ8 + x), (6*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 43:GROUP3D230_INDEX_CASES(43, (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 44:GROUP3D230_INDEX_CASES(44, (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 45:GROUP3D230_INDEX_CASES(45, (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 46:GROUP3D230_INDEX_CASES(46, (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 47:GROUP3D230_INDEX_CASES(47, (6*NXYZ8 + z), (6*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 48:GROUP3D230_INDEX_CASES(48, (4*NXYZ8 + x), (4*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 49:GROUP3D230_INDEX_CASES(49, (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 50:GROUP3D230_INDEX_CASES(50, (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 51:GROUP3D230_INDEX_CASES(51, (x), (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 52:GROUP3D230_INDEX_CASES(52, (4*NXYZ8 + y), (4*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 53:GROUP3D230_INDEX_CASES(53, (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (x));
                                            case 54:GROUP3D230_INDEX_CASES(54, (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 55:GROUP3D230_INDEX_CASES(55, (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 56:GROUP3D230_INDEX_CASES(56, (4*NXYZ8 + z), (4*NXYZ8 + x), (4*NXYZ8 + y));
                                            case 57:GROUP3D230_INDEX_CASES(57, (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 58:GROUP3D230_INDEX_CASES(58, (-1 + 4*NXYZ8 - z), (x), (-1 + 8*NXYZ8 - y));
                                            case 59:GROUP3D230_INDEX_CASES(59, (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 60:GROUP3D230_INDEX_CASES(60, (2*NXYZ8 + y), (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 61:GROUP3D230_INDEX_CASES(61, (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 62:GROUP3D230_INDEX_CASES(62, (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 63:GROUP3D230_INDEX_CASES(63, (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 64:GROUP3D230_INDEX_CASES(64, (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 65:GROUP3D230_INDEX_CASES(65, (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 66:GROUP3D230_INDEX_CASES(66, (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 67:GROUP3D230_INDEX_CASES(67, (2*NXYZ8 + x), (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 68:GROUP3D230_INDEX_CASES(68, (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 69:GROUP3D230_INDEX_CASES(69, (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 70:GROUP3D230_INDEX_CASES(70, (2*NXYZ8 + z), (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 71:GROUP3D230_INDEX_CASES(71, (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 72:GROUP3D230_INDEX_CASES(72, (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 73:GROUP3D230_INDEX_CASES(73, (x), (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 74:GROUP3D230_INDEX_CASES(74, (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 75:GROUP3D230_INDEX_CASES(75, (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 76:GROUP3D230_INDEX_CASES(76, (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 77:GROUP3D230_INDEX_CASES(77, (8*NXYZ8 + y), (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 78:GROUP3D230_INDEX_CASES(78, (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (x));
                                            case 79:GROUP3D230_INDEX_CASES(79, (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 80:GROUP3D230_INDEX_CASES(80, (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 81:GROUP3D230_INDEX_CASES(81, (8*NXYZ8 + z), (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y));
                                            case 82:GROUP3D230_INDEX_CASES(82, (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 83:GROUP3D230_INDEX_CASES(83, (-1 + 8*NXYZ8 - z), (x), (4*NXYZ8 + y));
                                            case 84:GROUP3D230_INDEX_CASES(84, (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 85:GROUP3D230_INDEX_CASES(85, (2*NXYZ8 + y), (2*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 86:GROUP3D230_INDEX_CASES(86, (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 87:GROUP3D230_INDEX_CASES(87, (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 88:GROUP3D230_INDEX_CASES(88, (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 89:GROUP3D230_INDEX_CASES(89, (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 90:GROUP3D230_INDEX_CASES(90, (2*NXYZ8 + x), (2*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 91:GROUP3D230_INDEX_CASES(91, (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 92:GROUP3D230_INDEX_CASES(92, (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 93:GROUP3D230_INDEX_CASES(93, (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 94:GROUP3D230_INDEX_CASES(94, (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 95:GROUP3D230_INDEX_CASES(95, (2*NXYZ8 + z), (2*NXYZ8 + y), (2*NXYZ8 + x));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 24:GROUP3D230_INDEX_CASES(24, (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 25:GROUP3D230_INDEX_CASES(25, (4*NXYZ8 + x), (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 26:GROUP3D230_INDEX_CASES(26, (x), (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 27:GROUP3D230_INDEX_CASES(27, (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 28:GROUP3D230_INDEX_CASES(28, (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 29:GROUP3D230_INDEX_CASES(29, (4*NXYZ8 + y), (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 30:GROUP3D230_INDEX_CASES(30, (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 31:GROUP3D230_INDEX_CASES(31, (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z), (x));
                                            case 32:GROUP3D230_INDEX_CASES(32, (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 33:GROUP3D230_INDEX_CASES(33, (4*NXYZ8 + z), (x), (-1 + 4*NXYZ8 - y));
                                            case 34:GROUP3D230_INDEX_CASES(34, (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 35:GROUP3D230_INDEX_CASES(35, (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x), (8*NXYZ8 + y));
                                            case 36:GROUP3D230_INDEX_CASES(36, (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 37:GROUP3D230_INDEX_CASES(37, (6*NXYZ8 + y), (6*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 38:GROUP3D230_INDEX_CASES(38, (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 39:GROUP3D230_INDEX_CASES(39, (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 40:GROUP3D230_INDEX_CASES(40, (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 41:GROUP3D230_INDEX_CASES(41, (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 42:GROUP3D230_INDEX_CASES(42, (6*NXYZ8 + x), (6*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 43:GROUP3D230_INDEX_CASES(43, (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 44:GROUP3D230_INDEX_CASES(44, (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 45:GROUP3D230_INDEX_CASES(45, (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 46:GROUP3D230_INDEX_CASES(46, (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 47:GROUP3D230_INDEX_CASES(47, (6*NXYZ8 + z), (6*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 48:GROUP3D230_INDEX_CASES(48, (4*NXYZ8 + x), (4*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 49:GROUP3D230_INDEX_CASES(49, (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 50:GROUP3D230_INDEX_CASES(50, (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 51:GROUP3D230_INDEX_CASES(51, (x), (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 52:GROUP3D230_INDEX_CASES(52, (4*NXYZ8 + y), (4*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 53:GROUP3D230_INDEX_CASES(53, (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (x));
                                            case 54:GROUP3D230_INDEX_CASES(54, (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 55:GROUP3D230_INDEX_CASES(55, (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 56:GROUP3D230_INDEX_CASES(56, (4*NXYZ8 + z), (4*NXYZ8 + x), (4*NXYZ8 + y));
                                            case 57:GROUP3D230_INDEX_CASES(57, (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 58:GROUP3D230_INDEX_CASES(58, (-1 + 4*NXYZ8 - z), (x), (-1 + 8*NXYZ8 - y));
                                            case 59:GROUP3D230_INDEX_CASES(59, (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 60:GROUP3D230_INDEX_CASES(60, (2*NXYZ8 + y), (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 61:GROUP3D230_INDEX_CASES(61, (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 62:GROUP3D230_INDEX_CASES(62, (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 63:GROUP3D230_INDEX_CASES(63, (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 64:GROUP3D230_INDEX_CASES(64, (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 65:GROUP3D230_INDEX_CASES(65, (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 66:GROUP3D230_INDEX_CASES(66, (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 67:GROUP3D230_INDEX_CASES(67, (2*NXYZ8 + x), (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 68:GROUP3D230_INDEX_CASES(68, (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 69:GROUP3D230_INDEX_CASES(69, (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 70:GROUP3D230_INDEX_CASES(70, (2*NXYZ8 + z), (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 71:GROUP3D230_INDEX_CASES(71, (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 72:GROUP3D230_INDEX_CASES(72, (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 73:GROUP3D230_INDEX_CASES(73, (x), (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 74:GROUP3D230_INDEX_CASES(74, (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 75:GROUP3D230_INDEX_CASES(75, (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 76:GROUP3D230_INDEX_CASES(76, (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 77:GROUP3D230_INDEX_CASES(77, (8*NXYZ8 + y), (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 78:GROUP3D230_INDEX_CASES(78, (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (x));
                                            case 79:GROUP3D230_INDEX_CASES(79, (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 80:GROUP3D230_INDEX_CASES(80, (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 81:GROUP3D230_INDEX_CASES(81, (8*NXYZ8 + z), (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y));
                                            case 82:GROUP3D230_INDEX_CASES(82, (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 83:GROUP3D230_INDEX_CASES(83, (-1 + 8*NXYZ8 - z), (x), (4*NXYZ8 + y));
                                            case 84:GROUP3D230_INDEX_CASES(84, (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 85:GROUP3D230_INDEX_CASES(85, (2*NXYZ8 + y), (2*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 86:GROUP3D230_INDEX_CASES(86, (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 87:GROUP3D230_INDEX_CASES(87, (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 88:GROUP3D230_INDEX_CASES(88, (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 89:GROUP3D230_INDEX_CASES(89, (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 90:GROUP3D230_INDEX_CASES(90, (2*NXYZ8 + x), (2*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 91:GROUP3D230_INDEX_CASES(91, (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 92:GROUP3D230_INDEX_CASES(92, (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 93:GROUP3D230_INDEX_CASES(93, (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 94:GROUP3D230_INDEX_CASES(94, (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 95:GROUP3D230_INDEX_CASES(95, (2*NXYZ8 + z), (2*NXYZ8 + y), (2*NXYZ8 + x));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 24:GROUP3D230_INDEX_CASES(24, (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 25:GROUP3D230_INDEX_CASES(25, (4*NXYZ8 + x), (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 26:GROUP3D230_INDEX_CASES(26, (x), (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 27:GROUP3D230_INDEX_CASES(27, (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 28:GROUP3D230_INDEX_CASES(28, (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 29:GROUP3D230_INDEX_CASES(29, (4*NXYZ8 + y), (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 30:GROUP3D230_INDEX_CASES(30, (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 31:GROUP3D230_INDEX_CASES(31, (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z), (x));
                                            case 32:GROUP3D230_INDEX_CASES(32, (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 33:GROUP3D230_INDEX_CASES(33, (4*NXYZ8 + z), (x), (-1 + 4*NXYZ8 - y));
                                            case 34:GROUP3D230_INDEX_CASES(34, (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 35:GROUP3D230_INDEX_CASES(35, (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x), (8*NXYZ8 + y));
                                            case 36:GROUP3D230_INDEX_CASES(36, (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 37:GROUP3D230_INDEX_CASES(37, (6*NXYZ8 + y), (6*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 38:GROUP3D230_INDEX_CASES(38, (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 39:GROUP3D230_INDEX_CASES(39, (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 40:GROUP3D230_INDEX_CASES(40, (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 41:GROUP3D230_INDEX_CASES(41, (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 42:GROUP3D230_INDEX_CASES(42, (6*NXYZ8 + x), (6*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 43:GROUP3D230_INDEX_CASES(43, (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 44:GROUP3D230_INDEX_CASES(44, (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 45:GROUP3D230_INDEX_CASES(45, (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 46:GROUP3D230_INDEX_CASES(46, (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 47:GROUP3D230_INDEX_CASES(47, (6*NXYZ8 + z), (6*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 48:GROUP3D230_INDEX_CASES(48, (4*NXYZ8 + x), (4*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 49:GROUP3D230_INDEX_CASES(49, (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 50:GROUP3D230_INDEX_CASES(50, (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 51:GROUP3D230_INDEX_CASES(51, (x), (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 52:GROUP3D230_INDEX_CASES(52, (4*NXYZ8 + y), (4*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 53:GROUP3D230_INDEX_CASES(53, (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (x));
                                            case 54:GROUP3D230_INDEX_CASES(54, (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 55:GROUP3D230_INDEX_CASES(55, (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 56:GROUP3D230_INDEX_CASES(56, (4*NXYZ8 + z), (4*NXYZ8 + x), (4*NXYZ8 + y));
                                            case 57:GROUP3D230_INDEX_CASES(57, (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 58:GROUP3D230_INDEX_CASES(58, (-1 + 4*NXYZ8 - z), (x), (-1 + 8*NXYZ8 - y));
                                            case 59:GROUP3D230_INDEX_CASES(59, (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 60:GROUP3D230_INDEX_CASES(60, (2*NXYZ8 + y), (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 61:GROUP3D230_INDEX_CASES(61, (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 62:GROUP3D230_INDEX_CASES(62, (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 63:GROUP3D230_INDEX_CASES(63, (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 64:GROUP3D230_INDEX_CASES(64, (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 65:GROUP3D230_INDEX_CASES(65, (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 66:GROUP3D230_INDEX_CASES(66, (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 67:GROUP3D230_INDEX_CASES(67, (2*NXYZ8 + x), (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 68:GROUP3D230_INDEX_CASES(68, (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 69:GROUP3D230_INDEX_CASES(69, (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 70:GROUP3D230_INDEX_CASES(70, (2*NXYZ8 + z), (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 71:GROUP3D230_INDEX_CASES(71, (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 72:GROUP3D230_INDEX_CASES(72, (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 73:GROUP3D230_INDEX_CASES(73, (x), (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 74:GROUP3D230_INDEX_CASES(74, (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 75:GROUP3D230_INDEX_CASES(75, (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 76:GROUP3D230_INDEX_CASES(76, (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 77:GROUP3D230_INDEX_CASES(77, (8*NXYZ8 + y), (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 78:GROUP3D230_INDEX_CASES(78, (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (x));
                                            case 79:GROUP3D230_INDEX_CASES(79, (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 80:GROUP3D230_INDEX_CASES(80, (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 81:GROUP3D230_INDEX_CASES(81, (8*NXYZ8 + z), (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y));
                                            case 82:GROUP3D230_INDEX_CASES(82, (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 83:GROUP3D230_INDEX_CASES(83, (-1 + 8*NXYZ8 - z), (x), (4*NXYZ8 + y));
                                            case 84:GROUP3D230_INDEX_CASES(84, (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 85:GROUP3D230_INDEX_CASES(85, (2*NXYZ8 + y), (2*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 86:GROUP3D230_INDEX_CASES(86, (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 87:GROUP3D230_INDEX_CASES(87, (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 88:GROUP3D230_INDEX_CASES(88, (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 89:GROUP3D230_INDEX_CASES(89, (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 90:GROUP3D230_INDEX_CASES(90, (2*NXYZ8 + x), (2*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 91:GROUP3D230_INDEX_CASES(91, (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 92:GROUP3D230_INDEX_CASES(92, (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 93:GROUP3D230_INDEX_CASES(93, (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 94:GROUP3D230_INDEX_CASES(94, (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 95:GROUP3D230_INDEX_CASES(95, (2*NXYZ8 + z), (2*NXYZ8 + y), (2*NXYZ8 + x));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 24:GROUP3D230_INDEX_CASES(24, (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 25:GROUP3D230_INDEX_CASES(25, (4*NXYZ8 + x), (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 26:GROUP3D230_INDEX_CASES(26, (x), (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 27:GROUP3D230_INDEX_CASES(27, (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 28:GROUP3D230_INDEX_CASES(28, (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 29:GROUP3D230_INDEX_CASES(29, (4*NXYZ8 + y), (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 30:GROUP3D230_INDEX_CASES(30, (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 31:GROUP3D230_INDEX_CASES(31, (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z), (x));
                                            case 32:GROUP3D230_INDEX_CASES(32, (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 33:GROUP3D230_INDEX_CASES(33, (4*NXYZ8 + z), (x), (-1 + 4*NXYZ8 - y));
                                            case 34:GROUP3D230_INDEX_CASES(34, (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 35:GROUP3D230_INDEX_CASES(35, (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x), (8*NXYZ8 + y));
                                            case 36:GROUP3D230_INDEX_CASES(36, (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 37:GROUP3D230_INDEX_CASES(37, (6*NXYZ8 + y), (6*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 38:GROUP3D230_INDEX_CASES(38, (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 39:GROUP3D230_INDEX_CASES(39, (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 40:GROUP3D230_INDEX_CASES(40, (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 41:GROUP3D230_INDEX_CASES(41, (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 42:GROUP3D230_INDEX_CASES(42, (6*NXYZ8 + x), (6*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 43:GROUP3D230_INDEX_CASES(43, (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 44:GROUP3D230_INDEX_CASES(44, (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 45:GROUP3D230_INDEX_CASES(45, (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 46:GROUP3D230_INDEX_CASES(46, (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 47:GROUP3D230_INDEX_CASES(47, (6*NXYZ8 + z), (6*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 48:GROUP3D230_INDEX_CASES(48, (4*NXYZ8 + x), (4*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 49:GROUP3D230_INDEX_CASES(49, (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 50:GROUP3D230_INDEX_CASES(50, (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 51:GROUP3D230_INDEX_CASES(51, (x), (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 52:GROUP3D230_INDEX_CASES(52, (4*NXYZ8 + y), (4*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 53:GROUP3D230_INDEX_CASES(53, (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (x));
                                            case 54:GROUP3D230_INDEX_CASES(54, (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 55:GROUP3D230_INDEX_CASES(55, (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 56:GROUP3D230_INDEX_CASES(56, (4*NXYZ8 + z), (4*NXYZ8 + x), (4*NXYZ8 + y));
                                            case 57:GROUP3D230_INDEX_CASES(57, (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 58:GROUP3D230_INDEX_CASES(58, (-1 + 4*NXYZ8 - z), (x), (-1 + 8*NXYZ8 - y));
                                            case 59:GROUP3D230_INDEX_CASES(59, (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 60:GROUP3D230_INDEX_CASES(60, (2*NXYZ8 + y), (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 61:GROUP3D230_INDEX_CASES(61, (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 62:GROUP3D230_INDEX_CASES(62, (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 63:GROUP3D230_INDEX_CASES(63, (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 64:GROUP3D230_INDEX_CASES(64, (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 65:GROUP3D230_INDEX_CASES(65, (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 66:GROUP3D230_INDEX_CASES(66, (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 67:GROUP3D230_INDEX_CASES(67, (2*NXYZ8 + x), (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 68:GROUP3D230_INDEX_CASES(68, (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 69:GROUP3D230_INDEX_CASES(69, (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 70:GROUP3D230_INDEX_CASES(70, (2*NXYZ8 + z), (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 71:GROUP3D230_INDEX_CASES(71, (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 72:GROUP3D230_INDEX_CASES(72, (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 73:GROUP3D230_INDEX_CASES(73, (x), (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 74:GROUP3D230_INDEX_CASES(74, (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 75:GROUP3D230_INDEX_CASES(75, (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 76:GROUP3D230_INDEX_CASES(76, (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 77:GROUP3D230_INDEX_CASES(77, (8*NXYZ8 + y), (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 78:GROUP3D230_INDEX_CASES(78, (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (x));
                                            case 79:GROUP3D230_INDEX_CASES(79, (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 80:GROUP3D230_INDEX_CASES(80, (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 81:GROUP3D230_INDEX_CASES(81, (8*NXYZ8 + z), (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y));
                                            case 82:GROUP3D230_INDEX_CASES(82, (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 83:GROUP3D230_INDEX_CASES(83, (-1 + 8*NXYZ8 - z), (x), (4*NXYZ8 + y));
                                            case 84:GROUP3D230_INDEX_CASES(84, (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 85:GROUP3D230_INDEX_CASES(85, (2*NXYZ8 + y), (2*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 86:GROUP3D230_INDEX_CASES(86, (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 87:GROUP3D230_INDEX_CASES(87, (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 88:GROUP3D230_INDEX_CASES(88, (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 89:GROUP3D230_INDEX_CASES(89, (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 90:GROUP3D230_INDEX_CASES(90, (2*NXYZ8 + x), (2*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 91:GROUP3D230_INDEX_CASES(91, (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 92:GROUP3D230_INDEX_CASES(92, (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 93:GROUP3D230_INDEX_CASES(93, (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 94:GROUP3D230_INDEX_CASES(94, (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 95:GROUP3D230_INDEX_CASES(95, (2*NXYZ8 + z), (2*NXYZ8 + y), (2*NXYZ8 + x));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 24:GROUP3D230_INDEX_CASES(24, (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 25:GROUP3D230_INDEX_CASES(25, (4*NXYZ8 + x), (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 26:GROUP3D230_INDEX_CASES(26, (x), (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 27:GROUP3D230_INDEX_CASES(27, (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 28:GROUP3D230_INDEX_CASES(28, (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 29:GROUP3D230_INDEX_CASES(29, (4*NXYZ8 + y), (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 30:GROUP3D230_INDEX_CASES(30, (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 31:GROUP3D230_INDEX_CASES(31, (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z), (x));
                                            case 32:GROUP3D230_INDEX_CASES(32, (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 33:GROUP3D230_INDEX_CASES(33, (4*NXYZ8 + z), (x), (-1 + 4*NXYZ8 - y));
                                            case 34:GROUP3D230_INDEX_CASES(34, (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 35:GROUP3D230_INDEX_CASES(35, (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x), (8*NXYZ8 + y));
                                            case 36:GROUP3D230_INDEX_CASES(36, (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 37:GROUP3D230_INDEX_CASES(37, (6*NXYZ8 + y), (6*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 38:GROUP3D230_INDEX_CASES(38, (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 39:GROUP3D230_INDEX_CASES(39, (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 40:GROUP3D230_INDEX_CASES(40, (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 41:GROUP3D230_INDEX_CASES(41, (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 42:GROUP3D230_INDEX_CASES(42, (6*NXYZ8 + x), (6*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 43:GROUP3D230_INDEX_CASES(43, (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 44:GROUP3D230_INDEX_CASES(44, (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 45:GROUP3D230_INDEX_CASES(45, (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 46:GROUP3D230_INDEX_CASES(46, (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 47:GROUP3D230_INDEX_CASES(47, (6*NXYZ8 + z), (6*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 48:GROUP3D230_INDEX_CASES(48, (4*NXYZ8 + x), (4*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 49:GROUP3D230_INDEX_CASES(49, (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 50:GROUP3D230_INDEX_CASES(50, (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 51:GROUP3D230_INDEX_CASES(51, (x), (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 52:GROUP3D230_INDEX_CASES(52, (4*NXYZ8 + y), (4*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 53:GROUP3D230_INDEX_CASES(53, (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (x));
                                            case 54:GROUP3D230_INDEX_CASES(54, (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 55:GROUP3D230_INDEX_CASES(55, (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 56:GROUP3D230_INDEX_CASES(56, (4*NXYZ8 + z), (4*NXYZ8 + x), (4*NXYZ8 + y));
                                            case 57:GROUP3D230_INDEX_CASES(57, (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 58:GROUP3D230_INDEX_CASES(58, (-1 + 4*NXYZ8 - z), (x), (-1 + 8*NXYZ8 - y));
                                            case 59:GROUP3D230_INDEX_CASES(59, (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 60:GROUP3D230_INDEX_CASES(60, (2*NXYZ8 + y), (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 61:GROUP3D230_INDEX_CASES(61, (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 62:GROUP3D230_INDEX_CASES(62, (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 63:GROUP3D230_INDEX_CASES(63, (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 64:GROUP3D230_INDEX_CASES(64, (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 65:GROUP3D230_INDEX_CASES(65, (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 66:GROUP3D230_INDEX_CASES(66, (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 67:GROUP3D230_INDEX_CASES(67, (2*NXYZ8 + x), (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 68:GROUP3D230_INDEX_CASES(68, (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 69:GROUP3D230_INDEX_CASES(69, (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 70:GROUP3D230_INDEX_CASES(70, (2*NXYZ8 + z), (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 71:GROUP3D230_INDEX_CASES(71, (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 72:GROUP3D230_INDEX_CASES(72, (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 73:GROUP3D230_INDEX_CASES(73, (x), (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 74:GROUP3D230_INDEX_CASES(74, (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 75:GROUP3D230_INDEX_CASES(75, (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 76:GROUP3D230_INDEX_CASES(76, (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 77:GROUP3D230_INDEX_CASES(77, (8*NXYZ8 + y), (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 78:GROUP3D230_INDEX_CASES(78, (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (x));
                                            case 79:GROUP3D230_INDEX_CASES(79, (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 80:GROUP3D230_INDEX_CASES(80, (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 81:GROUP3D230_INDEX_CASES(81, (8*NXYZ8 + z), (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y));
                                            case 82:GROUP3D230_INDEX_CASES(82, (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 83:GROUP3D230_INDEX_CASES(83, (-1 + 8*NXYZ8 - z), (x), (4*NXYZ8 + y));
                                            case 84:GROUP3D230_INDEX_CASES(84, (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 85:GROUP3D230_INDEX_CASES(85, (2*NXYZ8 + y), (2*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 86:GROUP3D230_INDEX_CASES(86, (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 87:GROUP3D230_INDEX_CASES(87, (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 88:GROUP3D230_INDEX_CASES(88, (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 89:GROUP3D230_INDEX_CASES(89, (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 90:GROUP3D230_INDEX_CASES(90, (2*NXYZ8 + x), (2*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 91:GROUP3D230_INDEX_CASES(91, (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 92:GROUP3D230_INDEX_CASES(92, (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 93:GROUP3D230_INDEX_CASES(93, (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 94:GROUP3D230_INDEX_CASES(94, (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 95:GROUP3D230_INDEX_CASES(95, (2*NXYZ8 + z), (2*NXYZ8 + y), (2*NXYZ8 + x));
                                            default: return;
                                            }
                                        }
                                        if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
                                            switch (outer)
                                            {
                                            case 0:GROUP3D230_INDEX_CASES(0, (x), (8*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 1:GROUP3D230_INDEX_CASES(1, (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 2:GROUP3D230_INDEX_CASES(2, (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 3:GROUP3D230_INDEX_CASES(3, (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 4:GROUP3D230_INDEX_CASES(4, (8*NXYZ8 + y), (8*NXYZ8 + z), (x));
                                            case 5:GROUP3D230_INDEX_CASES(5, (-1 + 4*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 6:GROUP3D230_INDEX_CASES(6, (-1 + 8*NXYZ8 - y), (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 7:GROUP3D230_INDEX_CASES(7, (4*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 8:GROUP3D230_INDEX_CASES(8, (8*NXYZ8 + z), (x), (8*NXYZ8 + y));
                                            case 9:GROUP3D230_INDEX_CASES(9, (-1 + 4*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 10:GROUP3D230_INDEX_CASES(10, (-1 + 8*NXYZ8 - z), (4*NXYZ8 + x), (-1 + 4*NXYZ8 - y));
                                            case 11:GROUP3D230_INDEX_CASES(11, (4*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 12:GROUP3D230_INDEX_CASES(12, (6*NXYZ8 + y), (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 13:GROUP3D230_INDEX_CASES(13, (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 14:GROUP3D230_INDEX_CASES(14, (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 15:GROUP3D230_INDEX_CASES(15, (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 16:GROUP3D230_INDEX_CASES(16, (2*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 17:GROUP3D230_INDEX_CASES(17, (-1 + 2*NXYZ8 - x), (6*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 18:GROUP3D230_INDEX_CASES(18, (-1 + 6*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 19:GROUP3D230_INDEX_CASES(19, (6*NXYZ8 + x), (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 20:GROUP3D230_INDEX_CASES(20, (-1 + 2*NXYZ8 - z), (6*NXYZ8 + y), (2*NXYZ8 + x));
                                            case 21:GROUP3D230_INDEX_CASES(21, (2*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 22:GROUP3D230_INDEX_CASES(22, (6*NXYZ8 + z), (2*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 23:GROUP3D230_INDEX_CASES(23, (-1 + 6*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 24:GROUP3D230_INDEX_CASES(24, (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z));
                                            case 25:GROUP3D230_INDEX_CASES(25, (4*NXYZ8 + x), (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z));
                                            case 26:GROUP3D230_INDEX_CASES(26, (x), (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z));
                                            case 27:GROUP3D230_INDEX_CASES(27, (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y), (8*NXYZ8 + z));
                                            case 28:GROUP3D230_INDEX_CASES(28, (-1 + 8*NXYZ8 - y), (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x));
                                            case 29:GROUP3D230_INDEX_CASES(29, (4*NXYZ8 + y), (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x));
                                            case 30:GROUP3D230_INDEX_CASES(30, (8*NXYZ8 + y), (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x));
                                            case 31:GROUP3D230_INDEX_CASES(31, (-1 + 4*NXYZ8 - y), (4*NXYZ8 + z), (x));
                                            case 32:GROUP3D230_INDEX_CASES(32, (-1 + 8*NXYZ8 - z), (-1 + 8*NXYZ8 - x), (-1 + 8*NXYZ8 - y));
                                            case 33:GROUP3D230_INDEX_CASES(33, (4*NXYZ8 + z), (x), (-1 + 4*NXYZ8 - y));
                                            case 34:GROUP3D230_INDEX_CASES(34, (8*NXYZ8 + z), (-1 + 4*NXYZ8 - x), (4*NXYZ8 + y));
                                            case 35:GROUP3D230_INDEX_CASES(35, (-1 + 4*NXYZ8 - z), (4*NXYZ8 + x), (8*NXYZ8 + y));
                                            case 36:GROUP3D230_INDEX_CASES(36, (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 37:GROUP3D230_INDEX_CASES(37, (6*NXYZ8 + y), (6*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 38:GROUP3D230_INDEX_CASES(38, (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 39:GROUP3D230_INDEX_CASES(39, (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 40:GROUP3D230_INDEX_CASES(40, (-1 + 2*NXYZ8 - x), (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 41:GROUP3D230_INDEX_CASES(41, (2*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 42:GROUP3D230_INDEX_CASES(42, (6*NXYZ8 + x), (6*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 43:GROUP3D230_INDEX_CASES(43, (-1 + 6*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 44:GROUP3D230_INDEX_CASES(44, (2*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 45:GROUP3D230_INDEX_CASES(45, (-1 + 2*NXYZ8 - z), (2*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 46:GROUP3D230_INDEX_CASES(46, (-1 + 6*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 47:GROUP3D230_INDEX_CASES(47, (6*NXYZ8 + z), (6*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 48:GROUP3D230_INDEX_CASES(48, (4*NXYZ8 + x), (4*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 49:GROUP3D230_INDEX_CASES(49, (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 50:GROUP3D230_INDEX_CASES(50, (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 51:GROUP3D230_INDEX_CASES(51, (x), (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 52:GROUP3D230_INDEX_CASES(52, (4*NXYZ8 + y), (4*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 53:GROUP3D230_INDEX_CASES(53, (-1 + 8*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (x));
                                            case 54:GROUP3D230_INDEX_CASES(54, (-1 + 4*NXYZ8 - y), (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 55:GROUP3D230_INDEX_CASES(55, (8*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 56:GROUP3D230_INDEX_CASES(56, (4*NXYZ8 + z), (4*NXYZ8 + x), (4*NXYZ8 + y));
                                            case 57:GROUP3D230_INDEX_CASES(57, (-1 + 8*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 58:GROUP3D230_INDEX_CASES(58, (-1 + 4*NXYZ8 - z), (x), (-1 + 8*NXYZ8 - y));
                                            case 59:GROUP3D230_INDEX_CASES(59, (8*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 60:GROUP3D230_INDEX_CASES(60, (2*NXYZ8 + y), (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z));
                                            case 61:GROUP3D230_INDEX_CASES(61, (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z));
                                            case 62:GROUP3D230_INDEX_CASES(62, (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x), (6*NXYZ8 + z));
                                            case 63:GROUP3D230_INDEX_CASES(63, (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z));
                                            case 64:GROUP3D230_INDEX_CASES(64, (6*NXYZ8 + x), (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y));
                                            case 65:GROUP3D230_INDEX_CASES(65, (-1 + 6*NXYZ8 - x), (2*NXYZ8 + z), (6*NXYZ8 + y));
                                            case 66:GROUP3D230_INDEX_CASES(66, (-1 + 2*NXYZ8 - x), (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y));
                                            case 67:GROUP3D230_INDEX_CASES(67, (2*NXYZ8 + x), (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y));
                                            case 68:GROUP3D230_INDEX_CASES(68, (-1 + 6*NXYZ8 - z), (2*NXYZ8 + y), (6*NXYZ8 + x));
                                            case 69:GROUP3D230_INDEX_CASES(69, (6*NXYZ8 + z), (-1 + 6*NXYZ8 - y), (2*NXYZ8 + x));
                                            case 70:GROUP3D230_INDEX_CASES(70, (2*NXYZ8 + z), (6*NXYZ8 + y), (-1 + 6*NXYZ8 - x));
                                            case 71:GROUP3D230_INDEX_CASES(71, (-1 + 2*NXYZ8 - z), (-1 + 2*NXYZ8 - y), (-1 + 2*NXYZ8 - x));
                                            case 72:GROUP3D230_INDEX_CASES(72, (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z));
                                            case 73:GROUP3D230_INDEX_CASES(73, (x), (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z));
                                            case 74:GROUP3D230_INDEX_CASES(74, (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z));
                                            case 75:GROUP3D230_INDEX_CASES(75, (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y), (4*NXYZ8 + z));
                                            case 76:GROUP3D230_INDEX_CASES(76, (-1 + 4*NXYZ8 - y), (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x));
                                            case 77:GROUP3D230_INDEX_CASES(77, (8*NXYZ8 + y), (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x));
                                            case 78:GROUP3D230_INDEX_CASES(78, (4*NXYZ8 + y), (-1 + 8*NXYZ8 - z), (x));
                                            case 79:GROUP3D230_INDEX_CASES(79, (-1 + 8*NXYZ8 - y), (8*NXYZ8 + z), (4*NXYZ8 + x));
                                            case 80:GROUP3D230_INDEX_CASES(80, (-1 + 4*NXYZ8 - z), (-1 + 4*NXYZ8 - x), (-1 + 4*NXYZ8 - y));
                                            case 81:GROUP3D230_INDEX_CASES(81, (8*NXYZ8 + z), (4*NXYZ8 + x), (-1 + 8*NXYZ8 - y));
                                            case 82:GROUP3D230_INDEX_CASES(82, (4*NXYZ8 + z), (-1 + 8*NXYZ8 - x), (8*NXYZ8 + y));
                                            case 83:GROUP3D230_INDEX_CASES(83, (-1 + 8*NXYZ8 - z), (x), (4*NXYZ8 + y));
                                            case 84:GROUP3D230_INDEX_CASES(84, (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z));
                                            case 85:GROUP3D230_INDEX_CASES(85, (2*NXYZ8 + y), (2*NXYZ8 + x), (2*NXYZ8 + z));
                                            case 86:GROUP3D230_INDEX_CASES(86, (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z));
                                            case 87:GROUP3D230_INDEX_CASES(87, (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z));
                                            case 88:GROUP3D230_INDEX_CASES(88, (-1 + 6*NXYZ8 - x), (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y));
                                            case 89:GROUP3D230_INDEX_CASES(89, (6*NXYZ8 + x), (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y));
                                            case 90:GROUP3D230_INDEX_CASES(90, (2*NXYZ8 + x), (2*NXYZ8 + z), (2*NXYZ8 + y));
                                            case 91:GROUP3D230_INDEX_CASES(91, (-1 + 2*NXYZ8 - x), (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y));
                                            case 92:GROUP3D230_INDEX_CASES(92, (6*NXYZ8 + z), (-1 + 2*NXYZ8 - y), (-1 + 6*NXYZ8 - x));
                                            case 93:GROUP3D230_INDEX_CASES(93, (-1 + 6*NXYZ8 - z), (6*NXYZ8 + y), (-1 + 2*NXYZ8 - x));
                                            case 94:GROUP3D230_INDEX_CASES(94, (-1 + 2*NXYZ8 - z), (-1 + 6*NXYZ8 - y), (6*NXYZ8 + x));
                                            case 95:GROUP3D230_INDEX_CASES(95, (2*NXYZ8 + z), (2*NXYZ8 + y), (2*NXYZ8 + x));
                                            default: return;
                                            }
                                        }

    #ifdef GROUP3D230_INDEX_CASES
    #undef GROUP3D230_INDEX_CASES
    #endif // GROUP3D230_INDEX_CASES
                                        if (element<NXYZ>(x, y, z)) {
                                            dst_real[real_index] = src_lebedev[lebdev_pos * NAsymunit + Asymunit_index];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /// <summary>
            /// SignedT must be a type of Signed interger!!!
            /// </summary>
            /// <typeparam name="SignedT">  must be a type of Signed interger </typeparam>
            template<size_t NXYZ, typename SignedT = int>
            struct GroupInfo {
                using KernelFuncT = void (*)(double*, const double*);

                constexpr static SignedT operation_number_ = detail::operation_number<SignedT>;
                constexpr static c_array<SignedT, 6> fields_required_table_ = detail::fields_required_table<SignedT>;
                constexpr static c_array<SignedT, 6> extract_num_table_ = detail::extract_num_table<SignedT>;
                constexpr static c_array<SignedT, 6> reconstruct_num_table_ = detail::reconstruct_num_table<SignedT>;
                constexpr static SignedT max_fields_required_ = detail::maximum(fields_required_table_);
                constexpr static SignedT N_x_ = NXYZ;
                constexpr static SignedT N_y_ = NXYZ;
                constexpr static SignedT N_z_ = NXYZ;
                constexpr static SignedT length_x_ = detail::lengthx<N_x_, SignedT>;
                constexpr static SignedT length_y_ = detail::lengthy<N_y_, SignedT>;
                constexpr static SignedT length_z_ = detail::lengthz<N_z_, SignedT>;
                constexpr static SignedT total_elements_ = detail::kernel::total_elements<NXYZ, SignedT>();
                
                template<LEBEDEV_POINT_TYPE PType>
                struct PointInfo {
                    constexpr static SignedT fields_required_ = detail::fields_required<PType, SignedT>;
                    constexpr static auto field_indexes_ = detail::field_indexes<PType, SignedT>;
                    constexpr static auto field_submultiplicities_ = detail::field_submultiplicities<PType, SignedT>;
                    constexpr static SignedT extract_num_ = detail::extract_num<PType, SignedT>;
                    constexpr static SignedT reconstruct_num_ = detail::reconstruct_num<PType, SignedT>;
                    constexpr static auto extract_fptr_ = detail::kernel::extract<NXYZ, PType, SignedT>;
                    constexpr static auto reconstruct_fptr_ = detail::kernel::reconstruct<NXYZ, PType, SignedT>;
                };

                constexpr static c_array<KernelFuncT, 6>  extract_fptr_table_ = {{
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN6_00C>::extract_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN12_0BB>::extract_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN8_AAA>::extract_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AAC>::extract_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AB0>::extract_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN48_ABC>::extract_fptr_
                }};
                constexpr static c_array<KernelFuncT, 6>  reconstruct_fptr_table_ = {{
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN6_00C>::reconstruct_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN12_0BB>::reconstruct_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN8_AAA>::reconstruct_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AAC>::reconstruct_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AB0>::reconstruct_fptr_,
                    PointInfo<LEBEDEV_POINT_TYPE::OPTRN48_ABC>::reconstruct_fptr_
                }};
                //note that this table is padded to let its elements to have the same length.
                //the valid length table is actually fields_required_table_
                constexpr static c_array<c_array<SignedT, max_fields_required_>, 6>  field_indexes_table_ = {{
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN6_00C>::field_indexes_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN12_0BB>::field_indexes_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN8_AAA>::field_indexes_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AAC>::field_indexes_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AB0>::field_indexes_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN48_ABC>::field_indexes_)
                }};
                //note that this table is padded to let its elements to have the same length.
                //the valid length table is actually fields_required_table_
                constexpr static c_array<c_array<SignedT, max_fields_required_>, 6>  field_submultiplicities_table_ = {{
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN6_00C>::field_submultiplicities_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN12_0BB>::field_submultiplicities_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN8_AAA>::field_submultiplicities_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AAC>::field_submultiplicities_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AB0>::field_submultiplicities_),
                    detail::add_padding<max_fields_required_>(PointInfo<LEBEDEV_POINT_TYPE::OPTRN48_ABC>::field_submultiplicities_)
                }};
            };
        }
    }

};