#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "../lebedev_info.h"

namespace lebedev
{
    namespace group2D009
    {
        using lebedev::c_array;
        using lebedev::LEBEDEV_POINT_TYPE;
        using lebedev::detail::point_type_multiplicity;

        namespace detail
        {
            template <typename SignedT = int>
            constexpr SignedT operation_number = 8;
            template <typename SignedT = int>
            constexpr c_array<SignedT, 6> fields_required_table = {{4, 5, 2, 6, 10, 12}};
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
            constexpr c_array<SignedT, 6> reconstruct_num_table = {{32, 40, 16, 48, 80, 96}};
            template <size_t NX, typename SignedT = int>
            constexpr SignedT lengthx = NX/4;
            template <size_t NY, typename SignedT = int>
            constexpr SignedT lengthy = NY/2;
            template <size_t NZ, typename SignedT = int>
            constexpr SignedT lengthz = NZ;
            template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
            constexpr SignedT fields_required = fields_required_table<SignedT>[static_cast<int>(PType)];
            template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
            constexpr SignedT extract_num = extract_num_table<SignedT>[static_cast<int>(PType)];
            template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
            constexpr SignedT reconstruct_num = reconstruct_num_table<SignedT>[static_cast<int>(PType)];
            template <LEBEDEV_POINT_TYPE PType, typename SignedT = int>
            constexpr c_array<SignedT, fields_required<PType, SignedT>> field_indexes;
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN6_00C, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN6_00C, SignedT> = {{0, 1, 2, 4}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN12_0BB, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN12_0BB, SignedT> = {{0, 2, 4, 5, 8}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN8_AAA, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN8_AAA, SignedT> = {{0, 2}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN24_AAC, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN24_AAC, SignedT> = {{0, 2, 4, 5, 8, 9}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN24_AB0, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN24_AB0, SignedT> = {{0, 4, 5, 8, 9, 12, 16, 17, 20, 21}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN48_ABC, SignedT>> field_indexes<LEBEDEV_POINT_TYPE::OPTRN48_ABC, SignedT> = {{0, 2, 4, 5, 8, 9, 12, 14, 16, 17, 20, 21}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN6_00C, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN6_00C, SignedT> = {{1, 1, 2, 2}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN12_0BB, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN12_0BB, SignedT> = {{2, 2, 2, 2, 4}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN8_AAA, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN8_AAA, SignedT> = {{4, 4}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN24_AAC, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN24_AAC, SignedT> = {{4, 4, 4, 4, 4, 4}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN24_AB0, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN24_AB0, SignedT> = {{4, 2, 2, 2, 2, 4, 2, 2, 2, 2}};
            template <typename SignedT>
            constexpr c_array<SignedT, fields_required<LEBEDEV_POINT_TYPE::OPTRN48_ABC, SignedT>> field_submultiplicities<LEBEDEV_POINT_TYPE::OPTRN48_ABC, SignedT> = {{4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};


            namespace kernel
            {
                /// <summary>
                /// This function determines whether an element is contained in the selected group2D009 asymmetric unit.
                /// The asymmetric unit is selected as 0 <= x <= 1/4 && 0 <= y <= 1/2 && 0 <= z <= 1.
                /// x falls in the range of [0, 1/4].
                /// y falls in the range of [0, 1/2].
                /// z falls in the range of [0, 1].
                /// Note that these restriction is never checked.
                /// </summary>
                template <size_t NX, size_t NY, size_t NZ, typename SignedT = int>
                __host__ __device__ __forceinline__ constexpr bool element(SignedT x, SignedT y, SignedT z)
                {
                    return true;
                }

                /// <summary>
                /// This function calculates the shift of a group2D009 asymmetric unit.
                /// The asymmetric unit is selected as 0 <= x <= 1/4 && 0 <= y <= 1/2 && 0 <= z <= 1.
                /// x falls in the range of [0, 1/4].
                /// y falls in the range of [0, 1/2].
                /// z falls in the range of [0, 1].
                /// Note that these restriction is never checked.
                /// </summary>
                template <size_t NX, size_t NY, size_t NZ, typename SignedT = int>
                __host__ __device__ __forceinline__ constexpr SignedT index_xyz(SignedT x, SignedT y, SignedT z)
                {
                    constexpr SignedT NX4 = (SignedT)(NX / 4);
                    constexpr SignedT NY2 = (SignedT)(NY / 2);
                    return (x * NY2 + y) * NZ + z;
                }

                /// <summary>
                /// This function calculates the number of the sample points in a group2D009 asymmetric unit.
                /// The asymmetric unit is selected as 0 <= x <= 1/4 && 0 <= y <= 1/2 && 0 <= z <= 1.
                /// </summary>
                template <size_t NX, size_t NY, size_t NZ, typename SignedT = int>
                __host__ __device__ __forceinline__ constexpr SignedT total_elements()
                {
                    constexpr SignedT NX4 = (SignedT)(NX / 4);
                    constexpr SignedT NY2 = (SignedT)(NY / 2);
                    return index_xyz<NX, NY, NZ, SignedT>(NX4 - 1, NY2 - 1, NZ - 1 + 1/*one-past-last*/);
                }

                /// <summary>
                /// This function return the relative positions of the symmetry operations in lebedev grids
                /// In detail, the data at lebedev point i in spatial-compressed format, 
                /// is equivalent to the data at extract_field_pos<i,...> in orientational-compressed format
                /// Ideally this should be consteval function which is not currently supported by CUDA 11.x
                /// </summary>
                template <decltype(dim3::x) LebedevID, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
                __host__ __device__ __forceinline__ constexpr SignedT extract_field_pos()
                {
                    if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
						constexpr SignedT list[6] =
                        {
                              0,   1,   2,   2,   3,   3
                        };
						return list[LebedevID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
						constexpr SignedT list[12] =
						{
						      0,   0,   1,   1,   2,   3,   2,   3,   4,   4,   4,   4
						};
						return list[LebedevID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
						constexpr SignedT list[8] =
						{
						      0,   0,   1,   1,   1,   1,   0,   0
						};
						return list[LebedevID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
						constexpr SignedT list[24] =
						{
						      0,   0,   1,   1,   2,   3,   2,   3,   4,   5,   5,   4, 
						      1,   1,   0,   0,   5,   4,   5,   4,   3,   2,   2,   3
						};
						return list[LebedevID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
						constexpr SignedT list[24] =
						{
						      0,   0,   0,   0,   1,   2,   1,   2,   3,   4,   4,   3, 
						      5,   5,   5,   5,   6,   7,   6,   7,   8,   9,   9,   8
						};
						return list[LebedevID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
						constexpr SignedT list[48] = 
                        {
                              0,   0,   1,   1,   2,   3,   2,   3,   4,   5,   5,   4, 
                              6,   6,   7,   7,   8,   9,   8,   9,  10,  11,  11,  10, 
                              1,   1,   0,   0,   3,   2,   3,   2,   5,   4,   4,   5, 
                              7,   7,   6,   6,   9,   8,   9,   8,  11,  10,  10,  11
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
                template <decltype(dim3::x) OptrnID, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
                __host__ __device__ __forceinline__ constexpr SignedT reconstruct_lebedev_point_pos()
                {
                    if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
						constexpr SignedT list[32] =
                        {
                              0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1, 
                              1,   1,   1,   1,   2,   3,   3,   2,   2,   3,   3,   2, 
                              4,   5,   4,   5,   4,   5,   4,   5
                        };
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
						constexpr SignedT list[40] =
						{
						      0,   1,   0,   1,   0,   1,   0,   1,   2,   3,   2,   3, 
						      2,   3,   2,   3,   4,   6,   6,   4,   4,   6,   6,   4, 
						      5,   7,   7,   5,   5,   7,   7,   5,   8,  11,   9,  10, 
						      8,  11,   9,  10
						};
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
						constexpr SignedT list[16] =
						{
						      0,   1,   7,   6,   0,   1,   7,   6,   2,   3,   4,   5, 
						      2,   3,   4,   5
						};
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
						constexpr SignedT list[48] =
						{
						      0,   1,  15,  14,   0,   1,  15,  14,   2,   3,  12,  13, 
						      2,   3,  12,  13,   4,   6,  22,  21,   4,   6,  22,  21, 
						      5,   7,  23,  20,   5,   7,  23,  20,   8,  11,  17,  19, 
						      8,  11,  17,  19,   9,  10,  16,  18,   9,  10,  16,  18
						};
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
						constexpr SignedT list[80] =
						{
						      0,   1,   2,   3,   0,   1,   2,   3,   4,   6,   4,   6, 
						      4,   6,   4,   6,   5,   7,   5,   7,   5,   7,   5,   7, 
						      8,  11,  11,   8,   8,  11,  11,   8,   9,  10,  10,   9, 
						      9,  10,  10,   9,  12,  13,  15,  14,  12,  13,  15,  14, 
						     16,  18,  18,  16,  16,  18,  18,  16,  17,  19,  19,  17, 
						     17,  19,  19,  17,  20,  23,  20,  23,  20,  23,  20,  23, 
						     21,  22,  21,  22,  21,  22,  21,  22
						};
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
						constexpr SignedT list[96] = 
                        {
                              0,   1,  27,  26,   0,   1,  27,  26,   2,   3,  25,  24, 
                              2,   3,  25,  24,   4,   6,  29,  31,   4,   6,  29,  31, 
                              5,   7,  28,  30,   5,   7,  28,  30,   8,  11,  34,  33, 
                              8,  11,  34,  33,   9,  10,  35,  32,   9,  10,  35,  32, 
                             12,  13,  38,  39,  12,  13,  38,  39,  14,  15,  36,  37, 
                             14,  15,  36,  37,  16,  18,  43,  41,  16,  18,  43,  41, 
                             17,  19,  42,  40,  17,  19,  42,  40,  20,  23,  45,  46, 
                             20,  23,  45,  46,  21,  22,  44,  47,  21,  22,  44,  47
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
                template <decltype(dim3::x) OptrnID, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
                __host__ __device__ __forceinline__ constexpr SignedT reconstruct_field_pos()
                {
                    if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
						constexpr SignedT list[32] =
                        {
                              0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1, 
                              1,   1,   1,   1,   2,   2,   2,   2,   2,   2,   2,   2, 
                              3,   3,   3,   3,   3,   3,   3,   3
                        };
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
						constexpr SignedT list[40] =
						{
						      0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1, 
						      1,   1,   1,   1,   2,   2,   2,   2,   2,   2,   2,   2, 
						      3,   3,   3,   3,   3,   3,   3,   3,   4,   4,   4,   4, 
						      4,   4,   4,   4
						};
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
						constexpr SignedT list[16] =
						{
						      0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1, 
						      1,   1,   1,   1
						};
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
						constexpr SignedT list[48] =
						{
						      0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1, 
						      1,   1,   1,   1,   2,   2,   2,   2,   2,   2,   2,   2, 
						      3,   3,   3,   3,   3,   3,   3,   3,   4,   4,   4,   4, 
						      4,   4,   4,   4,   5,   5,   5,   5,   5,   5,   5,   5
						};
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
						constexpr SignedT list[80] =
						{
						      0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1, 
						      1,   1,   1,   1,   2,   2,   2,   2,   2,   2,   2,   2, 
						      3,   3,   3,   3,   3,   3,   3,   3,   4,   4,   4,   4, 
						      4,   4,   4,   4,   5,   5,   5,   5,   5,   5,   5,   5, 
						      6,   6,   6,   6,   6,   6,   6,   6,   7,   7,   7,   7, 
						      7,   7,   7,   7,   8,   8,   8,   8,   8,   8,   8,   8, 
						      9,   9,   9,   9,   9,   9,   9,   9
						};
						return list[OptrnID];
					}
					if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
						constexpr SignedT list[96] = 
                        {
                              0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1, 
                              1,   1,   1,   1,   2,   2,   2,   2,   2,   2,   2,   2, 
                              3,   3,   3,   3,   3,   3,   3,   3,   4,   4,   4,   4, 
                              4,   4,   4,   4,   5,   5,   5,   5,   5,   5,   5,   5, 
                              6,   6,   6,   6,   6,   6,   6,   6,   7,   7,   7,   7, 
                              7,   7,   7,   7,   8,   8,   8,   8,   8,   8,   8,   8, 
                              9,   9,   9,   9,   9,   9,   9,   9,  10,  10,  10,  10, 
                             10,  10,  10,  10,  11,  11,  11,  11,  11,  11,  11,  11
                        };
						return list[OptrnID];
					}
					return -1;
                }

                /// <summary>
				/// Note that this kernel requires that blockIdx.x should be the number of lebedev points.
				/// However, this is not checked in the kernel. Any block exceeds this constraint ***will return without doing anything***
				/// </summary>
				template<size_t NX, size_t NY, size_t NZ, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
				__global__ void extract(double* dst_lebedev, const double* src_real) {
					constexpr SignedT Nx = NX;
					constexpr SignedT Ny = NY;
					constexpr SignedT Nz = NZ;

                    constexpr SignedT NX4 = NX / 4;
                    constexpr SignedT NY2 = NY / 2;
                    constexpr SignedT NZ1 = NZ / 1;
                    static_assert(NX == NX4 * 4, "template parameter NX must be divided by 4");
                    static_assert(NY == NY2 * 2, "template parameter NY must be divided by 2");
                    static_assert(NZ == NZ1 * 1, "template parameter NZ must be divided by 1");

					constexpr SignedT Sx = 0;
					constexpr SignedT Sy = 0;
					constexpr SignedT Sz = 0;

					constexpr SignedT NAsymunit = total_elements<NX, NY, NZ, SignedT>();

                    for(SignedT x = Sx + (SignedT)(blockIdx.y); x < Sx + lengthx<Nx, SignedT>; x += (SignedT)(gridDim.y))
                    {
                        for(SignedT y = Sy + (SignedT)(threadIdx.x); y < Sy + lengthy<Ny, SignedT>; y += (SignedT)(blockDim.x))
                        {
                            for(SignedT z = Sz + (SignedT)(threadIdx.y); z < Sz + lengthz<Nz, SignedT>; z += (SignedT)(blockDim.y))
                            {                               
                                SignedT Asymunit_index = index_xyz<NX, NY, NZ, SignedT>(x, y, z);
					            SignedT real_index = 0;
                                
#define GROUP2D009_INDEX_CASES(LebedevID, x, y, z) \
				                {\
                                    constexpr SignedT field_shift = extract_field_pos<LebedevID, PType, SignedT>() * Nx * Ny * Nz;\
                                    real_index = field_shift + (((x) % Nx) * Ny + ((y) % Ny)) * Nz + ((z) % Nz); break;\
                                }

                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (x), (y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (x), (y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (x), (y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (x), (y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (-1 + 4*NX4 - x), (y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (x), (-1 + 2*NY2 - y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (-1 + 4*NX4 - x), (y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (x), (-1 + 2*NY2 - y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (x), (-1 + 2*NY2 - y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (-1 + 4*NX4 - x), (y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (x), (y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (x), (y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (x), (y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (-1 + 4*NX4 - x), (y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (x), (-1 + 2*NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (x), (-1 + 2*NY2 - y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (-1 + 4*NX4 - x), (y), (z));
                                    case 16:GROUP2D009_INDEX_CASES(16, (-1 + 4*NX4 - x), (y), (z));
                                    case 17:GROUP2D009_INDEX_CASES(17, (-1 + 4*NX4 - x), (y), (z));
                                    case 18:GROUP2D009_INDEX_CASES(18, (x), (-1 + 2*NY2 - y), (z));
                                    case 19:GROUP2D009_INDEX_CASES(19, (x), (-1 + 2*NY2 - y), (z));
                                    case 20:GROUP2D009_INDEX_CASES(20, (x), (-1 + 2*NY2 - y), (z));
                                    case 21:GROUP2D009_INDEX_CASES(21, (x), (-1 + 2*NY2 - y), (z));
                                    case 22:GROUP2D009_INDEX_CASES(22, (-1 + 4*NX4 - x), (y), (z));
                                    case 23:GROUP2D009_INDEX_CASES(23, (-1 + 4*NX4 - x), (y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (-1 + 4*NX4 - x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (x), (y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (x), (y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (x), (y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (x), (y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (x), (-1 + 2*NY2 - y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (-1 + 4*NX4 - x), (y), (z));
                                    case 16:GROUP2D009_INDEX_CASES(16, (x), (y), (z));
                                    case 17:GROUP2D009_INDEX_CASES(17, (x), (y), (z));
                                    case 18:GROUP2D009_INDEX_CASES(18, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 19:GROUP2D009_INDEX_CASES(19, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 20:GROUP2D009_INDEX_CASES(20, (x), (y), (z));
                                    case 21:GROUP2D009_INDEX_CASES(21, (x), (y), (z));
                                    case 22:GROUP2D009_INDEX_CASES(22, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 23:GROUP2D009_INDEX_CASES(23, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (x), (y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (x), (y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (x), (y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (x), (y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (x), (y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 16:GROUP2D009_INDEX_CASES(16, (x), (y), (z));
                                    case 17:GROUP2D009_INDEX_CASES(17, (x), (y), (z));
                                    case 18:GROUP2D009_INDEX_CASES(18, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 19:GROUP2D009_INDEX_CASES(19, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 20:GROUP2D009_INDEX_CASES(20, (x), (y), (z));
                                    case 21:GROUP2D009_INDEX_CASES(21, (x), (y), (z));
                                    case 22:GROUP2D009_INDEX_CASES(22, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 23:GROUP2D009_INDEX_CASES(23, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 24:GROUP2D009_INDEX_CASES(24, (x), (-1 + 2*NY2 - y), (z));
                                    case 25:GROUP2D009_INDEX_CASES(25, (-1 + 4*NX4 - x), (y), (z));
                                    case 26:GROUP2D009_INDEX_CASES(26, (x), (-1 + 2*NY2 - y), (z));
                                    case 27:GROUP2D009_INDEX_CASES(27, (-1 + 4*NX4 - x), (y), (z));
                                    case 28:GROUP2D009_INDEX_CASES(28, (-1 + 4*NX4 - x), (y), (z));
                                    case 29:GROUP2D009_INDEX_CASES(29, (-1 + 4*NX4 - x), (y), (z));
                                    case 30:GROUP2D009_INDEX_CASES(30, (x), (-1 + 2*NY2 - y), (z));
                                    case 31:GROUP2D009_INDEX_CASES(31, (x), (-1 + 2*NY2 - y), (z));
                                    case 32:GROUP2D009_INDEX_CASES(32, (x), (-1 + 2*NY2 - y), (z));
                                    case 33:GROUP2D009_INDEX_CASES(33, (x), (-1 + 2*NY2 - y), (z));
                                    case 34:GROUP2D009_INDEX_CASES(34, (-1 + 4*NX4 - x), (y), (z));
                                    case 35:GROUP2D009_INDEX_CASES(35, (-1 + 4*NX4 - x), (y), (z));
                                    case 36:GROUP2D009_INDEX_CASES(36, (-1 + 4*NX4 - x), (y), (z));
                                    case 37:GROUP2D009_INDEX_CASES(37, (x), (-1 + 2*NY2 - y), (z));
                                    case 38:GROUP2D009_INDEX_CASES(38, (-1 + 4*NX4 - x), (y), (z));
                                    case 39:GROUP2D009_INDEX_CASES(39, (x), (-1 + 2*NY2 - y), (z));
                                    case 40:GROUP2D009_INDEX_CASES(40, (x), (-1 + 2*NY2 - y), (z));
                                    case 41:GROUP2D009_INDEX_CASES(41, (x), (-1 + 2*NY2 - y), (z));
                                    case 42:GROUP2D009_INDEX_CASES(42, (-1 + 4*NX4 - x), (y), (z));
                                    case 43:GROUP2D009_INDEX_CASES(43, (-1 + 4*NX4 - x), (y), (z));
                                    case 44:GROUP2D009_INDEX_CASES(44, (-1 + 4*NX4 - x), (y), (z));
                                    case 45:GROUP2D009_INDEX_CASES(45, (-1 + 4*NX4 - x), (y), (z));
                                    case 46:GROUP2D009_INDEX_CASES(46, (x), (-1 + 2*NY2 - y), (z));
                                    case 47:GROUP2D009_INDEX_CASES(47, (x), (-1 + 2*NY2 - y), (z));
                                    default: return;
                                    }
                                }
#ifdef GROUP2D009_INDEX_CASES
#undef GROUP2D009_INDEX_CASES
#endif // GROUP2D009_INDEX_CASES
                                if (element<NX, NY, NZ>(x, y, z)) {
						            dst_lebedev[blockIdx.x * NAsymunit + Asymunit_index] = src_real[real_index];
                                }
                            }
                        }
                    }

				}
				
                /// <summary>
				/// Note that this kernel requires that blockIdx.x should be the number of operations multiples by number of leading lebedev points.
				/// However, this is not checked in the kernel. Any block exceeds this constraint ***will return without doing anything***
				/// </summary>
				template<size_t NX, size_t NY, size_t NZ, LEBEDEV_POINT_TYPE PType = LEBEDEV_POINT_TYPE::OPTRN0_EMPTY, typename SignedT = int>
				__global__ void reconstruct(double* dst_real, const double* src_lebedev) {
					constexpr SignedT Nx = NX;
					constexpr SignedT Ny = NY;
					constexpr SignedT Nz = NZ;

                    constexpr SignedT NX4 = NX / 4;
                    constexpr SignedT NY2 = NY / 2;
                    constexpr SignedT NZ1 = NZ / 1;
                    static_assert(NX == NX4 * 4, "template parameter NX must be divided by 4");
                    static_assert(NY == NY2 * 2, "template parameter NY must be divided by 2");
                    static_assert(NZ == NZ1 * 1, "template parameter NZ must be divided by 1");

					constexpr SignedT Sx = 0;
					constexpr SignedT Sy = 0;
					constexpr SignedT Sz = 0;

					constexpr SignedT NAsymunit = total_elements<NX, NY, NZ, SignedT>();

                    for(SignedT x = Sx + (SignedT)(blockIdx.y); x < Sx + lengthx<Nx, SignedT>; x += (SignedT)(gridDim.y))
                    {
                        for(SignedT y = Sy + (SignedT)(threadIdx.x); y < Sy + lengthy<Ny, SignedT>; y += (SignedT)(blockDim.x))
                        {
                            for(SignedT z = Sz + (SignedT)(threadIdx.y); z < Sz + lengthz<Nz, SignedT>; z += (SignedT)(blockDim.y))
                            {                               
                                SignedT Asymunit_index = index_xyz<NX, NY, NZ, SignedT>(x, y, z);
                                SignedT lebdev_pos = 0;
					            SignedT real_index = 0;

#define GROUP2D009_INDEX_CASES(OptrnID, x, y, z) \
				                {\
                                    constexpr SignedT temppos = reconstruct_lebedev_point_pos<OptrnID, PType, SignedT>();\
                                    lebdev_pos = temppos;\
                                    constexpr SignedT field_shift = reconstruct_field_pos<OptrnID, PType, SignedT>() * Nx * Ny * Nz;\
                                    real_index = field_shift + (((x) % Nx) * Ny + ((y) % Ny)) * Nz + ((z) % Nz); break;\
                                }

                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN6_00C) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (-1 + 4*NX4 - x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (2*NX4 + x), (NY2 + y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (2*NX4 + x), (NY2 + y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 16:GROUP2D009_INDEX_CASES(16, (x), (y), (z));
                                    case 17:GROUP2D009_INDEX_CASES(17, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 18:GROUP2D009_INDEX_CASES(18, (-1 + 4*NX4 - x), (y), (z));
                                    case 19:GROUP2D009_INDEX_CASES(19, (x), (-1 + 2*NY2 - y), (z));
                                    case 20:GROUP2D009_INDEX_CASES(20, (2*NX4 + x), (NY2 + y), (z));
                                    case 21:GROUP2D009_INDEX_CASES(21, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 22:GROUP2D009_INDEX_CASES(22, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 23:GROUP2D009_INDEX_CASES(23, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 24:GROUP2D009_INDEX_CASES(24, (x), (y), (z));
                                    case 25:GROUP2D009_INDEX_CASES(25, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 26:GROUP2D009_INDEX_CASES(26, (-1 + 4*NX4 - x), (y), (z));
                                    case 27:GROUP2D009_INDEX_CASES(27, (x), (-1 + 2*NY2 - y), (z));
                                    case 28:GROUP2D009_INDEX_CASES(28, (2*NX4 + x), (NY2 + y), (z));
                                    case 29:GROUP2D009_INDEX_CASES(29, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 30:GROUP2D009_INDEX_CASES(30, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 31:GROUP2D009_INDEX_CASES(31, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN12_0BB) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (-1 + 4*NX4 - x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (2*NX4 + x), (NY2 + y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (2*NX4 + x), (NY2 + y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 16:GROUP2D009_INDEX_CASES(16, (x), (y), (z));
                                    case 17:GROUP2D009_INDEX_CASES(17, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 18:GROUP2D009_INDEX_CASES(18, (-1 + 4*NX4 - x), (y), (z));
                                    case 19:GROUP2D009_INDEX_CASES(19, (x), (-1 + 2*NY2 - y), (z));
                                    case 20:GROUP2D009_INDEX_CASES(20, (2*NX4 + x), (NY2 + y), (z));
                                    case 21:GROUP2D009_INDEX_CASES(21, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 22:GROUP2D009_INDEX_CASES(22, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 23:GROUP2D009_INDEX_CASES(23, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 24:GROUP2D009_INDEX_CASES(24, (x), (y), (z));
                                    case 25:GROUP2D009_INDEX_CASES(25, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 26:GROUP2D009_INDEX_CASES(26, (-1 + 4*NX4 - x), (y), (z));
                                    case 27:GROUP2D009_INDEX_CASES(27, (x), (-1 + 2*NY2 - y), (z));
                                    case 28:GROUP2D009_INDEX_CASES(28, (2*NX4 + x), (NY2 + y), (z));
                                    case 29:GROUP2D009_INDEX_CASES(29, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 30:GROUP2D009_INDEX_CASES(30, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 31:GROUP2D009_INDEX_CASES(31, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 32:GROUP2D009_INDEX_CASES(32, (x), (y), (z));
                                    case 33:GROUP2D009_INDEX_CASES(33, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 34:GROUP2D009_INDEX_CASES(34, (-1 + 4*NX4 - x), (y), (z));
                                    case 35:GROUP2D009_INDEX_CASES(35, (x), (-1 + 2*NY2 - y), (z));
                                    case 36:GROUP2D009_INDEX_CASES(36, (2*NX4 + x), (NY2 + y), (z));
                                    case 37:GROUP2D009_INDEX_CASES(37, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 38:GROUP2D009_INDEX_CASES(38, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 39:GROUP2D009_INDEX_CASES(39, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN8_AAA) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (-1 + 4*NX4 - x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (2*NX4 + x), (NY2 + y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (2*NX4 + x), (NY2 + y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AAC) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (-1 + 4*NX4 - x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (2*NX4 + x), (NY2 + y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (2*NX4 + x), (NY2 + y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 16:GROUP2D009_INDEX_CASES(16, (x), (y), (z));
                                    case 17:GROUP2D009_INDEX_CASES(17, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 18:GROUP2D009_INDEX_CASES(18, (-1 + 4*NX4 - x), (y), (z));
                                    case 19:GROUP2D009_INDEX_CASES(19, (x), (-1 + 2*NY2 - y), (z));
                                    case 20:GROUP2D009_INDEX_CASES(20, (2*NX4 + x), (NY2 + y), (z));
                                    case 21:GROUP2D009_INDEX_CASES(21, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 22:GROUP2D009_INDEX_CASES(22, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 23:GROUP2D009_INDEX_CASES(23, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 24:GROUP2D009_INDEX_CASES(24, (x), (y), (z));
                                    case 25:GROUP2D009_INDEX_CASES(25, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 26:GROUP2D009_INDEX_CASES(26, (-1 + 4*NX4 - x), (y), (z));
                                    case 27:GROUP2D009_INDEX_CASES(27, (x), (-1 + 2*NY2 - y), (z));
                                    case 28:GROUP2D009_INDEX_CASES(28, (2*NX4 + x), (NY2 + y), (z));
                                    case 29:GROUP2D009_INDEX_CASES(29, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 30:GROUP2D009_INDEX_CASES(30, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 31:GROUP2D009_INDEX_CASES(31, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 32:GROUP2D009_INDEX_CASES(32, (x), (y), (z));
                                    case 33:GROUP2D009_INDEX_CASES(33, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 34:GROUP2D009_INDEX_CASES(34, (-1 + 4*NX4 - x), (y), (z));
                                    case 35:GROUP2D009_INDEX_CASES(35, (x), (-1 + 2*NY2 - y), (z));
                                    case 36:GROUP2D009_INDEX_CASES(36, (2*NX4 + x), (NY2 + y), (z));
                                    case 37:GROUP2D009_INDEX_CASES(37, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 38:GROUP2D009_INDEX_CASES(38, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 39:GROUP2D009_INDEX_CASES(39, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 40:GROUP2D009_INDEX_CASES(40, (x), (y), (z));
                                    case 41:GROUP2D009_INDEX_CASES(41, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 42:GROUP2D009_INDEX_CASES(42, (-1 + 4*NX4 - x), (y), (z));
                                    case 43:GROUP2D009_INDEX_CASES(43, (x), (-1 + 2*NY2 - y), (z));
                                    case 44:GROUP2D009_INDEX_CASES(44, (2*NX4 + x), (NY2 + y), (z));
                                    case 45:GROUP2D009_INDEX_CASES(45, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 46:GROUP2D009_INDEX_CASES(46, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 47:GROUP2D009_INDEX_CASES(47, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN24_AB0) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (-1 + 4*NX4 - x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (2*NX4 + x), (NY2 + y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (2*NX4 + x), (NY2 + y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 16:GROUP2D009_INDEX_CASES(16, (x), (y), (z));
                                    case 17:GROUP2D009_INDEX_CASES(17, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 18:GROUP2D009_INDEX_CASES(18, (-1 + 4*NX4 - x), (y), (z));
                                    case 19:GROUP2D009_INDEX_CASES(19, (x), (-1 + 2*NY2 - y), (z));
                                    case 20:GROUP2D009_INDEX_CASES(20, (2*NX4 + x), (NY2 + y), (z));
                                    case 21:GROUP2D009_INDEX_CASES(21, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 22:GROUP2D009_INDEX_CASES(22, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 23:GROUP2D009_INDEX_CASES(23, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 24:GROUP2D009_INDEX_CASES(24, (x), (y), (z));
                                    case 25:GROUP2D009_INDEX_CASES(25, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 26:GROUP2D009_INDEX_CASES(26, (-1 + 4*NX4 - x), (y), (z));
                                    case 27:GROUP2D009_INDEX_CASES(27, (x), (-1 + 2*NY2 - y), (z));
                                    case 28:GROUP2D009_INDEX_CASES(28, (2*NX4 + x), (NY2 + y), (z));
                                    case 29:GROUP2D009_INDEX_CASES(29, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 30:GROUP2D009_INDEX_CASES(30, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 31:GROUP2D009_INDEX_CASES(31, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 32:GROUP2D009_INDEX_CASES(32, (x), (y), (z));
                                    case 33:GROUP2D009_INDEX_CASES(33, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 34:GROUP2D009_INDEX_CASES(34, (-1 + 4*NX4 - x), (y), (z));
                                    case 35:GROUP2D009_INDEX_CASES(35, (x), (-1 + 2*NY2 - y), (z));
                                    case 36:GROUP2D009_INDEX_CASES(36, (2*NX4 + x), (NY2 + y), (z));
                                    case 37:GROUP2D009_INDEX_CASES(37, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 38:GROUP2D009_INDEX_CASES(38, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 39:GROUP2D009_INDEX_CASES(39, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 40:GROUP2D009_INDEX_CASES(40, (x), (y), (z));
                                    case 41:GROUP2D009_INDEX_CASES(41, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 42:GROUP2D009_INDEX_CASES(42, (-1 + 4*NX4 - x), (y), (z));
                                    case 43:GROUP2D009_INDEX_CASES(43, (x), (-1 + 2*NY2 - y), (z));
                                    case 44:GROUP2D009_INDEX_CASES(44, (2*NX4 + x), (NY2 + y), (z));
                                    case 45:GROUP2D009_INDEX_CASES(45, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 46:GROUP2D009_INDEX_CASES(46, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 47:GROUP2D009_INDEX_CASES(47, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 48:GROUP2D009_INDEX_CASES(48, (x), (y), (z));
                                    case 49:GROUP2D009_INDEX_CASES(49, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 50:GROUP2D009_INDEX_CASES(50, (-1 + 4*NX4 - x), (y), (z));
                                    case 51:GROUP2D009_INDEX_CASES(51, (x), (-1 + 2*NY2 - y), (z));
                                    case 52:GROUP2D009_INDEX_CASES(52, (2*NX4 + x), (NY2 + y), (z));
                                    case 53:GROUP2D009_INDEX_CASES(53, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 54:GROUP2D009_INDEX_CASES(54, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 55:GROUP2D009_INDEX_CASES(55, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 56:GROUP2D009_INDEX_CASES(56, (x), (y), (z));
                                    case 57:GROUP2D009_INDEX_CASES(57, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 58:GROUP2D009_INDEX_CASES(58, (-1 + 4*NX4 - x), (y), (z));
                                    case 59:GROUP2D009_INDEX_CASES(59, (x), (-1 + 2*NY2 - y), (z));
                                    case 60:GROUP2D009_INDEX_CASES(60, (2*NX4 + x), (NY2 + y), (z));
                                    case 61:GROUP2D009_INDEX_CASES(61, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 62:GROUP2D009_INDEX_CASES(62, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 63:GROUP2D009_INDEX_CASES(63, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 64:GROUP2D009_INDEX_CASES(64, (x), (y), (z));
                                    case 65:GROUP2D009_INDEX_CASES(65, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 66:GROUP2D009_INDEX_CASES(66, (-1 + 4*NX4 - x), (y), (z));
                                    case 67:GROUP2D009_INDEX_CASES(67, (x), (-1 + 2*NY2 - y), (z));
                                    case 68:GROUP2D009_INDEX_CASES(68, (2*NX4 + x), (NY2 + y), (z));
                                    case 69:GROUP2D009_INDEX_CASES(69, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 70:GROUP2D009_INDEX_CASES(70, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 71:GROUP2D009_INDEX_CASES(71, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 72:GROUP2D009_INDEX_CASES(72, (x), (y), (z));
                                    case 73:GROUP2D009_INDEX_CASES(73, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 74:GROUP2D009_INDEX_CASES(74, (-1 + 4*NX4 - x), (y), (z));
                                    case 75:GROUP2D009_INDEX_CASES(75, (x), (-1 + 2*NY2 - y), (z));
                                    case 76:GROUP2D009_INDEX_CASES(76, (2*NX4 + x), (NY2 + y), (z));
                                    case 77:GROUP2D009_INDEX_CASES(77, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 78:GROUP2D009_INDEX_CASES(78, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 79:GROUP2D009_INDEX_CASES(79, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    default: return;
                                    }
                                }
                                if constexpr (PType == LEBEDEV_POINT_TYPE::OPTRN48_ABC) {
                                    switch (blockIdx.x)
                                    {
                                    case 0:GROUP2D009_INDEX_CASES(0, (x), (y), (z));
                                    case 1:GROUP2D009_INDEX_CASES(1, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 2:GROUP2D009_INDEX_CASES(2, (-1 + 4*NX4 - x), (y), (z));
                                    case 3:GROUP2D009_INDEX_CASES(3, (x), (-1 + 2*NY2 - y), (z));
                                    case 4:GROUP2D009_INDEX_CASES(4, (2*NX4 + x), (NY2 + y), (z));
                                    case 5:GROUP2D009_INDEX_CASES(5, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 6:GROUP2D009_INDEX_CASES(6, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 7:GROUP2D009_INDEX_CASES(7, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 8:GROUP2D009_INDEX_CASES(8, (x), (y), (z));
                                    case 9:GROUP2D009_INDEX_CASES(9, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 10:GROUP2D009_INDEX_CASES(10, (-1 + 4*NX4 - x), (y), (z));
                                    case 11:GROUP2D009_INDEX_CASES(11, (x), (-1 + 2*NY2 - y), (z));
                                    case 12:GROUP2D009_INDEX_CASES(12, (2*NX4 + x), (NY2 + y), (z));
                                    case 13:GROUP2D009_INDEX_CASES(13, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 14:GROUP2D009_INDEX_CASES(14, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 15:GROUP2D009_INDEX_CASES(15, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 16:GROUP2D009_INDEX_CASES(16, (x), (y), (z));
                                    case 17:GROUP2D009_INDEX_CASES(17, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 18:GROUP2D009_INDEX_CASES(18, (-1 + 4*NX4 - x), (y), (z));
                                    case 19:GROUP2D009_INDEX_CASES(19, (x), (-1 + 2*NY2 - y), (z));
                                    case 20:GROUP2D009_INDEX_CASES(20, (2*NX4 + x), (NY2 + y), (z));
                                    case 21:GROUP2D009_INDEX_CASES(21, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 22:GROUP2D009_INDEX_CASES(22, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 23:GROUP2D009_INDEX_CASES(23, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 24:GROUP2D009_INDEX_CASES(24, (x), (y), (z));
                                    case 25:GROUP2D009_INDEX_CASES(25, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 26:GROUP2D009_INDEX_CASES(26, (-1 + 4*NX4 - x), (y), (z));
                                    case 27:GROUP2D009_INDEX_CASES(27, (x), (-1 + 2*NY2 - y), (z));
                                    case 28:GROUP2D009_INDEX_CASES(28, (2*NX4 + x), (NY2 + y), (z));
                                    case 29:GROUP2D009_INDEX_CASES(29, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 30:GROUP2D009_INDEX_CASES(30, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 31:GROUP2D009_INDEX_CASES(31, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 32:GROUP2D009_INDEX_CASES(32, (x), (y), (z));
                                    case 33:GROUP2D009_INDEX_CASES(33, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 34:GROUP2D009_INDEX_CASES(34, (-1 + 4*NX4 - x), (y), (z));
                                    case 35:GROUP2D009_INDEX_CASES(35, (x), (-1 + 2*NY2 - y), (z));
                                    case 36:GROUP2D009_INDEX_CASES(36, (2*NX4 + x), (NY2 + y), (z));
                                    case 37:GROUP2D009_INDEX_CASES(37, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 38:GROUP2D009_INDEX_CASES(38, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 39:GROUP2D009_INDEX_CASES(39, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 40:GROUP2D009_INDEX_CASES(40, (x), (y), (z));
                                    case 41:GROUP2D009_INDEX_CASES(41, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 42:GROUP2D009_INDEX_CASES(42, (-1 + 4*NX4 - x), (y), (z));
                                    case 43:GROUP2D009_INDEX_CASES(43, (x), (-1 + 2*NY2 - y), (z));
                                    case 44:GROUP2D009_INDEX_CASES(44, (2*NX4 + x), (NY2 + y), (z));
                                    case 45:GROUP2D009_INDEX_CASES(45, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 46:GROUP2D009_INDEX_CASES(46, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 47:GROUP2D009_INDEX_CASES(47, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 48:GROUP2D009_INDEX_CASES(48, (x), (y), (z));
                                    case 49:GROUP2D009_INDEX_CASES(49, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 50:GROUP2D009_INDEX_CASES(50, (-1 + 4*NX4 - x), (y), (z));
                                    case 51:GROUP2D009_INDEX_CASES(51, (x), (-1 + 2*NY2 - y), (z));
                                    case 52:GROUP2D009_INDEX_CASES(52, (2*NX4 + x), (NY2 + y), (z));
                                    case 53:GROUP2D009_INDEX_CASES(53, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 54:GROUP2D009_INDEX_CASES(54, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 55:GROUP2D009_INDEX_CASES(55, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 56:GROUP2D009_INDEX_CASES(56, (x), (y), (z));
                                    case 57:GROUP2D009_INDEX_CASES(57, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 58:GROUP2D009_INDEX_CASES(58, (-1 + 4*NX4 - x), (y), (z));
                                    case 59:GROUP2D009_INDEX_CASES(59, (x), (-1 + 2*NY2 - y), (z));
                                    case 60:GROUP2D009_INDEX_CASES(60, (2*NX4 + x), (NY2 + y), (z));
                                    case 61:GROUP2D009_INDEX_CASES(61, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 62:GROUP2D009_INDEX_CASES(62, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 63:GROUP2D009_INDEX_CASES(63, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 64:GROUP2D009_INDEX_CASES(64, (x), (y), (z));
                                    case 65:GROUP2D009_INDEX_CASES(65, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 66:GROUP2D009_INDEX_CASES(66, (-1 + 4*NX4 - x), (y), (z));
                                    case 67:GROUP2D009_INDEX_CASES(67, (x), (-1 + 2*NY2 - y), (z));
                                    case 68:GROUP2D009_INDEX_CASES(68, (2*NX4 + x), (NY2 + y), (z));
                                    case 69:GROUP2D009_INDEX_CASES(69, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 70:GROUP2D009_INDEX_CASES(70, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 71:GROUP2D009_INDEX_CASES(71, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 72:GROUP2D009_INDEX_CASES(72, (x), (y), (z));
                                    case 73:GROUP2D009_INDEX_CASES(73, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 74:GROUP2D009_INDEX_CASES(74, (-1 + 4*NX4 - x), (y), (z));
                                    case 75:GROUP2D009_INDEX_CASES(75, (x), (-1 + 2*NY2 - y), (z));
                                    case 76:GROUP2D009_INDEX_CASES(76, (2*NX4 + x), (NY2 + y), (z));
                                    case 77:GROUP2D009_INDEX_CASES(77, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 78:GROUP2D009_INDEX_CASES(78, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 79:GROUP2D009_INDEX_CASES(79, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 80:GROUP2D009_INDEX_CASES(80, (x), (y), (z));
                                    case 81:GROUP2D009_INDEX_CASES(81, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 82:GROUP2D009_INDEX_CASES(82, (-1 + 4*NX4 - x), (y), (z));
                                    case 83:GROUP2D009_INDEX_CASES(83, (x), (-1 + 2*NY2 - y), (z));
                                    case 84:GROUP2D009_INDEX_CASES(84, (2*NX4 + x), (NY2 + y), (z));
                                    case 85:GROUP2D009_INDEX_CASES(85, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 86:GROUP2D009_INDEX_CASES(86, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 87:GROUP2D009_INDEX_CASES(87, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    case 88:GROUP2D009_INDEX_CASES(88, (x), (y), (z));
                                    case 89:GROUP2D009_INDEX_CASES(89, (-1 + 4*NX4 - x), (-1 + 2*NY2 - y), (z));
                                    case 90:GROUP2D009_INDEX_CASES(90, (-1 + 4*NX4 - x), (y), (z));
                                    case 91:GROUP2D009_INDEX_CASES(91, (x), (-1 + 2*NY2 - y), (z));
                                    case 92:GROUP2D009_INDEX_CASES(92, (2*NX4 + x), (NY2 + y), (z));
                                    case 93:GROUP2D009_INDEX_CASES(93, (-1 + 2*NX4 - x), (-1 + NY2 - y), (z));
                                    case 94:GROUP2D009_INDEX_CASES(94, (-1 + 2*NX4 - x), (NY2 + y), (z));
                                    case 95:GROUP2D009_INDEX_CASES(95, (2*NX4 + x), (-1 + NY2 - y), (z));
                                    default: return;
                                    }
                                }

#ifdef GROUP2D009_INDEX_CASES
#undef GROUP2D009_INDEX_CASES
#endif // GROUP2D009_INDEX_CASES
                                if (element<NX, NY, NZ>(x, y, z)) {
						            dst_real[real_index] = src_lebedev[lebdev_pos * NAsymunit + Asymunit_index];
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
		template<size_t NX, size_t NY, size_t NZ, typename SignedT = int>
		struct GroupInfo {
			using KernelFuncT = void (*)(double*, const double*);

			constexpr static SignedT operation_number_ = detail::operation_number<SignedT>;
			constexpr static c_array<SignedT, 6> fields_required_table_ = detail::fields_required_table<SignedT>;
			constexpr static c_array<SignedT, 6> extract_num_table_ = detail::extract_num_table<SignedT>;
			constexpr static c_array<SignedT, 6> reconstruct_num_table_ = detail::reconstruct_num_table<SignedT>;
			constexpr static SignedT N_x_ = NX;
			constexpr static SignedT N_y_ = NY;
			constexpr static SignedT N_z_ = NZ;
			constexpr static SignedT length_x_ = detail::lengthx<N_x_, SignedT>;
			constexpr static SignedT length_y_ = detail::lengthy<N_y_, SignedT>;
			constexpr static SignedT length_z_ = detail::lengthz<N_z_, SignedT>;
			constexpr static SignedT total_elements_ = detail::kernel::total_elements<NX, NY, NZ, SignedT>();
			
			template<LEBEDEV_POINT_TYPE PType>
			struct PointInfo {
				constexpr static SignedT fields_required_ = detail::fields_required<PType, SignedT>;
				constexpr static auto field_indexes_ = detail::field_indexes<PType, SignedT>;
				constexpr static auto field_submultiplicities_ = detail::field_submultiplicities<PType, SignedT>;
				constexpr static SignedT extract_num_ = detail::extract_num<PType, SignedT>;
				constexpr static SignedT reconstruct_num_ = detail::reconstruct_num<PType, SignedT>;
				constexpr static auto extract_fptr = detail::kernel::extract<NX, NY, NZ, PType, SignedT>;
				constexpr static auto reconstruct_fptr = detail::kernel::reconstruct<NX, NY, NZ, PType, SignedT>;
			};

			constexpr static c_array<KernelFuncT, 6>  extract_fptr_table_ = {
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN6_00C>::extract_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN12_0BB>::extract_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN8_AAA>::extract_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AAC>::extract_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AB0>::extract_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN48_ABC>::extract_fptr
			};
			constexpr static c_array<KernelFuncT, 6>  reconstruct_fptr_table_ = {
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN6_00C>::reconstruct_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN12_0BB>::reconstruct_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN8_AAA>::reconstruct_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AAC>::reconstruct_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN24_AB0>::reconstruct_fptr,
				PointInfo<LEBEDEV_POINT_TYPE::OPTRN48_ABC>::reconstruct_fptr
			};
		};
    }

};