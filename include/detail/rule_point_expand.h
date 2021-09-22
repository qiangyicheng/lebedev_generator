#pragma once

#include "qutility/c_array.h"

#include "point_type.h"

namespace lebedev{
    namespace detail{
        using qutility::c_array::c_array;

        template<LEBEDEV_POINT_TYPE type> 
        inline constexpr c_array<LEBEDEV_POINT_TYPE, point_type_multiplicity(type)> type_expand() {
            constexpr size_t size = point_type_multiplicity(type);
            c_array<LEBEDEV_POINT_TYPE, size> ans;

            for(size_t itr=0;itr<size;++itr){
                ans[itr]=type;
            }

            return ans;
        }

        template<LEBEDEV_POINT_TYPE type> 
        inline constexpr c_array<double, point_type_multiplicity(type)> weight_expand(double const & w ) {
            constexpr size_t size = point_type_multiplicity(type);
            c_array<double, size> ans;

            for(size_t itr=0;itr<size;++itr){
                ans[itr]=w;
            }

            return ans;        
        }

        template<LEBEDEV_POINT_TYPE type> 
        inline constexpr c_array<c_array<double,3>, point_type_multiplicity(type)> point_expand(c_array<double, 3> const & p ) {return {};}


        template<>
        inline constexpr c_array<c_array<double,3>, 6> point_expand<LEBEDEV_POINT_TYPE::OPTRN6_00C>(c_array<double, 3> const & p ){
            c_array<c_array<double,3>, 6> ans{};

            ans[0][0] = 0;
            ans[0][1] = 0;
            ans[0][2] = p[2];

            ans[1][0] = 0;
            ans[1][1] = 0;
            ans[1][2] = -p[2];

            ans[2][0] = p[2];
            ans[2][1] = 0;
            ans[2][2] = 0;

            ans[3][0] = -p[2];
            ans[3][1] = 0;
            ans[3][2] = 0;

            ans[4][0] = 0;
            ans[4][1] = p[2];
            ans[4][2] = 0;

            ans[5][0] = 0;
            ans[5][1] = -p[2];
            ans[5][2] = 0;

            return ans;
        }


        template<>
        inline constexpr c_array<c_array<double,3>, 12> point_expand<LEBEDEV_POINT_TYPE::OPTRN12_0BB>(c_array<double, 3> const & p ){
            c_array<c_array<double,3>, 12> ans{};

            ans[0][0] = 0;
            ans[0][1] = p[1];
            ans[0][2] = p[1];

            ans[1][0] = 0;
            ans[1][1] = -p[1];
            ans[1][2] = p[1];

            ans[2][0] = 0;
            ans[2][1] = p[1];
            ans[2][2] = -p[1];

            ans[3][0] = 0;
            ans[3][1] = -p[1];
            ans[3][2] = -p[1];

            ans[4][0] = p[1];
            ans[4][1] = 0;
            ans[4][2] = p[1];

            ans[5][0] = p[1];
            ans[5][1] = 0;
            ans[5][2] = -p[1];

            ans[6][0] = -p[1];
            ans[6][1] = 0;
            ans[6][2] = p[1];

            ans[7][0] = -p[1];
            ans[7][1] = 0;
            ans[7][2] = -p[1];

            ans[8][0] = p[1];
            ans[8][1] = p[1];
            ans[8][2] = 0;

            ans[9][0] = -p[1];
            ans[9][1] = p[1];
            ans[9][2] = 0;

            ans[10][0] = p[1];
            ans[10][1] = -p[1];
            ans[10][2] = 0;

            ans[11][0] = -p[1];
            ans[11][1] = -p[1];
            ans[11][2] = 0;

            return ans;
        }


        template<>
        inline constexpr c_array<c_array<double,3>, 8> point_expand<LEBEDEV_POINT_TYPE::OPTRN8_AAA>(c_array<double, 3> const & p ){
            c_array<c_array<double,3>, 8> ans{};

            ans[0][0] = p[0];
            ans[0][1] = p[0];
            ans[0][2] = p[0];

            ans[1][0] = -p[0];
            ans[1][1] = -p[0];
            ans[1][2] = p[0];

            ans[2][0] = -p[0];
            ans[2][1] = p[0];
            ans[2][2] = -p[0];

            ans[3][0] = p[0];
            ans[3][1] = -p[0];
            ans[3][2] = -p[0];

            ans[4][0] = p[0];
            ans[4][1] = p[0];
            ans[4][2] = -p[0];

            ans[5][0] = -p[0];
            ans[5][1] = -p[0];
            ans[5][2] = -p[0];

            ans[6][0] = p[0];
            ans[6][1] = -p[0];
            ans[6][2] = p[0];

            ans[7][0] = -p[0];
            ans[7][1] = p[0];
            ans[7][2] = p[0];

            return ans;
        }


        template<>
        inline constexpr c_array<c_array<double,3>, 24> point_expand<LEBEDEV_POINT_TYPE::OPTRN24_AAC>(c_array<double, 3> const & p ){
            c_array<c_array<double,3>, 24> ans{};

            ans[0][0] = p[0];
            ans[0][1] = p[0];
            ans[0][2] = p[2];

            ans[1][0] = -p[0];
            ans[1][1] = -p[0];
            ans[1][2] = p[2];

            ans[2][0] = -p[0];
            ans[2][1] = p[0];
            ans[2][2] = -p[2];

            ans[3][0] = p[0];
            ans[3][1] = -p[0];
            ans[3][2] = -p[2];

            ans[4][0] = p[2];
            ans[4][1] = p[0];
            ans[4][2] = p[0];

            ans[5][0] = p[2];
            ans[5][1] = -p[0];
            ans[5][2] = -p[0];

            ans[6][0] = -p[2];
            ans[6][1] = -p[0];
            ans[6][2] = p[0];

            ans[7][0] = -p[2];
            ans[7][1] = p[0];
            ans[7][2] = -p[0];

            ans[8][0] = p[0];
            ans[8][1] = p[2];
            ans[8][2] = p[0];

            ans[9][0] = -p[0];
            ans[9][1] = p[2];
            ans[9][2] = -p[0];

            ans[10][0] = p[0];
            ans[10][1] = -p[2];
            ans[10][2] = -p[0];

            ans[11][0] = -p[0];
            ans[11][1] = -p[2];
            ans[11][2] = p[0];

            ans[12][0] = p[0];
            ans[12][1] = p[0];
            ans[12][2] = -p[2];

            ans[13][0] = -p[0];
            ans[13][1] = -p[0];
            ans[13][2] = -p[2];

            ans[14][0] = p[0];
            ans[14][1] = -p[0];
            ans[14][2] = p[2];

            ans[15][0] = -p[0];
            ans[15][1] = p[0];
            ans[15][2] = p[2];

            ans[16][0] = p[0];
            ans[16][1] = p[2];
            ans[16][2] = -p[0];

            ans[17][0] = -p[0];
            ans[17][1] = p[2];
            ans[17][2] = p[0];

            ans[18][0] = -p[0];
            ans[18][1] = -p[2];
            ans[18][2] = -p[0];

            ans[19][0] = p[0];
            ans[19][1] = -p[2];
            ans[19][2] = p[0];

            ans[20][0] = p[2];
            ans[20][1] = p[0];
            ans[20][2] = -p[0];

            ans[21][0] = p[2];
            ans[21][1] = -p[0];
            ans[21][2] = p[0];

            ans[22][0] = -p[2];
            ans[22][1] = p[0];
            ans[22][2] = p[0];

            ans[23][0] = -p[2];
            ans[23][1] = -p[0];
            ans[23][2] = -p[0];

            return ans;
        }


        template<>
        inline constexpr c_array<c_array<double,3>, 24> point_expand<LEBEDEV_POINT_TYPE::OPTRN24_AB0>(c_array<double, 3> const & p ){
            c_array<c_array<double,3>, 24> ans{};

            ans[0][0] = p[0];
            ans[0][1] = p[1];
            ans[0][2] = 0;

            ans[1][0] = -p[0];
            ans[1][1] = -p[1];
            ans[1][2] = 0;

            ans[2][0] = -p[0];
            ans[2][1] = p[1];
            ans[2][2] = 0;

            ans[3][0] = p[0];
            ans[3][1] = -p[1];
            ans[3][2] = 0;

            ans[4][0] = 0;
            ans[4][1] = p[0];
            ans[4][2] = p[1];

            ans[5][0] = 0;
            ans[5][1] = -p[0];
            ans[5][2] = -p[1];

            ans[6][0] = 0;
            ans[6][1] = -p[0];
            ans[6][2] = p[1];

            ans[7][0] = 0;
            ans[7][1] = p[0];
            ans[7][2] = -p[1];

            ans[8][0] = p[1];
            ans[8][1] = 0;
            ans[8][2] = p[0];

            ans[9][0] = -p[1];
            ans[9][1] = 0;
            ans[9][2] = -p[0];

            ans[10][0] = p[1];
            ans[10][1] = 0;
            ans[10][2] = -p[0];

            ans[11][0] = -p[1];
            ans[11][1] = 0;
            ans[11][2] = p[0];

            ans[12][0] = p[1];
            ans[12][1] = p[0];
            ans[12][2] = 0;

            ans[13][0] = -p[1];
            ans[13][1] = -p[0];
            ans[13][2] = 0;

            ans[14][0] = p[1];
            ans[14][1] = -p[0];
            ans[14][2] = 0;

            ans[15][0] = -p[1];
            ans[15][1] = p[0];
            ans[15][2] = 0;

            ans[16][0] = p[0];
            ans[16][1] = 0;
            ans[16][2] = -p[1];

            ans[17][0] = -p[0];
            ans[17][1] = 0;
            ans[17][2] = p[1];

            ans[18][0] = -p[0];
            ans[18][1] = 0;
            ans[18][2] = -p[1];

            ans[19][0] = p[0];
            ans[19][1] = 0;
            ans[19][2] = p[1];

            ans[20][0] = 0;
            ans[20][1] = p[1];
            ans[20][2] = -p[0];

            ans[21][0] = 0;
            ans[21][1] = -p[1];
            ans[21][2] = p[0];

            ans[22][0] = 0;
            ans[22][1] = p[1];
            ans[22][2] = p[0];

            ans[23][0] = 0;
            ans[23][1] = -p[1];
            ans[23][2] = -p[0];

            return ans;
        }


        template<>
        inline constexpr c_array<c_array<double,3>, 48> point_expand<LEBEDEV_POINT_TYPE::OPTRN48_ABC>(c_array<double, 3> const & p ){
            c_array<c_array<double,3>, 48> ans{};

            ans[0][0] = p[0];
            ans[0][1] = p[1];
            ans[0][2] = p[2];

            ans[1][0] = -p[0];
            ans[1][1] = -p[1];
            ans[1][2] = p[2];

            ans[2][0] = -p[0];
            ans[2][1] = p[1];
            ans[2][2] = -p[2];

            ans[3][0] = p[0];
            ans[3][1] = -p[1];
            ans[3][2] = -p[2];

            ans[4][0] = p[2];
            ans[4][1] = p[0];
            ans[4][2] = p[1];

            ans[5][0] = p[2];
            ans[5][1] = -p[0];
            ans[5][2] = -p[1];

            ans[6][0] = -p[2];
            ans[6][1] = -p[0];
            ans[6][2] = p[1];

            ans[7][0] = -p[2];
            ans[7][1] = p[0];
            ans[7][2] = -p[1];

            ans[8][0] = p[1];
            ans[8][1] = p[2];
            ans[8][2] = p[0];

            ans[9][0] = -p[1];
            ans[9][1] = p[2];
            ans[9][2] = -p[0];

            ans[10][0] = p[1];
            ans[10][1] = -p[2];
            ans[10][2] = -p[0];

            ans[11][0] = -p[1];
            ans[11][1] = -p[2];
            ans[11][2] = p[0];

            ans[12][0] = p[1];
            ans[12][1] = p[0];
            ans[12][2] = -p[2];

            ans[13][0] = -p[1];
            ans[13][1] = -p[0];
            ans[13][2] = -p[2];

            ans[14][0] = p[1];
            ans[14][1] = -p[0];
            ans[14][2] = p[2];

            ans[15][0] = -p[1];
            ans[15][1] = p[0];
            ans[15][2] = p[2];

            ans[16][0] = p[0];
            ans[16][1] = p[2];
            ans[16][2] = -p[1];

            ans[17][0] = -p[0];
            ans[17][1] = p[2];
            ans[17][2] = p[1];

            ans[18][0] = -p[0];
            ans[18][1] = -p[2];
            ans[18][2] = -p[1];

            ans[19][0] = p[0];
            ans[19][1] = -p[2];
            ans[19][2] = p[1];

            ans[20][0] = p[2];
            ans[20][1] = p[1];
            ans[20][2] = -p[0];

            ans[21][0] = p[2];
            ans[21][1] = -p[1];
            ans[21][2] = p[0];

            ans[22][0] = -p[2];
            ans[22][1] = p[1];
            ans[22][2] = p[0];

            ans[23][0] = -p[2];
            ans[23][1] = -p[1];
            ans[23][2] = -p[0];

            ans[24][0] = -p[0];
            ans[24][1] = -p[1];
            ans[24][2] = -p[2];

            ans[25][0] = p[0];
            ans[25][1] = p[1];
            ans[25][2] = -p[2];

            ans[26][0] = p[0];
            ans[26][1] = -p[1];
            ans[26][2] = p[2];

            ans[27][0] = -p[0];
            ans[27][1] = p[1];
            ans[27][2] = p[2];

            ans[28][0] = -p[2];
            ans[28][1] = -p[0];
            ans[28][2] = -p[1];

            ans[29][0] = -p[2];
            ans[29][1] = p[0];
            ans[29][2] = p[1];

            ans[30][0] = p[2];
            ans[30][1] = p[0];
            ans[30][2] = -p[1];

            ans[31][0] = p[2];
            ans[31][1] = -p[0];
            ans[31][2] = p[1];

            ans[32][0] = -p[1];
            ans[32][1] = -p[2];
            ans[32][2] = -p[0];

            ans[33][0] = p[1];
            ans[33][1] = -p[2];
            ans[33][2] = p[0];

            ans[34][0] = -p[1];
            ans[34][1] = p[2];
            ans[34][2] = p[0];

            ans[35][0] = p[1];
            ans[35][1] = p[2];
            ans[35][2] = -p[0];

            ans[36][0] = -p[1];
            ans[36][1] = -p[0];
            ans[36][2] = p[2];

            ans[37][0] = p[1];
            ans[37][1] = p[0];
            ans[37][2] = p[2];

            ans[38][0] = -p[1];
            ans[38][1] = p[0];
            ans[38][2] = -p[2];

            ans[39][0] = p[1];
            ans[39][1] = -p[0];
            ans[39][2] = -p[2];

            ans[40][0] = -p[0];
            ans[40][1] = -p[2];
            ans[40][2] = p[1];

            ans[41][0] = p[0];
            ans[41][1] = -p[2];
            ans[41][2] = -p[1];

            ans[42][0] = p[0];
            ans[42][1] = p[2];
            ans[42][2] = p[1];

            ans[43][0] = -p[0];
            ans[43][1] = p[2];
            ans[43][2] = -p[1];

            ans[44][0] = -p[2];
            ans[44][1] = -p[1];
            ans[44][2] = p[0];

            ans[45][0] = -p[2];
            ans[45][1] = p[1];
            ans[45][2] = -p[0];

            ans[46][0] = p[2];
            ans[46][1] = -p[1];
            ans[46][2] = -p[0];

            ans[47][0] = p[2];
            ans[47][1] = p[1];
            ans[47][2] = p[0];

            return ans;
        }
    }
}