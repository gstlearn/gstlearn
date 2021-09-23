/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#define CONCATENATE(x,y) x ## y
#define EXPAND(x) x

// Inspired by :
// - https://stackoverflow.com/q/14732803
// - https://stackoverflow.com/questions/1872220/is-it-possible-to-iterate-over-arguments-in-variadic-macros/1872506#1872506
#define NARG(...) NARG_(__VA_ARGS__, RSEQ_N())
#define NARG_(...) EXPAND(ARG_N(__VA_ARGS__))
#define ARG_N( _1,  _2,  _3,  _4,  _5,  _6,  _7,  _8,\
               _9, _10, _11, _12, _13, _14, _15, _16,\
              _17, _18, _19, _20, _21, _22, _23, _24,\
              _25, _26, _27, _28, _29, _30, _31, _32,\
              N, ...) N
#define RSEQ_N() 32, 31, 30, 29, 28, 27, 26, 25,\
                 24, 23, 22, 21, 20, 19, 18, 17,\
                 16, 15, 14, 13, 12, 11, 10,  9,\
                  8,  7,  6,  5,  4,  3,  2,  1,  0

// REPEAT macro (repeat 'what' instruction for each remaining SINGLE argument)
#define REPEAT_1(what, x) what(x)
#define REPEAT_2(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_1(what, __VA_ARGS__))
#define REPEAT_3(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_2(what, __VA_ARGS__))
#define REPEAT_4(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_3(what, __VA_ARGS__))
#define REPEAT_5(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_4(what, __VA_ARGS__))
#define REPEAT_6(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_5(what, __VA_ARGS__))
#define REPEAT_7(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_6(what, __VA_ARGS__))
#define REPEAT_8(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_7(what, __VA_ARGS__))
#define REPEAT_9(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_8(what, __VA_ARGS__))
#define REPEAT_10(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_9(what, __VA_ARGS__))
#define REPEAT_11(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_10(what, __VA_ARGS__))
#define REPEAT_12(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_11(what, __VA_ARGS__))
#define REPEAT_13(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_12(what, __VA_ARGS__))
#define REPEAT_14(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_13(what, __VA_ARGS__))
#define REPEAT_15(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_14(what, __VA_ARGS__))
#define REPEAT_16(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_15(what, __VA_ARGS__))
#define REPEAT_17(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_16(what, __VA_ARGS__))
#define REPEAT_18(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_17(what, __VA_ARGS__))
#define REPEAT_19(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_18(what, __VA_ARGS__))
#define REPEAT_20(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_19(what, __VA_ARGS__))
#define REPEAT_21(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_20(what, __VA_ARGS__))
#define REPEAT_22(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_21(what, __VA_ARGS__))
#define REPEAT_23(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_22(what, __VA_ARGS__))
#define REPEAT_24(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_23(what, __VA_ARGS__))
#define REPEAT_25(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_24(what, __VA_ARGS__))
#define REPEAT_26(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_25(what, __VA_ARGS__))
#define REPEAT_27(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_26(what, __VA_ARGS__))
#define REPEAT_28(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_27(what, __VA_ARGS__))
#define REPEAT_29(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_28(what, __VA_ARGS__))
#define REPEAT_30(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_29(what, __VA_ARGS__))
#define REPEAT_31(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_30(what, __VA_ARGS__))
#define REPEAT_32(what, x, ...)\
  what(x)\
  EXPAND(REPEAT_31(what, __VA_ARGS__))
  
#define REPEAT_(N, what, ...) EXPAND(CONCATENATE(REPEAT_, N)(what, __VA_ARGS__))
#define REPEAT(what, ...) REPEAT_(NARG(__VA_ARGS__), what, __VA_ARGS__)


// REPEAT macro (repeat 'what' instruction for each remaining PAIRS of arguments)
#define NARG2(...) NARG2_(__VA_ARGS__, RSEQ2_N())
#define NARG2_(...) EXPAND(ARG2_N(__VA_ARGS__))
#define ARG2_N( _1a,  _1b,  _2a,  _2b,  _3a,  _3b,  _4a,  _4b,\
                _5a,  _5b,  _6a,  _6b,  _7a,  _7b,  _8a,  _8b,\
                _9a,  _9b, _10a, _10b, _11a, _11b, _12a, _12b,\
               _13a, _13b, _14a, _14b, _15a, _15b, _16a, _16b,\
               _17a, _17b, _18a, _18b, _19a, _19b, _20a, _20b,\
               _21a, _21b, _22a, _22b, _23a, _23b, _24a, _24b,\
               _25a, _25b, _26a, _26b, _27a, _27b, _28a, _28b,\
               _29a, _29b, _30a, _30b, _31a, _31b, _32a, _32b,\
               N, ...) N
#define RSEQ2_N() 32, 32, 31, 31, 30, 30, 29, 29,\
                  28, 28, 27, 27, 26, 26, 25, 25,\
                  24, 24, 23, 23, 22, 22, 21, 21,\
                  20, 20, 19, 19, 18, 18, 17, 17,\
                  16, 16, 15, 15, 14, 14, 13, 13,\
                  12, 12, 11, 11, 10, 10,  9,  9,\
                   8,  8,  7,  7,  6,  6,  5,  5,\
                   4,  4,  3,  3,  2,  2,  1,  1,  0,  0

#define REPEAT2_1(what, x, y) what(x, y)
#define REPEAT2_2(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_1(what, __VA_ARGS__))
#define REPEAT2_3(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_2(what, __VA_ARGS__))
#define REPEAT2_4(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_3(what, __VA_ARGS__))
#define REPEAT2_5(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_4(what, __VA_ARGS__))
#define REPEAT2_6(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_5(what, __VA_ARGS__))
#define REPEAT2_7(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_6(what, __VA_ARGS__))
#define REPEAT2_8(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_7(what, __VA_ARGS__))
#define REPEAT2_9(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_8(what, __VA_ARGS__))
#define REPEAT2_10(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_9(what, __VA_ARGS__))
#define REPEAT2_11(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_10(what, __VA_ARGS__))
#define REPEAT2_12(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_11(what, __VA_ARGS__))
#define REPEAT2_13(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_12(what, __VA_ARGS__))
#define REPEAT2_14(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_13(what, __VA_ARGS__))
#define REPEAT2_15(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_14(what, __VA_ARGS__))
#define REPEAT2_16(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_15(what, __VA_ARGS__))
#define REPEAT2_17(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_16(what, __VA_ARGS__))
#define REPEAT2_18(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_17(what, __VA_ARGS__))
#define REPEAT2_19(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_18(what, __VA_ARGS__))
#define REPEAT2_20(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_19(what, __VA_ARGS__))
#define REPEAT2_21(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_20(what, __VA_ARGS__))
#define REPEAT2_22(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_21(what, __VA_ARGS__))
#define REPEAT2_23(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_22(what, __VA_ARGS__))
#define REPEAT2_24(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_23(what, __VA_ARGS__))
#define REPEAT2_25(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_24(what, __VA_ARGS__))
#define REPEAT2_26(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_25(what, __VA_ARGS__))
#define REPEAT2_27(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_26(what, __VA_ARGS__))
#define REPEAT2_28(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_27(what, __VA_ARGS__))
#define REPEAT2_29(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_28(what, __VA_ARGS__))
#define REPEAT2_30(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_29(what, __VA_ARGS__))
#define REPEAT2_31(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_30(what, __VA_ARGS__))
#define REPEAT2_32(what, x, y, ...)\
  what(x, y)\
  EXPAND(REPEAT2_31(what, __VA_ARGS__))

#define REPEAT2_(N, what, ...) EXPAND(CONCATENATE(REPEAT2_, N)(what, __VA_ARGS__))
#define REPEAT2(what, ...) REPEAT2_(NARG2(__VA_ARGS__), what, __VA_ARGS__)

// REPEAT macro (repeat 'what' instruction for each remaining TRIPLET of arguments + 1 fixed argument named a)
#define NARG3(...) NARG3_(__VA_ARGS__, RSEQ3_N())
#define NARG3_(...) EXPAND(ARG3_N(__VA_ARGS__))
#define ARG3_N( _1a,  _1b,  _1c,  _2a,  _2b,  _2c,  _3a,  _3b,  _3c,\
                _4a,  _4b,  _4c,  _5a,  _5b,  _5c,  _6a,  _6b,  _6c,\
                _7a,  _7b,  _7c,  _8a,  _8b,  _8c,  _9a,  _9b,  _9c,\
               _10a, _10b, _10c, _11a, _11b, _11c, _12a, _12b, _12c,\
               _13a, _13b, _13c, _14a, _14b, _14c, _15a, _15b, _15c,\
               _16a, _16b, _16c, _17a, _17b, _17c, _18a, _18b, _18c,\
               _19a, _19b, _19c, _20a, _20b, _20c, _21a, _21b, _21c,\
               _22a, _22b, _22c, _23a, _23b, _23c, _24a, _24b, _24c,\
               _25a, _25b, _25c, _26a, _26b, _26c, _27a, _27b, _27c,\
               _28a, _28b, _28c, _29a, _29b, _29c, _30a, _30b, _30c,\
               _31a, _31b, _31c, _32a, _32b, _32c, N, ...) N
#define RSEQ3_N() 32, 32, 32, 31, 31, 31, 30, 30, 30, 29, 29, 29,\
                  28, 28 ,28, 27, 27, 27, 26, 26, 26, 25, 25, 25,\
                  24, 24, 24, 23, 23, 23, 22, 22, 22, 21, 21, 21,\
                  20, 20, 20, 19, 19, 19, 18, 18, 18, 17, 17, 17,\
                  16, 16, 16, 15, 15, 15, 14, 14, 14, 13, 13, 13,\
                  12, 12, 12, 11, 11, 11, 10, 10, 10,  9,  9,  9,\
                   8,  8,  8,  7,  7,  7,  6,  6,  6,  5,  5,  5,\
                   4,  4,  4,  3,  3,  3,  2,  2,  2,  1,  1,  1,  0,  0,  0

#define REPEAT3_1(what, a, x, y, z) what(a, x, y, z)
#define REPEAT3_2(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_1(what, a, __VA_ARGS__))
#define REPEAT3_3(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_2(what, a, __VA_ARGS__))
#define REPEAT3_4(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_3(what, a, __VA_ARGS__))
#define REPEAT3_5(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_4(what, a, __VA_ARGS__))
#define REPEAT3_6(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_5(what, a, __VA_ARGS__))
#define REPEAT3_7(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_6(what, a, __VA_ARGS__))
#define REPEAT3_8(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_7(what, a, __VA_ARGS__))
#define REPEAT3_9(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_8(what, a, __VA_ARGS__))
#define REPEAT3_10(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_9(what, a, __VA_ARGS__))
#define REPEAT3_11(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_10(what, a, __VA_ARGS__))
#define REPEAT3_12(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_11(what, a, __VA_ARGS__))
#define REPEAT3_13(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_12(what, a, __VA_ARGS__))
#define REPEAT3_14(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_13(what, a, __VA_ARGS__))
#define REPEAT3_15(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_14(what, a, __VA_ARGS__))
#define REPEAT3_16(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_15(what, a, __VA_ARGS__))
#define REPEAT3_17(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_16(what, a, __VA_ARGS__))
#define REPEAT3_18(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_17(what, a, __VA_ARGS__))
#define REPEAT3_19(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_18(what, a, __VA_ARGS__))
#define REPEAT3_20(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_19(what, a, __VA_ARGS__))
#define REPEAT3_21(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_20(what, a, __VA_ARGS__))
#define REPEAT3_22(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_21(what, a, __VA_ARGS__))
#define REPEAT3_23(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_22(what, a, __VA_ARGS__))
#define REPEAT3_24(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_23(what, a, __VA_ARGS__))
#define REPEAT3_25(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_24(what, a, __VA_ARGS__))
#define REPEAT3_26(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_25(what, a, __VA_ARGS__))
#define REPEAT3_27(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_26(what, a, __VA_ARGS__))
#define REPEAT3_28(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_27(what, a, __VA_ARGS__))
#define REPEAT3_29(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_28(what, a, __VA_ARGS__))
#define REPEAT3_30(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_29(what, a, __VA_ARGS__))
#define REPEAT3_31(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_30(what, a, __VA_ARGS__))
#define REPEAT3_32(what, a, x, y, z, ...)\
  what(a, x, y, z)\
  EXPAND(REPEAT3_31(what, a, __VA_ARGS__))

#define REPEAT3_(N, what, a, x, y, z, ...) EXPAND(CONCATENATE(REPEAT3_, N)(what, a, x, y, z, __VA_ARGS__))
#define REPEAT3(what, a, x, y, z, ...) REPEAT3_(NARG3(x, y, z, __VA_ARGS__), what, a, x, y, z, __VA_ARGS__)

