#pragma once

#define TEST_NUMERIC_TYPES(FUN)\
  FUN<uint8_t>();\
  FUN<uint16_t>();\
  FUN<uint32_t>();\
  FUN<uint64_t>();\
  FUN<int8_t>();\
  FUN<int16_t>();\
  FUN<int32_t>();\
  FUN<int64_t>();\
  FUN<float>();\
  FUN<double>()
