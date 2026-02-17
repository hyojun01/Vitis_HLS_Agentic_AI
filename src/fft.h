/**
 * @file fft.h
 * @brief 256-point Radix-2 DIT FFT IP Header
 *
 * Defines data types and interfaces for the FFT IP core.
 * - 64-bit AXI-Stream data: upper 32 bits = real, lower 32 bits = imaginary
 * - 32-bit signed fixed-point format: ap_fixed<32,16>
 * - AXI-Lite control interface for start/stop
 */

#ifndef FFT_H
#define FFT_H

#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>

/* ---------- Constants ---------- */
#define FFT_SIZE      256
#define LOG2_FFT_SIZE 8

/* ---------- Type Definitions ---------- */
// 32-bit signed fixed-point: 16 integer bits (including sign), 16 fractional bits
typedef ap_fixed<32, 16> fixed_t;

// 64-bit unsigned integer for packed complex data
typedef ap_uint<64> data_t;

// AXI-Stream packet type (64-bit data, no TUSER/TID/TDEST)
typedef ap_axiu<64, 0, 0, 0> axis_t;

/* ---------- Helper Functions ---------- */

/**
 * @brief Pack real and imaginary fixed-point values into a 64-bit word
 * @param real_val 32-bit signed fixed-point real part
 * @param imag_val 32-bit signed fixed-point imaginary part
 * @return 64-bit packed data (real in upper 32 bits, imag in lower 32 bits)
 */
inline data_t pack_data(fixed_t real_val, fixed_t imag_val) {
    data_t packed;
    packed.range(63, 32) = real_val.range();
    packed.range(31, 0)  = imag_val.range();
    return packed;
}

/**
 * @brief Unpack a 64-bit word into real and imaginary fixed-point values
 * @param packed 64-bit packed data
 * @param real_val Output: 32-bit signed fixed-point real part
 * @param imag_val Output: 32-bit signed fixed-point imaginary part
 */
inline void unpack_data(data_t packed, fixed_t &real_val, fixed_t &imag_val) {
    ap_uint<32> real_bits = packed.range(63, 32);
    ap_uint<32> imag_bits = packed.range(31, 0);
    real_val.range() = real_bits;
    imag_val.range() = imag_bits;
}

/* ---------- Top-Level Function ---------- */

/**
 * @brief 256-point Radix-2 DIT FFT
 *
 * Reads 256 complex samples from AXI-Stream input (TLAST=1 on last),
 * computes the FFT using bit-reversal and 8 butterfly stages,
 * and writes 256 complex results to AXI-Stream output (TLAST=1 on last).
 *
 * @param in_stream  AXI-Stream input (64-bit packed complex data)
 * @param out_stream AXI-Stream output (64-bit packed complex data)
 */
void fft(hls::stream<axis_t> &in_stream, hls::stream<axis_t> &out_stream);

#endif // FFT_H
