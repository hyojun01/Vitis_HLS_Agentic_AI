/**
 * @file fft_tb.cpp
 * @brief Testbench for 256-point Radix-2 DIT FFT IP
 *
 * Generates a dual-tone sine wave (10 Hz + 30 Hz within 256 samples),
 * feeds it through the FFT, and verifies:
 *   1. Spectral peaks appear at the expected frequency bins (10 and 30)
 *   2. TLAST is asserted on the last output sample
 *
 * Returns 0 on PASS, 1 on FAIL (compatible with Vitis HLS csim/cosim).
 */

#include "fft.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    /* ---- Configuration ---- */
    const int TONE1 = 10;   // First tone: 10 cycles in 256 samples
    const int TONE2 = 30;   // Second tone: 30 cycles in 256 samples

    hls::stream<axis_t> in_stream("in_stream");
    hls::stream<axis_t> out_stream("out_stream");

    /* ================================================================
     *  Generate Dual-Tone Sine Wave Input
     * ================================================================ */
    printf("=== 256-Point FFT Testbench ===\n");
    printf("Input: dual-tone sine wave (f1=%d, f2=%d)\n\n", TONE1, TONE2);

    for (int n = 0; n < FFT_SIZE; n++) {
        // Dual-tone: x[n] = sin(2*pi*f1*n/N) + sin(2*pi*f2*n/N)
        double val = sin(2.0 * M_PI * TONE1 * n / FFT_SIZE)
                   + sin(2.0 * M_PI * TONE2 * n / FFT_SIZE);

        fixed_t real_val = (fixed_t)val;
        fixed_t imag_val = (fixed_t)0.0;

        axis_t sample;
        sample.data = pack_data(real_val, imag_val);
        sample.keep = -1;
        sample.strb = -1;
        sample.last = (n == FFT_SIZE - 1) ? 1 : 0;

        in_stream.write(sample);
    }

    /* ================================================================
     *  Run FFT
     * ================================================================ */
    printf("Running FFT ...\n");
    fft(in_stream, out_stream);
    printf("FFT completed.\n\n");

    /* ================================================================
     *  Read Output and Compute Magnitude Spectrum
     * ================================================================ */
    float magnitude[FFT_SIZE];
    float real_out[FFT_SIZE], imag_out[FFT_SIZE];
    int   tlast_error = 0;

    for (int i = 0; i < FFT_SIZE; i++) {
        axis_t result = out_stream.read();

        fixed_t r, im;
        unpack_data(result.data, r, im);

        real_out[i]  = (float)r;
        imag_out[i]  = (float)im;
        magnitude[i] = sqrtf(real_out[i] * real_out[i] + imag_out[i] * imag_out[i]);

        // Check TLAST: must be 1 only on last sample
        if (i == FFT_SIZE - 1) {
            if (result.last != 1) {
                printf("ERROR: TLAST not asserted on last sample (index %d)!\n", i);
                tlast_error = 1;
            }
        } else {
            if (result.last != 0) {
                printf("ERROR: TLAST asserted on non-last sample (index %d)!\n", i);
                tlast_error = 1;
            }
        }
    }

    /* ================================================================
     *  Print Magnitude Spectrum
     * ================================================================ */
    printf("Magnitude Spectrum (first half):\n");
    printf("%-6s  %-14s\n", "Bin", "Magnitude");
    printf("------  --------------\n");
    for (int i = 0; i < FFT_SIZE / 2; i++) {
        printf("[%3d]   %12.4f", i, magnitude[i]);
        if (i == TONE1 || i == TONE2) {
            printf("  <-- expected peak (f=%d)", i);
        }
        printf("\n");
    }

    /* ================================================================
     *  Verify Results
     * ================================================================ */
    // Find maximum magnitude for threshold calculation
    float max_mag = 0.0f;
    for (int i = 0; i < FFT_SIZE; i++) {
        if (magnitude[i] > max_mag) max_mag = magnitude[i];
    }

    float threshold = max_mag * 0.3f;
    printf("\nMax magnitude: %.4f\n", max_mag);
    printf("Detection threshold (30%%): %.4f\n\n", threshold);

    // Check that peaks exist at expected bins
    bool peak_f1 = (magnitude[TONE1] > threshold);
    bool peak_f2 = (magnitude[TONE2] > threshold);

    // Check mirror peaks (N-f)
    bool peak_f1_mirror = (magnitude[FFT_SIZE - TONE1] > threshold);
    bool peak_f2_mirror = (magnitude[FFT_SIZE - TONE2] > threshold);

    // Count unexpected peaks in first half (excluding DC)
    int unexpected_peaks = 0;
    for (int i = 1; i < FFT_SIZE / 2; i++) {
        if (i != TONE1 && i != TONE2 && magnitude[i] > threshold) {
            printf("WARNING: Unexpected peak at bin %d (magnitude=%.4f)\n", i, magnitude[i]);
            unexpected_peaks++;
        }
    }

    /* ================================================================
     *  Report Results
     * ================================================================ */
    printf("\n=== Verification Results ===\n");
    printf("Peak at bin %d (f1):         %s (mag=%.4f)\n", TONE1, peak_f1 ? "FOUND" : "MISSING", magnitude[TONE1]);
    printf("Peak at bin %d (f2):         %s (mag=%.4f)\n", TONE2, peak_f2 ? "FOUND" : "MISSING", magnitude[TONE2]);
    printf("Peak at bin %d (f1 mirror):  %s (mag=%.4f)\n", FFT_SIZE - TONE1, peak_f1_mirror ? "FOUND" : "MISSING", magnitude[FFT_SIZE - TONE1]);
    printf("Peak at bin %d (f2 mirror):  %s (mag=%.4f)\n", FFT_SIZE - TONE2, peak_f2_mirror ? "FOUND" : "MISSING", magnitude[FFT_SIZE - TONE2]);
    printf("TLAST check:                 %s\n", tlast_error ? "FAIL" : "PASS");
    printf("Unexpected peaks:            %d\n", unexpected_peaks);

    // Final pass/fail
    int result = 0;
    if (!peak_f1 || !peak_f2) {
        printf("\nFAIL: Expected peaks not found at frequency bins %d and %d.\n", TONE1, TONE2);
        result = 1;
    }
    if (tlast_error) {
        printf("\nFAIL: TLAST signal error.\n");
        result = 1;
    }
    if (result == 0) {
        printf("\n=== TEST PASSED ===\n");
    } else {
        printf("\n=== TEST FAILED ===\n");
    }

    return result;
}
