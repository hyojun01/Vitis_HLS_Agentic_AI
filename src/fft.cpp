/**
 * @file fft.cpp
 * @brief 256-point Radix-2 DIT FFT Implementation
 *
 * Implements a decimation-in-time FFT with:
 *   - AXI-Stream input/output with TLAST
 *   - AXI-Lite control interface
 *   - DATAFLOW optimization (read -> bit_reverse -> fft_stages -> write)
 *   - Precomputed twiddle factors stored in internal ROM
 */

#include "fft.h"

/* ================================================================
 *  Precomputed Twiddle Factor Tables (stored in internal memory)
 *  W_N^k = cos(-2*pi*k/N) + j*sin(-2*pi*k/N), N=256, k=0..127
 * ================================================================ */
static const fixed_t TW_REAL[FFT_SIZE / 2] = {
     1.0000000000,  0.9996988187,  0.9987954562,  0.9972904567,
     0.9951847267,  0.9924795346,  0.9891765100,  0.9852776424,
     0.9807852804,  0.9757021300,  0.9700312532,  0.9637760658,
     0.9569403357,  0.9495281806,  0.9415440652,  0.9329927988,
     0.9238795325,  0.9142097557,  0.9039892931,  0.8932243012,
     0.8819212643,  0.8700869911,  0.8577286100,  0.8448535652,
     0.8314696123,  0.8175848132,  0.8032075315,  0.7883464276,
     0.7730104534,  0.7572088465,  0.7409511254,  0.7242470830,
     0.7071067812,  0.6895405447,  0.6715589548,  0.6531728430,
     0.6343932842,  0.6152315906,  0.5956993045,  0.5758081914,
     0.5555702330,  0.5349976199,  0.5141027442,  0.4928981922,
     0.4713967368,  0.4496113297,  0.4275550934,  0.4052413140,
     0.3826834324,  0.3598950365,  0.3368898534,  0.3136817404,
     0.2902846773,  0.2667127575,  0.2429801799,  0.2191012402,
     0.1950903220,  0.1709618888,  0.1467304745,  0.1224106752,
     0.0980171403,  0.0735645636,  0.0490676743,  0.0245412285,
     0.0000000000, -0.0245412285, -0.0490676743, -0.0735645636,
    -0.0980171403, -0.1224106752, -0.1467304745, -0.1709618888,
    -0.1950903220, -0.2191012402, -0.2429801799, -0.2667127575,
    -0.2902846773, -0.3136817404, -0.3368898534, -0.3598950365,
    -0.3826834324, -0.4052413140, -0.4275550934, -0.4496113297,
    -0.4713967368, -0.4928981922, -0.5141027442, -0.5349976199,
    -0.5555702330, -0.5758081914, -0.5956993045, -0.6152315906,
    -0.6343932842, -0.6531728430, -0.6715589548, -0.6895405447,
    -0.7071067812, -0.7242470830, -0.7409511254, -0.7572088465,
    -0.7730104534, -0.7883464276, -0.8032075315, -0.8175848132,
    -0.8314696123, -0.8448535652, -0.8577286100, -0.8700869911,
    -0.8819212643, -0.8932243012, -0.9039892931, -0.9142097557,
    -0.9238795325, -0.9329927988, -0.9415440652, -0.9495281806,
    -0.9569403357, -0.9637760658, -0.9700312532, -0.9757021300,
    -0.9807852804, -0.9852776424, -0.9891765100, -0.9924795346,
    -0.9951847267, -0.9972904567, -0.9987954562, -0.9996988187
};

static const fixed_t TW_IMAG[FFT_SIZE / 2] = {
     0.0000000000, -0.0245412285, -0.0490676743, -0.0735645636,
    -0.0980171403, -0.1224106752, -0.1467304745, -0.1709618888,
    -0.1950903220, -0.2191012402, -0.2429801799, -0.2667127575,
    -0.2902846773, -0.3136817404, -0.3368898534, -0.3598950365,
    -0.3826834324, -0.4052413140, -0.4275550934, -0.4496113297,
    -0.4713967368, -0.4928981922, -0.5141027442, -0.5349976199,
    -0.5555702330, -0.5758081914, -0.5956993045, -0.6152315906,
    -0.6343932842, -0.6531728430, -0.6715589548, -0.6895405447,
    -0.7071067812, -0.7242470830, -0.7409511254, -0.7572088465,
    -0.7730104534, -0.7883464276, -0.8032075315, -0.8175848132,
    -0.8314696123, -0.8448535652, -0.8577286100, -0.8700869911,
    -0.8819212643, -0.8932243012, -0.9039892931, -0.9142097557,
    -0.9238795325, -0.9329927988, -0.9415440652, -0.9495281806,
    -0.9569403357, -0.9637760658, -0.9700312532, -0.9757021300,
    -0.9807852804, -0.9852776424, -0.9891765100, -0.9924795346,
    -0.9951847267, -0.9972904567, -0.9987954562, -0.9996988187,
    -1.0000000000, -0.9996988187, -0.9987954562, -0.9972904567,
    -0.9951847267, -0.9924795346, -0.9891765100, -0.9852776424,
    -0.9807852804, -0.9757021300, -0.9700312532, -0.9637760658,
    -0.9569403357, -0.9495281806, -0.9415440652, -0.9329927988,
    -0.9238795325, -0.9142097557, -0.9039892931, -0.8932243012,
    -0.8819212643, -0.8700869911, -0.8577286100, -0.8448535652,
    -0.8314696123, -0.8175848132, -0.8032075315, -0.7883464276,
    -0.7730104534, -0.7572088465, -0.7409511254, -0.7242470830,
    -0.7071067812, -0.6895405447, -0.6715589548, -0.6531728430,
    -0.6343932842, -0.6152315906, -0.5956993045, -0.5758081914,
    -0.5555702330, -0.5349976199, -0.5141027442, -0.4928981922,
    -0.4713967368, -0.4496113297, -0.4275550934, -0.4052413140,
    -0.3826834324, -0.3598950365, -0.3368898534, -0.3136817404,
    -0.2902846773, -0.2667127575, -0.2429801799, -0.2191012402,
    -0.1950903220, -0.1709618888, -0.1467304745, -0.1224106752,
    -0.0980171403, -0.0735645636, -0.0490676743, -0.0245412285
};

/* ================================================================
 *  Internal Sub-Functions (DATAFLOW tasks)
 * ================================================================ */

/**
 * @brief Compute bit-reversed index for LOG2_FFT_SIZE-bit value
 */
static ap_uint<LOG2_FFT_SIZE> bit_reverse_idx(ap_uint<LOG2_FFT_SIZE> idx) {
    #pragma HLS INLINE
    ap_uint<LOG2_FFT_SIZE> rev = 0;
    for (int i = 0; i < LOG2_FFT_SIZE; i++) {
        #pragma HLS UNROLL
        rev[LOG2_FFT_SIZE - 1 - i] = idx[i];
    }
    return rev;
}

/**
 * @brief Read 256 complex samples from AXI-Stream into buffers
 */
static void read_input(hls::stream<axis_t> &in_stream,
                       fixed_t real_buf[FFT_SIZE],
                       fixed_t imag_buf[FFT_SIZE]) {
    READ_LOOP: for (int i = 0; i < FFT_SIZE; i++) {
        #pragma HLS PIPELINE II=1
        axis_t val = in_stream.read();
        unpack_data(val.data, real_buf[i], imag_buf[i]);
    }
}

/**
 * @brief Apply bit-reversal permutation on 256 complex samples
 *
 * In DIT FFT, input is rearranged so that element at position i
 * is stored at position bit_rev(i) before the butterfly stages.
 */
static void bit_reverse(fixed_t in_real[FFT_SIZE], fixed_t in_imag[FFT_SIZE],
                        fixed_t out_real[FFT_SIZE], fixed_t out_imag[FFT_SIZE]) {
    BIT_REV_LOOP: for (int i = 0; i < FFT_SIZE; i++) {
        #pragma HLS PIPELINE II=1
        ap_uint<LOG2_FFT_SIZE> rev = bit_reverse_idx(i);
        out_real[rev] = in_real[i];
        out_imag[rev] = in_imag[i];
    }
}

/**
 * @brief Perform 8 butterfly stages of the Radix-2 DIT FFT
 *
 * Uses precomputed twiddle factors from internal ROM.
 * Each stage has N/2 = 128 butterfly operations.
 */
static void fft_stages(fixed_t in_real[FFT_SIZE], fixed_t in_imag[FFT_SIZE],
                       fixed_t out_real[FFT_SIZE], fixed_t out_imag[FFT_SIZE]) {

    // Working buffers for in-place computation
    fixed_t wr[FFT_SIZE], wi[FFT_SIZE];

    // Copy input to working buffers
    COPY_IN: for (int i = 0; i < FFT_SIZE; i++) {
        #pragma HLS PIPELINE II=1
        wr[i] = in_real[i];
        wi[i] = in_imag[i];
    }

    // Process 8 butterfly stages (log2(256) = 8)
    STAGE_LOOP: for (int stage = 0; stage < LOG2_FFT_SIZE; stage++) {
        int half      = 1 << stage;               // butterflies per group
        int tw_stride = FFT_SIZE >> (stage + 1);   // twiddle index stride

        BUTTERFLY_LOOP: for (int k = 0; k < FFT_SIZE / 2; k++) {
            #pragma HLS PIPELINE II=1
            #pragma HLS DEPENDENCE variable=wr inter false
            #pragma HLS DEPENDENCE variable=wi inter false

            // Compute indices
            int j       = k & (half - 1);                  // position within group
            int group   = (k >> stage) << (stage + 1);     // group start index
            int idx_top = group + j;
            int idx_bot = group + j + half;
            int tw_idx  = j * tw_stride;

            // Load twiddle factor from ROM
            fixed_t tw_r = TW_REAL[tw_idx];
            fixed_t tw_i = TW_IMAG[tw_idx];

            // Load butterfly operands
            fixed_t ar = wr[idx_top];
            fixed_t ai = wi[idx_top];
            fixed_t br = wr[idx_bot];
            fixed_t bi = wi[idx_bot];

            // Complex multiply: (br + j*bi) * (tw_r + j*tw_i)
            fixed_t tr = (fixed_t)(br * tw_r) - (fixed_t)(bi * tw_i);
            fixed_t ti = (fixed_t)(br * tw_i) + (fixed_t)(bi * tw_r);

            // Butterfly output
            wr[idx_top] = ar + tr;
            wi[idx_top] = ai + ti;
            wr[idx_bot] = ar - tr;
            wi[idx_bot] = ai - ti;
        }
    }

    // Copy results to output
    COPY_OUT: for (int i = 0; i < FFT_SIZE; i++) {
        #pragma HLS PIPELINE II=1
        out_real[i] = wr[i];
        out_imag[i] = wi[i];
    }
}

/**
 * @brief Write 256 complex samples from buffers to AXI-Stream
 *
 * Sets TLAST=1 on the last (256th) sample.
 */
static void write_output(fixed_t real_buf[FFT_SIZE], fixed_t imag_buf[FFT_SIZE],
                         hls::stream<axis_t> &out_stream) {
    WRITE_LOOP: for (int i = 0; i < FFT_SIZE; i++) {
        #pragma HLS PIPELINE II=1
        axis_t val;
        val.data = pack_data(real_buf[i], imag_buf[i]);
        val.keep = -1;   // all bytes valid
        val.strb = -1;
        val.last = (i == FFT_SIZE - 1) ? 1 : 0;
        out_stream.write(val);
    }
}

/* ================================================================
 *  Top-Level Function
 * ================================================================ */
void fft(hls::stream<axis_t> &in_stream, hls::stream<axis_t> &out_stream) {
    // ---- Interface Pragmas ----
    #pragma HLS INTERFACE axis port=in_stream
    #pragma HLS INTERFACE axis port=out_stream
    #pragma HLS INTERFACE s_axilite port=return

    // ---- DATAFLOW: overlap I/O with computation ----
    #pragma HLS DATAFLOW

    // Intermediate buffers (become PIPO buffers in DATAFLOW)
    fixed_t in_real[FFT_SIZE],  in_imag[FFT_SIZE];
    fixed_t rev_real[FFT_SIZE], rev_imag[FFT_SIZE];
    fixed_t out_real[FFT_SIZE], out_imag[FFT_SIZE];

    // Stage 1: Read input stream
    read_input(in_stream, in_real, in_imag);

    // Stage 2: Bit-reversal permutation
    bit_reverse(in_real, in_imag, rev_real, rev_imag);

    // Stage 3: FFT butterfly stages
    fft_stages(rev_real, rev_imag, out_real, out_imag);

    // Stage 4: Write output stream
    write_output(out_real, out_imag, out_stream);
}
