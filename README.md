# Vitis_HLS_Agentic_AI

256-point Radix-2 DIT FFT IP core implemented in HLS C++ for the PYNQ-Z2 (xc7z020clg400-1), targeting Vitis Unified IDE 2025.1.

## Project Structure

```
Vitis_HLS_Agentic_AI/
├── src/
│   ├── fft.h          # Header: type definitions, AXI-Stream types, pack/unpack helpers
│   ├── fft.cpp        # FFT source: 256-point Radix-2 DIT with DATAFLOW optimization
│   └── fft_tb.cpp     # Testbench: dual-tone sine wave verification
├── run_hls.tcl        # TCL script: full Vitis HLS flow (csim → synth → cosim → export)
└── readme             # Original design specification
```

## Design Overview

| Parameter            | Value                                      |
|----------------------|--------------------------------------------|
| FFT Size             | 256-point                                  |
| Algorithm            | Radix-2 Decimation-In-Time (DIT)           |
| Data Format          | 64-bit unsigned int (32-bit real + 32-bit imag, `ap_fixed<32,16>`) |
| Input/Output Ports   | AXI-Stream (AXIS) with TLAST on last sample |
| Control Interface    | AXI-Lite (`s_axilite`) for start/stop      |
| Optimization         | DATAFLOW (pipelined read → bit-reverse → FFT stages → write) |
| Twiddle Factors      | Precomputed ROM (128 complex entries)      |
| Target FPGA          | xc7z020clg400-1 (PYNQ-Z2)                 |
| Clock Frequency      | 100 MHz (10 ns period)                     |
| Top Function         | `fft`                                      |

## Data Flow

```
AXI-Stream In → Read Input → Bit-Reverse → 8 Butterfly Stages → Write Output → AXI-Stream Out
                  (256 samples)                                     (TLAST on last)
```

## How to Run

### Prerequisites
- Vitis Unified IDE 2025.1 installed
- Environment sourced: `source <Vitis_install>/settings64.sh`

### Execute Full Flow
```bash
cd Vitis_HLS_Agentic_AI
vitis-run --tcl run_hls.tcl
```

This runs all four steps sequentially:
1. **C Simulation** — Validates functional correctness
2. **C Synthesis** — Generates RTL from the HLS C++ design
3. **C/RTL Co-Simulation** — Verifies RTL behavior matches C simulation
4. **IP Package** — Exports the design as a Vivado IP catalog entry

The packaged IP will be available under `FFT/solution1/impl/ip/`.
