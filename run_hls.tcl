# ==============================================================
# run_hls.tcl - Vitis HLS 2025.1 Project TCL Script
# ==============================================================
# Creates a 256-point Radix-2 DIT FFT HLS Component named "FFT"
#
# Usage: vitis_hls -f run_hls.tcl
#
# Flow: C Simulation -> C Synthesis -> C/RTL Co-Simulation -> IP Package
# ==============================================================

# Create HLS project
open_project FFT

# Set top-level function
set_top fft

# Add source files
add_files src/fft.cpp -cflags "-I./src"

# Add testbench files
add_files -tb src/fft_tb.cpp -cflags "-I./src"

# Create solution targeting Vivado flow
open_solution "solution1" -flow_target vivado

# Set target FPGA part (PYNQ-Z2: xc7z020clg400-1)
set_part {xc7z020clg400-1}

# Set clock period to 10ns (100 MHz)
create_clock -period 10 -name default

# Set clock uncertainty
set_clock_uncertainty 12.5%

# ==============================================================
# Step 1: C Simulation
# ==============================================================
puts "=========================================="
puts " Step 1: Running C Simulation"
puts "=========================================="
csim_design

# ==============================================================
# Step 2: C Synthesis
# ==============================================================
puts "=========================================="
puts " Step 2: Running C Synthesis"
puts "=========================================="
csynth_design

# ==============================================================
# Step 3: C/RTL Co-Simulation
# ==============================================================
puts "=========================================="
puts " Step 3: Running C/RTL Co-Simulation"
puts "=========================================="
cosim_design

# ==============================================================
# Step 4: Export IP (Package)
# ==============================================================
puts "=========================================="
puts " Step 4: Exporting IP Package"
puts "=========================================="
export_design -format ip_catalog

# Close project
close_project

puts "=========================================="
puts " All steps completed successfully!"
puts "=========================================="

exit
