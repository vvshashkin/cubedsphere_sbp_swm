#!/bin/bash

source $(dirname "$0")"/gen_namelist.sh"

EXE=$1/TS2_MAIN
Nprocs=$2

NAMELIST_TEMPLATE="
&domain\n
    N = %%%N,\n
    Nz = 1,\n
    topology_type    = 'cube',\n
    staggering_type  = 'Ah',\n
    metric_type      = 'ecs'\n
/\n
&metric\n
    rotation_axis(1) = 1.0,
    rotation_axis(2) = 0.0,
    rotation_axis(3) = 1.0,
/\n
&shallow_water_model\n
    swm_op_type       = 'vector_invariant'\n
    v_components_type = 'covariant'\n
    div_op_name       = '%%%divergence',\n
    grad_op_name      = '%%%gradient',\n
    coriolis_op_name  = 'coriolis_colocated',\n
    curl_op_name      = '%%%curl',\n
    KE_op_name        = 'KE_colocated',\n
    co2contra_op_name = 'co2contra_colocated',\n
    massflux_op_name  = 'massflux_colocated',\n
    quadrature_name   = '%%%quadrature',\n
    diff_time_scheme  = 'explicit_Eul1'\n
    uv_diff_coeff     =  0.03,\n
    hordiff_uv_name   = '%%%uv_diff',\n
    h_diff_coeff      =  0.01,\n
    hordiff_h_name    = '%%%h_diff',\n
    dt=%%%dt,\n
    tau_write = 86400.0,\n
    tau_diagnostics = 3600.0\n
    simulation_time_days  = 10.0,\n
    simulation_time_hours = 0.0,\n
    simulation_time_min   = 0.0,\n
    simulation_time_sec   = 0.0,\n
/"

run_ts2(){
    echo $4
	gen_namelist $1 $2 $3 "$4"  > namelist_swm
	mpirun -n $Nprocs $EXE &> swm_N$1_dt$2_Ah$3.out
    grep "l2_h" swm_N$1_dt$2_Ah$3.out > errors_N$1_dt$2_Ah$3.txt
    mv h.dat h_N$1_Ah$3.dat
}
run_ts2 020 800 21 "$NAMELIST_TEMPLATE"
run_ts2 040 400 21 "$NAMELIST_TEMPLATE"
run_ts2 080 200 21 "$NAMELIST_TEMPLATE"
run_ts2 160 100 21 "$NAMELIST_TEMPLATE"
run_ts2 020 800 42 "$NAMELIST_TEMPLATE"
run_ts2 040 400 42 "$NAMELIST_TEMPLATE"
run_ts2 080 200 42 "$NAMELIST_TEMPLATE"
run_ts2 160 100 42 "$NAMELIST_TEMPLATE"
run_ts2 020 800 43 "$NAMELIST_TEMPLATE"
run_ts2 040 400 43 "$NAMELIST_TEMPLATE"
run_ts2 080 200 43 "$NAMELIST_TEMPLATE"
run_ts2 160 100 43 "$NAMELIST_TEMPLATE"
run_ts2 020 800 63 "$NAMELIST_TEMPLATE"
run_ts2 040 400 63 "$NAMELIST_TEMPLATE"
run_ts2 080 200 63 "$NAMELIST_TEMPLATE"
run_ts2 160 100 63 "$NAMELIST_TEMPLATE"
