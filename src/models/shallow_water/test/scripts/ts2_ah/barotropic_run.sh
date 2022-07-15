#!/bin/bash

EXE=$1/BAROTROPIC_INST_MAIN
Nprocs=$2
source $(dirname "$0")"/gen_namelist.sh"

NAMELIST_TEMPLATE="
&domain\n
    N = %%%N,\n
    Nz = 1,\n
    topology_type    = 'cube',\n
    staggering_type  = 'Ah',\n
    metric_type      = 'ecs'\n
/\n
&metric\n
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
    simulation_time_days  = 30.0,\n
    simulation_time_hours = 0.0,\n
    simulation_time_min   = 0.0,\n
    simulation_time_sec   = 0.0,\n
/"

run_barotropic(){
	gen_namelist $1 $2 $3 "$4"  > namelist_swm
	mpirun -n $Nprocs $EXE &> swm_N$1_dt$2_Ah$3.out
    mv curl.dat curl_N$1_Ah$3.dat
    #grep "l2_h" swm_N$1_dt$2_Ah$3.out > errors_N$1_dt$2_Ah$3.txt
}
run_barotropic 096 200 21 "$NAMELIST_TEMPLATE"
run_barotropic 128 150 21 "$NAMELIST_TEMPLATE"
run_barotropic 192 100 21 "$NAMELIST_TEMPLATE"
run_barotropic 256 075 21 "$NAMELIST_TEMPLATE"
run_barotropic 096 200 42 "$NAMELIST_TEMPLATE"
run_barotropic 128 150 42 "$NAMELIST_TEMPLATE"
run_barotropic 192 100 42 "$NAMELIST_TEMPLATE"
run_barotropic 256 075 42 "$NAMELIST_TEMPLATE"
run_barotropic 096 200 43 "$NAMELIST_TEMPLATE"
run_barotropic 128 150 43 "$NAMELIST_TEMPLATE"
run_barotropic 192 100 43 "$NAMELIST_TEMPLATE"
run_barotropic 256 075 43 "$NAMELIST_TEMPLATE"
run_barotropic 096 200 63 "$NAMELIST_TEMPLATE"
run_barotropic 128 150 63 "$NAMELIST_TEMPLATE"
run_barotropic 192 100 63 "$NAMELIST_TEMPLATE"
run_barotropic 256 075 63 "$NAMELIST_TEMPLATE"
