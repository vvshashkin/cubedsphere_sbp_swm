gen_namelist(){
    N=$1
    DT=$2
    if [[ $3 -eq "21" ]]; then
        div="divergence_ah2"
        grad="gradient_ah21_sbp_ecs"
        curl="curl_"$div
        quad='SBP_Ah21_quadrature'
        uv_diff='hordiff_vec_xyz_Ah_sbp_21_narrow'
        h_diff='hordiff_scalar_Ah_sbp_21_narrow'
    elif  [[ $3 -eq "42" ]]; then
        div="divergence_ah42_sbp"
        grad="gradient_ah42_sbp_ecs"
        curl="curl_"$div
        quad='SBP_Ah42_quadrature'
        uv_diff='hordiff_vec_xyz_Ah_sbp_42_narrow'
        h_diff='hordiff_scalar_Ah_sbp_42_narrow'
    elif  [[ $3 -eq "43" ]]; then
        div="divergence_ah43_sbp"
        grad="gradient_ah43_sbp_ecs"
        curl="curl_"$div
        quad='SBP_Ah63_quadrature'
        uv_diff='hordiff_vec_xyz_Ah_sbp_42_narrow'
        h_diff='hordiff_scalar_Ah_sbp_42_narrow'
    elif  [[ $3 -eq "63" ]]; then
        div="divergence_ah63_sbp"
        grad="gradient_ah63_sbp_ecs"
        curl="curl_"$div
        quad='SBP_Ah63_quadrature'
        uv_diff='hordiff_vec_xyz_Ah_sbp_63_narrow'
        h_diff='hordiff_scalar_Ah_sbp_63_narrow'
	fi
    NAMELIST_TEMP=$4
	echo -e $NAMELIST_TEMP |
             sed "s/%%%divergence/$div/"  |
             sed "s/%%%gradient/$grad/"   |
             sed "s/%%%curl/$curl/"       |
             sed "s/%%%uv_diff/$uv_diff/" |
             sed "s/%%%h_diff/$h_diff/"   |
             sed "s/%%%quadrature/$quad/" |
             sed "s/%%%N/$N/" |
             sed "s/%%%dt/$DT/"
}

