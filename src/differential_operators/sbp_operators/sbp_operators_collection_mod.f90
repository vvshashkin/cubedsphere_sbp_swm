module sbp_operators_collection_mod

implicit none

!Colocated grid differentiation

!Ah 2-1 scheme second derivative
real(kind=8), parameter :: D2_21_in(3) = [1, -2, 1]
real(kind=8), parameter :: D2_21_shift = -1

real(kind=8),    parameter :: D2_21_edge(3,1) = reshape([1, -2, 1], [3,1])
integer(kind=4), parameter :: lastnonzeroD2_21_edge(1) = [3]
real(kind=8),    parameter :: D2_21_boundary_proj(3) = [3.0_8/2,-2.0_8, 1.0_8/2]
!Ah 2-1 scheme written as SBP
real(kind=8), parameter :: Q21(2,1) = reshape( [-1._8, 1._8],  [2,1])
integer(kind=4), parameter :: lastnonzeroQ21(1) =[2]
real(kind=8), parameter :: Q21_A(2) = [0.5_8, 1.0_8]
!2-th order diff non-staggered inner stencil and shift
real(kind=8),    parameter :: Da2_in(3) = [-0.5_8,0._8,0.5_8]
integer(kind=4), parameter :: Da2_inshift = -1

!boundary block of SBP diff matrix (inv(H)*Q)
real(kind=8), parameter :: Q42(6,4) = reshape( &
[-1.4117647058823528_8,   1.7352941176470589_8, -0.2352941176470588_8, -0.08823529411764705_8, 0.0_8,                 0.0_8, &
 -0.5_8,                  0.0_8,                0.5_8,                  0.0_8,                 0.0_8,                 0.0_8, &
  0.09302325581395347_8, -0.686046511627907_8,  0.0_8,                  0.686046511627907_8,  -0.09302325581395347_8, 0.0_8, &
  0.030612244897959186_8,-0.0_8,               -0.6020408163265307_8,   0.0_8,                 0.653061224489796_8,  -0.0816326530612245_8], &
  [6,4])
integer(kind=4), parameter :: lastnonzeroQ42(4) =[4,3,5,6]
real(kind=8), parameter :: Q42_A(4) = [17._8/48._8, 59._8/48._8, 43._8/48._8, 49._8/48._8]

!boundary block of SBP diff matrix (inv(H)*Q)
real(kind=8), parameter :: Q63_A(6) = [0.13649d5 / 0.43200d5, 0.12013d5 / 0.8640d4, &
                                        0.2711d4 / 0.4320d4,  0.5359d4  / 0.4320d4, &
                                        0.7877d4 / 0.8640d4,  0.43801d5 / 0.43200d5 ]

real(kind=8), parameter :: q56 =5591070156686698065364559.d0/7931626489314500743872000.d0

real(kind=8), parameter :: q11 = -1.d0/2.d0
real(kind=8), parameter :: q12 = -0.953d3 / 0.16200d5 + q56
real(kind=8), parameter :: q13 = 0.715489d6 / 0.259200d6 - (4.d0 * q56)
real(kind=8), parameter :: q14 = -0.62639d5 / 0.14400d5 + (6.d0 * q56)
real(kind=8), parameter :: q15 = 0.147127d6 / 0.51840d5 - (4.d0 * q56)
real(kind=8), parameter :: q16 = -0.89387d5 / 0.129600d6 + q56;
real(kind=8), parameter :: q23 = -0.57139d5 / 0.8640d4 + (10.d0 * q56)
real(kind=8), parameter :: q24 = 0.745733d6 / 0.51840d5 - (20.d0 * q56)
real(kind=8), parameter :: q25 = -0.18343d5 / 0.1728d4 + (15.d0 * q56)
real(kind=8), parameter :: q26 = 0.240569d6 / 0.86400d5 - (4.d0 * q56)
real(kind=8), parameter :: q34 = -0.176839d6 / 0.12960d5 + (20.d0 * q56)
real(kind=8), parameter :: q35 = 0.242111d6 / 0.17280d5 - (20.d0 * q56)
real(kind=8), parameter :: q36 = -0.182261d6 / 0.43200d5 + (6.d0 * q56)
real(kind=8), parameter :: q45 = -0.165041d6 / 0.25920d5 + (10.d0 * q56)
real(kind=8), parameter :: q46 = 0.710473d6 / 0.259200d6 - (4.d0 * q56)
real(kind=8), parameter :: q47 = 1.d0/6.d1
real(kind=8), parameter :: q57 = -3.D0/2.d1
real(kind=8), parameter :: q58 = 1.d0/6.d1
real(kind=8), parameter :: q67 = 3.d0/4.d0
real(kind=8), parameter :: q68 = -3.d0/2.d1
real(kind=8), parameter :: q69 = 1.d0/6.d1;

real(kind=8), parameter :: Q63(9,6) = reshape( &
[ q11/Q63_A(1),  q12/Q63_A(1), q13/Q63_A(1), q14/Q63_A(1), q15/Q63_A(1), q16/Q63_A(1),        0.0_8,        0.0_8,        0.0_8, &
 -q12/Q63_A(2),         0.0_8, q23/Q63_A(2), q24/Q63_A(2), q25/Q63_A(2), q26/Q63_A(2),        0.0_8,        0.0_8,        0.0_8, &
 -q13/Q63_A(3), -q23/Q63_A(3),        0.0_8, q34/Q63_A(3), q35/Q63_A(3), q36/Q63_A(3),        0.0_8,        0.0_8,        0.0_8, &
 -q14/Q63_A(4), -q24/Q63_A(4), -q34/Q63_A(4),       0.0_8, q45/Q63_A(4), q46/Q63_A(4), q47/Q63_A(4),        0.0_8,        0.0_8, &
 -q15/Q63_A(5), -q25/Q63_A(5), -q35/Q63_A(5), -q45/Q63_A(5),      0.0_8, q56/Q63_A(5), q57/Q63_A(5), q58/Q63_A(5),        0.0_8, &
 -q16/Q63_A(6), -q26/Q63_A(6), -q36/Q63_A(6), -q46/Q63_A(6), -q56/Q63_A(6),     0.0_8, q67/Q63_A(6), q68/Q63_A(6), q69/Q63_A(6)  ], &
  [9,6])
integer(kind=4), parameter :: lastnonzeroQ63(6) =[6,6,6,7,8,9]
!6-th order diff non-staggered inner stencil and shift
real(kind=8),    parameter :: Da6_in(7) = [-1.d0/6.d1, 3.d0/2.d1, -3.d0/4.d0, 0.d0, 3.d0/4.d0, -3.d0/2.d1, 1.d0/6.d1];
integer(kind=4), parameter :: Da6_inshift = -3


real(kind=8), parameter :: Q43(6,4) = reshape( &
[-1.8280236842718398_8,   2.978054584697369_8,  -1.4653148294044123_8,  0.3078538227474199_8,  0.008136925288119662_8, -0.0007068190566565205_8, &
- 0.37814159126691327_8, -0.31480210500706657_8, 0.707290999364066_8,   0.04835554461933453_8,-0.06866771096803422_8,   0.005964863258613564_8,  &
  0.1162478575921292_8,  -0.8028500807354999_8,  0.21558841368737428_8, 0.5078566674295846_8, -0.03231754093993847_8,  -0.004525317033649745_8,  &
 -0.010287398573128293_8, 0.12592344801071614_8,-0.7341531396449149_8,  0.04979271660173103_8, 0.6506171865540596_8,  -0.08189281294846368_8],   &
                                              [6,4])
integer(kind=4), parameter :: lastnonzeroQ43(4) =[6,6,6,6]

!4-th order diff non-staggered inner stencil and shift
real(kind=8),    parameter :: Da4_in(5) = [1._8/12._8,-2._8/3.,0._8,2._8/3._8,-1._8/12._8]
integer(kind=4), parameter :: Da4_inshift = -2

!Staggered 21 interpolation
!cell centers to interfaces
real(kind=8), parameter :: W21_staggered_c2i(1,1) = reshape([1.0_8],[1,1])
integer(kind=4), parameter :: W21_staggered_c2i_last_nonzero(1) = [1]
real(kind=8), parameter :: W21_staggered_c2i_in_shift = -1
!interfaces to cell centers
real(kind=8), parameter :: W21_staggered_i2c(2,1) = reshape([0.5_8,0.5_8],[2,1])
integer(kind=4), parameter :: W21_staggered_i2c_last_nonzero(1) = [2]
real(kind=8), parameter :: W21_staggered_i2c_in_shift = 0
real(kind=8), parameter :: W21_staggered_in(2) = [0.5_8, 0.5_8]

!Stagered grid interpolation
!Interpolation weights
!interpolation from cell interfaces to cell centers
!non-optimized option (minimum l2 of coefficients)
real(kind=8), parameter :: W42_staggered_i2c_noopt(6,4) = reshape( &
[0.3585442775006853_8,  0.6119832845967544_8,  0.20040059830444001_8, -0.1709281604018784_8,  0.0_8,    0.0_8,    &
 0.16724022138832148_8, 0.36635673412335484_8, 0.2655658675883291_8,   0.20083717689999595_8, 0.0_8,    0.0_8,    &
-0.09295505106784016_8, 0.1751328387308906_8,  0.36859947574175345_8,  0.609222736595201_8,  -0.06_8,   0.0_8,    &
-0.04904109392262206_8,-0.04097407434909711_8, 0.16657143046607598_8,  0.42344373780564915_8, 0.5625_8,-0.0625_8],&
   [6,4])
!Optimized vector->scalar x^2 interpolation
real(kind=8), parameter ::  W42_staggered_i2c_opt1(6,4) = reshape( &
[ 0.4892964885611128_8,    0.48123107353632727_8, 0.0696483872440123_8, -0.040175949341450676_8, 0.0_8,     0.0_8, &
-0.07913099957454726_8,   0.6127279550862234_8,  0.5119370885511979_8, -0.045534044062872786_8, 0.0_8,     0.0_8, &
-0.0869982983587545_8,    0.1691760860218048_8,  0.3626427230326682_8,  0.6151794893042863_8,  -0.06_8,    0.0_8, &
0.018680545032460714_8, -0.10869571330417979_8, 0.09884979151099299_8, 0.4911653767607321_8,   0.5625_8, -0.0625_8],&
    [6,4])
!Optimized vector <-> scalar interp in x^2
real(kind=8), parameter :: W42_staggered_i2c_opt2(6,4) = reshape( &
 [ 0.5396575759839645_8, 0.4475591035840635_8, -0.014090935120015191_8, 0.026874255551989024_8, 0.0_8, 0.0_8, &
  -0.13589939229291748_8, 0.6508224599464276_8, 0.6060532569859003_8, -0.12097632463940905_8, 0.0_8, 0.0_8,   &
  -0.14875399135119022_8, 0.2102338641077255_8, 0.46579424583813406_8, 0.5327258814053357_8, -0.06_8, 0.0_8,  &
   0.07812389082006573_8, -0.14831895644807316_8, -0.00023375956403533802_8, 0.5704288251920487_8, 0.5625_8, -0.0625_8],&
     [6,4])
integer(kind=4), parameter :: W42_staggered_i2c_last_nonzero(4) = [4,4,5,6]

!interpolation from cell centers to interfaces
real(kind=8), parameter :: W42_staggered_c2i_noopt(5,4) = reshape( &
[0.9988019158947662_8,  0.37629049812372334_8, -0.24898674393171474_8, -0.12610567008674245_8,  0.0_8, &
 0.5893172370190968_8,  0.2849441265403871_8,   0.16216003586193575_8, -0.03642139942141965_8,  0.0_8, &
 0.21710064816314334_8, 0.23237013413978796_8,  0.3839577872309932_8,   0.16657143046607598_8,  0.0_8, &
-0.18778023255417622_8, 0.17820763584084146_8,  0.6435451442907053_8,   0.42940773411277094_8, -0.06338028169014084_8],&
   [5,4])
!optimized only vector to scalar interp
real(kind=8), parameter ::  W42_staggered_c2i_opt1(5,4) = reshape( &
[1.3630402181345285_8, -0.17804474904273135_8, -0.23303115631809243_8,  0.048035687226327554_8, 0.0_8, &
0.46340770044238916_8, 0.4765661872892848_8,   0.15664452409426372_8, -0.09661841182593758_8,  0.0_8, &
0.07545241951434666_8, 0.44794495248229815_8,  0.37775283649236274_8,  0.09884979151099299_8,  0.0_8, &
-0.04413695843145285_8,-0.04040344754874627_8,  0.649837488701711_8,    0.4980831989686297_8,  -0.06338028169014084_8],&
    [5,4])
!optimize both vector to scalar and scalar to vector interp
real(kind=8), parameter ::  W42_staggered_c2i_opt2(5,4) = reshape( &
[1.5033318188124725_8, -0.3057736326590643_8, -0.3984481911192596_8, 0.2008900049658833_8, 0.0_8, &
 0.43098284048835733_8, 0.5061952466249993_8, 0.194660985284931_8, -0.13183907239828724_8, 0.0_8, &
-0.01526517971334979_8, 0.5302965998626628_8, 0.485202339414723_8, -0.00023375956403533802_8, 0.0_8, &
 0.029523830043030195_8, -0.1073451894687714_8, 0.5627386071183124_8, 0.5784630339975705_8, -0.06338028169014084_8],&
   [5,4])
integer(kind=4), parameter :: W42_staggered_c2i_last_nonzero(4) = [4,4,4,5]
!inner interpolation stencil
real(kind=8), parameter :: W42_staggered_in(4) = [-1._8/16._8, 9._8/16._8, 9._8/16._8, -1._8/16._8]
real(kind=8), parameter :: W42_staggered_i2c_in_shift = -1
real(kind=8), parameter :: W42_staggered_c2i_in_shift = -2


!SBP staggered differences
real(kind=8), parameter :: D21_staggered_c2i(2,1) = reshape( [-1.0_8, 1.0_8], [2,1])
integer(kind=4), parameter :: D21_staggered_c2i_last_nonzero(1) = [2]
real(kind=8), parameter :: D21_staggered_i2c(2,1) = reshape( [-1.0_8, 1.0_8], [2,1])
integer(kind=4), parameter :: D21_staggered_i2c_last_nonzero(1) = [2]
real(kind=8), parameter :: D21_staggered_in(2) = [-1.0_8, 1.0_8]
real(kind=8), parameter :: D21_staggered_c2i_in_shift = -1
real(kind=8), parameter :: D21_staggered_i2c_in_shift =  0
!Projection operator for left side os stencil
real(kind=8), parameter :: D21_boundary_proj(2) = [1.5_8,-0.5_8]
!Mass matrices for cell-center and cell-interface variables (left sude of stencil)
real(kind=8), parameter :: D21_A_interfaces(2) = [0.5_8, 1._8]
real(kind=8), parameter :: D21_A_centers(2) = [1._8, 1._8]

real(kind=8), parameter :: D42_staggered_c2i(5,4) = reshape( &
                           [-2.0_8,       3.0_8,     -1.0_8,       0.0_8,       0.0_8,       &
                            -1.0_8,       1.0_8,      0.0_8,       0.0_8,       0.0_8,       &
                             1._8/24._8, -9._8/8._8,  9._8/8._8,  -1._8/24._8,  0.0_8,       &
                            -1._8/71._8,  6._8/71._8,-83._8/71._8, 81._8/71._8,-3._8/71._8], &
                            [5,4])

integer(kind=4), parameter :: D42_staggered_c2i_last_nonzero(4) = [3,2,4,5]

real(kind=8), parameter :: D42_staggered_i2c(5,3) = reshape( &
                           [-79.0_8/78.0_8, 27.0_8/26.0_8, -1.0_8/26.0_8, 1.0_8/78.0_8, 0.0_8,       &
                             2.0_8/21.0_8, -9.0_8/7.0_8, 9.0_8/7.0_8, -2.0_8/21.0_8, 0.0_8,      &
                             1.0_8/75.0_8, 0.0_8, -27.0_8/25.0_8, 83.0_8/75.0_8, -1.0_8/25.0_8], &
                                              [5,3])

integer(kind=4), parameter :: D42_staggered_i2c_last_nonzero(3) = [4,4,5]

real(kind=8), parameter :: D42_staggered_in(4) = [1.0_8/24.0_8, -9.0_8/8.0_8, 9.0_8/8.0_8, -1.0_8/24.0_8]
real(kind=8), parameter :: D42_staggered_c2i_in_shift = -2
real(kind=8), parameter :: D42_staggered_i2c_in_shift = -1
!Projection operator for left side os stencil
real(kind=8), parameter :: D42_boundary_proj(3) = [15.0_8/8.0_8, -5.0_8/4.0_8, 3.0_8/8.0_8]

!Mass matrices for cell-center and cell-interface variables (left sude of stencil)
real(kind=8), parameter :: D42_A_interfaces(5) = [7._8/18._8, 9._8/8._8, 1._8, 71._8/72._8, 1._8]
real(kind=8), parameter :: D42_A_centers(4) = [13._8/12._8, 7._8/8._8, 25._8/24._8, 1._8]

end module sbp_operators_collection_mod
