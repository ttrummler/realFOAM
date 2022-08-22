# Define colors
lc_0 = "lc rgb 'black' "
lc_1 = "lc rgb '#A50034' "
lc_2 = "lc rgb '#ED6E00' "
lc_3 = "lc rgb '#9D9D9C' "

# Define LineStyles
st_wl = " with lines "
st_lp = " with lp "
st_wp = " with points pt 7"

st_l0 = " lt  1 lw 2" # solid

# Define data sets
refCP   = '../refData/ref_data.csv'
refPR   = '../refData/PR.csv'
refSRK  = '../refData/SRK.csv'
simPR   = '../simData/PengRobinson.txt'
simSRK  = '../simData/SoaveRedlichKwong.txt'
refRKPR = '../refData/RKPR.csv'
simRKPR = '../simData/RKPR.txt'

# Plot 1
set terminal png size 800,500 enhanced font "Helvetica,12"
set output 'Density.png'
set title 'Density'

set ylabel "Density [kg/m^3]"
set xlabel "Temperature [K]"

set grid  

plot refPR   u 3:4 @st_wl @st_l0  @lc_1 title 'PR   python tool '  ,\
     refSRK  u 3:4 @st_wl @st_l0  @lc_2 title 'SRK  python tool '  ,\
     refRKPR u 3:4 @st_wl @st_l0  @lc_3 title 'RKPR python tool '  ,\
     simPR   u 2:3 @st_wp @st_l0  @lc_1 title 'PR   OpenFOAM    '  ,\
     simSRK  u 2:3 @st_wp @st_l0  @lc_2 title 'SRK  OpenFOAM    '  ,\
     simRKPR u 2:3 @st_wp @st_l0  @lc_3 title 'RKPR OpenFOAM    '  ,\
     refCP   u 3:4 @st_wl @st_l0  @lc_0 title 'Refdata(CoolProp)'
