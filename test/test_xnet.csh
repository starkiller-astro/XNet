#!/bin/csh -f
# Regression Testing for XNet
# Selection of test problem controlled by argument.
# All serial test problems = 0
# Thermonuclear SN     with alpha (= 1),   47 species (= 2) or 160 species (= 3)
# Tidally-Crushed WD   with alpha (= 4),  160 species (= 5)
# Nova        with 16 species CNO (= 7),  169 species (= 8)
# X-ray burst with 16 species CNO (= 9),  304 species (=10)
# Core Collapse SN     with alpha (=11),  160 species (=12) or 1072 network(=13) 
# Neutron Star Wind    with alpha (=15), 7852 species(=16)
# All parallel test problems = 30
# 4 different SN       with alpha (=31),  160 species (=32)

set xnet = ../source/xnet
echo 'Testing ' $xnet
if (! -d Test_Results) then
  mkdir Test_Results
endif

# TN SN tracer, from Ed Brown, with alpha network 
if ($argv[1] == 0 || $argv[1] == 1) then
  cat test_settings_small Test_Problems/setup_tnsn_alpha >! control
  echo 'Test: Thermonuclear SN with alpha network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_tnsn_alpha
endif

# TN SN tracer, from Ed Brown, with 160 species network 
if ($argv[1] == 0 || $argv[1] == 2) then
  cat test_settings Test_Problems/setup_tnsn_torch47 >! control
  echo 'Test: Thermonuclear SN with 47 species network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_tnsn_torch47
endif

# TN SN tracer, from Ed Brown, with 160 species network 
if ($argv[1] == 0 || $argv[1] == 3) then
  cat test_settings Test_Problems/setup_tnsn_sn160 >! control
  echo 'Test: Thermonuclear SN with 160 species network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_tnsn_sn160
endif

# TI SN tracer, from Stephan Rosswog, with alpha network 
if ($argv[1] == 0 || $argv[1] == 4) then
  cat test_settings_small Test_Problems/setup_tisn_alpha >! control
  echo 'Test: Tidally Induced SN with alpha network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_tisn_alpha
endif

# TI SN tracer, from Stephan Rosswog, with 160 species network 
if ($argv[1] == 0 || $argv[1] == 5) then
  cat test_settings Test_Problems/setup_tisn_sn160 >! control
  echo 'Test: Tidally Induced SN with 160 species network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_tisn_sn160
endif

# Nova zone, from Sumner Starrfield, with 16 species CNO network 
if ($argv[1] == 0 || $argv[1] == 7) then
  cat test_settings_small Test_Problems/setup_nova_cno >! control
  echo 'Test: Nova with minimal CNO network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_nova_cno
endif

# Nova zone, from Sumner Starrfield, with 189 species network from Chritian Iliadis
if ($argv[1] == 0 || $argv[1] == 8) then
  cat test_settings Test_Problems/setup_nova_Iliadis >! control
  echo 'Test: Nova with 189 species network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_nova_Iliadis
endif

# XRB tracer, from Jacob Fisker, with minimal CNO network 
if ($argv[1] == 0 || $argv[1] == 9) then
  cat test_settings_small Test_Problems/setup_xrb_cno >! control
  echo 'Test: X-ray burst with CNO network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_xrb_cno
endif

# XRB tracer, from Jacob Fisker, with 304 species rp-process network 
if ($argv[1] == 0 || $argv[1] == 10) then
  cat test_settings Test_Problems/setup_xrb_fisker >! control
  echo 'Test: X-ray burst with Fiskers network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_xrb_fisker
endif

# CC SN zone, from Carla Froehlich, with alpha network 
if ($argv[1] == 0 || $argv[1] == 11) then
  cat test_settings_small Test_Problems/setup_ccsn_alpha >! control
  echo 'Test: Core-Collapse SN with alpha network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_ccsn_alpha
endif

#  CC SN zone, from Carla Froehlich with 160 species network
if ($argv[1] == 0 || $argv[1] == 12) then
  cat test_settings Test_Problems/setup_ccsn_sn160 >! control
  echo 'Test: Core-Collapse SN with 160 species network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_ccsn_sn160
endif

#  CC SN zone, from Carla Froehlich with nu p-process network 
if ($argv[1] == 0 || $argv[1] == 13) then
  cat test_settings Test_Problems/setup_ccsn_nup >! control
  echo 'Test: Core-Collapse SN with nu-p process network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_ccsn_nup
endif

# Neutino driven wind example, from Josh Beun, with alpha network 
if ($argv[1] == 0 || $argv[1] == 15) then
  cat test_settings_small Test_Problems/setup_nuwind_alpha >! control
  echo 'Test: Neutrino-driven wind with alpha network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_nuwind_alpha
endif

# Neutino driven wind (Meyer & Brown) with full JINA REACLIB network
if ($argv[1] == 0 || $argv[1] == 16) then
  cat test_settings Test_Problems/setup_nuwind_rprocess >! control
  echo 'Test: Neutrino-driven wind with 7852 species network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_nuwind_rprocess
endif

# Parallel test, runs 4 different zones using SN160 
if ($argv[1] == 30 || $argv[1] == 31) then
  cat test_settings_parallel Test_Problems/setup_parallel_sn160 >! control
  echo 'Test: 4 different zones in parallel'
  mpirun -n 4 ../source/xnet_mpi
  mv -f net_diag01 Test_Results/net_diag_parallel_sn160_1
  mv -f net_diag11 Test_Results/net_diag_parallel_sn160_2
  mv -f net_diag21 Test_Results/net_diag_parallel_sn160_3
  mv -f net_diag31 Test_Results/net_diag_parallel_sn160_4
endif

# Parallel test, runs 4 different zones using alpha
if ($argv[1] == 30 || $argv[1] == 32) then
  cat test_settings_parallel Test_Problems/setup_parallel_alpha >! control
  echo 'Test: 4 different zones in parallel'
  mpirun -n 4 ../source/xnet_mpi
  mv -f net_diag01 Test_Results/net_diag_parallel_alpha_1
  mv -f net_diag11 Test_Results/net_diag_parallel_alpha_2
  mv -f net_diag21 Test_Results/net_diag_parallel_alpha_3
  mv -f net_diag31 Test_Results/net_diag_parallel_alpha_4
endif

# NSE initial abundance test
if ($argv[1] == 40 || $argv[1] == 41) then
  cat test_settings_nse Test_Problems/setup_nse_nup >! control
  echo 'Test NSE: Core-Collapse SN with nu-p process network'
  ../source/xnet_nse
  mv -f net_diag01 Test_Results/net_diag_nse_ccsn_nup
endif

# Self-heating test using alpha (explosive burning of degenerate C/O)
if ($argv[1] == 50 || $argv[1] == 51) then
  cat test_settings_heat Test_Problems/setup_heat_alpha >! control
  echo 'Test: Self-heating from explosive burning of degenerate C/O with alpha network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_heat_alpha
endif

# Self-heating test using torch47 (explosive burning of degenerate C/O)
if ($argv[1] == 50 || $argv[1] == 52) then
  cat test_settings_heat Test_Problems/setup_heat_torch47 >! control
  echo 'Test: Self-heating from explosive burning of degenerate C/O with torch47 network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_heat_torch47
endif

# Self-heating test using SN160 (explosive burning of degenerate C/O)
if ($argv[1] == 50 || $argv[1] == 53) then
  cat test_settings_heat Test_Problems/setup_heat_sn160 >! control
  echo 'Test: Self-heating from explosive burning of degenerate C/O with SN160 network'
  $xnet
  mv -f net_diag01 Test_Results/net_diag_heat_sn160
endif
