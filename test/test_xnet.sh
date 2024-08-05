#!/bin/bash -f
###############################################################################
## Regression Testing for XNet
## This script will execute pre-defined test problems and diff the resulting
## net_diag files with previously accepted output.
##
## Usage:
##   ./test_xnet.sh <xnet executables to test> <test problem ID #'s>
##
## Description:
##   <xnet executables to test> (Optional)      default: ../source/xnetp
##     Paths to XNet exeuctables to test.
##     This is useful for testing custom builds of that do not correspond to a specific Makefile target.
##     Currently, this only works with serial tests.
##
##   <test problem ID #'s>                      default: none
##     Any combination of the test problem numbers given in the table below.
##       =====================================================================
##       Description                    Network (# species)             ID #
##       =====================================================================
##       All serial test problems                                         0
##
##       Thermonuclear SN                 alpha (14)                      1
##                                        torch47 (47)                    2
##                                        SN160 (160)                     3
##
##       Tidally-Crushed WD               alpha (14)                      4
##                                        SN160 (160)                     5
##
##       Nova                             CNO (16)                        7
##                                        nova169 (169)                   8
##
##       X-ray burst                      CNO (16)                        9
##                                        Fisker (304)                    10
##
##       Core Collapse SN                 alpha (16)                      11
##                                        SN160 (160)                     12
##                                        nup (1072)                      13
##
##       Neutron Star Wind                alpha (16)                      15
##                                        Reaclib20180621 (7852)          16
##
##       Self-heating C/O burn                                            50
##                                        alpha (16)                      51
##                                        torch47 (47)                    52
##                                        SN160 (160)                     53
##                                        SN160 (160) (BDF integrator)    54
##                                        SN231 (231) (log-ft rates)      55
##       =====================================================================
##       All parallel test problems                                       30
##       4 different SN                   alpha (16)                      31
##                                        SN160 (160)                     32
##       =====================================================================
##       All NSE initialization problems                                  40
##
##       CCSN                             nup (1072)                      41
##
###############################################################################

xnetd=../source/xnetd
xnetm=../source/xnetm
xnetp=../source/xnetp
xnet_nse=../source/xnet_nse
xnetd_mpi=../source/xnetd_mpi
xnetm_mpi=../source/xnetm_mpi
xnetp_mpi=../source/xnetp_mpi
xnet_nse_mpi=../source/xnet_nse_mpi
xnse=../source/xnse

xnet_mpi=$xnetp_mpi

function test_diff {
  # Remove timers from files for diff
  sed -e '/^Timers Summary:/,+14d' $1 >| diff1.txt
  sed -e '/^Timers Summary:/,+14d' $2 >| diff2.txt
  if ! diff -q diff1.txt diff2.txt >/dev/null 2>&1; then
    echo "Warning: $1 differs from $2"
    echo "File diff_$3 contains diff output"
    diff diff1.txt diff2.txt >| diff_$3
    #read -rsp $'Press any key to continue...\n' -n1 key
  fi
  rm diff1.txt diff2.txt
}

function do_test {
  cat test_settings Test_Problems/setup_$2 >| control
  $1
  mv -f net_diag01 Test_Results/net_diag_$2
  test_diff Test_Results/net_diag_$2 Test_Problems/Results/net_diag_$2 $2
}

function do_test_small {
  cat test_settings_small Test_Problems/setup_$2 >| control
  $1
  mv -f net_diag01 Test_Results/net_diag_$2
  test_diff Test_Results/net_diag_$2 Test_Problems/Results/net_diag_$2 $2
}

function do_test_heat {
  cat test_settings_heat Test_Problems/setup_$2 >| control
  $1
  mv -f net_diag01 Test_Results/net_diag_$2
  test_diff Test_Results/net_diag_$2 Test_Problems/Results/net_diag_$2 $2
}

function do_test_bdf {
  cat test_settings_bdf Test_Problems/setup_$2 >| control
  $1
  mv -f net_diag01 Test_Results/net_diag_$2
  test_diff Test_Results/net_diag_$2 Test_Problems/Results/net_diag_$2 $2
}

function do_test_logft {
  cat test_settings_logft Test_Problems/setup_$2 >| control
  $1
  mv -f net_diag01 Test_Results/net_diag_$2
  test_diff Test_Results/net_diag_$2 Test_Problems/Results/net_diag_$2 $2
}

function do_test_batch {
  cat test_settings_batch Test_Problems/setup_$2 >| control
  $1
  mv -f net_diag01 Test_Results/net_diag_$2
  test_diff Test_Results/net_diag_$2 Test_Problems/Results/net_diag_$2 $2
}

function do_test_parallel {
  cat test_settings_parallel Test_Problems/setup_$2 >| control
  mpirun -n 4 $1
  mv -f net_diag0* Test_Results/net_diag_${2}_1
  mv -f net_diag1* Test_Results/net_diag_${2}_2
  mv -f net_diag2* Test_Results/net_diag_${2}_3
  mv -f net_diag3* Test_Results/net_diag_${2}_4
  test_diff Test_Results/net_diag_${2}_1 Test_Problems/Results/net_diag_${2}_1 ${2}_1
  test_diff Test_Results/net_diag_${2}_2 Test_Problems/Results/net_diag_${2}_2 ${2}_2
  test_diff Test_Results/net_diag_${2}_3 Test_Problems/Results/net_diag_${2}_3 ${2}_3
  test_diff Test_Results/net_diag_${2}_4 Test_Problems/Results/net_diag_${2}_4 ${2}_4
}

function do_test_nse {
  cat test_settings_nse Test_Problems/setup_nse_$2 >| control
  $1
  mv -f net_diag01 Test_Results/net_diag_nse_$2
  test_diff Test_Results/net_diag_nse_$2 Test_Problems/Results/net_diag_nse_$2 nse_$2
}

function do_test_xnse {
  cat test_settings_xnse Test_Problems/setup_xnse_$2 >| control
  $1 < Test_Problems/input_xnse_test > xnse_$2.out
  mv -f nse_diag01 Test_Results/nse_diag_$2
  if [ -f Test_Problems/Results/nse_diag_$2 ]; then
    test_diff Test_Results/nse_diag_$2 Test_Problems/Results/nse_diag_$2 nse_$2
  else
    cp -v Test_Results/nse_diag_$2 Test_Problems/Results/nse_diag_$2
  fi
}

xnet_list=()
# Use user-supplied executable if provided as argument
for arg in $*; do
  if [ -f $arg ]; then
    xnet_list+=($arg)
  elif [ "$arg" -eq "$arg" ] 2>/dev/null; then
    test_list+=($arg)
  fi
done
# If no executable specified, use serial PARDISO
if [ ${#xnet_list[@]} -lt 1 ]; then
  xnet_list+=($xnetp)
fi

mkdir -pv Test_Results
mkdir -pv Test_Problems/Results

n=1
for xnet in ${xnet_list[@]}; do

  echo "Testing $xnet"

  for itest in ${test_list[@]}; do

    # TN SN tracer, from Ed Brown, with alpha network 
    if [ $itest -eq 0 -o $itest -eq 1 ]; then
      echo "Test: Thermonuclear SN with alpha network"
      test_th="tnsn"; test_net="alpha"; test_name=${test_th}_${test_net}
      do_test_small $xnet $test_name
    fi

    # TN SN tracer, from Ed Brown, with 47 species network 
    if [ $itest -eq 0 -o $itest -eq 2 ]; then
      echo "Test: Thermonuclear SN with 47 species network"
      test_th="tnsn"; test_net="torch47"; test_name=${test_th}_${test_net}
      do_test $xnet $test_name
    fi

    # TN SN tracer, from Ed Brown, with 160 species network 
    if [ $itest -eq 0 -o $itest -eq 3 ]; then
      echo "Test: Thermonuclear SN with 160 species network"
      test_th="tnsn"; test_net="sn160"; test_name=${test_th}_${test_net}
      do_test $xnet $test_name
    fi

    # TI SN tracer, from Stephan Rosswog, with alpha network 
    if [ $itest -eq 0 -o $itest -eq 4 ]; then
      echo "Test: Tidally Induced SN with alpha network"
      test_th="tisn"; test_net="alpha"; test_name=${test_th}_${test_net}
      do_test_small $xnet $test_name
    fi

    # TI SN tracer, from Stephan Rosswog, with 160 species network 
    if [ $itest -eq 0 -o $itest -eq 5 ]; then
      echo "Test: Tidally Induced SN with 160 species network"
      test_th="tisn"; test_net="sn160"; test_name=${test_th}_${test_net}
      do_test $xnet $test_name
    fi

    # Nova zone, from Sumner Starrfield, with 16 species CNO network 
    if [ $itest -eq 0 -o $itest -eq 7 ]; then
      echo "Test: Nova with minimal CNO network"
      test_th="nova"; test_net="cno"; test_name=${test_th}_${test_net}
      do_test_small $xnet $test_name
    fi

    # Nova zone, from Sumner Starrfield, with 169 species network from Chritian Iliadis
    if [ $itest -eq 0 -o $itest -eq 8 ]; then
      echo "Test: Nova with 169 species network"
      test_th="nova"; test_net="Iliadis"; test_name=${test_th}_${test_net}
      do_test $xnet $test_name
    fi

    # XRB tracer, from Jacob Fisker, with minimal CNO network 
    if [ $itest -eq 0 -o $itest -eq 9 ]; then
      echo "Test: X-ray burst with CNO network"
      test_th="xrb"; test_net="cno"; test_name=${test_th}_${test_net}
      do_test_small $xnet $test_name
    fi

    # XRB tracer, from Jacob Fisker, with 304 species rp-process network 
    if [ $itest -eq 0 -o $itest -eq 10 ]; then
      echo "Test: X-ray burst with Fiskers network"
      test_th="xrb"; test_net="fisker"; test_name=${test_th}_${test_net}
      do_test $xnet $test_name
    fi

    # CC SN zone, from Carla Froehlich, with alpha network 
    if [ $itest -eq 0 -o $itest -eq 11 ]; then
      echo "Test: Core-Collapse SN with alpha network"
      test_th="ccsn"; test_net="alpha"; test_name=${test_th}_${test_net}
      do_test_small $xnet $test_name
    fi

    #  CC SN zone, from Carla Froehlich with 160 species network
    if [ $itest -eq 0 -o $itest -eq 12 ]; then
      echo "Test: Core-Collapse SN with 160 species network"
      test_th="ccsn"; test_net="sn160"; test_name=${test_th}_${test_net}
      do_test $xnet $test_name
    fi

    #  CC SN zone, from Carla Froehlich with nu p-process network 
    if [ $itest -eq 0 -o $itest -eq 13 ]; then
      echo "Test: Core-Collapse SN with nu-p process network"
      test_th="ccsn"; test_net="nup"; test_name=${test_th}_${test_net}
      do_test $xnet $test_name
    fi

    # Neutino driven wind example, from Josh Beun, with alpha network 
    if [ $itest -eq 0 -o $itest -eq 15 ]; then
      echo "Test: Neutrino-driven wind with alpha network"
      test_th="nuwind"; test_net="alpha"; test_name=${test_th}_${test_net}
      do_test_small $xnet $test_name
    fi

    # Neutino driven wind (Meyer & Brown) with full JINA REACLIB network
    if [ $itest -eq 0 -o $itest -eq 16 ]; then
      echo "Test: Neutrino-driven wind with 7852 species network"
      test_th="nuwind"; test_net="rprocess"; test_name=${test_th}_${test_net}
      do_test $xnet $test_name
    fi

    # Self-heating test using alpha (explosive burning of degenerate C/O)
    if [ $itest -eq 50 -o $itest -eq 51 ]; then
      echo "Test: Self-heating from explosive burning of degenerate C/O with alpha network"
      test_th="heat"; test_net="alpha"; test_name=${test_th}_${test_net}
      do_test_heat $xnet $test_name
    fi

    # Self-heating test using torch47 (explosive burning of degenerate C/O)
    if [ $itest -eq 50 -o $itest -eq 52 ]; then
      echo "Test: Self-heating from explosive burning of degenerate C/O with torch47 network"
      test_th="heat"; test_net="torch47"; test_name=${test_th}_${test_net}
      do_test_heat $xnet $test_name
    fi

    # Self-heating test using SN160 (explosive burning of degenerate C/O)
    if [ $itest -eq 50 -o $itest -eq 53 ]; then
      echo "Test: Self-heating from explosive burning of degenerate C/O with SN160 network"
      test_th="heat"; test_net="sn160"; test_name=${test_th}_${test_net}
      do_test_heat $xnet $test_name
    fi

    # Self-heating test using SN160 (explosive burning of degenerate C/O) (with BDF integrator)
    if [ $itest -eq 50 -o $itest -eq 54 ]; then
      echo "Test: Self-heating from explosive burning of degenerate C/O with SN160 network and BDF integrator"
      test_th="bdf"; test_net="sn160"; test_name=${test_th}_${test_net}
      do_test_bdf $xnet $test_name
    fi

    # Self-heating test using SN160 (explosive burning of degenerate C/O) (with BDF integrator)
    if [ $itest -eq 50 -o $itest -eq 55 ]; then
      echo "Test: Self-heating from explosive burning of degenerate C/O with SN231 network and log(ft) rates"
      test_th="logft"; test_net="sn231"; test_name=${test_th}_${test_net}
      do_test_logft $xnet $test_name
    fi

    # Zone batching test using alpha (explosive burning of degenerate C/O)
    if [ $itest -eq 60 -o $itest -eq 61 ]; then
      echo "Test: Zone batching w/ self-heating from explosive burning of degenerate C/O with alpha network"
      test_th="batch"; test_net="alpha"; test_name=${test_th}_${test_net}
      do_test_batch $xnet $test_name
    fi

    # Zone batching test using torch47 (explosive burning of degenerate C/O)
    if [ $itest -eq 60 -o $itest -eq 62 ]; then
      echo "Test: Zone batching w/ self-heating from explosive burning of degenerate C/O with torch47 network"
      test_th="batch"; test_net="torch47"; test_name=${test_th}_${test_net}
      do_test_batch $xnet $test_name
    fi

  done

done

if [ -f $xnet_mpi ]; then

  for itest in ${test_list[@]}; do

    # Parallel test, runs 4 different zones using SN160 
    if [ $itest -eq 30 -o $itest -eq 31 ]; then
      echo "Test: 4 different zones in parallel"
      test_th="parallel"; test_net="sn160"; test_name=${test_th}_${test_net}
      do_test_parallel $xnet_mpi $test_name
    fi

    # Parallel test, runs 4 different zones using alpha
    if [ $itest -eq 30 -o $itest -eq 32 ]; then
      echo "Test: 4 different zones in parallel"
      test_th="parallel"; test_net="alpha"; test_name=${test_th}_${test_net}
      do_test_parallel $xnet_mpi $test_name
    fi

  done

fi

if [ -f $xnet_nse ]; then

  for itest in ${test_list[@]}; do

    # NSE initial abundance test
    if [ $itest -eq 40 -o $itest -eq 41 ]; then
      echo "Test NSE: Core-Collapse SN with nu-p process network"
      test_th="ccsn"; test_net="nup"; test_name=${test_th}_${test_net}
      do_test_nse $xnet_nse $test_name
    fi

  done

fi

if [ -f $xnse ]; then

  for itest in ${test_list[@]}; do

    # NSE initial abundance test
    if [ $itest -eq 80 -o $itest -eq 81 ]; then
      echo "Test NSE: SN160 network"
      test_net="sn160"; test_name=${test_net}
      do_test_xnse $xnse $test_name
    fi

    # NSE initial abundance test
    if [ $itest -eq 80 -o $itest -eq 82 ]; then
      echo "Test NSE: SN231 network"
      test_net="sn231"; test_name=${test_net}
      do_test_xnse $xnse $test_name
    fi

  done

fi
