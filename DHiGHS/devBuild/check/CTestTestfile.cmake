# CMake generated Testfile for 
# Source directory: /home/edeakins/LP/EQLPSolver/DHiGHS/check
# Build directory: /home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/check
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(unit-test-build "/usr/bin/cmake" "--build" "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild" "--target" "unit_tests")
set_tests_properties(unit-test-build PROPERTIES  RESOURCE_LOCK "unittestbin" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;65;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(unit_tests_all "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/unit_tests" "--success")
set_tests_properties(unit_tests_all PROPERTIES  DEPENDS "unit-test-build" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;78;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(25fv47--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/25fv47.mps")
set_tests_properties(25fv47--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 2888
Objective value     :  5.501846" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(25fv47--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/25fv47.mps")
set_tests_properties(25fv47--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  5.501846" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(80bau3b--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/80bau3b.mps")
set_tests_properties(80bau3b--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 3760
Objective value     :  9.872242" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(80bau3b--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/80bau3b.mps")
set_tests_properties(80bau3b--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  9.872242" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(adlittle--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/adlittle.mps")
set_tests_properties(adlittle--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 81
Objective value     :  2.254950" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(adlittle--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/adlittle.mps")
set_tests_properties(adlittle--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  2.254950" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(afiro--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/afiro.mps")
set_tests_properties(afiro--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 21
Objective value     : -4.647531" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(afiro--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/afiro.mps")
set_tests_properties(afiro--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -4.647531" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(etamacro--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/etamacro.mps")
set_tests_properties(etamacro--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 507
Objective value     : -7.557152" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(etamacro--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/etamacro.mps")
set_tests_properties(etamacro--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.557152" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(greenbea--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/greenbea.mps")
set_tests_properties(greenbea--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 5249
Objective value     : -7.255525" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(greenbea--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/greenbea.mps")
set_tests_properties(greenbea--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -7.255525" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(pilot87--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/pilot87.mps")
set_tests_properties(pilot87--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 9919
Objective value     :  3.017103" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(pilot87--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/pilot87.mps")
set_tests_properties(pilot87--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  3.017103" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(shell--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/shell.mps")
set_tests_properties(shell--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 624
Objective value     :  1.208825" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(shell--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/shell.mps")
set_tests_properties(shell--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.208825" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(stair--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/stair.mps")
set_tests_properties(stair--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 555
Objective value     : -2.512670" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(stair--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/stair.mps")
set_tests_properties(stair--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     : -2.512670" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(standata--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/standata.mps")
set_tests_properties(standata--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 74
Objective value     :  1.257699" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(standata--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/standata.mps")
set_tests_properties(standata--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.257699" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(standgub--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/standgub.mps")
set_tests_properties(standgub--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 67
Objective value     :  1.257699" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(standgub--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/standgub.mps")
set_tests_properties(standgub--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.257699" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(standmps--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/standmps.mps")
set_tests_properties(standmps--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Simplex   iterations: 215
Objective value     :  1.406017" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(standmps--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/standmps.mps")
set_tests_properties(standmps--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Objective value     :  1.406017" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;222;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(bgetam--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/bgetam.mps")
set_tests_properties(bgetam--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(bgetam--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/bgetam.mps")
set_tests_properties(bgetam--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(box1--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/box1.mps")
set_tests_properties(box1--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(box1--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/box1.mps")
set_tests_properties(box1--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(ex72a--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/ex72a.mps")
set_tests_properties(ex72a--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(ex72a--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/ex72a.mps")
set_tests_properties(ex72a--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(forest6--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/forest6.mps")
set_tests_properties(forest6--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(forest6--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/forest6.mps")
set_tests_properties(forest6--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(galenet--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/galenet.mps")
set_tests_properties(galenet--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(galenet--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/galenet.mps")
set_tests_properties(galenet--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(gams10am--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/gams10am.mps")
set_tests_properties(gams10am--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(gams10am--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/gams10am.mps")
set_tests_properties(gams10am--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(klein1--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/klein1.mps")
set_tests_properties(klein1--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(klein1--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/klein1.mps")
set_tests_properties(klein1--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(woodinfe--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/woodinfe.mps")
set_tests_properties(woodinfe--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(woodinfe--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/woodinfe.mps")
set_tests_properties(woodinfe--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;224;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(gas11--presolve=off "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=off" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/gas11.mps")
set_tests_properties(gas11--presolve=off PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Unbounded" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;225;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(gas11--presolve=on "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/gas11.mps")
set_tests_properties(gas11--presolve=on PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Unbounded" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;197;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;225;add_instancetests;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(presolve-adlittle "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/adlittle.mps")
set_tests_properties(presolve-adlittle PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Optimal" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;235;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(presolve-scrs8 "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/scrs8.mps")
set_tests_properties(presolve-scrs8 PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Optimal" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;240;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(presolve-woodinfe "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/woodinfe.mps")
set_tests_properties(presolve-woodinfe PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Infeasible" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;246;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")
add_test(presolve-gas11 "/home/edeakins/LP/EQLPSolver/DHiGHS/devBuild/bin/highs" "--presolve=on" "/home/edeakins/LP/EQLPSolver/DHiGHS/check/instances/gas11.mps")
set_tests_properties(presolve-gas11 PROPERTIES  DEPENDS "unit_tests_all" PASS_REGULAR_EXPRESSION "Model   status      : Unbounded" _BACKTRACE_TRIPLES "/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;252;add_test;/home/edeakins/LP/EQLPSolver/DHiGHS/check/CMakeLists.txt;0;")