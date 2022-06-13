integer, parameter :: TS_RESTO=1
integer, parameter :: TS_DYNAMICS=2
integer, parameter :: TS_PHYSICS=3
integer, parameter :: TS_RK_RESTO=4
integer, parameter :: TS_RK_ADV=5
integer, parameter :: TS_RK_ADVMON=6
character(len=8), parameter :: names(6) = (/&
"RESTO   ", &
"DYNAMICS", &
"PHYSICS ", &
"RK_RESTO", &
"RK_ADV  ", &
"RK_ADVMO" /)
