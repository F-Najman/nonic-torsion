SetColumns(0);SetQuitOnError(true);
load "common/data/dcm_list_all_min_orders_X1_10k.m";
a:=dcm_list_all_min_orders_X1[133];assert a[1] eq 133 and a[2] eq 36;
printf "NONIC133_CM_MINIMUM level=133 minimum_degree=36 exceeds9=true complete=true\n";quit;
