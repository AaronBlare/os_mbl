use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

$period = 0.1;

$Nc = 8;
$num_periods = 100000;
$init_state_id = 49;
$num_dumps = 251;
$dump_type = 1;
$num_periods_in_trans_proc = 0;
$num_omp_threads = 16;
$num_trajectories = 10000;
$rnd_max = 2000000;
$rnd_cur = 0;
$calc_characteristics = 1;
$dump_rho = 0;

$start = 0;
$finish = 1;
%exp = ();
$i = 0;
$i = $start;

$dt = 1;
$Ns = 70;
$W = 10.0;
$U = 1.0;
$J = 1.0;
$g = 0.1;

for($seed = 1; $seed <= 100; $seed+=1)
{
	$input_file_name = sprintf("input_data_period%0.2f_dt%d_Ns%d_W%0.2f_U%0.2f_J%0.2f_g%0.2f_seed%d.bin", $period, $dt, $Ns, $W, $U, $J, $g, $seed);
	$aux_file_name = sprintf("aux_data_period%0.2f_dt%d_Ns%d_W%0.2f_U%0.2f_J%0.2f_g%0.2f_seed%d.bin", $period, $dt, $Ns, $W, $U, $J, $g, $seed);

	$period_str = sprintf("%.2f", $period);
	$W_str = sprintf("%.2f", $W);
	$U_str = sprintf("%.2f", $U);
	$J_str = sprintf("%.2f", $J);
	$g_str = sprintf("%.2f", $g);

	sub ForderName{
		$key_str = $_[0];
		return  "results/period_${period_str}/dt_${dt}/Ns_${Ns}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/is_${init_state_id}/seed_${seed}/dump_${dump_type}/rnd_${key_str}";
	}

	mkdir "results";
	mkdir "results/period_${period_str}";
	mkdir "results/period_${period_str}/dt_${dt}";
	mkdir "results/period_${period_str}/dt_${dt}/Ns_${Ns}";
	mkdir "results/period_${period_str}/dt_${dt}/Ns_${Ns}/W_${W_str}";
	mkdir "results/period_${period_str}/dt_${dt}/Ns_${Ns}/W_${W_str}/U_${U_str}";
	mkdir "results/period_${period_str}/dt_${dt}/Ns_${Ns}/W_${W_str}/U_${U_str}/J_${J_str}";
	mkdir "results/period_${period_str}/dt_${dt}/Ns_${Ns}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}";
	mkdir "results/period_${period_str}/dt_${dt}/Ns_${Ns}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/is_${init_state_id}";
	mkdir "results/period_${period_str}/dt_${dt}/Ns_${Ns}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/is_${init_state_id}/seed_${seed}";
	mkdir "results/period_${period_str}/dt_${dt}/Ns_${Ns}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/is_${init_state_id}/seed_${seed}/dump_${dump_type}";

	for($val = $start; $val < $finish; $val+=1)
	{
		$exp{ForderName($i)} = $val;
		$i++;
	}


	for($i = $start; $i < $finish; $i++)
	{
		$key = ForderName($i);    
		mkdir "$key";
		
		$rnd = $rnd_cur + $exp{$key} * $num_trajectories;
		
		open( WF,">$key/config.txt");
		print WF "Nc = $Nc \n"; 
		print WF "num_periods = $num_periods \n"; 
		print WF "init_state_id = $init_state_id \n";
		print WF "num_dumps = $num_dumps \n";
		print WF "dump_type = $dump_type \n";
		print WF "num_periods_in_trans_proc = $num_periods_in_trans_proc \n";
		print WF "num_omp_threads = $num_omp_threads \n";
		print WF "num_trajectories = $num_trajectories \n"; 
		print WF "rnd_max = $rnd_max \n";
		print WF "rnd_cur = $rnd \n";
		print WF "calc_characteristics = $calc_characteristics \n";
		print WF "dump_rho = $dump_rho";
		close WF;

		print $key."\n";
		
		$test_file = "characteristics.txt";
		
		unless (-e "$key/$test_file")
		{	
			print "qsub -wd $dir script_multi_thread.sh $key $dir $input_file_name $aux_file_name \n";
			system "qsub -wd $dir script_multi_thread.sh $key $dir $input_file_name $aux_file_name ";
		}
	
	}
}

