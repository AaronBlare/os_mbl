use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

$W = 10;
$U = 1;
$J = 1;
$g = 0.1;

$Nc = 8;

$diss_type = 1;

$init_state_type = 1;
$init_state_id = 50;

$num_decade_dumps = 50;
$begin_decade = -1;
$end_decade = 4;

$start = 6;
$num_seeds = 1;

for($g = 0.1; $g <= 1.0; $g += 0.1)
{

	$W_str = sprintf("%.2f", $W);
	$U_str = sprintf("%.2f", $U);
	$J_str = sprintf("%.2f", $J);
	$g_str = sprintf("%.4f", $g);

	sub ForderName{
		$key_str = $_[0];
		return "results/dt_${diss_type}/Ns_${Nc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/init_type_${init_state_type}/init_id_${init_state_id}/seed_${key_str}";
	}

	mkdir "results";
	mkdir "results/dt_${diss_type}";
	mkdir "results/dt_${diss_type}/Ns_${Nc}";
	mkdir "results/dt_${diss_type}/Ns_${Nc}/W_${W_str}";
	mkdir "results/dt_${diss_type}/Ns_${Nc}/W_${W_str}/U_${U_str}";
	mkdir "results/dt_${diss_type}/Ns_${Nc}/W_${W_str}/U_${U_str}/J_${J_str}";
	mkdir "results/dt_${diss_type}/Ns_${Nc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}";
	mkdir "results/dt_${diss_type}/Ns_${Nc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/init_type_${init_state_type}";
	mkdir "results/dt_${diss_type}/Ns_${Nc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/init_type_${init_state_type}/init_id_${init_state_id}";


	for($val = $start; $val < $start + $num_seeds; $val+=1)
	{
		$exp{ForderName($i)} = $val;
		$i++;
	}


	for($seed = $start; $seed < $start + $num_seeds; $seed++)
	{	
		$key = ForderName($seed);    
		mkdir "$key";
		
		open( WF,">$key/config.txt");
		print WF "$W\n"; 
		print WF "$U\n";  
		print WF "$J\n"; 
		print WF "$g\n";
		print WF "$Nc\n";
		print WF "$diss_type\n";
		print WF "$seed\n";
		print WF "$init_state_type\n";
		print WF "$init_state_id\n"; 
		print WF "$num_decade_dumps\n";
		print WF "$begin_decade\n";
		print WF "$end_decade";
		close WF;
		
		$time_file = sprintf("times_Nc(%d)_dt(%d)_W(%.2f)_U(%.2f)_J(%.2f)_gamma(%.2f)_it(%d)_is(%d).txt", $Nc, $diss_type, $W, $U, $J, $g, $init_state_type, $init_state_id);
		
		unless (-e "$key/$time_file") 
		{
			print "qsub -wd $dir script.sh $key\n";
			system "qsub -wd $dir script.sh $key";
		}
	}

}

