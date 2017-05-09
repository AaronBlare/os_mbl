use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
$dir = getcwd;

$W = 8;
$U = 1;
$J = 1;
$g = 0.1;

$Nc = 14;

$diss_type = 1;
$alpha = 0;
$energy_type = 0;

$start = 1;
$num_seeds = 1000;

$dump_aux = 1;

$W_str = sprintf("%.4f", $W);
$U_str = sprintf("%.4f", $U);
$J_str = sprintf("%.4f", $J);
$g_str = sprintf("%.4f", $g);
$alpha_str = sprintf("%.4f", $alpha);


sub ForderName{
	$key_str = $_[0];
	return "results/dt_${diss_type}/alpha_${alpha_str}/et_${energy_type}/Ns_${Nc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/seed_${key_str}";
}

mkdir "results";
mkdir "results/dt_${diss_type}";
mkdir "results/dt_${diss_type}/alpha_${alpha_str}";
mkdir "results/dt_${diss_type}/alpha_${alpha_str}/et_${energy_type}";
mkdir "results/dt_${diss_type}/alpha_${alpha_str}/et_${energy_type}/Ns_${Nc}";
mkdir "results/dt_${diss_type}/alpha_${alpha_str}/et_${energy_type}/Ns_${Nc}/W_${W_str}";
mkdir "results/dt_${diss_type}/alpha_${alpha_str}/et_${energy_type}/Ns_${Nc}/W_${W_str}/U_${U_str}";
mkdir "results/dt_${diss_type}/alpha_${alpha_str}/et_${energy_type}/Ns_${Nc}/W_${W_str}/U_${U_str}/J_${J_str}";
mkdir "results/dt_${diss_type}/alpha_${alpha_str}/et_${energy_type}/Ns_${Nc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}";



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
	print WF "$alpha\n";
	print WF "$energy_type\n";
	print WF "$seed\n";
	print WF "$dump_aux\n";
	close WF;
	
	$test_file = sprintf('random_energies_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_seed(%d).txt', $Nc, $diss_type, $alpha, $energy_type, $W, $U, $J, $g, $seed);

	
	unless (-e "$key/$test_file") 
	{
		print "qsub -wd $dir script.sh $key\n";
		system "qsub -wd $dir script.sh $key";
	}
}

