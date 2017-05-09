use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
$dir = getcwd;

$Nc = 14;
$W = 1.0;
$U = 1.0;
$J = 1.0;
$g = 0.1;
$dt = 1;
$alpha = 0;
$et = 0;
$init = 0;

$start_seed = 1;
$num_seeds = 100;

$max_num_seeds = 1000000;

$dump = 0;

$W_str = sprintf("%.4f", $W);
$U_str = sprintf("%.4f", $U);
$J_str = sprintf("%.4f", $J);
$g_str = sprintf("%.4f", $g);
$alpha_str = sprintf("%.4f", $alpha);


sub ForderName{
	$key_str = $_[0];
	return "results/Nc_${Nc}/dt_${dt}/alpha_${alpha_str}/et_${et}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/max_num_seeds_${max_num_seeds}/seed_${key_str}";
}

mkdir "results";
mkdir "results/Nc_${Nc}";
mkdir "results/Nc_${Nc}/dt_${dt}";
mkdir "results/Nc_${Nc}/dt_${dt}/alpha_${alpha_str}";
mkdir "results/Nc_${Nc}/dt_${dt}/alpha_${alpha_str}/et_${et}";
mkdir "results/Nc_${Nc}/dt_${dt}/alpha_${alpha_str}/et_${et}/W_${W_str}";
mkdir "results/Nc_${Nc}/dt_${dt}/alpha_${alpha_str}/et_${et}/W_${W_str}/U_${U_str}";
mkdir "results/Nc_${Nc}/dt_${dt}/alpha_${alpha_str}/et_${et}/W_${W_str}/U_${U_str}/J_${J_str}";
mkdir "results/Nc_${Nc}/dt_${dt}/alpha_${alpha_str}/et_${et}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}";
mkdir "results/Nc_${Nc}/dt_${dt}/alpha_${alpha_str}/et_${et}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/max_num_seeds_${max_num_seeds}";



for($val = $start_seed; $val < $start_seed + $num_seeds; $val+=1)
{
	$exp{ForderName($i)} = $val;
	$i++;
}


for($seed = $start_seed; $seed < $start_seed + $num_seeds; $seed++)
{	
	$key = ForderName($seed);    
	mkdir "$key";
	
	open( WF,">$key/config.txt");

	print WF "Nc $Nc\n";
	print WF "W $W\n";
	print WF "U $U\n";
	print WF "J $J\n";
	print WF "g $g\n";
	print WF "dt $dt\n";
	print WF "alpha $alpha\n";
	print WF "et $et\n";
	print WF "init aaa.txt \n";
	print WF "seed $seed\n";
	print WF "max_num_seeds $max_num_seeds\n";
	print WF "dump $dump\n";

	close WF;
	
	$test_file = sprintf('characteristics_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', $Nc, $dt, $alpha, $et, $W, $U, $J, $g, $max_num_seeds, $seed);
	
	if ($dump == 1)
	{
		$test_file = sprintf('hamiltonian_eg_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_max_num_seeds(%d)_seed(%d).txt', $Nc, $dt, $alpha, $et, $W, $U, $J, $g, $max_num_seeds, $seed);
	}
	
	
	unless (-e "$key/$test_file") 
	{
		print "qsub -wd $dir script.sh $key $dir \n";
		system "qsub -wd $dir script.sh $key $dir";
	}
}

