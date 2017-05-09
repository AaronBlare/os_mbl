use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
$dir = getcwd;

$data_path = "/data/biophys/yusipov/mbl_zero_super_vector";

$Nc                  = 8;
$dissipator_type     = 1;
$alpha               = 0.0;		
$energy_type         = 0;
$border_conditions   = 1;
$W                   = 4;
$U                   = 1;
$J                   = 1;
$g                   = 0.1;
$seed_start          = 1;
$seed_num            = 100;
$save_type           = 0;

$num_seed_start = 1;

$alpha_start = 0.0;
$alpha_shift = pi/10.0;
$alpha_num = 1;

$W_start = 1.0;
$W_shift = 1.0;
$W_num = 1;

for($W = $W_start; $W < $W_start + $W_shift * $W_num; $W += $W_shift)
{

	for($alpha = $alpha_start; $alpha < $alpha_start + $alpha_shift * $alpha_num; $alpha += $alpha_shift)
	{

		$alpha_str = sprintf("%.4f", $alpha);
		$W_str = sprintf("%.4f", $W);
		$U_str = sprintf("%.4f", $U);
		$J_str = sprintf("%.4f", $J);
		$g_str = sprintf("%.4f", $g);

		sub ForderName{
			$key_str = $_[0];
			return "$data_path/results/Nc_${Nc}/dt_${dissipator_type}/alpha_${alpha_str}/et_${energy_type}/bc_${border_conditions}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/seed_start_${key_str}";
		}

		mkdir "$data_path/results";
		mkdir "$data_path/results/Nc_${Nc}";
		mkdir "$data_path/results/Nc_${Nc}/dt_${dissipator_type}";
		mkdir "$data_path/results/Nc_${Nc}/dt_${dissipator_type}/alpha_${alpha_str}";
		mkdir "$data_path/results/Nc_${Nc}/dt_${dissipator_type}/alpha_${alpha_str}/et_${energy_type}";
		mkdir "$data_path/results/Nc_${Nc}/dt_${dissipator_type}/alpha_${alpha_str}/et_${energy_type}/bc_${border_conditions}";
		mkdir "$data_path/results/Nc_${Nc}/dt_${dissipator_type}/alpha_${alpha_str}/et_${energy_type}/bc_${border_conditions}/W_${W_str}";
		mkdir "$data_path/results/Nc_${Nc}/dt_${dissipator_type}/alpha_${alpha_str}/et_${energy_type}/bc_${border_conditions}/W_${W_str}/U_${U_str}";
		mkdir "$data_path/results/Nc_${Nc}/dt_${dissipator_type}/alpha_${alpha_str}/et_${energy_type}/bc_${border_conditions}/W_${W_str}/U_${U_str}/J_${J_str}";
		mkdir "$data_path/results/Nc_${Nc}/dt_${dissipator_type}/alpha_${alpha_str}/et_${energy_type}/bc_${border_conditions}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}";


		for($val = $seed_start; $val < $seed_start + $num_seed_start * $seed_num; $val += $seed_num)
		{
			$exp{ForderName($i)} = $val;
			$i++;
		}


		for($val = $seed_start; $val < $seed_start + $num_seed_start * $seed_num; $val += $seed_num)
		{	
			$key = ForderName($val);    
			mkdir "$key";
			
			open( WF,">$key/config.txt");
			print WF "$Nc\n"; 
			print WF "$dissipator_type\n";  
			print WF "$alpha\n"; 
			print WF "$energy_type\n";
			print WF "$border_conditions\n";
			print WF "$W\n";
			print WF "$U\n";
			print WF "$J\n";
			print WF "$g\n";
			print WF "$val\n";
			print WF "$seed_num\n";
			print WF "$save_type\n";
			close WF;	
			
			$test_file = sprintf('%s/characteristics_Nc(%d)_dt(%d)_alpha(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_ss(%d)_sn(%d)_seed(%d).txt.txt', $key, $Nc, $dissipator_type, $alpha, $energy_type, $border_conditions, $W, $U, $J, $g, $val, $seed_num, $val);
			
			unless (-e "$test_file") 
			{
				print "qsub -wd $dir run.sh $key\n";
				system "qsub -wd $dir run.sh $key";
			}
		}

	}
}

