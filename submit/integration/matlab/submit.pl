use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

$data_path = "/data/biophys/yusipov/os_mbl/int/matlab";
$prefix = "characteristics_pdf";


$Nc = 8;
$dissipator_type = 0;
$alpha = 0.0;
$energy_type = 0;
$bc = 0;
$W = 10.0;
$U = 1.0;
$J = 1.0;
$g = 0.1;
$seed = 0;
$init_state_type = 0;
$init_state_id = 50;
$dump_type = 1;
$begin_dump = 0.1;
$end_dump = 10000;
$num_dumps = 250;
$save_type = 0;
$file_system_type = 0;

$seed_begin = 1;
$seed_num = 100;

$alpha_str = sprintf("%0.4f", $alpha);
$W_str = sprintf("%.4f", $W);
$U_str = sprintf("%.4f", $U);
$J_str = sprintf("%.4f", $J);
$g_str = sprintf("%.4f", $g);

sub ForderName{
	$key_str = $_[0];
	
	return "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}/bc_${bc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/is_${init_state_type}_${init_state_id}/dump_${dump_type}/seed_${key_str}";
}

mkdir "$data_path/$prefix";
mkdir "$data_path/$prefix/Nc_${Nc}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}/bc_${bc}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}/bc_${bc}/W_${W_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}/bc_${bc}/W_${W_str}/U_${U_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}/bc_${bc}/W_${W_str}/U_${U_str}/J_${J_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}/bc_${bc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}/bc_${bc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/is_${init_state_type}_${init_state_id}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${dissipator_type}/dp_${alpha_str}/et_${energy_type}/bc_${bc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/is_${init_state_type}_${init_state_id}/dump_${dump_type}";


for($val = $seed_begin; $val < $seed_begin + $seed_num; $val+=1)
{
	$exp{ForderName($i)} = $val;
	$i++;
}

for($seed = $seed_begin; $seed < $seed_begin + $seed_num; $seed++)
{	
	$key = ForderName($seed);    
	mkdir "$key";
	
	open( WF,">$key/config.txt");
	print WF "$Nc\n"; 
	print WF "$dissipator_type\n";  
	print WF "$alpha\n"; 
	print WF "$energy_type\n";
	print WF "$bc\n";
	print WF "$W\n";
	print WF "$U\n";
	print WF "$J\n";
	print WF "$g\n"; 
	print WF "$seed\n";
	print WF "$init_state_type\n";
	print WF "$init_state_id\n";
	print WF "$dump_type\n";
	print WF "$begin_dump\n";
	print WF "$end_dump\n";
	print WF "$num_dumps\n";
	print WF "$save_type\n";
	print WF "$file_system_type\n";
	close WF;
	
	$test_file = sprintf("random_energies_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_gamma(%0.4f)_ist(%d)_iss(%d)_dump_type(%d)_seed(%d).txt", $Nc, $dissipator_type, $alpha, $energy_type, $bc, $W, $U, $J, $g, $init_state_type, $init_state_id, $dump_type, $seed);
	
	unless (-e "$key/$test_file") 
	{
		print "qsub -wd $dir run.sh $key\n";
		system "qsub -wd $dir run.sh $key";
	}
}

