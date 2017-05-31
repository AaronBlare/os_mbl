use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

$data_path = "/data/biophys/yusipov/os_mbl/int/matlab";
$prefix = "characteristics";

$Nc = 8;
$diss_type = 1;
$diss_phase = 0.0;
$energy_type = 0;
$periodic_bc = 0;
$W = 10.0;
$U = 1.0;
$J = 1.0;
$g = 0.1;
$seed_start = 0;
$seed_num = 100;
$is_int = 0;
$int_ist = 0;
$int_isi = 50;
$int_dt = 1;
$int_db = 0.1;
$int_de = 10000;
$int_dn = 250;
$is_zev = 1;
$is_zev_check = 1;
$is_eg_stat = 1;
$is_save_vec = 1;
$is_save_mtx = 0
$fs_type = 0;

$seed_start_begin = 1;
$seed_start_num = 1000;

$diss_phase_str = sprintf("%0.4f", $diss_phase);
$W_str = sprintf("%.4f", $W);
$U_str = sprintf("%.4f", $U);
$J_str = sprintf("%.4f", $J);
$g_str = sprintf("%.4f", $g);

sub ForderName{
	$key_str = $_[0];
	
	return "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}/dp_${diss_phase_str}/et_${energy_type}/bc_${periodic_bc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}/seed_start${key_str}";
}

mkdir "$data_path/$prefix";
mkdir "$data_path/$prefix/Nc_${Nc}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}/dp_${diss_phase_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}/dp_${diss_phase_str}/et_${energy_type}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}/dp_${diss_phase_str}/et_${energy_type}/bc_${periodic_bc}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}/dp_${diss_phase_str}/et_${energy_type}/bc_${periodic_bc}/W_${W_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}/dp_${diss_phase_str}/et_${energy_type}/bc_${periodic_bc}/W_${W_str}/U_${U_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}/dp_${diss_phase_str}/et_${energy_type}/bc_${periodic_bc}/W_${W_str}/U_${U_str}/J_${J_str}";
mkdir "$data_path/$prefix/Nc_${Nc}/dt_${diss_type}/dp_${diss_phase_str}/et_${energy_type}/bc_${periodic_bc}/W_${W_str}/U_${U_str}/J_${J_str}/g_${g_str}";

for($val = $seed_start_begin; $val < $seed_start_begin + $seed_start_num * $seed_num; $val+=$seed_num)
{
	$exp{ForderName($i)} = $val;
	$i++;
}

for($seed_start = $seed_begin; $seed_start < $seed_start_begin + $seed_start_num * $seed_num; $seed_start+=$seed_num)
{	
	$key = ForderName($seed);    
	mkdir "$key";
	
	open( WF,">$key/config.txt");
	print WF "$Nc\n"; 
	print WF "$diss_type\n";  
	print WF "$diss_phase\n"; 
	print WF "$energy_type\n";
	print WF "$periodic_bc\n";
	print WF "$W\n";
	print WF "$U\n";
	print WF "$J\n";
	print WF "$g\n"; 
	print WF "$seed_start\n";
	print WF "$seed_num\n";
	print WF "$is_int\n";
	print WF "$int_ist\n";
	print WF "$int_isi\n";
	print WF "$int_dt\n";
	print WF "$int_db\n";
	print WF "$int_de\n";
	print WF "$int_dn\n";
	print WF "$is_zev\n";
	print WF "$is_zev_check\n";
	print WF "$is_eg_stat\n";
	print WF "$is_save_vec\n";
	print WF "$is_save_mtx\n";
	print WF "$fs_type\n";
	close WF;
	
	$test_file = sprintf("random_energies_Nc(%d)_dt(%d)_dp(%0.4f)_et(%d)_bc(%d)_W(%0.4f)_U(%0.4f)_J(%0.4f)_g(%0.4f)_seed(%d).txt", $Nc, $diss_type, $diss_phase, $energy_type, $periodic_bc, $W, $U, $J, $g, $seed_start + $seed_num - 1);
	
	unless (-e "$key/$test_file") 
	{
		print "qsub -wd $dir run.sh $key\n";
		system "qsub -wd $dir run.sh $key";
	}
}

