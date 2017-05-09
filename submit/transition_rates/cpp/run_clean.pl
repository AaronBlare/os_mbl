use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
$dir = getcwd;

print "qsub -wd $dir clean.sh \n";
system "qsub -wd $dir clean.sh";


