#!/usr/bin/perl
use Cwd;

@dirs = glob("*");
$cwd = Cwd::getcwd();
foreach $dir (@dirs) {
    if (-d $dir) {
        print("$dir\n");
        chdir($dir);
        chdir("Q_val_dist");
        &smear (0.2);
        # chdir("../..");
        chdir($cwd);
    }
}

sub smear {
    my ($a) = @_;
    my (@files,@name,$j,@e,@P);
    my ($x,$n,$f,$log,$pi);

    $pi = 3.141592653589793;
    @name = @files;
    @files = glob("Q*.dat");
    foreach (@files) {
        /(.*)\.dat/;
        open(in,"$1.dat");
        $j = 0;
        while (<in>) {
            /[\s\d]\S+/;
            @e[$j] = $&;
            @P[$j] = $';
            ++$j;
        }
        close(in);        
        $n = @e;
        open(out,">smear_$1.dat");
        $j = 0;
        $x = -2.0;
        $dx = 0.01;
        while ($x <= @e[$n-1]) {
           $f = 0.0 ;
           for($j=0; $j<=$n-1; ++$j) {
#              $f += $a*@P[$j] / (($x-@e[$j])**2 + $a*$a)/$pi;
               $f += @P[$j] / (sqrt(2.0*$pi) * $a)
                   * exp(-($x-@e[$j])**2/(2.0*$a*$a));
           }
           print out "$x  $f\n";
           $x += $dx;
           ++$j;
       }
       close(out);
   }
}

