#!/gsc/bin/perl

# P R A G M A S ###############################################################
use warnings;
use strict;

# M O D U L E S ###############################################################
use POSIX;
use YAML::Syck;
use Sys::Hostname;
use Time::HiRes qw(gettimeofday tv_interval);
use Path::Class;
use File::Basename;
use Getopt::Long;

# M A I N #####################################################################

$| = 1; # enable AUTOFLUSH mode

# Setup Default options
my $opt_fastq =
  '/gscmnt/sata141/techd/twylie/miRNA_20100507/seq/s_1_sequence.fastq';
my @opt_mismatches;
my $opt_adaptor = 'ATCTCGTATGCCGTCTTCTGCTTGC.*$';
my $opt_no_trim  = 0;
my $opt_no_plot  = 0;
my $opt_read_length_histogram;
my $opt_read_mismatch_count_histogram;
my $fqgrep = '/gscuser/idas/git/fqgrep/fqgrep';

# parse options
my $result = GetOptions(
    "fastq=s"               => \$opt_fastq,
    "mismatches=i{,}"        => \@opt_mismatches,
    "adaptor=s"             => \$opt_adaptor,
    "no-trim"               => \$opt_no_trim,
    "no-plot"               => \$opt_no_plot,
    "read-length-histogram" => \$opt_read_length_histogram,
    "read-mismatch-count-histogram"  => \$opt_read_mismatch_count_histogram
);

# ensure input fastq/a file existance
unless ($opt_fastq) {
    die "[err] Need to pass in a fastq/a file to process via --fastq !\n";
}
my $fastq = Path::Class::File->new($opt_fastq);
die "[err] Could not find fastq $fastq \n" unless(-e $fastq);

# ensure adaptor sequence
unless($opt_adaptor) {
    die "[err] Need to pass an adaptor sequence to trim for via --adaptor!\n";
}

# ensure a mismatch level
unless (@opt_mismatches) {
    @opt_mismatches = qw(0);
}

# ensure output stat dump file names
unless ($opt_read_length_histogram) {
    $opt_read_length_histogram = $fastq->basename . '.rlh.dat';
}

unless ($opt_read_mismatch_count_histogram) {
    $opt_read_mismatch_count_histogram = $fastq->basename . '.rcmh.dat';
}

# ensure that we are running on a 64-bit machine (for fqgrep)
check_architecture();

print "Processing original fastq file: $fastq \n";

my $stats = trim(
    input      => $fastq,
    mismatches => \@opt_mismatches,
    adaptor    => $opt_adaptor,
    no_trim    => $opt_no_trim,
);

my $rl_hist_file = Path::Class::File->new($opt_read_length_histogram);
my $rmc_hist_file = Path::Class::File->new($opt_read_mismatch_count_histogram);
my $overall_filter_count_file = Path::Class::File->new($fastq->basename . '.cnt');

dump_stats(
    stats                         => $stats,
    mismatches                    => \@opt_mismatches,
    read_length_histogram         => $rl_hist_file,
#    read_mismatch_count_histogram => $rmc_hist_file,
    overall_filter_count          => $overall_filter_count_file,
    fastq                         => $fastq
);

unless ($opt_no_plot) {
    plot(
        mismatches                    => \@opt_mismatches,
        read_length_histogram         => $rl_hist_file,
        overall_filter_count          => $overall_filter_count_file,
    );
}

exit(0);

# S U B R O U T I N E S #######################################################
sub trim {
    my %args = @_;
    my ($input, $mismatches, $adaptor, $no_trim) =
      @args{'input', 'mismatches', 'adaptor', 'no_trim'};
    
    my %all_stats;
    for my $mismatch (@{$mismatches}) {
        print "Working on mismatch : $mismatch", "\n";
        my $stats = run_fqgrep(
            mismatch => $mismatch,
            pattern  => $adaptor,
            fastq    => $input,
            no_trim  => $no_trim
        );
        $all_stats{$mismatch} = $stats;
    }

    return \%all_stats;
}

sub run_fqgrep {
    my %args = @_;
    my ($mismatch, $pattern, $fastq, $no_trim) =
      @args{'mismatch', 'pattern', 'fastq', 'no_trim'};

    my $cmd;
    if ($mismatch == 0) {
        $cmd = qq($fqgrep -e -r -p '$pattern' $fastq);
    }
    else {
        $cmd = qq($fqgrep -r -m $mismatch -p '$pattern' $fastq);
    }

    print $cmd, "\n";

    my %read_length_histogram;
    my %mismatch_count_histogram;
    my $match_counter = 0;

    my ($trim_file, $trim_fh, $omit_file, $omit_fh, $utrim_file);
    my %trimmed_reads;
    unless ($no_trim) {
#        my $dir = $fastq->dir;
        my ($file_name, $dir_path, $extension) =
          File::Basename::fileparse($fastq->basename, qr/\.[^.]*/);

        $trim_file = 
            Path::Class::File->new($file_name . "${extension}.${mismatch}.trim");
        $omit_file = 
            Path::Class::File->new($file_name . "${extension}.${mismatch}.omit");
        $utrim_file = 
            Path::Class::File->new($file_name . "${extension}.${mismatch}.utrim");

        $trim_fh  = $trim_file->openw;
        $omit_fh  = $omit_file->openw;
    }

    my $t0 = [gettimeofday];

    open(my $fh, '-|', $cmd)
      or die "[err] Cannot start fqgrep : $!\n";

    while (my $line = <$fh>) {
        next if ($match_counter == 0 && $line =~ /^read name/);
        $match_counter++;
        chomp($line);
        my @cols = split(/\t/, $line);
        my $num_mismatches  = $cols[1];
        my $match_start_pos = $cols[5];
#        my $match_string    = $cols[7];
        $mismatch_count_histogram{$num_mismatches}++;
        $read_length_histogram{$match_start_pos}++;

        unless($no_trim) {
            my $read_name = $cols[0];
            $trimmed_reads{$read_name} = 1;
            if ($match_start_pos == 0) {
                print $omit_fh '@', $read_name, "\n";
                print $omit_fh $cols[8], "\n";
                print $omit_fh '+', "\n";
                print $omit_fh $cols[9], "\n";
            }
            else {
                print $trim_fh '@', $read_name, "\n";
                print $trim_fh substr($cols[8], 0, $match_start_pos), "\n";
                print $trim_fh '+', "\n";
                print $trim_fh substr($cols[9], 0, $match_start_pos), "\n";
            }
        }
    }

    my $elapsed = tv_interval ( $t0 );
    print "Elapsed fqgrep processing time: $elapsed secs\n";
    close($fh);

    unless ($no_trim) {
        print 'Dumping out the un-trimmed reads into ', $utrim_file->basename,
              "\n";
        $t0 = [gettimeofday];
        dump_untrimmed_reads(
            input      => $fastq,
            output     => $utrim_file,
            trim_index => \%trimmed_reads
        );
        $elapsed = tv_interval($t0);
        print "Elapsed un-trimmed reads dumping processing time: ",
              $elapsed, " secs\n";
          
        close($trim_fh);
        close($omit_fh);
    }

    my %stats = (
        match_counter            => $match_counter,
        read_length_histogram    => \%read_length_histogram,
        mismatch_count_histogram => \%mismatch_count_histogram
    );

    return \%stats;
}

sub dump_untrimmed_reads {
    my %args = @_;
    my ($input, $output, $trim_index) = @args{'input', 'output', 'trim_index'};

    my $ifh = $input->openr;
    my $ofh = $output->openw;

    while (my $name = <$ifh>) {
        my ($read_name) = $name =~ /^@(\S+)/;
        my $sequence = <$ifh>;
        my $comment  = <$ifh>;
        my $quality  = <$ifh>;

        next if exists $trim_index->{$read_name};

        print $ofh $name;
        print $ofh $sequence;
        print $ofh $comment;
        print $ofh $quality;
    }
    
    close($ofh);
}

sub dump_stats {
    my %args = @_;
    my ($stats, $mismatches, $rlh_file, $rmch_file, $overall_cnt_file, $fastq)
      = @args{
        'stats',                 'mismatches',
        'read_length_histogram', 'read_mismatch_count_histogram',
        'overall_filter_count',  'fastq'
      };

    my $total_read_count = `wc -l $fastq`;
    chomp($total_read_count);
    ($total_read_count) = $total_read_count =~ /(\d+)\s*(\S+)/;
    $total_read_count = int($total_read_count/4);

    # Dump out overall_filter_count file stats ($fastq.cnt)
    my $fh = $overall_cnt_file->openw;
    print $fh join("\t", '#Mismatch', 'cumulative_filtered_read_count', 'filtered_reads', 'total_reads', 'pct'), "\n";
    for my $mismatch (sort { $a <=> $b } @{$mismatches}) {
        my $cumulative_filter_count = $stats->{$mismatch}->{'match_counter'};
        my $filter_pct =
          $total_read_count
          ? ($cumulative_filter_count / $total_read_count) * 100
          : 'na';

        my $filter_count =
          $stats->{$mismatch}->{'mismatch_count_histogram'}->{$mismatch};
        $filter_count = 0 unless $filter_count;

        print $fh join("\t",
            $mismatch, 
            $cumulative_filter_count,
            $filter_count,
            $total_read_count, 
            sprintf("%.2f", $filter_pct)), "\n";
    }
    close($fh);

#    # Dump out read mismatch count histogram file stats ($fastq.rmch.dat)
#    my $max_mismatch = (sort { $b <=> $a } @{$mismatches})[0];
#    my $rmc_hist = $stats->{$max_mismatch}->{'mismatch_count_histogram'};
#    $fh = $rmch_file->openw;
#    print $fh join("\t", '#Mismatch', 'filtered read count'), "\n";
#    for my $m (0 .. $max_mismatch) {
#        my $count = exists $rmc_hist->{$m} ? $rmc_hist->{$m} : 0;
#        print $fh join("\t", $m, $count), "\n";
#    }
#    close($fh);

    # Dump out read length count histogram file stats ($fastq.rlh.dat)
    my $max_length = max_filtered_read_length($stats);
    $fh = $rlh_file->openw;
    for my $m (sort { $a <=> $b } @{$mismatches}) {
        print $fh "# Mismatch : $m", "\n";
        my $rlh = $stats->{$m}->{'read_length_histogram'};
        print $fh join("\t", 'Read Length', 'Read Count'), "\n";
        for my $l (0 .. $max_length) {
            my $count = exists $rlh->{$l} ? $rlh->{$l} : 0;
            print $fh join("\t", $l, $count), "\n";
        }
        print $fh "\n\n"; # separate out the data set
    }
    close($fh);

    return 1;
}

sub max_filtered_read_length {
    my $rl_hist = shift;

    my $overall_max_read_length = 0;
    for my $mismatch (sort {$a <=> $b} keys %{$rl_hist}) {
        my $max_set_read_length =
          (sort { $b <=> $a }
              keys %{$rl_hist->{$mismatch}->{'read_length_histogram'}})[0];
        if ($overall_max_read_length < $max_set_read_length) {
            $overall_max_read_length = $max_set_read_length;
        }
    }

    return $overall_max_read_length;
}

sub plot {
}

sub check_architecture {
    my $arch = POSIX::uname;
    unless ($arch =~ /64/) {
        my $host = Sys::Hostname::hostname;
        die "[err] Sorry pard'ner but you got to be on a ",
            "64-bit machine to confirm this pse!", "\n",
            "($host) is a $arch machine!\n";
    }
    return 1;
}

__END__

# S U B R O U T I N E S #######################################################
=head1 NAME

fqgrep-trim.pl - a "trimmer" wrapped around fqgrep

=head1 SYNOPSIS

fqgrep-trim.pl [OPTIONS] --adaptor [adaptor regex] --fastq [fastq/a file]

=head1 OPTIONS

=over 4

=item B<--fastq=FASTQ>

The fastq/a file of interest to process.

=item B<--adaptor=ADAPTOR>

The adpator of interest to trim for.

=item B<--no-trim>

Skip producing the actual trimming files and just produce statistics files

=item B<--no-plot>

Skip producing statistics plot file.

=item B<--mismatches=MISMATCH-VALUES>

The number of allowable mismatch levels to filter adaptors for.

=item B<--read-length-histogram=FILENAME>

Name of the read length histogram file (default: <fastq-file>.rlh.dat)

=item B<--read-mismatch-count-histogram=FILENAME>

Name of the read mismatch count histogram file (default: <fastq-file>.rcmh.dat)

=back

=head1 DESCRIPTION

This script a wrapper around the fqgrep tool to perform "trimming".

=head1 EXAMPLES

=head1 NOTES

=head1 DEPENDENCIES

    POSIX
    YAML::Syck
    Sys::Hostname
    Time::HiRes
    Path::Class
    File::Basename
    Getopt::Long

    gnuplot

=head1 ORIGINAL AUTHOR

Indraniel Das <idas@wustl.edu> 2010-11-06

=head1 BUGS

Please report bugs to LIMS <lims@genome.wustl.edu>

=cut

gnuplot> plot "./read-lengths.dat2" index 0:11 using 1:2 w linespoints
gnuplot> plot "./read-lengths.dat2" index 0:11 using 1:3 w linespoints
gnuplot> plot "./mismatch-cnt.dat2" index 0 using 1:3 w linespoints

