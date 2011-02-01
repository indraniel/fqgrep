#!/usr/bin/env perl

# L I C E N S E ###############################################################
#    Copyright (C) 2010, 2011 Indraniel Das <indraniel@gmail.com>
#                             and Washington University in St. Louis
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, see <http://www.gnu.org/licenses/>

# P R A G M A S ###############################################################
use warnings;
use strict;

# M O D U L E S ###############################################################
#use YAML::Syck;
use Time::HiRes qw(gettimeofday tv_interval);
use Path::Class;
use File::Basename;
use Pod::Usage;
use Getopt::Long;

# M A I N #####################################################################

$| = 1; # enable AUTOFLUSH mode

# Setup Default options
my ($opt_fastq, @opt_mismatches, $opt_adaptor, $opt_help);
my $opt_read_length_histogram;
my $opt_read_count_histogram;
my $fqgrep = qx(which fqgrep) || undef; # assuming fqgrep is somewhere in $PATH
my $opt_trim = "left";
my $opt_format = 'FASTQ';

# parse options
my $result = GetOptions(
    "input=s"                        => \$opt_fastq,
    "mismatches=i{,}"                => \@opt_mismatches,
    "adaptor=s"                      => \$opt_adaptor,
    "read-length-histogram"          => \$opt_read_length_histogram,
    "read-count-histogram"           => \$opt_read_count_histogram,
    "fqgrep=s"                       => \$fqgrep,
    "trim"                           => \$opt_trim,
    "help"                           => \$opt_help
);

if ($opt_help) {
    pod2usage(-exitval => 0, -verbose => 0);
}

# ensure that fqgrep is found and runnable
unless ($fqgrep) {
    die "[err] Could not find the location to fqgrep!\n";
}
chomp($fqgrep);
$fqgrep = Path::Class::File->new($fqgrep);

unless (-e $fqgrep) {
    die "[err] Did not find fqgrep at : $fqgrep !\n";
}

unless (-x $fqgrep) {
    die "[err] fqgrep ($fqgrep) is not executable!\n";
}

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

# ensure that there is a trim type
unless ($opt_trim =~/(left|right)/i) {
    die "[err] Need to provide a correct --trim type: 'left' or 'right' !\n";
}

# ensure a mismatch level
unless (@opt_mismatches) {
    @opt_mismatches = qw(0);
}

# ensure there is an appropriate output format
unless ($opt_format =~ /^FAST(A|Q)$/) {
    die "[err] The format option ($opt_format) is not valid. ",
        "It should be either 'FASTA' or 'FASTQ'!\n";
}

# ensure output stat dump file names
unless ($opt_read_length_histogram) {
    $opt_read_length_histogram = $fastq->basename . '.read_length_histogram.dat';
}

unless ($opt_read_count_histogram) {
    $opt_read_count_histogram =
      $fastq->basename . '.read_count_histogram.dat';
}

print "Processing original fastq file: $fastq \n";

# Performing trimming on the file and collect relevant statistics
my $stats = trim_file(
    input      => $fastq,
    mismatches => \@opt_mismatches,
    adaptor    => $opt_adaptor,
    trim_type  => $opt_trim,
    format     => $opt_format
);

# setup the files that the collected statistics will be filled into
my $rl_hist_file = Path::Class::File->new($opt_read_length_histogram);
my $rmc_hist_file =
  Path::Class::File->new($opt_read_count_histogram);
my $overall_filter_count_file =
  Path::Class::File->new($fastq->basename . '.cnt');

dump_stats(
    stats                         => $stats,
    mismatches                    => \@opt_mismatches,
    read_length_histogram         => $rl_hist_file,
#    read_count_histogram          => $rmc_hist_file,
    overall_filter_count          => $overall_filter_count_file,
    fastq                         => $fastq
);

exit(0);

# S U B R O U T I N E S #######################################################
sub trim_file {
    my %args = @_;
    my ($input, $mismatches, $adaptor, $trim_type, $format) =
      @args{'input', 'mismatches', 'adaptor', 'trim_type', 'format'};

    my $adaptor_pattern = adjust_adaptor_regex($adaptor, $trim_type);
    
    my %all_stats;
    for my $mismatch (@{$mismatches}) {
        print "Working on mismatch : $mismatch", "\n";
        my $stats = run_fqgrep_and_trim_reads(
            mismatch  => $mismatch,
            pattern   => $adaptor_pattern,
            fastq     => $input,
            trim_type => $trim_type,
            format    => $format
        );
        $all_stats{$mismatch} = $stats;
    }

    return \%all_stats;
}

# ensure that the '^' or '$' anchor is placed for respective
# 'left' or 'right' end trimming.
sub adjust_adaptor_regex {
    my ($adaptor, $trim_type) = @_;

    my $adaptor_pattern = $adaptor;

    # handle left-side "start" trimming
    if ($trim_type =~ /left/i) {
        unless ($adaptor =~ /^\^/) {
            $adaptor_pattern = '^' . $adaptor;
        }
        return $adaptor_pattern;
    }

    # handle right-side "end" trimming
    unless ($adaptor =~ /\$$/) {
        $adaptor_pattern = $adaptor . '$';
    }
    return $adaptor_pattern;
}

sub run_fqgrep_and_trim_reads {
    my %args = @_;
    my ($mismatch, $pattern, $fastq, $trim_type, $format) =
      @args{'mismatch', 'pattern', 'fastq', 'trim_type', 'format'};

    my $cmd;
    if ($mismatch == 0) {
        $cmd = qq($fqgrep -e -r -p '$pattern' $fastq);
    }
    else {
        $cmd = qq($fqgrep -r -m $mismatch -p '$pattern' $fastq);
    }

    print $cmd, "\n";

    # trim statistics recording variables
    my %trimmed_reads;
    my %read_length_histogram;
    my %mismatch_count_histogram;
    my ($match_counter, $omit_counter, $trim_counter) = (0, 0, 0);

    # setup output trimming files
    my ($file_name, $dir_path, $extension) =
      File::Basename::fileparse($fastq->basename, qr/\.[^.]*/);

    my $trim_file =
        Path::Class::File->new($file_name . "${extension}.${mismatch}.trim");
    my $omit_file =
        Path::Class::File->new($file_name . "${extension}.${mismatch}.omit");
    my $utrim_file =
        Path::Class::File->new($file_name . "${extension}.${mismatch}.utrim");

    my $trim_fh  = $trim_file->openw;
    my $omit_fh  = $omit_file->openw;

    # invoke the fqgrep command and read its output line-by-line
    open(my $fh, '-|', $cmd)
      or die "[err] Cannot start fqgrep : $!\n";

    # go through each of the fqgrep outputs and attempt to
    # append to the trim and/or omit files
    my $t0 = [gettimeofday]; # start stopwatch
    while (my $line = <$fh>) {
        next if ($match_counter == 0 && $line =~ /^read name/);
        $match_counter++;
        chomp($line);

        # parse out the fqgrep output line for the relevant info
        my @cols = split(/\t/, $line);
        my $read_name       = $cols[0];
        my $num_mismatches  = $cols[1];
        my $match_start_pos = $cols[5];
        my $match_end_pos   = $cols[6];
#        my $match_string    = $cols[7];
        my $sequence        = $cols[8];
        my $quality         = $cols[9] ? $cols[9] : '';
        $mismatch_count_histogram{$num_mismatches}++;

        my ($trim_or_omit_flag, $trimmed_read_length) = trim_or_omit_read(
            trim_type       => $trim_type,
            read_name       => $read_name,
            sequence        => $sequence,
            quality         => $quality,
            match_start_pos => $match_start_pos,
            match_end_pos   => $match_end_pos,
            omit_fh         => $omit_fh,
            trim_fh         => $trim_fh,
            format          => $format
        );

        $read_length_histogram{$trimmed_read_length}++;
        $trimmed_reads{$read_name} = 1;
        $omit_counter++ if ($trim_or_omit_flag eq 'omit');
        $trim_counter++ if ($trim_or_omit_flag eq 'trim');
    }
    my $elapsed = tv_interval ( $t0 ); # end stopwatch
    print "Elapsed fqgrep processing time: $elapsed secs\n";

    close($fh);

    # all the remaining reads unaccounted for by fqgrep are keepable
    # reads.  Place them into the un-trimmed file.
    print 'Dumping out the un-trimmed reads into ', $utrim_file->basename,
          "\n";
    $t0 = [gettimeofday];
    dump_untrimmed_reads(
        input      => $fastq,
        output     => $utrim_file,
        trim_index => \%trimmed_reads,
        format     => $format
    );
    $elapsed = tv_interval($t0);
    print "Elapsed un-trimmed reads dumping processing time: ",
          $elapsed, " secs\n";

    close($trim_fh);
    close($omit_fh);

    my %stats = (
        match_counter            => $match_counter,
        omit_counter             => $omit_counter,
        trim_counter             => $trim_counter,
        read_length_histogram    => \%read_length_histogram,
        mismatch_count_histogram => \%mismatch_count_histogram
    );

    return \%stats;
}

sub trim_or_omit_read {
    my %args = @_;

    if ($args{'type'} =~ /left/i) {
        return trim_left_adaptor_on_read(%args);
    }
    else {
        return trim_right_adaptor_on_read(%args);
    }
}

# perform 'left' side trim (adaptor removal)
# <adaptor> followed by <read>
sub trim_left_adaptor_on_read {
    my %args = @_;

    my ($format,        $read_name, $sequence, $quality,
        $match_end_pos, $trim_type, $omit_fh,  $trim_fh)
      = @args{
        'format',        'read_name', 'sequence', 'quality',
        'match_end_pos', 'trim_type', 'omit_fh',  'trim_fh',
      };

    # Found the adaptor at the end of the read. That is bad!
    my $trimmed_read_length = 0;
    if ($match_end_pos == (length($sequence) - 1) ) {
        print_read(
            fh        => $omit_fh,
            read_name => $read_name,
            sequence  => $sequence,
            quality   => $quality,
            format    => $format
        );
        return ('omit', $trimmed_read_length);
    }

    # keep only the read sequence after the match end position
    my $trimmed_read    = substr($sequence, $match_end_pos);
    my $trimmed_quality = $quality ? substr($quality,  $match_end_pos) : '';
    $trimmed_read_length = length($trimmed_read);

    print_read(
        fh        => $trim_fh,
        read_name => $read_name,
        sequence  => $trimmed_read,
        quality   => $trimmed_quality,
        format    => $format
    );

    return ('trim', $trimmed_read_length);
}

# perform 'right' side trim (adaptor removal)
# <read> followed by <adaptor>
sub trim_right_adaptor_on_read {
    my %args = @_;

    my ($format,           $read_name, $sequence, $quality,
        $match_start_pos,  $trim_type, $omit_fh, $trim_fh)
      = @args{
        'format',          'read_name', 'sequence', 'quality',
        'match_start_pos', 'trim_type', 'omit_fh',  'trim_fh'
      };

    # Found the adaptor at the beginning of the read. That is bad!
    my $trimmed_read_length = 0;
    if ($match_start_pos == 0) {
        print_read(
            omit_fh   => $omit_fh,
            read_name => $read_name,
            sequence  => $sequence,
            quality   => $quality,
            format    => $format
        );
        return ('omit', $trimmed_read_length);
    }

    # keep only the read sequence before the match start position
    my $trimmed_read = substr($sequence, 0, $match_start_pos);
    my $trimmed_quality =
      $quality ? substr($quality, 0, $match_start_pos) : '';
    $trimmed_read_length = length($trimmed_read);

    print_read(
        fh        => $trim_fh,
        read_name => $read_name,
        sequence  => $trimmed_read,
        quality   => $trimmed_quality,
        format    => $format
    );

    return ('trim', $trimmed_read_length);
}

sub print_read {
    my %args = @_;

    my ($read_name, $sequence, $quality, $fh, $format) =
        @args{'read_name', 'sequence', 'quality', 'fh', 'format'};

    if ($format eq 'FASTQ') {
        print $fh '@', $read_name, "\n";
        print $fh $sequence, "\n";
        print $fh '+', "\n";
        print $fh $quality, "\n";
    }
    else {
        print $fh '>', $read_name, "\n";
        print $fh $sequence, "\n";
    }
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
        'read_length_histogram', 'read_count_histogram',
        'overall_filter_count',  'fastq'
      };

    my $total_read_count = `wc -l $fastq`;
    chomp($total_read_count);
    ($total_read_count) = $total_read_count =~ /(\d+)\s*(\S+)/;

    # ascertain if we are working with an input FASTQ or FASTA file
    my $fh = $fastq->openr;
    my $header = <$fh>;
    chomp($header);
    close($fh);
    my $file_type = 'FASTA';
    if ($header = /^@/) {
        $file_type = 'FASTQ';
    }

    if ($file_type eq 'FASTQ') {
        $total_read_count = int($total_read_count/4);
    }
    else {
        $total_read_count = int($total_read_count/2);
    }

    # Dump out overall_filter_count file stats ($fastq.cnt)
    $fh = $overall_cnt_file->openw;
    print $fh join("\t", 'Mismatch', 'cumulative_filtered_read_count', 'filtered_reads', 'omitted_reads', 'total_reads', 'pct'), "\n";
    for my $mismatch (sort { $a <=> $b } @{$mismatches}) {
        my $cumulative_filter_count = $stats->{$mismatch}->{'match_counter'};
        my $omit_count = $stats->{$mismatch}->{'omit_counter'};
        my $filter_pct =
          $total_read_count
          ? ($cumulative_filter_count / $total_read_count) * 100
          : 'NA';

        my $filter_count =
          $stats->{$mismatch}->{'mismatch_count_histogram'}->{$mismatch};
        $filter_count = 0 unless $filter_count;

        print $fh join("\t",
            $mismatch, 
            $cumulative_filter_count,
            $filter_count,
            $omit_count,
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
        my $rlh = $stats->{$m}->{'read_length_histogram'};
        print $fh join("\t", 'Mismatch', 'Read Length', 'Read Count'), "\n";
        for my $l (0 .. $max_length) {
            my $count = exists $rlh->{$l} ? $rlh->{$l} : 0;
            print $fh join("\t", $m, $l, $count), "\n";
        }
        print $fh "\n";
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

__END__

# D O C U M E N T A T I O N ###################################################
=head1 NAME

fqgrep-trim.pl - an adaptive FASTQ/FASTA "trimmer" based upon fqgrep

=head1 SYNOPSIS

fqgrep-trim.pl [OPTIONS] --adaptor [adaptor regex] --input [fastq file]

  Options:
    --help                  brief help message
    --input                 input FASTQ or FASTA file
    --adaptor               the adaptor pattern to match for
                            (can be a regular expression)
    --mismatches            adaptor pattern mismatch level threshold
    --read-length-histogram filename for the read length histogram stats
    --read-count-histogram  filename for the read count histogram stats
    --trim                  adaptor trimming direction "left" or "right"
    --format                format of the untrimmed, trimmed, and omit files
                            ("FASTQ" or "FASTA").

=head1 OPTIONS

=over 4

=item B<--input=FASTQ|FASTA>

The FASTQ or FASTA file of interest to process B<(required parameter)>.

=item B<--adaptor=ADAPTOR>

The adpator of interest to trim out B<(required parameter)>.

=item B<--mismatches=MISMATCH-VALUES>

The number of allowable mismatches to filter adaptors for (default
0). This can be multiple integer values separated by a single space.

=item B<--read-length-histogram=FILENAME>

Name of the read length histogram statistics file (default:
<fastq-file>.rlh.dat)

=item B<--read-count-histogram=FILENAME>

Name of the read mismatch count histogram statistics file (default:
<fastq-file>.rcmh.dat)

=item B<--trim="left" or "right">

Performs "left"-end adaptor trimming (e.g. when the desired adaptor
is found towards the beginning of the read), or "right"-end adaptor
trimming (e.g. when the desired adaptor is found towards the end
of the read).  B<(default: "left")>

=item B<--format=FASTQ|FASTA>

The format type to output the omit, trimmed, and untrimmed sequence
read data. Can be either 'FASTQ' or 'FASTA'. (default: 'FASTQ')

=back

=head1 DESCRIPTION

This script a wrapper around the fqgrep tool to perform either
"left"-end or "right"-end adaptive adaptor trimming, and notes relevant
statistics.

Upon invocation, the script will split the reads from the input file
into 3 (FASTA or FASTQ formatted) file categories in the current working
directory. These 3 categories are:

  * an omit file -- <input-filename>.omit

    These represent reads that were improper, and can be omitted. By
    improper, the adaptor was found at the opposite end of the read from
    where it was originally intended.

  * a trimmed file -- <input-filename>.trim

    Reads from the original input-file that had its adaptor identified
    and excised.

  * an untrimmed file -- <input-filename>.utrim

    Reads from the original input-file that did not have have an identifiable
    adaptor, and have remained unaltered from the orignal file.

In addition 2 trimming statistic files will be produced in the
current working directory as well:

  * <input-filename>.read_length_histogram.dat

    Contains a table of the read length counts for a given mismatch level.

  * <input-filename>.read_count_histogram.dat

    Contains a table of the number of counts of reads trimmed and omitted
    for a given mismatch level.

=head1 EXAMPLES

=head1 NOTES

This script performs a simplistic form of read trimming based on each
individual read in a given FASTQ and/or FASTA file. It uses fqgrep to
identify a start or end adaptor for each indivdual read in the supplied
input file, and subsequently excises the adaptor from the read. It
does not make any prior assumptions about the input FASTQ/FASTA file
and relate the consequence of a particular read trimming to other reads
in the input file. For example, when trimming a read from a FASTQ that
contains Illumina Paired-End sequence reads, the script does not account
for and appropriately segregate the the opposite read 'mate-pair' that
may also be affected by the trim.

Feel free to use this script as a starting point for more advanced
trimming thoughts.

=head1 DEPENDENCIES

The following perl modules are dependent upon this script:

    YAML::Syck
    Time::HiRes
    Path::Class
    File::Basename
    Getopt::Long

=head1 ORIGINAL AUTHOR

Indraniel Das <indraniel at gmail dot com> 2010-01-30

=cut
