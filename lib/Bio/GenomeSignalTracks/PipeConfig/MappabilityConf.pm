package Bio::GenomeSignalTracks::PipeConfig::MappabilityConf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{ $self->SUPER::pipeline_wide_parameters },
        samtools         => $self->o('samtools'),
        bowtie           => $self->o('bowtie'),
        bedtools         => $self->o('bedtools'),
        fasta_suffix     => $self->o('fasta_suffix'),
        bedGraphToBigWig => $self->o('bedGraphToBigWig'),

    };
}

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() }
        ,    # inherit other stuff from the base class

        # name used by the beekeeper to prefix job names on the farm
        'pipeline_name' => 'mappability',

        # suffix used when looking for fasta files
        'fasta_suffix' => 'fa',

        # maximum number of kmers to write into a batch for bowtie
        'kmer_split_limit' => 50000000,

        # gzip kmer files
        'gzip_kmer_files' => 1,
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {
            -logic_name => 'start',
            -module     => 'Bio::GenomeSignalTracks::Process::MappabilityPreFlightChecks',
            -flow_into => {
                '1' => {
                    'kmer_factory' => {
                        output_dir => '#output_dir#',
                        kmer_sizes => '#kmer_sizes#',
                        fasta_dir  => '#fasta_dir#',
                        index_dir  => '#index_dir#',
                        index_name => '#index_name#',
                        chrom_list => '#chrom_list#',
                    },
                },
            },
        },
        {
            -logic_name => 'kmer_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                inputlist    => '#expr( [eval(#kmer_sizes#)] )expr#',
                column_names => ['kmer_size'],
            },
            -flow_into => {
                '2' => {
                    'make_output_dir' => {
                        output_dir => '#output_dir#',
                        kmer_size  => '#kmer_size#',
                        fasta_dir  => '#fasta_dir#',
                        index_dir  => '#index_dir#',
                        index_name => '#index_name#',
                        chrom_list => '#chrom_list#',
                    },
                },
            },
        },
        {
            -logic_name => 'make_output_dir',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                kmer_out_dir => '#output_dir#/#index_name#/k#kmer_size#',
                cmd          => 'mkdir -p #kmer_out_dir#',
            },
            -flow_into => {
                '1->A' => {
                    'fasta_factory' => {
                        fasta_dir  => '#fasta_dir#',
                        output_dir => '#kmer_out_dir#',
                        kmer_size  => '#kmer_size#',
                        index_dir  => '#index_dir#',
                        index_name => '#index_name#',
                    }
                },
                'A->1' => {
                    'bam_merge' => {
                        output_dir => '#output_dir#',
                        name       => '#index_name#.k#kmer_size#',
                        kmer_size  => '#kmer_size#',
                        chrom_list => '#chrom_list#',
                    }
                }
            }
        },
        {
            -logic_name => 'fasta_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                inputcmd => 'find #fasta_dir# -type f -name "*.#fasta_suffix#"',
                column_names => ['filename'],
            },
            -flow_into => {
                '2' => {
                    'fasta_kmers_factory' => {
                        fasta_file => '#filename#',
                        fasta_dir  => '#fasta_dir#',
                        output_dir => '#output_dir#',
                        kmer_size  => '#kmer_size#',
                        index_dir  => '#index_dir#',
                        index_name => '#index_name#',
                    }
                },
            }
        },
        {
            -logic_name => 'fasta_kmers_factory',
            -module => 'Bio::GenomeSignalTracks::Process::FastaKmerSplitter',
            -parameters => {
                gzip        => $self->o('gzip_kmer_files'),
                split_limit => $self->o('kmer_split_limit'),
            },
            -rc_name   => '2Gb_job',
            -flow_into => {
                2 => {
                    'bowtie' => {
                        'kmer_file'  => '#kmer_file#',
                        'index_dir'  => '#index_dir#',
                        'index_name' => '#index_name#',
                    }
                }
              }

        },
        {
            -logic_name => 'bowtie',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                bam => '#kmer_file#.bam',
                cmd =>
'gunzip -c #kmer_file# | #bowtie# -v 0 -k 1 -m 1 -r -S #index_dir#/#index_name# - | #samtools# view -SbF4 - > #bam#',
            },
            -rc_name   => '2Gb_job',
            -flow_into => {
                1 => {
                    'sort_bam' => { 'bam'  => '#bam#' },
                    'rm_file'  => { 'file' => '#kmer_file#' }
                }
            },
        },
        {
            -logic_name => 'sort_bam',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                sorted_bam => '#bam#.sorted.bam',
                cmd =>
                  '#samtools# sort -O bam -o #sorted_bam# -T #bam#.tmp #bam#',
            },
            -rc_name   => '2Gb_job',
            -flow_into => {
                1 => {
                    ':////accu?sorted_bam=[]' => {},
                    'rm_file'                 => { file => '#bam#' }
                }
            },
        },
        {
            -logic_name => 'rm_file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { cmd => 'rm #file#', },
        },

        {
            -logic_name => 'bam_merge',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name    => '2Gb_job',
            -parameters => {
                bam => '#output_dir#/#name#.bam',
                cmd =>
'#samtools# merge -f #bam# #expr(join(" ",@{#sorted_bam#}))expr#',
            },
            -flow_into => {
                1 => {
                    'bam2bg' => {
                        bam        => '#bam#',
                        kmer_size  => '#kmer_size#',
                        name       => '#name#',
                        output_dir => '#output_dir#',
                        chrom_list => '#chrom_list#',
                    },
                    'rm_file' =>
                      { 'file' => '#expr(join(" ",@{#sorted_bam#}))expr#' }
                }
            }
        },
        {
            -logic_name => 'bam2bg',
            -rc_name    => '3Gb_job',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                scaling_factor => '#expr( 1/#kmer_size# )expr#',
                bedgraph       => '#output_dir#/#name#.bg',
                cmd =>
'#bedtools# genomecov -ibam #bam# -bg -scale #scaling_factor# > #bedgraph#',
            },
            -flow_into => {
                1 => {
                    'rm_file'          => { file => '#bam#' },
                    'bedGraphToBigWig' => {
                        name       => '#name#',
                        bedgraph   => '#bedgraph#',
                        output_dir => '#output_dir#',
                        chrom_list => '#chrom_list#',
                    }
                }
            }
        },
        {
            -logic_name => 'bedGraphToBigWig',
            -rc_name    => '2Gb_job',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                bigwig => '#output_dir#/#name#.bw',
                cmd => '#bedGraphToBigWig# #bedgraph# #chrom_list# #bigwig# ',
            },
            -flow_into => { 1 => { 'rm_file' => { file => '#bedgraph#' } } }
        },
    ];
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{ $self->SUPER::resource_classes }
        ,    # inherit 'default' from the parent class

        'default' =>
          { 'LSF' => '-M100   -R"select[mem>100]   rusage[mem=100]"' }
        ,    # to make sure it fails similarly on both farms
        '200Mb_job' =>
          { 'LSF' => '-M200   -R"select[mem>200]   rusage[mem=200]"' },
        '400Mb_job' =>
          { 'LSF' => '-M400   -R"select[mem>400]   rusage[mem=400]"' },
        '1Gb_job' =>
          { 'LSF' => '-M1000  -R"select[mem>1000]  rusage[mem=1000]"' },
        '2Gb_job' =>
          { 'LSF' => '-M2000  -R"select[mem>2000]  rusage[mem=2000]"' },
        '3Gb_job' =>
          { 'LSF' => '-M3000  -R"select[mem>3000]  rusage[mem=3000]"' },
        '5Gb_job' =>
          { 'LSF' => '-M1000  -R"select[mem>5000]  rusage[mem=5000]"' },
    };
}
1;
