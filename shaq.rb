require "parse_fasta"
require_relative "graph"

KSIZE = 9
graph = Graph.new KSIZE

SeqFile.open(ARGV.first).each_record_fast do |head, seq|
  graph.add head, seq
end

graph.make_contigs

outdir  = File.dirname(ARGV.first)
ext     = File.extname(ARGV.first)
base    = File.basename ARGV.first, ext

fasta_f = File.join outdir, "shaq.#{base}.contigs.fa"
graph_f = File.join outdir, "shaq.#{base}.contig_graph.txt"

graph.print_fasta fasta_f
warn "LOG -- contigs: #{fasta_f}"

graph.print_contig_connections graph_f
warn "LOG -- graph:   #{graph_f}"
