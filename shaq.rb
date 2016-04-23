require "parse_fasta"
require_relative "graph"

KSIZE = ARGV[0].to_i
inf = ARGV[1]
graph = Graph.new KSIZE

SeqFile.open(inf).each_record_fast do |head, seq|
  graph.add head, seq.downcase
end

graph.make_contigs

outdir  = File.dirname(inf)
ext     = File.extname(inf)
base    = File.basename inf, ext

fasta_f = File.join outdir, "shaq.#{base}.contigs.fa"
graph_f = File.join outdir, "shaq.#{base}.contig_graph.txt"

graph.print_fasta fasta_f
warn "LOG -- contigs: #{fasta_f}"

graph.print_contig_connections graph_f
warn "LOG -- graph:   #{graph_f}"
