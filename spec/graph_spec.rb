# Copyright 2016 Ryan Moore
# Contact: moorer@udel.edu
#
# This file is part of Shaq.
#
# Shaq is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Shaq is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Shaq.  If not, see <http://www.gnu.org/licenses/>.

require_relative "../graph"

describe Graph do
  let(:outdir) { File.join File.dirname(__FILE__),
                           "test_files",
                           "output" }
  let(:fasta_f) { File.join outdir, "output.fa" }
  let(:contig_connections_f) { File.join outdir, "connections.txt" }

  let(:ksize) { 3 }
  let(:graph) { Graph.new ksize }

  let(:header) { "read 1" }
  let(:read) { "aaataaac" }

  let(:db_graph) {
    { "aaa" => Set.new(["aat", "aac"]),
      "aat" => Set.new(["ata"]),
      "ata" => Set.new(["taa"]),
      "taa" => Set.new(["aaa"])}
  }

  let(:kmer_to_read) {
    { "aaa" => [["read 1", 0], ["read 1", 4]],
      "aat" => [["read 1", 1]],
      "ata" => [["read 1", 2]],
      "taa" => [["read 1", 3]],
      "aac" => [["read 1", 5]], }
  }

  let(:kmer_abundance) {
    { "aaa" => 2,
      "aat" => 1,
      "ata" => 1,
      "taa" => 1,
      "aac" => 1,}
  }

  let(:contigs) {
    ["gaaa", "aaa", "aataaa", "aacac"]
  }

  let(:fasta_contents) {
    %w[>contig_1 >contig_2 >contig_3 >contig_4].zip(contigs).flatten.join("\n")
  }

  let(:contig_connections) {
    { "aaa" => Set.new([0, 1, 2, 3]) }
  }

  let(:cc_contents) {
    [
      ["contig_0", "contig_1", "aaa"].join("\t"),
      ["contig_0", "contig_2", "aaa"].join("\t"),
      ["contig_0", "contig_3", "aaa"].join("\t"),
      ["contig_1", "contig_2", "aaa"].join("\t"),
      ["contig_1", "contig_3", "aaa"].join("\t"),
      ["contig_2", "contig_3", "aaa"].join("\t"),
    ].join "\n"
  }


  describe "#add" do
    it "adds the read to the graph" do
      graph.add header, read

      expect(graph.graph).to eq db_graph
    end

    it "adds the header and pos to the kmer_to_read hash" do
      graph.add header, read

      expect(graph.kmer_to_read).to eq kmer_to_read
    end
  end

  describe "#kmer_abundance" do
    it "returns a coll with kmer abundance" do
      graph.add header, read

      expect(graph.kmer_abundance).to eq kmer_abundance
    end
  end

  describe "#contigs" do
    it "returns the contigs" do
      graph.add header, "gaaaataaaacac"
      graph.make_contigs

      # NOTE if you just let this walk in order, it will spit back the
      # whole, thing, however, when you have actual reads, the order
      # isn't guaranteed to be correct, and you can get an incorrect
      # contig if you don't break at every branch point when the order
      # isn't perfect and the repeats are longer than ksize
      expect(graph.contigs).to match_array contigs
    end

    # TODO I have a bad feeling about the code for this...
    it "tracks the contigs connected with branch point kmers" do
      graph.add header, "gaaaataaaacac"
      graph.make_contigs

      expect(graph.contig_connections).to eq contig_connections
    end
  end

  describe "#print_fasta" do
    it "prints contigs to fasta file" do
      graph.add header, "gaaaataaaacac"
      graph.make_contigs
      graph.print_fasta fasta_f

      data = File.open(fasta_f).read.chomp

      expect(data).to eq fasta_contents

      FileUtils.rm fasta_f
      FileUtils.rmdir outdir
    end
  end

  describe "#print_contig_connections" do
    it "prints contigs connections to file" do
      graph.add header, "gaaaataaaacac"
      graph.make_contigs
      graph.print_contig_connections contig_connections_f

      data = File.open(contig_connections_f).read.chomp

      expect(data).to eq cc_contents

      FileUtils.rm contig_connections_f
      FileUtils.rmdir outdir
    end
  end
end
